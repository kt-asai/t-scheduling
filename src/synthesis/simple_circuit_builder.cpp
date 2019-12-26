#include <random>
#include <chrono>

#include "simple_circuit_builder.hpp"

#include "../decomposer/gaussian_decomposer.hpp"

namespace tskd {

bool SimpleCircuitBuilder::Init(const std::vector<util::xor_func>& in,
                                const std::vector<util::xor_func>& out)
{
    bool is_io_different = true;

    bits_ = std::vector<util::xor_func>(qubit_num_);
    preparation_ = std::vector<util::xor_func>(qubit_num_);
    restoration_ = std::vector<util::xor_func>(qubit_num_);

    for (int i = 0; i < qubit_num_; i++)
    {
        is_io_different &= (in[i] == out[i]);
        preparation_[i] = util::xor_func(qubit_num_ + 1, 0);
        restoration_[i] = util::xor_func(qubit_num_ + 1, 0);
        preparation_[i].set(i);
        restoration_[i].set(i);
    }

    return is_io_different;
}

void SimpleCircuitBuilder::InitBits(const std::set<int>& phase_exponent_index_set,
                                    std::unordered_map<int, int>& target_phase_map)
{
    std::set<int>::iterator ti;
    int counter = 0;
    target_phase_map.clear();
    for (ti = phase_exponent_index_set.begin(), counter = 0; counter < qubit_num_; counter++)
    {
        if (counter < static_cast<int>(phase_exponent_index_set.size()))
        {
            bits_[counter] = phase_exponent_[*ti].second;
            target_phase_map.emplace(counter, *ti);
            ti++;
        }
        else
        {
            bits_[counter] = util::xor_func(dimension_ + 1, 0);
        }
    }
}

void SimpleCircuitBuilder::Prepare(std::list<Gate>& gate_list,
                                   const std::vector<util::xor_func>& in,
                                   const int num_partition,
                                   std::unordered_map<int, int>& target_phase_map)
{
    util::ToUpperEchelon(num_partition, dimension_, bits_, &restoration_, std::vector<std::string>());
    util::FixBasis(qubit_num_, dimension_, num_partition, in, bits_, &restoration_, std::vector<std::string>());
    util::Compose(qubit_num_, preparation_, restoration_);

//    std::cout << "-- before preparation_" << std::endl;
//    for (auto&& e : preparation_)
//    {
//        std::cout << e << std::endl;
//    }
//    std::cout << "---- restoration_" << std::endl;
//    for (auto&& e : restoration_)
//    {
//        std::cout << e << std::endl;
//    }

    if (option_.change_row_order())
    {
        ChangeRowOrder(target_phase_map, preparation_);
    }

//    std::cout << "-- after preparation_" << std::endl;
//    for (auto&& e : preparation_)
//    {
//        std::cout << e << std::endl;
//    }

    gate_list.splice(gate_list.end(), (*decomposer_)(qubit_num_, 0, preparation_, qubit_names_));
}

void SimpleCircuitBuilder::ApplyPhaseGates(std::list<Gate>& gate_list,
                                           const std::unordered_map<int, int>& target_phase_map)
{
    for (auto&& tp : target_phase_map)
    {
        const int target = tp.first;
        const int phase_exponent_index = tp.second;

        if (phase_exponent_[phase_exponent_index].first <= 4)
        {
            if (phase_exponent_[phase_exponent_index].first / 4 == 1)
            {
                gate_list.emplace_back("Z", qubit_names_[target]);
            }
            if (phase_exponent_[phase_exponent_index].first / 2 == 1)
            {
                gate_list.emplace_back("P", qubit_names_[target]);
            }
            if (phase_exponent_[phase_exponent_index].first % 2 == 1)
            {
                gate_list.emplace_back("T", qubit_names_[target]);
            }
        }
        else
        {
            if (phase_exponent_[phase_exponent_index].first == 5 || phase_exponent_[phase_exponent_index].first == 6)
            {
                gate_list.emplace_back("P*", qubit_names_[target]);
            }
            if (phase_exponent_[phase_exponent_index].first % 2 == 1)
            {
                gate_list.emplace_back("T*", qubit_names_[target]);
            }
        }
    }
}

void SimpleCircuitBuilder::UnPrepare()
{
    preparation_ = std::move(restoration_);
    restoration_ = std::vector<util::xor_func>(qubit_num_);
    // re-initialize
    for (int i = 0; i < qubit_num_; i++)
    {
        restoration_[i] = util::xor_func(qubit_num_ + 1, 0);
        restoration_[i].set(i);
    }
}

void SimpleCircuitBuilder::PrepareLastPart(std::list<Gate>& gate_list,
                                           const std::vector<util::xor_func>& in,
                                           const std::vector<util::xor_func>& out)
{
    for (int i = 0; i < qubit_num_; i++)
    {
        bits_[i] = out[i];
    }
    util::ToUpperEchelon(qubit_num_, dimension_, bits_, &restoration_, std::vector<std::string>());
    util::FixBasis(qubit_num_, dimension_, qubit_num_, in, bits_, &restoration_, std::vector<std::string>());
    util::Compose(qubit_num_, preparation_, restoration_);


    gate_list.splice(gate_list.end(), (*decomposer_)(qubit_num_, 0, preparation_, qubit_names_));
}

static int EvaluateMatrix(const int n,
                          std::vector<util::xor_func> matrix)
{
    constexpr int not_cost = 1;
    constexpr int cnot_cost = 1;
    constexpr int swap_cost = cnot_cost * 3;

    int result = 0;

    bool flg = false;
    bool is_create = false;

    for (int j = 0; j < n; j++)
    {
        if (matrix[j].test(n))
        {
            matrix[j].reset(n);

            if (!is_create)
            {
                result += not_cost;
                is_create = true;
            }
        }
    }

    // Make triangular
    for (int i = 0; i < n; i++)
    {
        flg = false;
        is_create = false;
        for (int j = i; j < n; j++)
        {
            if (matrix[j].test(i))
            {
                if (!flg)
                {
                    if (j != i)
                    {
                        swap(matrix[i], matrix[j]);
                        result += swap_cost;
                    }
                    flg = true;
                }
                else
                {
                    matrix[j] ^= matrix[i];

                    if (!is_create)
                    {
                        result += cnot_cost;
                        is_create = true;
                    }
                }
            }
        }
        if (!flg)
        {
            std::cerr << "ERROR: not full rank" << std::endl;

            exit(1);
        }
    }

    // Finish the job
    for (int i = n - 1; i > 0; i--)
    {
        is_create = false;
        for (int j = i - 1; j >= 0; j--)
        {
            if (matrix[j].test(i))
            {
                matrix[j] ^= matrix[i];

                if (!is_create)
                {
                    result += cnot_cost;
                    is_create = true;
                }
            }
        }
    }

    return result;
}

void SimpleCircuitBuilder::ChangeRowOrder(std::unordered_map<int, int>& target_phase_map,
                                          std::vector<util::xor_func>& matrix)
{
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::uniform_int_distribution<> dist_index(0, matrix.size() - 1);

    // SA parameters
    const auto req_time = std::chrono::milliseconds(100);
    const auto start = std::chrono::system_clock::now();
    const auto end = std::chrono::system_clock::now() + req_time;
    auto current_time = start;

    constexpr auto rate = 10000;
    std::uniform_int_distribution<> dist(1, rate);

    int best_eval = EvaluateMatrix(qubit_num_, matrix);
    int current_eval = best_eval;
    std::vector<util::xor_func> current_matrix = matrix;
    std::unordered_map<int, int> current_target_phase_map = target_phase_map;

    /*
     * implement SA
     */
    while (true)
    {
        current_time = std::chrono::system_clock::now();
        if (current_time > end)
        {
            break;
        }

        const int index_a = dist_index(engine);
        const int index_b = dist_index(engine);

        if (index_a == index_b)
        {
            continue;
        }

        std::unordered_map<int, int> next_target_phase_map(current_target_phase_map);
        std::vector<util::xor_func> next_matrix(current_matrix);

        std::swap(next_matrix[index_a], next_matrix[index_b]);
        int target_a = -1;
        int target_b = -1;
        if (next_target_phase_map.count(index_a))
        {
            target_a = next_target_phase_map[index_a];
            next_target_phase_map.erase(index_a);
        }
        if (next_target_phase_map.count(index_b))
        {
            target_b = next_target_phase_map[index_b];
            next_target_phase_map.erase(index_b);
        }

        if (target_a > -1)
        {
            next_target_phase_map.emplace(index_b, target_a);
        }
        if (target_b > -1)
        {
            next_target_phase_map.emplace(index_a, target_b);
        }

        // evaluate matrix
        const int next_eval = EvaluateMatrix(qubit_num_, next_matrix);
        const auto time = current_time - start;
        const bool force_next = (rate * (req_time - time)) > (req_time * dist(seed_gen));

        // update
        if (current_eval > next_eval || force_next)
        {
            current_eval = next_eval;
            current_matrix = next_matrix;
            current_target_phase_map = next_target_phase_map;
        }

        if (best_eval > current_eval)
        {
            best_eval = current_eval;
            matrix = current_matrix;
            target_phase_map = current_target_phase_map;
        }
    }
}

std::list<Gate> SimpleCircuitBuilder::Build(const tpar::partitioning& partition,
                                            std::vector<util::xor_func>& in,
                                            const std::vector<util::xor_func>& out)
{
    std::list<Gate> ret;

    std::unordered_map<int, int> target_phase_map;

    if (Init(in, out) && partition.empty())
    {
        return ret;
    }

    /*
    std::cout << "-->> build circuit" << std::endl;
    std::cout << "# partition size: " << static_cast<int>(partition.size()) << std::endl;
    for (auto&& e : partition)
    {
        std::cout << "# --" << std::endl;
        for (auto&& ee : e)
        {
            std::cout << "# " << phase_exponent_[ee].second << std::endl;
        }
    }
    std::cout << "# input " << std::endl;
    for (auto&& e : in)
    {
        std::cout << "# " << e << std::endl;
    }
    std::cout << "# output " << std::endl;
    for (auto&& e : out)
    {
        std::cout << "# " << e << std::endl;
    }
     */

    /*
     * Reduce in to echelon form to decide on a basis
     */
    util::ToUpperEchelon(qubit_num_, dimension_, in, &preparation_, std::vector<std::string>());
    /*
     * For each partition... Compute *it, apply T gates, uncompute
     */
    for (auto&& it : partition)
    {
        InitBits(it, target_phase_map);

        /*
         * Prepare the bits
         */
        Prepare(ret, in, static_cast<int>(it.size()), target_phase_map);

        /*
         * Apply the phase gates
         */
        ApplyPhaseGates(ret, target_phase_map);

        /*
         * Unprepare the bits
         */
        UnPrepare();
    }

    /*
     * Reduce out to the basis of in
     */
    PrepareLastPart(ret, in, out);

    return ret;
}

std::list<Gate> SimpleCircuitBuilder::BuildGlobalPhase(int qubit_num,
                                                       int phase,
                                                       const std::vector<std::string>& qubit_names)
{
    std::list<Gate> acc;
    int qubit = 0;

    if (phase % 2 == 1)
    {
        acc.splice(acc.end(), util::ComposeOM(qubit, qubit_names));
        qubit = (qubit + 1) % qubit_num;
    }
    for (int i = phase / 2; i > 0; i--)
    {
        acc.splice(acc.end(), util::ComposeImaginaryUnit(qubit, qubit_names));
        qubit = (qubit + 1) % qubit_num;
    }

    return acc;
}

}