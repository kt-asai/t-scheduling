#include <random>
#include <chrono>

#include "greedy_circuit_builder.hpp"

#include "../util/util.hpp"

#include "../circuit/gate.hpp"

#include "../tpar/partition.hpp"

namespace tskd {

bool GreedyCircuitBuilder::Init(const std::vector<util::xor_func>& in,
                                const std::vector<util::xor_func>& out)
{
    bool is_io_different = true;

//    bits_ = std::vector<util::xor_func>(qubit_num_);
//    preparation_ = std::vector<util::xor_func>(qubit_num_);
//    restoration_ = std::vector<util::xor_func>(qubit_num_);

    for (int i = 0; i < qubit_num_; i++)
    {
        is_io_different &= (in[i] == out[i]);
//        preparation_[i] = util::xor_func(qubit_num_ + 1, 0);
//        restoration_[i] = util::xor_func(qubit_num_ + 1, 0);
//        preparation_[i].set(i);
//        restoration_[i].set(i);
    }

    return is_io_different;
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

void GreedyCircuitBuilder::ChangeRowOrder(std::unordered_map<int, int>& target_phase_map,
                                          std::vector<util::xor_func>& matrix)
{
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::uniform_int_distribution<> dist(0, matrix.size() - 1);

    int eval = EvaluateMatrix(qubit_num_, matrix);

    std::vector<util::xor_func> current_matrix = matrix;
    std::unordered_map<int, int>& current_target_phase_map = target_phase_map;

    // SA parameters
    auto req_time = std::chrono::seconds(1);
    auto start = std::chrono::system_clock::now();
    auto end = std::chrono::system_clock::now() + req_time;

    // TODO: implement SA
    int cnt = 5;
    while (cnt > 0)
    {
        const int index_a = dist(engine);
        const int index_b = dist(engine);

        if (index_a == index_b)
        {
            continue;
        }

        current_matrix = matrix;
        current_target_phase_map = target_phase_map;

        std::swap(current_matrix[index_a], current_matrix[index_b]);
        int target_a = -1;
        int target_b = -1;
        if (current_target_phase_map.count(index_a))
        {
            target_a = current_target_phase_map[index_a];
            current_target_phase_map.erase(index_a);
        }
        if (current_target_phase_map.count(index_b))
        {
            target_b = current_target_phase_map[index_b];
            current_target_phase_map.erase(index_b);
        }

        if (target_a > -1)
        {
            current_target_phase_map.emplace(index_b, target_a);
        }
        if (target_b > -1)
        {
            current_target_phase_map.emplace(index_a, target_b);
        }

        // evaluate matrix
        const int current_eval = EvaluateMatrix(qubit_num_, current_matrix);
        if (current_eval > eval)
        {
            // update
            eval = current_eval;
            matrix = current_matrix;
            target_phase_map = current_target_phase_map;
        }

        cnt--;
    }
}

int GreedyCircuitBuilder::ComputeTimeStep(const std::list<Gate>& gate_list)
{
    return static_cast<int>(gate_list.size());
}

void GreedyCircuitBuilder::ApplyPhaseGates(std::list<Gate>& gate_list,
                                           const std::set<int>& phase_exponent_index_set)
{
    int index = 0;
    for (auto&& phase_exponent_index : phase_exponent_index_set)
    {
        if (phase_exponent_[phase_exponent_index].first <= 4)
        {
            if (phase_exponent_[phase_exponent_index].first / 4 == 1)
            {
                gate_list.emplace_back("Z", qubit_names_[index]);
            }
            if (phase_exponent_[phase_exponent_index].first / 2 == 1)
            {
                gate_list.emplace_back("P", qubit_names_[index]);
            }
            if (phase_exponent_[phase_exponent_index].first % 2 == 1)
            {
                gate_list.emplace_back("T", qubit_names_[index]);
            }
        }
        else
        {
            if (phase_exponent_[phase_exponent_index].first == 5 || phase_exponent_[phase_exponent_index].first == 6)
            {
                gate_list.emplace_back("P*", qubit_names_[index]);
            }
            if (phase_exponent_[phase_exponent_index].first % 2 == 1)
            {
                gate_list.emplace_back("T*", qubit_names_[index]);
            }
        }
        index++;
    }
}

void GreedyCircuitBuilder::ApplyPhaseGates(std::list<Gate>& gate_list,
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

int GreedyCircuitBuilder::CheckDimension(const Character& chr,
                                         std::vector <util::xor_func>& wires,
                                         int current_dimension)
{
    int new_dimension = 0;
    const int updated_dimension = util::ComputeRank(chr.num_qubit(), chr.num_data_qubit() + chr.num_hadamard(), wires);
    if (updated_dimension > current_dimension)
    {
        new_dimension = updated_dimension;
        oracle_.set_dim(new_dimension);
    }

    return new_dimension;
}

std::list<Gate> GreedyCircuitBuilder::Build(std::list<int>& index_list,
                                            std::list<int>& carry_index_list,
                                            std::vector<util::xor_func>& in,
                                            const std::vector<util::xor_func>& out)
{
    /*
    std::cout << std::endl;
    std::cout << "-->> greedy circuit builder" << std::endl;
    std::cout << "# input" << std::endl;
    std::cout << "-- index list" << std::endl;
    for (auto&& e : index_list)
    {
        std::cout << e << ":" << phase_exponent_[e].second << std::endl;
    }
    std::cout << "-- carry index list" << std::endl;
    for (auto&& e : carry_index_list)
    {
        std::cout << e << ":" << phase_exponent_[e].second << std::endl;
    }
     */

    std::list<Gate> ret;

    if (Init(in, out) && index_list.empty())
    {
        return ret;
    }

    std::vector<util::xor_func> bits(qubit_num_);
    std::vector<util::xor_func> preparation(qubit_num_);
    std::vector<util::xor_func> restoration(qubit_num_);
    for (int i = 0; i < qubit_num_; i++)
    {
        preparation[i] = util::xor_func(qubit_num_ + 1, 0);
        restoration[i] = util::xor_func(qubit_num_ + 1, 0);
        preparation[i].set(i);
        restoration[i].set(i);
    }
    std::vector<std::pair<int, int>> phase_target_list;

    /*
     * Reduce in to echelon form to decide on a basis
     */
    util::ToUpperEchelon(qubit_num_, dimension_, in, &preparation, std::vector<std::string>());

    while (!index_list.empty())
    {
//        std::cout << "-- index list size: " << index_list.size() << std::endl;
//        std::cout << "-- index list in loop" << std::endl;
//        for (auto&& e : index_list)
//        {
//            std::cout << e << ":" << phase_exponent_[e].second << std::endl;
//        }

        /**
         * first build sub-circuit
         */
        std::vector<util::xor_func> result_bits(qubit_num_);
        std::vector<util::xor_func> result_preparation(qubit_num_);
        std::vector<util::xor_func> result_restoration(qubit_num_);
        for (int i = 0; i < qubit_num_; i++)
        {
            result_preparation[i] = util::xor_func(qubit_num_ + 1, 0);
            result_restoration[i] = util::xor_func(qubit_num_ + 1, 0);
            result_preparation[i].set(i);
            result_restoration[i].set(i);
        }
        std::list<Gate> result_gate_list;
        std::set<int> result_sub_part;
        std::list<int> delete_index_list;
        std::unordered_map<int, int> result_target_phase_map;

        for (auto it = index_list.begin(); it != index_list.end();)
        {
            std::list<Gate> tmp_gate_list;
            std::set<int> tmp_sub_part = result_sub_part;
            tmp_sub_part.insert(*it);

            std::vector<util::xor_func> tmp_bits(qubit_num_);
            std::vector<util::xor_func> tmp_preparation = preparation;
            std::vector<util::xor_func> tmp_restoration = restoration;

            std::unordered_map<int, int> tmp_target_phase_map;

            if (oracle_(phase_exponent_, tmp_sub_part))
            {
//                std::cout << "oracle part num: " << tmp_sub_part.size() << std::endl;
                /**
                 * create bits matrix
                 */
                std::set<int>::iterator ti;
                int counter = 0;
                for (ti = tmp_sub_part.begin(), counter = 0; counter < qubit_num_; counter++)
                {
                    if (counter < static_cast<int>(tmp_sub_part.size()))
                    {
                        tmp_bits[counter] = phase_exponent_[*ti].second;
                        tmp_target_phase_map.emplace(counter, *ti);
                        ti++;
                    }
                    else
                    {
                        tmp_bits[counter] = util::xor_func(dimension_ + 1, 0);
                    }
                }

                /**
                 * prepare preparation matrix
                 */

                const int num_partition = static_cast<int>(tmp_sub_part.size());
                /*
                std::cout << std::endl;
                std::cout << "# prepare" << std::endl;
                std::cout << "## num qubit: " << qubit_num_ << std::endl;
                std::cout << "## dimension: " << dimension_ << std::endl;
                std::cout << "## num partition: " << num_partition << std::endl;
                std::cout << "## in" << std::endl;
                for (auto&& e : in)
                {
                    std::cout << e << std::endl;
                }
                std::cout << "## tmp_bits" << std::endl;
                for (auto&& e : tmp_bits)
                {
                    std::cout << e << std::endl;
                }
                std::cout << "## tmp_restoration" << std::endl;
                for (auto&& e : tmp_restoration)
                {
                    std::cout << e << std::endl;
                }
                 */
                util::ToUpperEchelon(num_partition, dimension_, tmp_bits, &tmp_restoration, std::vector<std::string>());
                util::FixBasis(qubit_num_, dimension_, num_partition, in, tmp_bits, &tmp_restoration,
                               std::vector<std::string>());
                util::Compose(qubit_num_, tmp_preparation, tmp_restoration);

                /**
                 * change row order in the matrix
                 */
                if (option_.change_row_order())
                {
                    // ChangeRowOrder();
                    ChangeRowOrder(tmp_target_phase_map, tmp_preparation);
                }

                /**
                 * create gate list
                 */
                tmp_gate_list.splice(tmp_gate_list.end(),
                                     (*decomposer_)(qubit_num_, 0, tmp_preparation, qubit_names_));

                /**
                 * check time step
                 */
                if (ComputeTimeStep(tmp_gate_list) < option_.distillation_step())
                {
                    result_bits = tmp_bits;
                    result_preparation = tmp_preparation;
                    result_restoration = tmp_restoration;
                    result_sub_part = tmp_sub_part;
                    result_gate_list = tmp_gate_list;
                    result_target_phase_map = tmp_target_phase_map;
                    it = index_list.erase(it);
                }
                else
                {
                    it++;
                }
            }
            else
            {
                it++;
            }
        }

//        index_list.splice(index_list.end(), std::move(delete_index_list));
//        index_list.unique();
//        std::cout << "after unique" << std::endl;
//        for (auto&& e : index_list)
//        {
//            std::cout << e << " ";
//        }
//        std::cout << std::endl;

        /**
         * second build sub-circuit
         */
        for (auto it = carry_index_list.begin(); it != carry_index_list.end();)
        {
            std::list<Gate> tmp_gate_list;
            std::set<int> tmp_sub_part = result_sub_part;
            tmp_sub_part.insert(*it);

            std::vector<util::xor_func> tmp_bits(qubit_num_);
            std::vector<util::xor_func> tmp_preparation = preparation;
            std::vector<util::xor_func> tmp_restoration = restoration;

            std::unordered_map<int, int> tmp_target_phase_map;

            if (oracle_(phase_exponent_, tmp_sub_part))
            {
                /**
                 * create bits matrix
                 */
                std::set<int>::iterator ti;
                int counter = 0;
                for (ti = tmp_sub_part.begin(), counter = 0; counter < qubit_num_; counter++)
                {
                    if (counter < static_cast<int>(tmp_sub_part.size()))
                    {
                        tmp_bits[counter] = phase_exponent_[*ti].second;
                        tmp_target_phase_map.emplace(counter, *ti);
                        ti++;
                    }
                    else
                    {
                        tmp_bits[counter] = util::xor_func(dimension_ + 1, 0);
                    }
                }

                /**
                 * prepare preparation matrix
                 */
                const int num_partition = static_cast<int>(tmp_sub_part.size());
                util::ToUpperEchelon(num_partition, dimension_, tmp_bits, &tmp_restoration, std::vector<std::string>());
                util::FixBasis(qubit_num_, dimension_, num_partition, in, tmp_bits, &tmp_restoration,
                               std::vector<std::string>());
                util::Compose(qubit_num_, tmp_preparation, tmp_restoration);

                /**
                 * change row order in the matrix
                 */
                if (option_.change_row_order())
                {
                     ChangeRowOrder(tmp_target_phase_map, tmp_preparation);
                }

                /**
                 * create gate list
                 */
                tmp_gate_list.splice(tmp_gate_list.end(),
                                     (*decomposer_)(qubit_num_, 0, tmp_preparation, qubit_names_));

                /**
                 * check time step
                 */
                if (ComputeTimeStep(tmp_gate_list) < option_.distillation_step())
                {
                    result_bits = tmp_bits;
                    result_preparation = tmp_preparation;
                    result_restoration = tmp_restoration;
                    result_sub_part = tmp_sub_part;
                    result_gate_list = tmp_gate_list;
                    result_target_phase_map = tmp_target_phase_map;
                    it = carry_index_list.erase(it);
                }
                else
                {
                    it++;
                }
            }
            else
            {
                it++;
            }
        }

        ret.splice(ret.end(), result_gate_list);

        /*
         * Apply the phase gates
         */
//        ApplyPhaseGates(ret, result_sub_part);
        ApplyPhaseGates(ret, result_target_phase_map);

        /*
         * Unprepare the bits
         */
        preparation = std::move(result_restoration);
        restoration = std::vector<util::xor_func>(qubit_num_);
        // re-initialize
        for (int i = 0; i < qubit_num_; i++)
        {
            restoration[i] = util::xor_func(qubit_num_ + 1, 0);
            restoration[i].set(i);
        }
    }

    /*
     * Reduce out to the basis of in
     */
    for (int i = 0; i < qubit_num_; i++)
    {
        bits[i] = out[i];
    }
    util::ToUpperEchelon(qubit_num_, dimension_, bits, &restoration, std::vector<std::string>());
    util::FixBasis(qubit_num_, dimension_, qubit_num_, in, bits, &restoration, std::vector<std::string>());
    util::Compose(qubit_num_, preparation, restoration);

    ret.splice(ret.end(), (*decomposer_)(qubit_num_, 0, preparation, qubit_names_));

    return ret;
}

std::list<Gate> GreedyCircuitBuilder::BuildGlobalPhase(int qubit_num,
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
