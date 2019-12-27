#include "matrix_reconstructor.hpp"

namespace tskd {

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

std::vector<util::xor_func> MatrixReconstructor::Execute(const int num_partition,
                                                         const std::vector<util::xor_func>& bits,
                                                         std::unordered_map<int, int>& target_phase_map)
{
    // random generator
    std::uniform_int_distribution<> dist_index(0, num_partition - 1);

    // SA parameters
    const auto start = std::chrono::system_clock::now();
    const auto end = std::chrono::system_clock::now() + req_time_;
    auto current_time = start;

    std::uniform_int_distribution<> dist(1, rate_);

    // current (temporary) parameters
    std::vector<util::xor_func> init_restoration(num_qubit_);
    for (int i = 0; i < num_qubit_; i++)
    {
        init_restoration[i] = util::xor_func(num_qubit_ + 1, 0);
        init_restoration[i].set(i);
    }
    std::vector<util::xor_func> bits_backup(bits);
    util::ToUpperEchelon(num_partition, dimension_, bits_backup, &init_restoration, std::vector<std::string>());
    util::FixBasis(num_qubit_, dimension_, num_partition, input_, bits_backup, &init_restoration, std::vector<std::string>());

    int current_eval = EvaluateMatrix(num_qubit_, init_restoration);
    std::vector<util::xor_func> current_bits = bits;
    std::unordered_map<int, int> current_target_phase_map = target_phase_map;

    // best parameter
    int best_eval = current_eval;
    std::vector<util::xor_func> best_bits = bits;

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

        const int target_a = dist_index(engine_);
        const int target_b = dist_index(engine_);

        if (target_a == target_b)
        {
            continue;
        }

        /*
        // swap prcess
//        std::swap(current_matrix[target_a], current_matrix[target_b]);
//        int phase_index_a = -1;
//        int phase_index_b = -1;
//        if (current_target_phase_map.count(target_a))
//        {
//            phase_index_a = current_target_phase_map[target_a];
//            current_target_phase_map.erase(target_a);
//        }
//        if (current_target_phase_map.count(target_b))
//        {
//            phase_index_b = current_target_phase_map[target_b];
//            current_target_phase_map.erase(target_b);
//        }
//
//        if (phase_index_a != -1)
//        {
//            current_target_phase_map.emplace(target_b, phase_index_a);
//        }
//        if (phase_index_b != -1)
//        {
//            current_target_phase_map.emplace(target_a, phase_index_b);
//        }
//
//        // evaluate matrix
//        const int next_eval = EvaluateMatrix(qubit_num_, current_matrix);
//        const auto time = current_time - start;
//        const bool force_current = false; //(rate * (req_time - time)) > (req_time * dist(seed_gen));
//
//        // update
//        if (current_eval > next_eval || force_current)
//        {
//            current_eval = next_eval;
//        }
//        else
//        {
//            // recover
//            std::swap(current_matrix[target_a], current_matrix[target_b]);
//            if (phase_index_a != -1)
//            {
//                current_target_phase_map.erase(target_b);
//                current_target_phase_map.emplace(target_a, phase_index_a);
//            }
//            if (phase_index_b != -1)
//            {
//                current_target_phase_map.erase(target_a);
//                current_target_phase_map.emplace(target_b, phase_index_b);
//            }
//        }
//
//        if (best_eval > current_eval)
//        {
//            best_eval = current_eval;
//            matrix = current_matrix;
//            target_phase_map = current_target_phase_map;
//        }
        */

        // swap process
        std::vector<util::xor_func> next_bits(current_bits);
        std::unordered_map<int, int> next_target_phase_map(current_target_phase_map);

        std::swap(next_bits[target_a], next_bits[target_b]);
        int phase_index_a = -1;
        int phase_index_b = -1;
        if (next_target_phase_map.count(target_a))
        {
            phase_index_a = next_target_phase_map[target_a];
            next_target_phase_map.erase(target_a);
        }
        if (next_target_phase_map.count(target_b))
        {
            phase_index_b = next_target_phase_map[target_b];
            next_target_phase_map.erase(target_b);
        }

        if (phase_index_a != -1)
        {
            next_target_phase_map.emplace(target_b, phase_index_a);
        }
        if (phase_index_b != -1)
        {
            next_target_phase_map.emplace(target_a, phase_index_b);
        }

        // evaluate matrix
        std::vector<util::xor_func> tmp_restoration(num_qubit_);
        for (int i = 0; i < num_qubit_; i++)
        {
            tmp_restoration[i] = util::xor_func(num_qubit_ + 1, 0);
            tmp_restoration[i].set(i);
        }
        std::vector<util::xor_func> next_bits_backup(next_bits);

        util::ToUpperEchelon(num_partition, dimension_, next_bits, &tmp_restoration, std::vector<std::string>());
        util::FixBasis(num_qubit_, dimension_, num_partition, input_, next_bits, &tmp_restoration, std::vector<std::string>());
        const int next_eval = EvaluateMatrix(num_qubit_, tmp_restoration);
        const auto time = current_time - start;
        const bool force_next = (rate_ * (req_time_ - time)) > (req_time_ * dist(seed_generator_));

        // update
        if (current_eval > next_eval || force_next)
        {
            current_eval = next_eval;
            current_bits = next_bits_backup;
            current_target_phase_map = next_target_phase_map;
        }

        if (best_eval > current_eval)
        {
            best_eval = current_eval;
            best_bits = current_bits;
            target_phase_map = current_target_phase_map;
        }
    }

    return best_bits;
}

}