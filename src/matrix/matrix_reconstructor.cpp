#include "matrix_reconstructor.hpp"

namespace tskd {

static int evaluate_matrix(const int n,
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

std::vector<util::xor_func> MatrixReconstructor::execute(const std::vector<util::xor_func>& identity,
                                                         const std::vector<util::xor_func>& preparation,
                                                         const std::vector<util::xor_func>& restoration,
                                                         std::unordered_map<int, int>& target_phase_map)
{
    // random generator
    const int matrix_size = static_cast<int>(preparation.size());
    std::uniform_int_distribution<> dist_index(0, matrix_size - 1);

    // SA parameters
    const auto start = std::chrono::system_clock::now();
    const auto end = std::chrono::system_clock::now() + req_time_;
    auto current_time = start;
    std::uniform_int_distribution<> dist(1, rate_);

    // current (temporary) parameters
    std::vector<util::xor_func> tmp_prep(preparation);
    std::vector<util::xor_func> current_rest(restoration);
    util::compose(num_qubit_, tmp_prep, current_rest);
    int current_eval = evaluate_matrix(num_qubit_, tmp_prep);
    std::unordered_map<int, int> current_target_phase_map = target_phase_map;

    // best parameters
    int best_eval = current_eval;
    std::vector<util::xor_func> best_rest = current_rest;

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

        // choice index randomly
        const int target_a = dist_index(engine_);
        const int target_b = dist_index(engine_);

        if (target_a == target_b)
        {
            continue;
        }

        // swap process
        std::vector<util::xor_func> next_prep(identity);
        util::compose(num_qubit_, next_prep, current_rest);
        std::unordered_map<int, int> next_target_phase_map(current_target_phase_map);

        std::swap(next_prep[target_a], next_prep[target_b]);
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
        std::vector<util::xor_func> next_rest(identity);
        util::compose(num_qubit_, next_rest, next_prep);
        tmp_prep = preparation;
        util::compose(num_qubit_, tmp_prep, next_rest);
        const int next_eval = evaluate_matrix(num_qubit_, tmp_prep);

        // sa update param
        const auto time = current_time - start;
        const bool force_next = (rate_ * (req_time_ - time)) > (req_time_ * dist(seed_generator_));

        // update
        if (current_eval > next_eval || force_next)
        {
            current_eval = next_eval;
            current_rest = next_rest;
            current_target_phase_map = next_target_phase_map;
        }

        if (best_eval > current_eval)
        {
            best_eval = current_eval;
            best_rest = current_rest;
            target_phase_map = current_target_phase_map;
        }

    }

    return best_rest;
}

}