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

void MatrixReconstructor::init()
{
    // initialize some variables
    engine_ = std::mt19937(seed_generator_());
    rate_ = 10000;
    req_time_ = std::chrono::milliseconds(1000);

    /*
     * construct identity matrix
     */
    identity_.resize(num_qubit_);
    for (int i = 0; i < num_qubit_; i++)
    {
        identity_[i] = util::xor_func(num_qubit_ + 1, 0);
        identity_[i].set(i);
    }

}

std::vector<util::xor_func> MatrixReconstructor::execute(const std::vector<util::xor_func>& preparation,
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

    // initial parameters
    std::vector<util::xor_func> init_prep(identity_);
    util::compose(num_qubit_, init_prep, restoration);

    // current (temporary) parameters
    std::vector<util::xor_func> current_prep(init_prep);
    std::vector<util::xor_func> tmp_prep(preparation);
    util::compose(num_qubit_, tmp_prep, restoration);
    int current_eval = evaluate_matrix(num_qubit_, tmp_prep);

    // best parameters
    int best_eval = current_eval;
    std::vector<util::xor_func> best_prep = current_prep;

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
        std::vector<util::xor_func> next_prep(current_prep);
        std::swap(next_prep[target_a], next_prep[target_b]);

        // evaluate matrix
        std::vector<util::xor_func> tmp_rest(identity_);
        tmp_prep = preparation;
        util::compose(num_qubit_, tmp_rest, next_prep);
        util::compose(num_qubit_, tmp_prep, tmp_rest);
        const int next_eval = evaluate_matrix(num_qubit_, tmp_prep);

        // sa update param
        const auto time = std::chrono::duration_cast<std::chrono::milliseconds>(current_time - start);
        const bool force_next = (rate_ * (req_time_.count() - time.count())) > (req_time_.count() * dist(seed_generator_));

        // update
        if (current_eval > next_eval || force_next)
        {
            current_eval = next_eval;
            current_prep = next_prep;
        }

        if (best_eval > current_eval)
        {
            best_eval = current_eval;
            best_prep = current_prep;
        }
    }

    // set result restoration
    std::vector<util::xor_func> result_rest(identity_);
    util::compose(num_qubit_, result_rest, best_prep);

    // update target phase
    std::unordered_map<int, int> result_target_phase_map;
    for (auto&& map : target_phase_map)
    {
        const int bit = map.first;
        const int phase_index = map.second;
        const util::xor_func func = init_prep[bit];
        for (size_t i = 0; i < best_prep.size(); i++)
        {
            if (func == best_prep[i])
            {
                result_target_phase_map.emplace(i, phase_index);
            }
        }
    }
    target_phase_map = result_target_phase_map;

    return result_rest;
}

}