#include <map>
#include <unordered_map>

#include "parallel_decomposer.hpp"

#include "../parallel/parallelization_oracle.hpp"

namespace tskd {

static void update_gate_set_list(const int sign,
                                 const int pivot,
                                 const std::vector<int>& one_array,
                                 std::unordered_map<std::string, int>& depth,
                                 std::vector<std::list<Gate>>& gate_set_list,
                                 const std::vector<std::string>& qubit_names,
                                 ParallelizationOracle& oracle)
{
    /*
     * Init depth of each qubit
     */
    auto l = [](int lhs, int rhs){return lhs < rhs;};
    auto g = [](int lhs, int rhs){return lhs > rhs;};
    auto comp = sign ? g : l;
    std::map<int, std::vector<int>> depth_bit_map; // <depth, qubit_index>
    depth_bit_map.emplace(depth[qubit_names[pivot]], std::vector<int>(1, pivot));
    for (auto&& o : one_array)
    {
        const int bit_depth = depth[qubit_names[o]];
        if (depth_bit_map.count(bit_depth))
        {
            depth_bit_map[bit_depth].push_back(o);
        }
        else
        {
            depth_bit_map.emplace(bit_depth, std::vector<int>(1, o));
        }
    }

    /*
     * generate cnot gate
     */
    std::vector<int> carry_bit_set;
    for (auto it = depth_bit_map.begin(); it != depth_bit_map.end(); it++)
    {
        const int bit_depth = it->first;
        std::vector<int> bit_set = it->second;

        // carry move process
        for (auto&& carry : carry_bit_set)
        {
            const std::string carry_qubit_name = qubit_names[carry];
            depth[carry_qubit_name] = bit_depth;
        }
        bit_set.insert(bit_set.end(), carry_bit_set.begin(), carry_bit_set.end());
        carry_bit_set.clear();


        if (bit_set.size() == 1)
        {
            carry_bit_set.push_back(bit_set.front());

            continue;
        }

        std::sort(bit_set.begin(), bit_set.end(), comp);

        /*
         * Search for gates that can be parallelized
         */
        for (size_t i = 0; i < bit_set.size(); i++)
        {
            if (bit_set[i] == -1)
            {
                continue;
            }

            const std::string control = qubit_names[bit_set[i]];
            std::vector<std::string> target_list;
            std::vector<Gate> candidates;

            for (size_t j = i + 1; j < bit_set.size(); j++)
            {
                if (bit_set[j] == -1)
                {
                    continue;
                }

                const std::string candidate_target = qubit_names[bit_set[j]];
                target_list.push_back(candidate_target);
                const Gate new_gate("tof", control, target_list);

                if (oracle.check(gate_set_list[bit_depth], new_gate))
                {
                    bit_set[j] = -1;
                    depth[candidate_target]++;
                    candidates.push_back(new_gate);
                }
                else
                {
                    target_list.pop_back();
                }
            }

            if (!candidates.empty())
            {
                gate_set_list[bit_depth].push_back(candidates.back());
            }
        }

        for (auto&& bit : bit_set)
        {
            if (bit != -1)
            {
                const std::string qubit_name = qubit_names[bit];
                depth[qubit_name]++;
                carry_bit_set.push_back(bit);
            }
        }

        if (carry_bit_set.size() >= 2 && std::next(it) == depth_bit_map.end())
        {
            const int next_bit_depth = bit_depth + 1;
            depth_bit_map.insert(std::make_pair(next_bit_depth, std::vector<int>()));
        }
    }
}

static std::list<Gate> generate_gate_list(const std::vector<std::list<Gate>>& gate_set_list)
{
    std::list<Gate> ret;

    if (gate_set_list.empty())
    {
        return ret;
    }

    for (auto&& gate_group : gate_set_list)
    {
        if (gate_group.empty()) break;

        if (gate_group.size() >= 2)
        {
            ret.emplace_front("cnot");
        }
        else
        {
            ret.emplace_front(gate_group.front());
        }
    }

    return ret;
}

std::list<Gate> ParallelDecomposer::execute(std::vector<util::xor_func>& matrix,
                                            std::vector<int>& func_map)
{
    constexpr int swap_step = 3;
    constexpr int cnot_step = 1;
    const int max_num_gate = (swap_step + cnot_step) * (2 * matrix.size());
    std::list<Gate> ret;
    std::unordered_map<std::string, int> depth;
    std::vector<std::list<Gate>> gate_set_list(max_num_gate);

    ParallelizationOracle oracle(layout());

    // init
    for (auto&& name : qubit_names())
    {
        depth.emplace(name, 0);
    }

    for (int j = 0; j < n(); j++)
    {
        if (matrix[j].test(n()))
        {
            matrix[j].reset(n());
            ret.splice(ret.begin(), util::compose_x(j, qubit_names()));
        }
    }

    // Make triangular
    std::vector<int> one_array(static_cast<int>(matrix.size()));
    for (int i = 0; i < n(); i++)
    {
        bool flg = false;
        one_array.clear();
        for (int j = i; j < n() + m(); j++)
        {

            if (matrix[j].test(i))
            {
                // If we haven't yet seen a vector with bit i set...
                if (!flg) {
                    // If it wasn't the first vector we tried, swap to the front
                    if (j != i)
                    {
                        swap(matrix[i], matrix[j]);
                        std::swap(func_map[i], func_map[j]);
                    }
                    flg = true;
                }
                else
                {
                    matrix[j] ^= matrix[i];
                    one_array.push_back(func_map[j]);
                }
            }
        }
        if (!flg)
        {
            std::cerr << "ERROR: not full rank" << std::endl;

            exit(1);
        }

        // generate candidate cnot list
        update_gate_set_list(0, func_map[i], one_array, depth, gate_set_list, qubit_names(), oracle);
    }

    //Finish the job
    for (int i = n() - 1; i > 0; i--)
    {
        one_array.clear();
        for (int j = i - 1; j >= 0; j--)
        {
            if (matrix[j].test(i))
            {
                matrix[j] ^= matrix[i];
                one_array.push_back(func_map[j]);
            }
        }

        // generate candidate cnot list
        update_gate_set_list(1, func_map[i], one_array, depth, gate_set_list, qubit_names(), oracle);
    }

    // add gate
    ret.splice(ret.end(), generate_gate_list(gate_set_list));

    return ret;
}

}