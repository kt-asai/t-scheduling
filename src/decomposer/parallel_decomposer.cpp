#include <map>
#include <unordered_map>

#include "parallel_decomposer.hpp"

#include "../circuit/gate.hpp"

namespace tskd {

static
bool CanParallelization(const std::list<Gate>& mappled_gate_list,
                        const Gate& new_gate)
{
    if (mappled_gate_list.empty())
    {
        return true;
    }

    return true;
}

static void UpdateUpperGateSetList(const int pivot,
                                   const std::vector<int>& one_array,
                                   std::unordered_map<std::string, int>& depth,
                                   std::vector<std::list<Gate>>& gate_set_list,
                                   std::vector<util::xor_func>& matrix,
                                   const std::vector<std::string>& qubit_names)
{
    std::map<int, std::vector<int>> depth_bit_map; // <depth, qubit_index>

    /*
     * Init
     */
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
     * check and generate cnot gate
     */
    std::vector<int> carry_bit_set;
    for (auto&& m : depth_bit_map)
    {
        const int bit_depth = m.first;
        std::vector<int> bit_set = m.second;

        bit_set.insert(bit_set.end(), carry_bit_set.begin(), carry_bit_set.end());
        carry_bit_set.clear();

        if (bit_set.size() == 1)
        {
            carry_bit_set.push_back(bit_set.front());

            continue;
        }

        std::sort(bit_set.begin(), bit_set.end());

        /*
         * Search for gates that can be parallelized
         */
        for (size_t i = 0; i < bit_set.size(); i++)
        {
            if (bit_set[i] == -1)
            {
                continue;
            }

            const std::string control = qubit_names[i];
            std::vector<std::string> target_list;
            std::vector<Gate> candidates;

            for (size_t j = i + 1; j < bit_set.size(); j++)
            {
                if (bit_set[j] == -1)
                {
                    continue;
                }

                target_list.push_back(qubit_names[j]);
                const Gate new_gate("tof", control, target_list);

                if (CanParallelization(gate_set_list[bit_depth], new_gate))
                {
                    bit_set[j] = -1;
                    depth[qubit_names[j]]++;
                    candidates.push_back(new_gate);
                }
                else
                {
                    target_list.pop_back();
                }
            }

            if (!candidates.empty())
            {
                depth[qubit_names[i]]++;
                gate_set_list[bit_depth].push_back(candidates.back());
            }
        }

        for (auto&& e : bit_set)
        {
            if (e != -1)
            {
                carry_bit_set.push_back(e);
            }
        }
    }

    // debug
//    std::cout << "# update gate set list" << std::endl;
//    int d = 0;
//    for (auto&& gate_group : gate_set_list)
//    {
//        std::cout << "## " << d << "part" << std::endl;
//        for (auto&& gate : gate_group)
//        {
//            gate.print();
//        }
//        d++;
//    }
}

static std::list<Gate> GenerateGateList(const std::vector<std::list<Gate>>& gate_set_list)
{
    std::list<Gate> ret;

    if (gate_set_list.empty())
    {
        return ret;
    }

    for (auto&& gate_group : gate_set_list)
    {
        if (gate_group.empty()) break;

        ret.emplace_back("block");
        for (auto&& gate : gate_group)
        {
            ret.push_back(gate);
        }
        ret.emplace_back("block");
    }

    return ret;
}

static void UpdateLowerGateSetList(const int pivot,
                                   const std::vector<int>& one_array,
                                   std::unordered_map<std::string, int>& depth,
                                   std::vector<std::list<Gate>>& gate_set_list,
                                   const std::vector<std::string>& qubit_names)
{
    std::map<int, std::vector<int>> depth_bit_map; // <depth, qubit_index>

    /*
     * Init
     */
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
     * check and generate cnot gate
     */
    std::vector<int> carry_bit_set;
    for (auto&& m : depth_bit_map)
    {
        const int bit_depth = m.first;
        std::vector<int> bit_set = m.second;

        bit_set.insert(bit_set.end(), carry_bit_set.begin(), carry_bit_set.end());
        carry_bit_set.clear();

        if (bit_set.size() == 1)
        {
            carry_bit_set.push_back(bit_set.front());

            continue;
        }

        std::sort(bit_set.begin(), bit_set.end(), std::greater<int>());

        /*
         * Search for gates that can be parallelized
         */
        for (size_t i = 0; i < bit_set.size(); i++)
        {
            if (bit_set[i] == -1)
            {
                continue;
            }

            const std::string control = qubit_names[i];
            std::vector<std::string> target_list;
            std::vector<Gate> candidates;

            for (size_t j = i + 1; j < bit_set.size(); j++)
            {
                if (bit_set[j] == -1)
                {
                    continue;
                }

                target_list.push_back(qubit_names[j]);
                const Gate new_gate("tof", control, target_list);

                if (CanParallelization(gate_set_list[bit_depth], new_gate))
                {
                    bit_set[j] = -1;
                    depth[qubit_names[j]]++;
                    candidates.push_back(new_gate);
                }
                else
                {
                    target_list.pop_back();
                }
            }

            if (!candidates.empty())
            {
                depth[qubit_names[i]]++;
                gate_set_list[bit_depth].push_back(candidates.back());
            }
        }

        for (auto&& e : bit_set)
        {
            if (e != -1)
            {
                carry_bit_set.push_back(e);
            }
        }
    }

    // debug
//    std::cout << "# update gate set list" << std::endl;
//    int d = 0;
//    for (auto&& gate_group : gate_set_list)
//    {
//        std::cout << "## " << d << "part" << std::endl;
//        for (auto&& gate : gate_group)
//        {
//            gate.print();
//        }
//        d++;
//    }
}

std::list<Gate> ParallelDecomposer::operator()(const int n,
                                               const int m,
                                               std::vector<util::xor_func>& matrix,
                                               const std::vector<std::string>& qubit_names)
{
    std::list<Gate> ret;
    std::vector<int> one_array(static_cast<int>(matrix.size()));
    std::unordered_map<std::string, int> depth;
    std::vector<std::list<Gate>> gate_set_list(4 * (2 * matrix.size()));

    // init
    for (auto&& name : qubit_names)
    {
        depth.emplace(name, 0);
    }

    for (int j = 0; j < n; j++)
    {
        if (matrix[j].test(n))
        {
            matrix[j].reset(n);
            ret.splice(ret.begin(), util::ComposeX(j, qubit_names));
        }
    }

//    std::cout << "--- init matrix decompose" << std::endl;
//    for (auto&& e : matrix)
//    {
//        std::cout << e << std::endl;
//    }

    // Make triangular
    for (int i = 0; i < n; i++)
    {
        bool flg = false;
        one_array.clear();
        for (int j = i; j < n + m; j++)
        {

            if (matrix[j].test(i))
            {
                // If we haven't yet seen a vector with bit i set...
                if (!flg) {
                    // If it wasn't the first vector we tried, swap to the front
                    if (j != i)
                    {
                        swap(matrix[i], matrix[j]);
                        ret.splice(ret.begin(), util::ComposeSwap(i, j, qubit_names));
                    }
                    flg = true;
                }
                else
                {
                    matrix[j] ^= matrix[i];
                    one_array.push_back(j);
                }
            }
        }
        if (!flg)
        {
            std::cerr << "ERROR: not full rank" << std::endl;

            exit(1);
        }

        // generate candidate cnot list
        UpdateUpperGateSetList(i, one_array, depth, gate_set_list, matrix, qubit_names);
    }

    // add gate
    ret.splice(ret.end(), GenerateGateList(gate_set_list));

//    std::cout << "--- upper after matrix decompose" << std::endl;
//    for (auto&& e : matrix)
//    {
//        std::cout << e << std::endl;
//    }

    //Finish the job
    std::vector<int> lower_one_array(static_cast<int>(matrix.size()));
    for (int i = n - 1; i > 0; i--)
    {
        lower_one_array.clear();
        for (int j = i - 1; j >= 0; j--)
        {
            if (matrix[j].test(i))
            {
                matrix[j] ^= matrix[i];
                lower_one_array.push_back(j);
                ret.splice(ret.begin(), util::ComposeCNOT(i, j, qubit_names));
            }
        }

        // generate candidate cnot list
//        UpdateLowerGateSetList(i, lower_one_array, depth, gate_set_list, qubit_names);
    }

//    // add gate
//    ret.splice(ret.end(), GenerateGateList(gate_set_list));

    return ret;
}

}