#include <unordered_map>

#include "parallel_decomposer.hpp"

namespace tskd {

static
std::vector<std::list<Gate>> UpdateGateSetList(const int pivot,
                                               const std::vector<int>& one_array,
                                               std::unordered_map<std::string, int>& depth,
                                               std::vector<std::list<Gate>>& gate_set_list,
                                               const std::vector<std::string>& qubit_names)
{
    
}

std::list<Gate> ParallelDecomposer::operator()(const int n,
                                               const int m,
                                               std::vector<util::xor_func>& matrix,
                                               const std::vector<std::string>& qubit_names)
{
    std::list<Gate> ret;
    std::vector<int> one_array(static_cast<int>(matrix.size()));
    std::unordered_map<std::string, int> depth;
    std::vector<std::list<Gate>> gate_set_list;

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
                    one_array.push_back(j);
//                    matrix[j] ^= matrix[i];
//                    ret.splice(ret.begin(), util::ComposeCNOT(i, j, qubit_names));
                }
            }
        }
        if (!flg)
        {
            std::cerr << "ERROR: not full rank" << std::endl;

            exit(1);
        }

        // generate candidate cnot list
        AddGate(i, one_array, depth, qubit_names)

    }

    //Finish the job
    for (int i = n - 1; i > 0; i--)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            if (matrix[j].test(i))
            {
                matrix[j] ^= matrix[i];
                ret.splice(ret.begin(), util::ComposeCNOT(i, j, qubit_names));
            }
        }
    }

    return ret;
}

}