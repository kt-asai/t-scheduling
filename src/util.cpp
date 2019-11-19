#include <map>

#include "util.hpp"

namespace tskd {
namespace util {

int ComputeRankDestructive(int num_qubit,
                           int num_qubit_and_hadamard,
                           std::vector<xor_func>& bits)
{
    int rank = 0;
    for (int row = 0; row < num_qubit_and_hadamard; ++row)
    {
        bool flag = false;
        for (int col = 0; col < num_qubit; ++col)
        {
            if (bits[col].test(row))
            {
                if (!flag)
                {
                    if (col != rank)
                    {
                        std::swap(bits[rank], bits[col]);
                    }
                    flag = true;
                }
                else
                {
                    bits[col] ^= bits[rank];
                }
            }
        }
        if (flag)
        {
            rank++;
        }
    }

    return rank;
}

int ComputeRank(int num_qubit,
                int num_qubit_and_hadamard,
                std::vector<xor_func>& bits)
{
    std::vector<xor_func> parity_matrix = bits;

    return ComputeRankDestructive(num_qubit, num_qubit_and_hadamard, bits);;
}

bool IsIndependentDestructive(int num_qubit,
                              const std::vector<xor_func>& bits,
                              xor_func& parity)
{
    std::map<int, int> pivots;
    for (int row = 0, col = 0; row < num_qubit && col < static_cast<int>(bits.size());)
    {
        if (bits[col].test(row))
        {
            pivots.insert(std::make_pair(row, col));
            row++;
            col++;
        }
        else
        {
            col++;
        }
    }

    for (int i = 0; i < num_qubit; ++i)
    {
        if (parity.test(i))
        {
            auto it = pivots.find(i);
            if (it == pivots.end())
            {
                return true;
            }
            else
            {
                parity ^= bits[it->second];
            }
        }
    }

    return false;
}

bool IsIndependent(int num_qubit,
                   const std::vector<xor_func>& bits,
                   const xor_func& parity)
{
    xor_func temp_parity = parity;

    return IsIndependentDestructive(num_qubit, bits, temp_parity);
}

}
}