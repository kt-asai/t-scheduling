#include <map>

#include "util.hpp"

namespace tskd {
namespace util {

bool IndependentOracle::operator()(const std::vector <phase_exponent>& expnts,
                                   const std::set<int>& lst) const
{
    if (lst.size() > num_)
    {
        return false;
    }

    if (lst.size() == 1 || (num_ - lst.size()) >= dim_)
    {
        return true;
    }

    std::set<int>::const_iterator it;
    int i, j, rank = 0;
    auto tmp = std::vector<xor_func>(lst.size());

    for (i = 0, it = lst.begin(); it != lst.end(); it++, i++)
    {
        tmp[i] = expnts[*it].second;
    }

    for (i = 0; i < length_; i++)
    {
        bool flg = false;
        for (j = rank; j < lst.size(); j++)
        {
            if (tmp[j].test(i))
            {
                if (!flg)
                {
                    if (j != rank) swap(tmp[rank], tmp[j]);
                    flg = true;
                }
                else
                {
                    tmp[j] ^= tmp[rank];
                }
            }
        }
        if (flg) rank++;
    }

    return (num_ - lst.size()) >= (dim_ - rank);
}

int IndependentOracle::retrieve_lin_dep(const std::vector<phase_exponent>& expnts,
                                        const std::set<int>& lst) const
{
    std::set<int>::const_iterator it;
    int i, j, rank = 0, tmpr;
    std::map<int, int> mp;
    auto tmp = std::vector<xor_func>(lst.size());

    for (i = 0, it = lst.begin(); it != lst.end(); it++, i++)
    {
        tmp[i] = expnts[*it].second;
        mp[i] = *it;
    }

    for (j = 0; j < lst.size(); j++)
    {
        if (tmp[j].test(length_)) tmp[j].reset(length_);
    }

    for (i = 0; i < length_; i++)
    {
        bool flg = false;
        for (j = rank; j < lst.size(); j++)
        {
            if (tmp[j].test(i))
            {
                if (!flg)
                {
                    if (j != rank)
                    {
                        swap(tmp[rank], tmp[j]);
                        tmpr = mp[rank];
                        mp[rank] = mp[j];
                        mp[j] = tmpr;
                    }
                    flg = true;
                }
                else
                {
                    tmp[j] ^= tmp[rank];
                    if (tmp[j].none()) return mp[j];
                }
            }
        }
        if (flg) rank++;
    }

    assert((num_ - lst.size()) >= (dim_ - rank));
    return -1;
}

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