#include <map>

#include "util.hpp"

#include "../circuit/gate.hpp"

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

int compute_rank_destructive(int num_qubit,
                           int num_qubit_and_hadamard,
                           std::vector<xor_func>& bits)
{
    int rank = 0;
    for (int row = 0; row < num_qubit_and_hadamard; ++row)
    {
        bool flag = false;
        for (int col = rank; col < num_qubit; ++col)
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

int compute_rank(int num_qubit,
                 int num_qubit_and_hadamard,
                 std::vector<xor_func>& bits)
{
    std::vector<xor_func> parity_matrix = bits;

    return compute_rank_destructive(num_qubit, num_qubit_and_hadamard, parity_matrix);;
}

bool is_independent_destructive(int num_qubit,
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

bool is_independent(int num_qubit,
                   const std::vector<xor_func>& bits,
                   const xor_func& parity)
{
    xor_func temp_parity = parity;

    return is_independent_destructive(num_qubit, bits, temp_parity);
}

std::list<Gate> to_upper_echelon(int m,
                                 int n,
                                 std::vector<xor_func>& bits,
                                 std::vector<xor_func> *mat,
                                 const std::vector<std::string>& qubit_names)
{
    std::list<Gate> acc;
    int rank = 0;
    for (int j = 0; j < m; j++)
    {
        if (bits[j].test(n))
        {
            bits[j].reset(n);
            if (mat == NULL)
            {
                acc.splice(acc.end(), compose_x(j, qubit_names));
            }
            else
            {
                (*mat)[j].set(m);
            }
        }
    }

    /*
     * Make triangular
     */
    for (int i = 0; i < n; i++)
    {
        bool flg = false;
        for (int j = rank; j < m; j++)
        {
            if (bits[j].test(i))
            {
                // If we haven't yet seen a vector with bit i set...
                if (!flg)
                {
                    // If it wasn't the first vector we tried, swap to the front
                    if (j != rank)
                    {
                        swap(bits[rank], bits[j]);
                        if (mat == nullptr)
                        {
                            acc.splice(acc.end(), compose_swap(rank, j, qubit_names));
                        }
                        else
                        {
                            swap((*mat)[rank], (*mat)[j]);
                        }
                    }
                    flg = true;
                }
                else
                {
                    bits[j] ^= bits[rank];
                    if (mat == nullptr)
                    {
                        acc.splice(acc.end(), compose_cnot(rank, j, qubit_names));
                    }
                    else
                    {
                        (*mat)[j] ^= (*mat)[rank];
                    }
                }
            }
        }
        if (flg) rank++;
    }

    return acc;
}

std::list<Gate> to_lower_echelon(int m,
                                 int n,
                                 std::vector<xor_func>& bits,
                                 std::vector<xor_func>* mat,
                                 const std::vector<std::string>& qubit_names)
{
    std::list<Gate> acc;
    int i, j;

    for (i = n - 1; i > 0; i--)
    {
        for (j = i - 1; j >= 0; j--)
        {
            if (bits[j].test(i))
            {
                bits[j] ^= bits[i];
                if (mat == NULL)
                {
                    acc.splice(acc.end(), compose_cnot(i, j, qubit_names));
                }
                else
                {
                    (*mat)[j] ^= (*mat)[i];
                }
            }
        }
    }

    return acc;
}

std::list<Gate> fix_basis(int m,
                          int n,
                          int k,
                          const std::vector<xor_func>& fst,
                          std::vector<xor_func>& snd,
                          std::vector<xor_func>* mat,
                          const std::vector<std::string>& qubit_names) {
    std::list<Gate> acc;
    int j = 0;
    bool flg = false;
    std::map<int, int> pivots;  // mapping from columns to rows that have that column as pivot
    for (int i = 0; i < n; i++)
    {
        pivots[i] = -1;
    }

    // First pass makes sure tmp has the same pivots as fst
    for (int i = 0; i < m; i++)
    {
        // Find the next pivot
        while (j < n && !fst[i].test(j)) j++;
        if (j < n)
        {
            pivots[j] = i;
            flg = false;
            for (int h = i; !flg && h < k; h++)
            {
                // We found a vector with the same pivot
                if (snd[h].test(j))
                {
                    flg = true;
                    if (h != i)
                    {
                        swap(snd[h], snd[i]);
                        if (mat == NULL)
                        {
                            acc.splice(acc.end(), compose_swap(h, i, qubit_names));
                        }
                        else
                        {
                            swap((*mat)[h], (*mat)[i]);
                        }
                    }
                }
            }
            // There was no vector with the same pivot
            if (!flg)
            {
                if (k >= m)
                {
                    std::cerr << "FATAL ERROR: second space not a subspace\n" << std::endl;
                    exit(1);
                }
                snd[k] = fst[i];
                if (k != i)
                {
                    swap(snd[k], snd[i]);
                    if (mat == NULL)
                    {
                        acc.splice(acc.end(), compose_swap(k, i, qubit_names));
                    }
                    else
                    {
                        swap((*mat)[k], (*mat)[i]);
                    }
                }
                k++;
            }
        }
    }

    // Second pass makes each row of tmp equal to that row of fst
    for (int i = 0; i < m; i++) {
        for (int j = i + 1; j < n; j++)
        {
            if (fst[i][j] != snd[i][j])
            {
                if (pivots[j] == -1)
                {
                    std::cerr << "FATAL ERROR: cannot fix basis\n" << std::endl;
                    exit(1);
                }
                else
                {
                    snd[i] ^= snd[pivots[j]];
                    if (mat == NULL)
                    {
                        acc.splice(acc.end(), compose_cnot(pivots[j], i, qubit_names));
                    }
                    else
                    {
                        (*mat)[i] ^= (*mat)[pivots[j]];
                    }
                }
            }
        }
        if (!(snd[i] == fst[i]))
        {
            std::cerr << "FATAL ERROR: basis differs\n" << std::endl;
            exit(1);
        }
    }

    return acc;
}

void compose(int num,
             std::vector<xor_func>& A,
             const std::vector<xor_func>& B)
{
    auto tmp = std::vector<xor_func>(num);
    for (int i = 0; i < num; i++) {
        tmp[i] = B[i];
    }
    to_upper_echelon(num, num, tmp, &A, std::vector<std::string>());
    to_lower_echelon(num, num, tmp, &A, std::vector<std::string>());
}

std::list<Gate> compose_x(int target,
                          const std::vector<std::string>& qubit_names)
{
    std::list<Gate> ret;

    ret.emplace_back("tof", qubit_names[target]);

    return ret;
}

std::list<Gate> compose_swap(int a,
                             int b,
                             const std::vector<std::string>& qubit_names)
{
    std::list<Gate> ret;

    ret.emplace_back("tof", qubit_names[a], qubit_names[b]);
    ret.emplace_back("tof", qubit_names[b], qubit_names[a]);
    ret.emplace_back("tof", qubit_names[a], qubit_names[b]);

    return ret;
}

std::list<Gate> compose_cnot(int target,
                             int control,
                             const std::vector<std::string>& qubit_names)
{
    std::list<Gate> ret;

    ret.emplace_back("tof", qubit_names[target], qubit_names[control]);

    return ret;
}

std::list<Gate> compose_om(int target,
                           const std::vector<std::string>& qubit_names)
{
    std::list<Gate> ret;

    ret.emplace_back("H", qubit_names[target]);
    ret.emplace_back("P", qubit_names[target]);
    ret.emplace_back("H", qubit_names[target]);
    ret.emplace_back("P", qubit_names[target]);
    ret.emplace_back("H", qubit_names[target]);
    ret.emplace_back("P", qubit_names[target]);

    return ret;
}

std::list<Gate> compose_imaginary_unit(int target,
                                       const std::vector<std::string>& qubit_names)
{
    std::list<Gate> ret;

    ret.emplace_back("tof", qubit_names[target]);
    ret.emplace_back("Z", qubit_names[target]);
    ret.emplace_back("Y", qubit_names[target]);

    return ret;
}

}
}