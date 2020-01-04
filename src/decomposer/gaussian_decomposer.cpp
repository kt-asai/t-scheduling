#include "gaussian_decomposer.hpp"

#include "../util/util.hpp"

namespace tskd {

std::list<Gate> GaussianDecomposer::execute(std::vector<util::xor_func>& matrix,
                                            std::vector<int>& func_map)
{
    std::list<Gate> lst;

    for (int j = 0; j < n(); j++)
    {
        if (matrix[j].test(n()))
        {
            matrix[j].reset(n());
            lst.splice(lst.begin(), util::compose_x(j, qubit_names()));
        }
    }

    // Make triangular
    for (int i = 0; i < n(); i++)
    {
        bool flg = false;
        for (int j = i; j < n() + m(); j++)
        {
            if (matrix[j].test(i))
            {
                // If we haven't yet seen a vector with bit i set...
                if (!flg)
                {
                    // If it wasn't the first vector we tried, swap to the front
                    if (j != i)
                    {
                        std::swap(matrix[i], matrix[j]);
                        std::swap(func_map[i], func_map[j]);
                    }
                    flg = true;
                }
                else
                {
                    matrix[j] ^= matrix[i];
                    lst.splice(lst.begin(), util::compose_cnot(func_map[i], func_map[j], qubit_names()));
                }
            }
        }
        if (!flg)
        {
            std::cerr << "ERROR: not full rank" << std::endl;

            exit(1);
        }
    }

    //Finish the job
    for (int i = n() - 1; i > 0; i--)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            if (matrix[j].test(i))
            {
                matrix[j] ^= matrix[i];
                lst.splice(lst.begin(), util::compose_cnot(func_map[i], func_map[j], qubit_names()));
            }
        }
    }

    return lst;
}

}