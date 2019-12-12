#include "parallel_decomposer.hpp"

namespace tskd {

std::list<Gate> ParallelDecomposer::operator()(const int n,
                                               const int m,
                                               std::vector<util::xor_func>& matrix,
                                               const std::vector<std::string>& qubit_names)
{
    std::list<Gate> ret;

    return ret;
}

}