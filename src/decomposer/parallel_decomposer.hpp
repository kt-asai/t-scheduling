#ifndef T_SCHEDULING_PARALLEL_DECOMPOSER_HPP
#define T_SCHEDULING_PARALLEL_DECOMPOSER_HPP

#include "matrix_decomposer.hpp"

#include "../layout/layout.hpp"

namespace tskd {

class ParallelDecomposer : public MatrixDecomposer
{
private:

public:
    ParallelDecomposer() = default;

    ~ParallelDecomposer() final = default;

    std::list<Gate> operator()(const Layout& layout,
                               const int n,
                               const int m,
                               std::vector<util::xor_func>& matrix,
                               const std::vector<std::string>& qubit_names) final;
};

}

#endif //T_SCHEDULING_PARALLEL_DECOMPOSER_HPP
