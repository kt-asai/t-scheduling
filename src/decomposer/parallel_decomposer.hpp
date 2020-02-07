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

    ParallelDecomposer(const Layout& layout,
                       const int n,
                       const int m,
                       const std::vector<std::string>& qubit_names)
            : MatrixDecomposer(layout, n, m, qubit_names) { }

    ~ParallelDecomposer() final = default;

    std::list<Gate> execute(std::vector<util::xor_func>& matrix,
                            std::vector<int>& func_map) final;
};

}

#endif //T_SCHEDULING_PARALLEL_DECOMPOSER_HPP
