#ifndef T_SCHEDULING_GAUSSIAN_DECOMPOSER_HPP
#define T_SCHEDULING_GAUSSIAN_DECOMPOSER_HPP

#include "matrix_decomposer.hpp"

namespace tskd {

class GaussianDecomposer : public MatrixDecomposer {
private:

public:
    GaussianDecomposer() = default;

    ~GaussianDecomposer() final = default;

    std::list<Gate> operator()(const int n,
                               const int m,
                               std::vector<util::xor_func>& matrix,
                               const std::vector<std::string>& qubit_names) final;
};

}

#endif //T_SCHEDULING_GAUSSIAN_DECOMPOSER_HPP
