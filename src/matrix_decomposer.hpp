#ifndef T_SCHEDULING_MATRIX_DECOMPOSER_HPP
#define T_SCHEDULING_MATRIX_DECOMPOSER_HPP

#include <list>
#include <vector>
#include <string>

#include "gate.hpp"
#include "util.hpp"

namespace tskd {

class MatrixDecomposer {
private:

public:
    MatrixDecomposer() = default;

    virtual ~MatrixDecomposer();

    virtual std::list<Gate> operator()(const int n,
                                       const int m,
                                       std::vector<util::xor_func> &matrix,
                                       const std::vector<std::string> &qubit_names);
};

}

#endif //T_SCHEDULING_MATRIX_DECOMPOSER_HPP
