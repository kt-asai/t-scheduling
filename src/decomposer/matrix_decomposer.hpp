#ifndef T_SCHEDULING_MATRIX_DECOMPOSER_HPP
#define T_SCHEDULING_MATRIX_DECOMPOSER_HPP

#include <list>
#include <vector>
#include <string>

#include "../util/util.hpp"

#include "../layout/layout.hpp"

#include "../circuit/gate.hpp"


namespace tskd {

class MatrixDecomposer
{
private:

public:
    MatrixDecomposer() = default;

    virtual ~MatrixDecomposer() = default;

    virtual std::list<Gate> operator()(const Layout& layout,
                                       const int n,
                                       const int m,
                                       std::vector<util::xor_func>& matrix,
                                       const std::vector<std::string>& qubit_names) = 0;
};

}

#endif //T_SCHEDULING_MATRIX_DECOMPOSER_HPP
