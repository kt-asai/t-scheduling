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
    Layout layout_;

    int n_;
    int m_;

    std::vector<std::string> qubit_names_;

public:
    MatrixDecomposer() = default;

    MatrixDecomposer(const Layout& layout,
                     const int n,
                     const int m,
                     const std::vector<std::string>& qubit_names)
        : layout_(layout),
          n_(n),
          m_(m),
          qubit_names_(qubit_names) { }

    virtual ~MatrixDecomposer() = default;

    Layout layout() const
    {
        return layout_;
    }

    int n() const
    {
        return n_;
    }

    int m() const
    {
        return m_;
    }

    std::vector<std::string> qubit_names() const
    {
        return qubit_names_;
    }

    virtual std::list<Gate> execute(std::vector<util::xor_func>& matrix,
                                    std::vector<int>& func_map) = 0;
};

}

#endif //T_SCHEDULING_MATRIX_DECOMPOSER_HPP
