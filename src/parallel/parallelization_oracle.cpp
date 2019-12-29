#include "z3++.h"

#include "parallelization_oracle.hpp"

namespace tskd {

void ParallelizationOracle::init()
{

}

bool ParallelizationOracle::check(const std::vector<Gate>& gate_list)
{
    const int num_gate = static_cast<int>(gate_list.size());
    constexpr int num_distillation = 4;



    z3::context c;
    z3::expr x = c.int_const("x");
    z3::expr y = c.int_const("y");
    z3::solver s(c);

    s.add(x >= 1);
    s.add(y < x + 3);
    std::cout << s.check() << "\n";

    z3::model m = s.get_model();
    std::cout << m << "\n";
    // traversing the model
    for (unsigned i = 0; i < m.size(); i++) {
        z3::func_decl v = m[i];
        // this problem contains only constants
        assert(v.arity() == 0);
        std::cout << v.name() << " = " << m.get_const_interp(v) << "\n";
    }
    // we can evaluate expressions in the model.
    std::cout << "x + y + 1 = " << m.eval(x + y + 1) << "\n";

    return false;
}

}