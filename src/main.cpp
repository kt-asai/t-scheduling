#include <iostream>

#include "util/option.hpp"
#include "io/circuit_reader.hpp"
#include "circuit/circuit.hpp"
#include "character/character.hpp"

#include "synthesis/synthesis_method_factory.hpp"
#include "synthesis/synthesis.hpp"

#include "parallel/parallelization_oracle.hpp"

//#include "z3++.h"
//
//using namespace z3;
//
//int main()
//{
//    std::cout << "unsat core example1\n";
//    context c;
//    // We use answer literals to track assertions.
//    // An answer literal is essentially a fresh Boolean marker
//    // that is used to track an assertion.
//    // For example, if we want to track assertion F, we
//    // create a fresh Boolean variable p and assert (p => F)
//    // Then we provide p as an argument for the check method.
//    expr p1 = c.bool_const("p1");
//    expr p2 = c.bool_const("p2");
//    expr p3 = c.bool_const("p3");
//    expr x  = c.int_const("x");
//    expr y  = c.int_const("y");
//    solver s(c);
//    s.add(implies(p1, x > 10));
//    s.add(implies(p1, y > x));
//    s.add(implies(p2, y < 5));
//    s.add(implies(p3, y > 0));
//    expr assumptions[3] = { p1, p2, p3 };
//    std::cout << s.check(3, assumptions) << "\n";
//    expr_vector core = s.unsat_core();
//    std::cout << core << "\n";
//    std::cout << "size: " << core.size() << "\n";
//    for (unsigned i = 0; i < core.size(); i++) {
//        std::cout << core[i] << "\n";
//    }
//    // Trying again without p2
//    expr assumptions2[2] = { p1, p3 };
//    std::cout << s.check(2, assumptions2) << "\n";
//}

int main(int argc, char** argv)
{
    std::cout << "T-scheduling" << std::endl;

    /*
     * Parse arguments
     */
    tskd::util::Option option;

    const std::string input = argv[1];
    const std::string path = "../benchmarks/" + input;

    const std::string dec_type = argv[2];
    if (dec_type == "gauss")
    {
        option.set_dec_type(DecompositionType::kgauss);
    }
    else if (dec_type == "par")
    {
        option.set_dec_type(DecompositionType::kparallel);
    }
    else
    {
        option.set_dec_type(DecompositionType::kgauss);
    }

    const std::string change = argv[3];
    if (change == "true")
    {
        option.set_change_row_order(true);
    }
    else
    {
        option.set_change_row_order(false);
    }

    option.set_input_path(path);
    option.set_num_distillation(4);
    option.set_distillation_step(1000);
    option.set_syn_method((SynthesisMethod::ktskd));
    option.show();

    std::cout << "-->> read file" << std::endl;
    tskd::CircuitReader reader(option.input_path());
    tskd::Circuit qc = reader.read();

    std::cout << "# ----------------" << std::endl;
    std::cout << "# Original circuit" << std::endl;
    qc.print();
//    qc.print_gate_list();

    std::cout << "-->> construct character" << std::endl;
    tskd::Character chr(qc);
    chr.parse();

    std::cout << "-->> synthsis" << std::endl;
    tskd::Synthesis* synthesis = tskd::SynthesisMethodFactory().create(option.syn_method(), chr, option);
    tskd::Circuit result = synthesis->execute();

    std::cout << "# ----------------" << std::endl;
    std::cout << "# Optimized circuit" << std::endl;
    result.print();
    result.print_gate_list();

    return 0;
}