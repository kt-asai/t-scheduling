#include <iostream>

#include "util/option.hpp"
#include "io/circuit_reader.hpp"
#include "circuit/circuit.hpp"
#include "character/character.hpp"

#include "synthesis/synthesis_method_factory.hpp"
#include "synthesis/synthesis.hpp"

#include "parallel/parallelization_oracle.hpp"
#include "layout/layout.hpp"

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
    option.set_num_distillation(2);
    option.set_distillation_step(1000);
    option.set_syn_method((SynthesisMethod::ktskd));
    option.show();


    tskd::CircuitReader reader(option.input_path());
    tskd::Circuit qc = reader.read();
    std::cout << "-->> read file complete" << std::endl;

    std::cout << "# ----------------" << std::endl;
    std::cout << "# Original circuit" << std::endl;
//    qc.print();
//    qc.print_gate_list();
    std::cout << "# ----------------" << std::endl;

    tskd::Layout layout(option, qc);
    layout.print();
//    layout.print_edge();
    std::cout << "-->> construct layout complete" << std::endl;

    // test z3
    tskd::ParallelizationOracle oracle(layout);
//    std::vector<std::string> targets = {"8", "9"};
    tskd::Gate gate1("cnot", "1", "3");
    tskd::Gate gate2("cnot", "2", "4");
    std::vector<tskd::Gate> gate_list = {gate1, gate2};
    bool result = oracle.check(gate_list);
    std::cout << "result:" << result << std::endl;

//    tskd::Character chr(qc);
//    chr.parse();
//    std::cout << "-->> construct character complete" << std::endl;
//
//    tskd::Synthesis* synthesis = tskd::SynthesisMethodFactory().create(option.syn_method(), option, layout, chr);
//    tskd::Circuit result = synthesis->execute();
//    std::cout << "-->> synthsis complete" << std::endl;
//
//    std::cout << "# ----------------" << std::endl;
//    std::cout << "# Optimized circuit" << std::endl;
//    result.print();
////    result.print_gate_list();
//    std::cout << "# ----------------" << std::endl;

    return 0;
}