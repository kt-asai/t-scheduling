#include <iostream>

#include "util/option.hpp"
#include "io/circuit_reader.hpp"
#include "circuit/circuit.hpp"
#include "character/character.hpp"

#include "synthesis/synthesis_method_factory.hpp"
#include "synthesis/synthesis.hpp"


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
    chr.Parse();

    std::cout << "-->> synthsis" << std::endl;
    tskd::Synthesis* synthesis = tskd::SynthesisMethodFactory().Create(option.syn_method(), chr, option);
    tskd::Circuit result = synthesis->Execute();

    std::cout << "# ----------------" << std::endl;
    std::cout << "# Optimized circuit" << std::endl;
    result.print();
    result.print_gate_list();

    return 0;
}