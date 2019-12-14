#include <iostream>

#include "util/option.hpp"
#include "io/circuit_reader.hpp"
#include "circuit/circuit.hpp"
#include "character.hpp"

#include "synthesis/synthesis_method_factory.hpp"
#include "synthesis/synthesis.hpp"
#include "synthesis/tpar_synthesis.hpp"


int main(int argc, char** argv)
{
    std::cout << "T-scheduling" << std::endl;

    /*
     * Parse arguments
     */
    tskd::util::Option option;
    option.set_input_path("../benchmarks/tof_2.qc");
    option.set_num_distillation(4);
    option.set_distillation_step(10);
    option.set_change_row_order(true);
    option.set_syn_method((SynthesisMethod::ktskd));
    option.set_part_type(PartitionType::kmatroid);
    option.set_dec_type(DecompositionType::kgauss);
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
//    result.print_gate_list();

    return 0;
}