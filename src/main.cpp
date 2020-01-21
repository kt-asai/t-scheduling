#include <iostream>
#include <memory>

#include "util/option.hpp"
#include "io/circuit_reader.hpp"
#include "circuit/circuit.hpp"
#include "character/character.hpp"

#include "synthesis/synthesis_method_factory.hpp"
#include "synthesis/synthesis.hpp"

#include "parallel/parallelization_oracle.hpp"
#include "layout/layout.hpp"
#include "simulator/simulator.hpp"


int main(int argc, char** argv)
{
    std::cout << "T-scheduling" << std::endl;

    /*
     * Parse arguments
     */
    tskd::util::Option option;

    const std::string circuit_file = argv[1];
    const std::string path = "../benchmarks/" + circuit_file;

    const std::string syn_type = argv[2];
    if (syn_type == "tpar")
    {
        option.set_syn_method((SynthesisMethod::ktpar));
    }
    else if (syn_type == "tskd")
    {
        option.set_syn_method((SynthesisMethod::ktskd));
    }
    else
    {
        option.set_syn_method((SynthesisMethod::ktpar));
    }

    const std::string dec_type = argv[3];
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

    // TODO: fix bug
    const std::string change = argv[4];
    if (change == "true")
    {
        option.set_change_row_order(true);
    }
    else
    {
        option.set_change_row_order(false);
    }

    const std::string s_nd = argv[5];
    const int num_distillation = std::stoi(s_nd);
    option.set_num_distillation(num_distillation);

    // set number of buffer row
    option.set_num_buffer_row(0);
    const std::string nbr = argv[6];
    const int num_row = std::stoi(nbr);
    option.set_num_buffer_row(num_row);


    option.set_input_path(path);
    option.set_distillation_step(10);
    option.set_num_buffer(0);


    tskd::CircuitReader reader(option.input_path());
    tskd::Circuit qc = reader.read();
    std::cout << "-->> read file complete" << std::endl;

    std::cout << "# ----------------" << std::endl;
    tskd::Layout layout(option, qc);
    layout.print();
//    layout.print_edge();
    std::cout << "-->> construct layout complete" << std::endl;

    // set buffer
    const int num_buffer = option.num_buffer_row() * layout.width();
    option.set_num_buffer(num_buffer);
    option.show();

//    qc.print_gate_list();
    std::cout << "# ----------------" << std::endl;
    std::cout << "# Original circuit" << std::endl;
    tskd::Simulator sim_init(option, qc);
    sim_init.print();
    std::cout << "# ----------------" << std::endl;

    // test z3
//    tskd::ParallelizationOracle oracle(layout);
//    std::vector<std::string> targets = {"8", "9"};
//    tskd::Gate gate1("cnot", "1", "3");
//    tskd::Gate gate2("cnot", "2", "4");
//    std::list<tskd::Gate> gate_list = {gate1, gate2};
//    bool result = oracle.check(gate_lisct);
//    std::cout << "result:" << result << std::endl;

    tskd::Character chr(qc);
    chr.parse();
    std::cout << "-->> construct character complete" << std::endl;

    const auto start = std::chrono::system_clock::now();
    std::shared_ptr<tskd::Synthesis> synthesis = tskd::SynthesisMethodFactory().create(option.syn_method(), option, layout, chr);
    tskd::Circuit result = synthesis->execute();
    const auto end = std::chrono::system_clock::now();
    std::cout << "-->> synthsis complete" << std::endl;

//    result.print_gate_list();
    std::cout << "# ----------------" << std::endl;
    std::cout << "# Optimized circuit" << std::endl;
    tskd::Simulator sim_result(option, result);
    sim_result.print();
    const auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << "# " << elapsed_time.count() << " sec" << std::endl;
    std::cout << "# ----------------" << std::endl;

    return 0;
}