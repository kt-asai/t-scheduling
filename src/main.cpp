#include <iostream>

#include "io/circuit_reader.hpp"
#include "circuit.hpp"
#include "character.hpp"
#include "synthesis.hpp"

int main()
{
    std::cout << "T-scheduling" << std::endl;

    const std::string path("../benchmarks/tof_3.qc");

    std::cout << "-->> read file" << std::endl;
    tskd::CircuitReader reader(path);
    tskd::Circuit qc = reader.read();

    std::cout << "# ----------------" << std::endl;
    std::cout << "# Original circuit" << std::endl;
    qc.print();
//    qc.print_gate_list();

    std::cout << "-->> construct character" << std::endl;
    tskd::Character character(qc);
    character.Parse();

    std::cout << "-->> synthsis" << std::endl;
    tskd::Synthesis synthesis(character);
    tskd::Circuit result = synthesis.Execute();

    std::cout << "# ----------------" << std::endl;
    std::cout << "# Optimized circuit" << std::endl;
    result.print();
//    result.print_gate_list();

    return 0;
}