#include <iostream>

#include "circuit_reader.hpp"
#include "circuit.hpp"
#include "character.hpp"

int main()
{
    std::cout << "T-scheduling" << std::endl;

    std::string path("../benchmarks/sample.qc");

    std::cout << "-->> read file" << std::endl;
    tskd::CircuitReader reader(path);
    tskd::Circuit qc = reader.read();

    std::cout << "# ----------------" << std::endl;
    std::cout << "# Original circuit" << std::endl;
    qc.print();
    qc.print_gate_list();

    std::cout << "-->> construct character" << std::endl;
    tskd::Character character(qc);
    character.Parse();

    return 0;
}