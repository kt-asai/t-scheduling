#include <iostream>

#include "circuit_reader.hpp"
#include "circuit.hpp"
#include "character.hpp"

int main()
{
    std::cout << "T-scheduling" << std::endl;

    std::string path("../benchmarks/tof_1.qc");

    CircuitReader reader(path);
    Circuit qc =  reader.read();

    std::cout << "-----------------------------" << std::endl;
    std::cout << "Original circuit" << std::endl;
    qc.print();

    Character character(qc);
    character.Parse();

    return 0;
}