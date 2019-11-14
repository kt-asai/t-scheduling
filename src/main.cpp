#include <iostream>

#include "circuit_reader.hpp"
#include "circuit.hpp"

int main()
{
    std::cout << "T-scheduling" << std::endl;

    std::string path("../benchmarks/tof_3.qc");

    CircuitReader reader(path);
    Circuit qc =  reader.read();

    std::cout << "-----------------------------" << std::endl;
    std::cout << "Original circuit" << std::endl;
    qc.printStatus();
//    qc.printGateList();

    std::cout << "-----------------------------" << std::endl;
    std::cout << "Decomposed circuit" << std::endl;
    qc.decomposeCZZ();
    qc.printStatus();
//    qc.printGateList();

    return 0;
}