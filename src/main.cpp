#include <iostream>

#include "circuit_reader.hpp"

int main()
{
    std::cout << "T-scheduling" << std::endl;

    std::string path("../benchmarks/tof_3.qc");

    CircuitReader reader(path);
    reader.read();

    return 0;
}