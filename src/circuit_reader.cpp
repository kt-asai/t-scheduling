#include <iostream>

#include "circuit_reader.hpp"

void CircuitReader::readQC_(Circuit& qc)
{
    std::cout << "start qc reader" << std::endl;
}

Circuit CircuitReader::read()
{
    Circuit qc = Circuit();

    std::string extension = getExtension_(path_);
    std::cout << "extension is " << extension << std::endl;

    if (extension == "qc")
    {
        readQC_(qc);
    }
    else
    {
        std::cout << "invalid format" << std::endl;

        exit(1);
    }

    return qc;
}