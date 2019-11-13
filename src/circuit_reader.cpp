#include <iostream>
#include <fstream>
#include <algorithm>

#include "circuit_reader.hpp"

void CircuitReader::readQC_(Circuit& qc)
{
    std::cout << "start qc reader" << std::endl;

    std::ifstream ifs(path_);
    std::string str;

    if (ifs.fail())
    {
        std::cerr << "failed to open " << path_ << std::endl;

        exit(1);
    }

    while (getline(ifs, str))
    {
        // remove space
        str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());

        // read parameters
        if (str.find(".v") != std::string::npos)
        {
            std::cout << ".v line parse" << std::endl;
            std::cout << str << std::endl;
        }
        else if (str.find(".i") != std::string::npos)
        {
            std::cout << ".i line parse" << std::endl;
            std::cout << str << std::endl;
        }
        else
        {
            std::cout << "gate parse" << std::endl;
            std::cout << str << std::endl;
        }
    }
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
        std::cerr << "invalid format" << std::endl;

        exit(1);
    }

    return qc;
}