#include <iostream>
#include <fstream>
#include <algorithm>

#include "circuit_reader.hpp"

void CircuitReader::ReadQC_(Circuit& qc)
{
    std::ifstream ifs(path_);
    std::string str;

    if (ifs.fail())
    {
        std::cerr << "failed to open " << path_ << std::endl;
        exit(1);
    }

    while (getline(ifs, str))
    {
        // skip empty line
        if (str.empty()) continue;

        // preparation
        std::vector<std::string> buf = split_(str, ' ');
        std::string id = buf.front();
        buf.erase(buf.begin());

        // read parameters
        if (id == ".v")
        {
            for (std::string s : buf)
            {
                qc.add_qubit(s);
            }
        }
        else if (id == ".i")
        {
            continue;
        }
        else
        {
            if (id == "H" || id == "X")
            {
                qc.add_gate(id, buf.front());
            }
            if (id == "Z")
            {
                std::string target = buf.back();
                buf.pop_back();
                qc.add_gate("czz", buf, target);
            }
            if (id == "tof")
            {
                qc.add_gate("cnot", buf.front(), buf.back());
            }
        }
    }
}

Circuit CircuitReader::read()
{
    Circuit qc = Circuit();
    const std::string extension = extension_(path_);

    if (extension == "qc")
    {
        ReadQC_(qc);
        qc.DecomposeCZZ();
    }
    else
    {
        std::cerr << "invalid format" << std::endl;
        exit(1);
    }

    return qc;
}