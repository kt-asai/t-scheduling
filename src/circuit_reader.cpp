#include <iostream>
#include <fstream>
#include <algorithm>

#include "circuit_reader.hpp"

namespace tskd {

void CircuitReader::ReadQC_(Circuit& circuit)
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
            for (const std::string& s : buf)
            {
                circuit.add_qubit(s);
                circuit.set_ancilla(s);
            }
        }
        else if (id == ".i")
        {
            for (const std::string& s : buf)
            {
                circuit.set_ancilla(s, false);
            }
        }
        else
        {
            if (id == "H" || id == "X")
            {
                circuit.add_gate(id, buf.front());
            }
            if (id == "Z")
            {
                std::string target = buf.back();
                buf.pop_back();
                circuit.add_gate("czz", buf, target);
            }
            if (id == "tof")
            {
                circuit.add_gate("cnot", buf.front(), buf.back());
            }
        }
    }
}

Circuit CircuitReader::read()
{
    Circuit circuit = Circuit();
    const std::string extension = extension_(path_);

    if (extension == "qc")
    {
        ReadQC_(circuit);
        circuit.DecomposeCZZ();
    }
    else
    {
        std::cerr << "invalid format" << std::endl;

        exit(1);
    }

    circuit.RemoveIdentities();

    return circuit;
}

}