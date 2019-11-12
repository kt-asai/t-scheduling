#include "circuit_reader.hpp"
#include "circuit.hpp"

Circuit CircuitReader::read()
{
    Circuit qc();

    return std::move(qc);
}