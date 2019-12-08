#ifndef T_SCHEDULING_CIRCUIT_BUILDER_HPP
#define T_SCHEDULING_CIRCUIT_BUILDER_HPP

#include "circuit.hpp"
#include "synthesis.hpp"

namespace tskd {

/**
 * this class build {CNOT, T} subcircuit for given partitions
 */
class CircuitBuilder
{
private:
    int n_;

public:
    CircuitBuilder() = default;

    Circuit build();
};

}

#endif //T_SCHEDULING_CIRCUIT_BUILDER_HPP
