#ifndef T_SCHEDULING_CIRCUIT_BUILDER_HPP
#define T_SCHEDULING_CIRCUIT_BUILDER_HPP

#include <vector>
#include <string>

#include "gate.hpp"
#include "circuit.hpp"
#include "synthesis.hpp"

namespace tskd {

/**
 * this class build {CNOT, T} subcircuit for given partitions
 */
class CircuitBuilder
{
private:
    int qubit_num_;
    int dimension_;

    std::vector<std::string> qubit_names_;
    std::vector<util::phase_exponent> phase_exponent_;

public:
    CircuitBuilder() = default;

    CircuitBuilder(int qubit_num,
                   int dimension,
                   const std::vector<std::string>& qubit_names,
                   const std::vector<util::phase_exponent>& phase_exponent)
        : qubit_num_(qubit_num),
          dimension_(dimension),
          qubit_names_(qubit_names),
          phase_exponent_(phase_exponent) { }

    std::vector<Gate> build(const tpar::partitioning& partition,
                            std::vector<util::xor_func>& in,
                            const std::vector<util::xor_func>& out);
};

}

#endif //T_SCHEDULING_CIRCUIT_BUILDER_HPP
