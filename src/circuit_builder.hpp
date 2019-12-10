#ifndef T_SCHEDULING_CIRCUIT_BUILDER_HPP
#define T_SCHEDULING_CIRCUIT_BUILDER_HPP

#include <vector>
#include <list>
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

    std::vector<util::xor_func> bits_;
    std::vector<util::xor_func> preparation_;
    std::vector<util::xor_func> restoration_;

    bool Init(const std::vector<util::xor_func>& in,
              const std::vector<util::xor_func>& out);

    void Prepare(std::list<Gate>& gate_list,
                 const std::vector<util::xor_func>& in,
                 const int num_partition);

    void ApplyPhaseGates(std::list<Gate>& gate_list,
                         const std::set<int>& phase_exponent_index_set);

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

    std::list<Gate> Build(const tpar::partitioning& partition,
                          std::vector<util::xor_func>& in,
                          const std::vector<util::xor_func>& out);

    static std::list<Gate> BuildGlobalPhase(int qubit_num,
                                            int phase,
                                            const std::vector<std::string>& qubit_names);
};

}

#endif //T_SCHEDULING_CIRCUIT_BUILDER_HPP
