#ifndef T_SCHEDULING_CHARACTER_HPP
#define T_SCHEDULING_CHARACTER_HPP

#include <vector>
#include <string>
#include <map>
#include <unordered_map>

#include "circuit.hpp"
#include "gate.hpp"
#include "util.hpp"

namespace tskd {

class Character
{
private:
    struct Hadamard
    {
        int target_;
        int previous_qubit_index_;
        std::vector<util::xor_func> input_wires_parity_;
        std::set<int> in_;

        Hadamard(int target,
                 int previous_qubit_index,
                 std::vector<util::xor_func>& input_wires_parity)
                : target_(target),
                  previous_qubit_index_(previous_qubit_index),
                  input_wires_parity_(input_wires_parity) { }
    };

    Circuit circuit_;

    int num_qubit_;
    int num_ancilla_;
    int num_hadamard_;
    std::vector<std::string> qubit_names_;
    std::vector<bool> ancilla_list_;
    std::map<int, int> value_map_;
    std::vector<util::phase_exponent> phase_exponents_;
    std::vector<util::xor_func> outputs_;
    std::vector<Hadamard> hadamards_;

    int InsertPhase_(int coefficient,
                     const util::xor_func& function);

public:
    /**
     * constructor
     */
    Character() = default;

    /**
     * constructor
     * @param qc circuit info
     */
    Character(const Circuit& circuit)
            : circuit_(circuit),
              num_qubit_(circuit_.qubit_num()),
              num_ancilla_(circuit_.ancilla_qubit_num()),
              num_hadamard_(circuit_.count_gate("H"))
    {
        qubit_names_.reserve(num_qubit_ + num_hadamard_);
        ancilla_list_.reserve(num_qubit_);
    }

    /**
     * parse circuit info and calculate parity
     */
    void Parse();
};

}

#endif //T_SCHEDULING_CHARACTER_HPP
