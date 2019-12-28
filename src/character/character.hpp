#ifndef T_SCHEDULING_CHARACTER_HPP
#define T_SCHEDULING_CHARACTER_HPP

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <set>

#include "../circuit/circuit.hpp"
#include "../circuit/gate.hpp"

#include "../util/util.hpp"

namespace tskd {

class Character
{
public:
    struct Hadamard
    {
        int target_;
        int previous_qubit_index_;
        std::vector<util::xor_func> input_wires_parity_;
        std::set<int> in_;

        Hadamard() { }

        Hadamard(int target,
                 int previous_qubit_index,
                 std::vector<util::xor_func>& input_wires_parity)
                : target_(target),
                  previous_qubit_index_(previous_qubit_index),
                  input_wires_parity_(input_wires_parity) { }
    };

private:
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

    int insert_phase(int coefficient,
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
        qubit_names_.resize(num_qubit_ + num_hadamard_);
        ancilla_list_.resize(num_qubit_);
    }

    /**
     * return number of qubit in circuit
     * @return number of qubit
     */
    int num_qubit() const
    {
        return num_qubit_;
    }

    /**
     * return number of data qubit in circuit
     * @return number of data qubit
     */
    int num_data_qubit() const
    {
        return num_qubit_ - num_ancilla_;
    }

    /**
     * return number of qubit in circuit
     * @return number of ancilla qubit
     */
    int num_ancilla_qubit() const
    {
        return num_ancilla_;
    }

    /**
     * return number of hadamard gate
     * @return number of hadamard gate
     */
    int num_hadamard() const
    {
        return num_hadamard_;
    }

    /**
     * return qubit names in the circuit
     * @return array of qubit name
     */
    std::vector<std::string> qubit_names() const
    {
        return qubit_names_;
    }

    /**
     * return ancilla list
     * @return ancilla list
     */
    std::vector<bool> ancilla_list() const
    {
        return ancilla_list_;
    }

    /**
     * return phase exponents
     * @return phase exponents
     */
    std::vector<util::phase_exponent> phase_exponents() const
    {
        return phase_exponents_;
    }

    /**
     * return hadamard gate list
     * @return hadmard gates
     */
    std::vector<Hadamard> hadamards() const
    {
        return hadamards_;
    }

    /**
     * return output parities of circuit
     * @return output parities of circuit
     */
    std::vector<util::xor_func> outputs() const
    {
        return outputs_;
    }

    /**
     * parse circuit info and calculate parity
     */
    void parse();
};

}

#endif //T_SCHEDULING_CHARACTER_HPP
