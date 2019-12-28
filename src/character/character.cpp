#include <iostream>
#include <string>

#include "character.hpp"

namespace tskd {

int Character::insert_phase(const int coefficient,
                            const util::xor_func& function)
{
    int index = 0;
    for (auto& phase_exponent : phase_exponents_)
    {
        if (phase_exponent.second == function)
        {
            phase_exponent.first = (phase_exponent.first + coefficient) % 8;

            return index;
        }
        index++;
    }

    phase_exponents_.emplace_back(std::make_pair(coefficient, util::xor_func(function)));

    return index;
}

void Character::parse()
{
    // gate value map
    const std::unordered_map<std::string, int> gate_lookup{
            {"T",  1},
            {"T*", 7},
            {"P",  2},
            {"P*", 6},
            {"Z",  4},
            {"Y",  4}
    };

    int name_max = 0;
    int value_max = 0;
    int gate_index = 0;
    std::vector<util::xor_func> parity_list;
    std::unordered_map<std::string, int> name_map;
    auto wires = std::vector<util::xor_func>(num_qubit_);

    /**
     * initialization
     */
    for (const std::string& name : circuit_.qubit_names())
    {
        name_map.insert(std::make_pair(name, name_max));
        qubit_names_[name_max] = name;
        ancilla_list_[name_max] = circuit_.is_ancilla_map().at(name);
        wires[name_max] = util::xor_func((num_qubit_ - num_ancilla_) + num_hadamard_ + 1, 0);

        if (!ancilla_list_[name_max])
        {
            wires[name_max].set(value_max);
            value_map_[value_max] = name_max;
            value_max++;
        }
        name_max++;
    }

    /**
     *  compute phase exponent and hadamard info
     */
    for (const Gate& gate : circuit_.gate_list())
    {
        if (gate.type() == "cnot")
        {
            wires[name_map.at(gate.target_list().front())] ^= wires[name_map[gate.control_list().front()]];
        }
        else if (gate.type() == "X")
        {
            wires[name_map.at(gate.target_list().front())].flip((num_qubit_ - num_ancilla_) + num_hadamard_);
        }
        else if (gate.type() == "Y")
        {
            gate_index = name_map.at(gate.target_list().front());
            insert_phase(gate_lookup.at(gate.type()), wires[gate_index]);
            wires[name_map.at(gate.target_list().front())].flip((num_qubit_ - num_ancilla_) + num_hadamard_);
        }
        else if (gate.type() == "T" || gate.type() == "T*"
                 || gate.type() == "P" || gate.type() == "P*"
                 || gate.type() == "Z")
        {
            gate_index = name_map.at(gate.target_list().front());
            insert_phase(gate_lookup.at(gate.type()), wires[gate_index]);
        }
        else if (gate.type() == "H")
        {
            // hadamard process
            Hadamard new_hadamard(name_map.at(gate.target_list().front()), value_max, wires);
            value_max++;

            // compute rank
            wires[new_hadamard.target_].reset();

            util::compute_rank_destructive(num_qubit_, (num_qubit_ - num_ancilla_) + num_hadamard_, wires);
            int index = 0;
            for (const auto& phase_exponent : phase_exponents_)
            {
                if (phase_exponent.first != 0)
                {
                    if (util::is_independent((num_qubit_ - num_ancilla_) + num_hadamard_, wires, phase_exponent.second))
                    {
                        new_hadamard.in_.insert(index);
                    }
                }
                index++;
            }

            // reset the current wire values
            wires = new_hadamard.input_wires_parity_;

            // done creating the new hadamard
            hadamards_.push_back(new_hadamard);

            // prepare the new value
            wires[new_hadamard.target_].reset();
            wires[new_hadamard.target_].set(new_hadamard.previous_qubit_index_);

            // give this value a name
            value_map_.insert(std::make_pair(new_hadamard.previous_qubit_index_, name_max));
            qubit_names_[name_max] = qubit_names_[new_hadamard.target_];
            qubit_names_[name_max].append(std::to_string(new_hadamard.previous_qubit_index_));
            name_max++;
        }
        else
        {
            std::cerr << "invalid gate" << std::endl;

            exit(1);
        }
    }

    outputs_ = std::move(wires);
}

}