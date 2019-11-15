#ifndef T_SCHEDULING_CIRCUIT_HPP
#define T_SCHEDULING_CIRCUIT_HPP

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "gate.hpp"

class Circuit
{
private:
    std::vector<std::string> qubit_names_;
    std::vector<Gate> gate_list_;

public:
    /**
     * constructor
     */
    Circuit() = default;

    /**
     * return number of qubit in circuit
     * @return number of qubit
     */
    int qubit_num() const
    {
        return static_cast<int>(qubit_names_.size());
    }

    /**
     * return number of gate in circuit
     * @return number of gate
     */
    int gate_num() const
    {
        return static_cast<int>(gate_list_.size());
    }

    /**
     * return qubit names in the circuit
     * @return array of qubit name
     */
    const std::vector<std::string> qubit_names() const
    {
        return qubit_names_;
    }

    /**
     * return gate list in the circuit
     * @return array of gate
     */
    const std::vector<Gate> gate_list() const
    {
        return gate_list_;
    }

    /**
     * add qubit to the circuit
     * @param qubit qubit name
     */
    void add_qubit(const std::string& qubit)
    {
        qubit_names_.push_back(qubit);
    }

    /**
     * add multi-control and multi-target gate to the circuit
     * @param type gate name
     * @param control_list control qubit names
     * @param target_list target qubit names
     */
    void add_qate(const std::string& type,
                 const std::vector<std::string>& control_list,
                 const std::vector<std::string>& target_list)
    {
        Gate gate(type, control_list, target_list);
        gate_list_.push_back(gate);
    }

    /**
     * add multi-control gate (e.g.toffoli) to the circuit
     * @param type gate name
     * @param control_list control qubit names
     * @param target target qubit name
     */
    void add_gate(const std::string& type,
                 const std::vector<std::string>& control_list,
                 const std::string& target)
    {
        Gate gate(type, control_list, target);
        gate_list_.push_back(gate);
    }

    /**
     * add muti-target gate (e.g.multi-target cnot) to the circuit
     * @param type gate name
     * @param control control qubit name
     * @param target_list target qubit names
     */
    void add_gate(const std::string& type,
                 const std::string& control,
                 const std::vector<std::string>& target_list)
    {
        Gate gate(type, control, target_list);
        gate_list_.push_back(gate);
    }

    /**
     * add 2 qubit gate (e.g.cnot, cz) to the circuit
     * @param type gate name
     * @param control control qubit name
     * @param target target qubit name
     */
    void add_gate(const std::string& type,
                 const std::string& control,
                 const std::string& target)
    {
        Gate gate(type, control, target);
        gate_list_.push_back(gate);
    }

    /**
     * add single qubit gate (e.g.T, P) to the circuit
     * @param type gate name
     * @param target target qubit name
     */
    void add_gate(const std::string& type,
                 const std::string& target)
    {
        Gate gate(type, target);
        gate_list_.push_back(gate);
    }

    /**
     * decompose a czz gate to {CNOT, T} gates
     */
    void DecomposeCZZ();

    /**
     * print circuit status
     */
    void print()
    {
        // counter
        int h_num = std::count_if(gate_list_.begin(), gate_list_.end(),
                [](Gate& gate) { return gate.type() == "H"; });
        int t_num = std::count_if(gate_list_.begin(), gate_list_.end(),
                [](Gate& gate) { return gate.type() == "T" || gate.type() == "T*"; });
        int x_num = std::count_if(gate_list_.begin(), gate_list_.end(),
                [](Gate& gate) { return gate.type() == "X"; });
        int toffoli_num = std::count_if(gate_list_.begin(), gate_list_.end(),
                [](Gate& gate) { return gate.type() == "toffoli"; });
        int cnot_num = std::count_if(gate_list_.begin(), gate_list_.end(),
                [](Gate& gate) { return gate.type() == "cnot"; });

        std::cout << "qubits: " << qubit_num() << std::endl;
        std::cout << "gates: " << gate_num() << std::endl;
        std::cout << "T: " << t_num << std::endl;
        std::cout << "H: " << h_num << std::endl;
        std::cout << "X: " << x_num << std::endl;
        std::cout << "CNOT: " << cnot_num << std::endl;
        std::cout << "Toffoli: " << toffoli_num << std::endl;
    }

    /**
     * print gate list in the circuit
     */
    void print_gate_list()
    {
        for (Gate& gate : gate_list_)
        {
            gate.print();
        }
    }
};

#endif //T_SCHEDULING_CIRCUIT_HPP
