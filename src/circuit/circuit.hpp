#ifndef T_SCHEDULING_CIRCUIT_HPP
#define T_SCHEDULING_CIRCUIT_HPP

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <unordered_map>
#include <algorithm>

#include "gate.hpp"

namespace tskd {

class Circuit
{
private:
    int num_qubit_;
    int num_ancilla_qubit_;
    int num_gate_;

    std::vector<std::string> qubit_names_;
    std::list<Gate> gate_list_;
    std::unordered_map<std::string, bool> is_ancilla_map_;

    bool equal_gate(const Gate& gate_a,
                    const Gate& gate_b);

public:
    /**
     * constructor
     */
    Circuit()
        : num_qubit_(0),
          num_ancilla_qubit_(0),
          num_gate_(0) { }

    /**
     * return number of qubit in circuit
     * @return number of qubit
     */
    int qubit_num() const
    {
        return num_qubit_;
    }

    int ancilla_qubit_num() const
    {
        return num_ancilla_qubit_;
    }

    /**
     * return number of gate in circuit
     * @return number of gate
     */
    int gate_num() const
    {
        return num_gate_;
    }

    /**
     * count specified gates
     * @param type gate name
     * @return number of specified gate
     */
    int count_gate(const std::string& type) const
    {
        const int gate_count = std::count_if(gate_list_.begin(), gate_list_.end(),
                                             [type](const Gate& gate)
                                             {
                                                 return type == gate.type();
                                             });
        return gate_count;
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
     * return info of whether each qubit is ancilla
     * @return ancilla map
     */
    std::unordered_map<std::string, bool> is_ancilla_map() const
    {
        return is_ancilla_map_;
    }

    /**
     * return gate list in the circuit
     * @return array of gate
     */
    std::list<Gate> gate_list() const
    {
        return gate_list_;
    }

    /**
     * set whether ancilla qubit
     * @param qubit qubit name
     * @param is_ancilla whether qubit is ancilla qubit
     */
    void set_ancilla(const std::string& qubit,
                     const bool is_ancilla = true)
    {
        const size_t exist = is_ancilla_map_.count(qubit);
        if (exist)
        {
            num_ancilla_qubit_ += static_cast<int>(is_ancilla) - static_cast<int>(is_ancilla_map_[qubit]);
            is_ancilla_map_[qubit] = is_ancilla;
        }
        else
        {
            is_ancilla_map_.emplace(qubit, is_ancilla);
            num_ancilla_qubit_ += static_cast<int>(is_ancilla);
        }
    }

    /**
     * add qubit to the circuit
     * @param qubit qubit name
     */
    void add_qubit(const std::string& qubit)
    {
        qubit_names_.push_back(qubit);
        num_qubit_++;
    }

    /**
     * add gate list
     * @param gate_list gate list
     */
    void add_gate_list(std::list<Gate>&& gate_list)
    {
        num_gate_ += static_cast<int>(gate_list.size());
        gate_list_.splice(gate_list_.end(), std::forward<std::list<Gate>>(gate_list));
    }

    /**
     * add Gate object to the circuit list
     * @param gate gate object
     */
    void add_gate(const Gate& gate)
    {
        gate_list_.push_back(gate);
        num_gate_++;
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
        num_gate_++;
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
        num_gate_++;
    }

    /**
     * add multi-target gate (e.g.multi-target cnot) to the circuit
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
        num_gate_++;
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
        num_gate_++;
    }

    /**
     * add single qubit gate (e.g.H, Z, T) to the circuit
     * @param type gate name
     * @param target target qubit name
     */
    void add_gate(const std::string& type,
                  const std::string& target)
    {
        Gate gate(type, target);
        gate_list_.push_back(gate);
        num_gate_++;
    }

    /**
     * remove identity gates
     */
    void remove_identities();

    /**
     * decompose a czz gate to {CNOT, T} gates
     */
    void decompose_ccz();

    /**
     * print circuit status
     */
    void print()
    {
        // counter
        int num_h = count_gate("H");
        int num_x = count_gate("X");
        int num_t = count_gate("T") + count_gate("T*");
        int num_p = count_gate("P") + count_gate("P*");
        int num_z = count_gate("Z");
        int num_toffoli = count_gate("ccz");
        int num_cnot = count_gate("cnot") + count_gate("tof");

        std::cout << "# qubits: " << num_qubit_ << std::endl;
        std::cout << "# ancilla: " << num_ancilla_qubit_ << std::endl;
        std::cout << "# gates: " << num_gate_ << std::endl;
        std::cout << "# H: " << num_h << std::endl;
        std::cout << "# CNOT: " << num_cnot << std::endl;
        std::cout << "# X: " << num_x << std::endl;
        std::cout << "# T: " << num_t << std::endl;
        std::cout << "# P: " << num_p << std::endl;
        std::cout << "# Z: " << num_z << std::endl;
        std::cout << "# Toffoli: " << num_toffoli << std::endl;
        std::cout << "# T-depth: " << count_t_depth() << std::endl;
    }

    /**
     * print gate list in the circuit
     */
    void print_gate_list()
    {
        for (auto&& gate : gate_list_)
        {
            gate.print();
        }
    }

    /**
     * count t depth
     */
    int count_t_depth()
    {
        std::map<std::string, int> t_map;
        int t_depth = 0;

        // init
        for (auto&& name : qubit_names_)
        {
            t_map.insert(std::make_pair(name, 0));
        }

        // count t gate applied each qubit
        for (auto&& gate : gate_list_)
        {
            if (gate.type() == "T" || gate.type() == "T*")
            {
                t_map[gate.target_list().front()]++;
            }
        }

        // count T depth
        for (auto&& e : t_map)
        {
            t_depth = std::max(t_depth, e.second);
        }

        return t_depth;
    }
};

}

#endif //T_SCHEDULING_CIRCUIT_HPP
