#include <algorithm>

#include "circuit.hpp"

namespace tskd {

bool Circuit::EqualBits_(const Gate& gate_a,
                         const Gate& gate_b)
{
    return true;
}

void Circuit::RemoveIdentities()
{
    std::unordered_map<std::string, std::string> identity_map{
            {"cnot", "cnot"},
            {"Z", "Z"},
            {"X", "X"},
            {"Y", "Y"},
            {"H", "H"},
            {"P", "P*"},
            {"P*", "P"},
            {"T", "T*"},
            {"T*", "T"}
    };

    bool update = true;
    bool equal_bits = false;

    /*
     * search the identity matrix until the gate is no longer removed
     */
    while (update)
    {
        update = false;
        for (auto it_a = gate_list_.begin(); it_a != gate_list_.end(); it_a++)
        {
            auto it_b = it_a + 1;
            for (; it_b != gate_list_.end(); it_b++)
            {

                equal_bits = EqualBits_(*it_a, *it_b);
                if (equal_bits && identity_map[it_a->type()] == it_b->type())
                {
                    gate_list_.erase(it_a);
                    gate_list_.erase(it_b);
                    update = true;
                    break;
                }
                if (equal_bits && identity_map[it_a->type()] != it_b->type())
                {
                    break;
                }
            }
        }
    }
}

void Circuit::DecomposeCZZ()
{
    for (auto it = gate_list_.begin(); it != gate_list_.end(); it++)
    {
        if (it->type() == "czz")
        {
            const std::string control_a = it->control_list().front();
            const std::string control_b = it->control_list().back();
            const std::string target = it->target_list().front();

            // remove czz
            gate_list_.erase(it);
            num_gate_--;

            /*
             *   {CNOT + T} tempalte of CZZ
             *   CNOT: 6
             *   T:    7
             *   T-depth: 5
             *   ---------+------------+-----+--T---+--
             *   --+------------+---------T--X--T*--X--
             *   --X--T*--X--T--X--T*--X--T------------
             */
            const std::vector<Gate> czz = {Gate("cnot", control_b, target),
                                           Gate("T*", target),
                                           Gate("cnot", control_a, target),
                                           Gate("T", target),
                                           Gate("cnot", control_b, target),
                                           Gate("T*", target),
                                           Gate("cnot", control_a, target),
                                           Gate("T", control_b),
                                           Gate("T", target),
                                           Gate("cnot", control_a, control_b),
                                           Gate("T", control_a),
                                           Gate("T*", control_b),
                                           Gate("cnot", control_a, control_b)};

            for (const Gate& gate : czz)
            {
                it = gate_list_.insert(it, gate);
                it++;
                num_gate_++;
            }
        }
    }
}

}