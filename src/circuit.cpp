#include "circuit.hpp"

void Circuit::decomposeCZZ()
{
    for (auto it = gate_list_.begin(); it != gate_list_.end(); it++)
    {
        if (it->type() == "czz")
        {
            std::string control_a = it->control_list().front();
            std::string control_b = it->control_list().back();
            std::string target = it->target_list().front();

            // remove czz
            gate_list_.erase(it);

            /*
             *   {CNOT + T} tempalte of CZZ
             *   CNOT: 6
             *   T:    7
             *   T-depth: 5
             *   ---------+------------+-----+--T---+--
             *   --+------------+---------T--X--T*--X--
             *   --X--T*--X--T--X--T*--X--T------------
             */
            std::vector<Gate> czz = {Gate("cnot", control_b, target),
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

            for (Gate& gate : czz)
            {
                it = gate_list_.insert(it, gate);
                it++;
            }
        }
    }
}