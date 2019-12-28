#include <algorithm>

#include "circuit.hpp"

namespace tskd {

bool Circuit::equal_gate(const Gate& gate_a,
                         const Gate& gate_b)
{
    if (gate_a.type() != gate_b.type())
    {
        return false;
    }

    bool equal_bits = true;
    const int num_control_bit = static_cast<int>(gate_a.control_list().size());
    const int num_target_bit = static_cast<int>(gate_a.target_list().size());

    for (int i = 0; i < num_control_bit; ++i)
    {
        if (gate_a.control_list()[i] != gate_b.control_list()[i])
        {
            equal_bits = false;
        }
    }

    for (int i = 0; i < num_target_bit; ++i)
    {
        if (gate_a.target_list()[i] != gate_b.target_list()[i])
        {
            equal_bits = false;
        }
    }

    return equal_bits;
}

void Circuit::remove_identities()
{
    std::unordered_map<std::string, std::string> identity_map{
            {"ccz", "ccz"},
            {"toffoli", "toffoli"},
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
    std::size_t is_identity = 0;
    std::size_t same_bit_counter = 0;
    std::unordered_map<std::string, std::list<Gate>::iterator> gate_map;
    std::string compared_bit;

    /*
     * search the identity matrix until the gate is no longer removed
     */
    while (update)
    {
        update = false;
        for (auto it = gate_list_.begin(); it != gate_list_.end(); ++it)
        {
            // check the same gate as the gate in the gate map
            for (const std::string& bit : it->control_list())
            {
                same_bit_counter = gate_map.count(bit);
                if (same_bit_counter)
                {
                    is_identity += same_bit_counter;
                    compared_bit = bit;
                }
            }
            for (const std::string& bit : it->target_list())
            {
                same_bit_counter = gate_map.count(bit);
                if (same_bit_counter)
                {
                    is_identity += same_bit_counter;
                    compared_bit = bit;
                }
            }

            // remove identity gate
            if (is_identity > 0 && equal_gate(*it, *gate_map[compared_bit]))
            {
                gate_list_.erase(it);
                gate_list_.erase(gate_map.at(compared_bit));
                num_gate_ -= 2;

                gate_map.clear();
                is_identity = 0;
                update = true;
                break;
            }

            // register the gate to the gate_map
            for (const std::string& bit : it->control_list())
            {
                gate_map[bit] = it;
            }
            for (const std::string& bit : it->target_list())
            {
                gate_map[bit] = it;
            }
        }
    }
}

void Circuit::decompose_ccz()
{
    for (auto it = gate_list_.begin(); it != gate_list_.end(); it++)
    {
        if (it->type() == "ccz")
        {
            const std::string control_a = it->control_list().front();
            const std::string control_b = it->control_list().back();
            const std::string target = it->target_list().front();

            // remove czz
            it = gate_list_.erase(it);
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
            std::list<Gate> czz = {Gate("cnot", control_b, target),
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

            num_gate_ += static_cast<int>(czz.size());
            gate_list_.splice(it, std::move(czz));
            it--;
        }
    }
}

}