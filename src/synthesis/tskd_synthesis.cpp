#include <iomanip>

#include "tskd_synthesis.hpp"

namespace tskd {

void TskdSynthesis::init(const Character& chr)
{
    global_phase_ = 0;
    index_list_.resize(2);
    carry_index_list_.resize(2);
    remaining_.resize(2);

    /*
     * initialize some stuff
     */
    mask_ = util::xor_func(chr.num_data_qubit() + chr.num_hadamard() + 1, 0);
    mask_.set(chr.num_data_qubit() + chr.num_hadamard());
    for (int i = 0, j = 0; i < chr.num_qubit(); i++)
    {
        circuit_.add_qubit(chr.qubit_names()[i]);
        circuit_.set_ancilla(chr.qubit_names()[i], chr.ancilla_list()[i]);
        wires_.emplace_back(chr.num_data_qubit() + chr.num_hadamard()+ 1, 0);
        if (!chr.ancilla_list()[i])
        {
            wires_[i].set(j);
            mask_.set(j);
            j++;
        }
    }

    /*
     * initialize the remaining list
     */
    int index = 0;
    for (auto&& phase_exponent : chr.phase_exponents())
    {
        if (phase_exponent.second == util::xor_func(chr.num_data_qubit() + chr.num_hadamard() + 1, 0))
        {
            global_phase_ = phase_exponent.first;
        }
        else if (phase_exponent.first % 2 == 1)
        {
            remaining_[0].push_back(index);
        }
        else if (phase_exponent.first != 0)
        {
            remaining_[1].push_back(index);
        }
        index++;
    }

    /**
     * sort index of phase exponents
     */
    for (int i = 0; i < 2; ++i)
    {
        remaining_[i].sort([&chr](int lhs, int rhs) -> bool
                           {
                               if (chr.phase_exponents()[lhs].second.count() == chr.phase_exponents()[rhs].second.count())
                               {
                                   return chr.phase_exponents()[lhs].second < chr.phase_exponents()[rhs].second;
                               }
                               else
                               {
                                   return chr.phase_exponents()[lhs].second.count() < chr.phase_exponents()[rhs].second.count();
                               }
                           });
    }
}

void TskdSynthesis::determine_apply_phase_set(Character::Hadamard& hadamard)
{
    std::vector<std::list<int>> tmp_index_list(2);
    std::vector<std::list<int>> tmp_carry_index_list(2);
    for (int i = 0; i < 2; ++i)
    {
        for (auto it = remaining_[i].begin(); it != remaining_[i].end();)
        {
            util::xor_func tmp = (~mask_) & (chr_.phase_exponents()[*it].second);
            if (tmp.none())
            {
                auto ti = hadamard.in_.find(*it);
                if (ti != hadamard.in_.end())
                {
                    tmp_index_list[i].push_back(*it);
                }
                else
                {
                    tmp_carry_index_list[i].push_back(*it);
                }
                it = remaining_[i].erase(it);
            }
            else
            {
                it++;
            }
        }
    }

    index_list_[0] = tmp_index_list[0];
    index_list_[1] = tmp_index_list[1];
    carry_index_list_[0] = tmp_carry_index_list[0];
    carry_index_list_[1] = tmp_carry_index_list[1];
}

void TskdSynthesis::construct_subcircuit(const Character::Hadamard& hadamard)
{
    circuit_.add_gate_list(builder_.build(index_list_[0], carry_index_list_[0], wires_, wires_));
    circuit_.add_gate_list(builder_.build(index_list_[1], carry_index_list_[1], wires_, hadamard.input_wires_parity_));
    remaining_[0].splice(remaining_[0].begin(), carry_index_list_[0]);
    remaining_[1].splice(remaining_[1].begin(), carry_index_list_[1]);
    for (int i = 0; i < chr_.num_qubit(); i++)
    {
        wires_[i] = hadamard.input_wires_parity_[i];
    }
}

void TskdSynthesis::apply_hadamard(const Character::Hadamard& hadamard)
{
    circuit_.add_gate("H", chr_.qubit_names()[hadamard.target_]);
    wires_[hadamard.target_].reset();
    wires_[hadamard.target_].set(hadamard.previous_qubit_index_);
    mask_.set(hadamard.previous_qubit_index_);
}

void TskdSynthesis::construct_final_subcircuit()
{
    std::list<int> none_list;
    std::vector<std::list<int>> final_index_list(2);
    for (int i = 0; i < 2; ++i)
    {
        for (auto it = remaining_[i].begin(); it != remaining_[i].end(); it++)
        {
            final_index_list[i].push_back(*it);
        }
    }

    circuit_.add_gate_list(builder_.build(final_index_list[0], none_list, wires_, wires_));
    circuit_.add_gate_list(builder_.build(final_index_list[1], none_list, wires_, chr_.outputs()));

    /*
     * Add the global phase
     */
    circuit_.add_gate_list(builder_.build_global_phase(chr_.num_qubit(), global_phase_, chr_.qubit_names()));
}

Circuit TskdSynthesis::execute()
{
    std::cout << "t-scheduling running..." << std::endl;

    int dimension = chr_.num_data_qubit();

    for (auto&& hadamard : chr_.hadamards())
    {
        /*
         * determine apply (carry) index list
         */
        determine_apply_phase_set(hadamard);

        /**
         * Construct sub-circuit
         */
        construct_subcircuit(hadamard);

        /*
         * Apply Hadamard gate
        */
        apply_hadamard(hadamard);

        /*
         * Check for increases in dimension
         */
        dimension = builder_.check_dimension(chr_, wires_, dimension);
    }

    /*
     * Construct the final {CNOT, T} subcircuit
     */
    construct_final_subcircuit();

    return circuit_;
}

}
