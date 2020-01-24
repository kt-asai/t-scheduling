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
    bit_map_.resize(chr.num_qubit());
    mask_ = util::xor_func(chr.num_data_qubit() + chr.num_hadamard() + 1, 0);
    mask_.set(chr.num_data_qubit() + chr.num_hadamard());
    for (int i = 0, j = 0; i < chr.num_qubit(); i++)
    {
        bit_map_[i] = i;
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
        else
        {
            remaining_.push_back(index);
        }
        index++;
    }
}

void TskdSynthesis::determine_apply_phase_set(Character::Hadamard& hadamard)
{
    /**
     * sort index of phase exponents
     */
    remaining_.sort([this](int lhs, int rhs) -> bool {
                        if (chr_.phase_exponents()[lhs].second.count() == chr_.phase_exponents()[rhs].second.count())
                        {
                            return chr_.phase_exponents()[lhs].second < chr_.phase_exponents()[rhs].second;
                        }
                        else
                        {
                            return chr_.phase_exponents()[lhs].second.count() < chr_.phase_exponents()[rhs].second.count();
                        }});

    std::list<int> tmp_index_list;
    std::list<int> tmp_carry_index_list;
    for (auto it = remaining_.begin(); it != remaining_.end();)
    {
        util::xor_func tmp = (~mask_) & (chr_.phase_exponents()[*it].second);
        if (tmp.none())
        {
            auto ti = hadamard.in_.find(*it);
            if (ti != hadamard.in_.end())
            {
                tmp_index_list.push_back(*it);
            }
            else
            {
                tmp_carry_index_list.push_back(*it);
            }
            it = remaining_.erase(it);
        }
        else
        {
            it++;
        }
    }

    index_list_ = tmp_index_list;
    carry_index_list_ = tmp_carry_index_list;
}

void TskdSynthesis::construct_subcircuit(Character::Hadamard& hadamard)
{
    std::vector<util::xor_func> original_hadamard_outputs = hadamard.input_wires_parity_;
    circuit_.add_gate_list(builder_.build(index_list_, carry_index_list_, wires_, hadamard.input_wires_parity_));
    update_bit_map(original_hadamard_outputs, hadamard.input_wires_parity_);

    remaining_.splice(remaining_.begin(), carry_index_list_);
    for (int i = 0; i < chr_.num_qubit(); i++)
    {
        wires_[i] = hadamard.input_wires_parity_[i];
    }
}

void TskdSynthesis::apply_hadamard(const Character::Hadamard& hadamard)
{
    const int hadamard_target = bit_map_[hadamard.target_];
    circuit_.add_gate("H", chr_.qubit_names()[hadamard_target]);
    wires_[hadamard_target].reset();
    wires_[hadamard_target].set(hadamard.previous_qubit_index_);
    mask_.set(hadamard.previous_qubit_index_);
}

void TskdSynthesis::construct_final_subcircuit()
{
    std::list<int> none_list;
    std::list<int> final_index_list;
    for (auto it = remaining_.begin(); it != remaining_.end(); it++)
    {
        final_index_list.push_back(*it);
    }

    std::vector<util::xor_func> outputs = chr_.outputs();
    circuit_.add_gate_list(builder_.build(final_index_list, none_list, wires_, outputs));

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
