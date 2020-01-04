#include "tpar_synthesis.hpp"

#include "../tpar/partition.hpp"
#include "../tpar/matroid.hpp"

#include "../character/character.hpp"

namespace tskd {

void TparSynthesis::init(const Character& chr)
{
    global_phase_ = 0;
    floats_.resize(2);
    frozen_.resize(2);
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
        wires_.emplace_back(chr.num_data_qubit() + chr.num_hadamard() + 1, 0);
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
    for (auto&& phase_exponent : chr_.phase_exponents())
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

void TparSynthesis::create_partition()
{
    for (auto it = remaining_.begin(); it != remaining_.end();)
    {
        util::xor_func tmp = (~mask_) & (chr_.phase_exponents()[*it].second);
        if (tmp.none())
        {
            tpar::add_to_partition(floats_, *it, chr_.phase_exponents(), oracle_);
            it = remaining_.erase(it);
        }
        else
        {
            it++;
        }
    }
}

void TparSynthesis::determine_apply_partition(Character::Hadamard& hadamard)
{
    // TODO: freeze_partisionsの第二引数のconst化
    frozen_ = tpar::freeze_partitions(floats_, hadamard.in_);
}

void TparSynthesis::construct_subcircuit(Character::Hadamard& hadamard)
{
    std::vector<util::xor_func> original_hadamard_outputs = hadamard.input_wires_parity_;
    circuit_.add_gate_list(builder_.build(frozen_, wires_, hadamard.input_wires_parity_, bit_map_));
    update_bit_map(original_hadamard_outputs, hadamard.input_wires_parity_);

    for (int i = 0; i < chr_.num_qubit(); i++)
    {
        wires_[i] = hadamard.input_wires_parity_[i];
    }
}

void TparSynthesis::apply_hadamard(const Character::Hadamard& hadamard)
{
    const int hadamard_target = bit_map_[hadamard.target_];
    circuit_.add_gate("H", chr_.qubit_names()[hadamard_target]);
    wires_[hadamard_target].reset();
    wires_[hadamard_target].set(hadamard.previous_qubit_index_);
    mask_.set(hadamard.previous_qubit_index_);
}

int TparSynthesis::check_dimension(int current_dimension)
{
    int new_dimension = 0;
    const int updated_dimension = util::compute_rank(chr_.num_qubit(), chr_.num_data_qubit() + chr_.num_hadamard(), wires_);
    if (updated_dimension > current_dimension)
    {
        new_dimension = updated_dimension;
        oracle_.set_dim(new_dimension);
        tpar::repartition(floats_, chr_.phase_exponents(), oracle_);
    }

    return new_dimension;
}

void TparSynthesis::construct_final_subcircuit()
{
    std::vector<util::xor_func> outputs = chr_.outputs();
    circuit_.add_gate_list(builder_.build(floats_, wires_, outputs, bit_map_));

    /*
     * Add the global phase
     */
    circuit_.add_gate_list(builder_.build_global_phase(chr_.num_qubit(), global_phase_, chr_.qubit_names()));
}


Circuit TparSynthesis::execute()
{
    std::cout << "t-par running..." << std::endl;

    int dimension = chr_.num_data_qubit();

    /*
     * create an initial partition
     */
    create_partition();

    /**
     * Synthesize the circuit by applying H-gate to greed
     * 1. freeze partitions that are not disjoint from the hadamard input
     * 2. construct CNOT+T circuit
     * 3. apply the hadamard gate
     * 4. add new functions to the partition
     */
    for (auto&& hadamard : chr_.hadamards())
    {
        /*
         * determine frozen partitions
         */
        determine_apply_partition(hadamard);

        /*
         * Construct {CNOT, T} subcircuit for the frozen partitions
         */
        construct_subcircuit(hadamard);

        /*
         * Apply Hadamard gate
         */
        apply_hadamard(hadamard);

        /*
         * Check for increases in dimension
         */
        dimension = check_dimension(dimension);

        /*
         * Add new functions to the partition
         */
        create_partition();
    }

    /*
     * Construct the final {CNOT, T} subcircuit
     */
    construct_final_subcircuit();

    return circuit_;
}

}