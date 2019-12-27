#include "tpar_synthesis.hpp"

#include "../tpar/partition.hpp"
#include "../tpar/matroid.hpp"

#include "../character/character.hpp"

namespace tskd {

void TparSynthesis::Init(const Character& chr)
{
    global_phase_ = 0;
    floats_.resize(2);
    frozen_.resize(2);
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
}

void TparSynthesis::CreatePartition()
{
    for (int j = 0; j < 2; j++)
    {
        for (auto it = remaining_[j].begin(); it != remaining_[j].end();)
        {
            util::xor_func tmp = (~mask_) & (chr_.phase_exponents()[*it].second);
            if (tmp.none())
            {
                tpar::add_to_partition(floats_[j], *it, chr_.phase_exponents(), oracle_);
                it = remaining_[j].erase(it);
            }
            else
            {
                it++;
            }
        }
    }
}

void TparSynthesis::DetermineApplyPartition(Character::Hadamard& hadamard)
{
    for (int i = 0; i < 2; i++)
    {
        // TODO: freeze_partisionsの第二引数のconst化
        frozen_[i] = tpar::freeze_partitions(floats_[i], hadamard.in_);
    }
}

void TparSynthesis::ConstructSubCircuit(const Character::Hadamard& hadamard)
{
    circuit_.add_gate_list(builder_.Build(frozen_[0], wires_, wires_));
    circuit_.add_gate_list(builder_.Build(frozen_[1], wires_, hadamard.input_wires_parity_));
    for (int i = 0; i < chr_.num_qubit(); i++)
    {
        wires_[i] = hadamard.input_wires_parity_[i];
    }
}

void TparSynthesis::ApplyHadamard(const Character::Hadamard& hadamard)
{
    circuit_.add_gate("H", chr_.qubit_names()[hadamard.target_]);
    wires_[hadamard.target_].reset();
    wires_[hadamard.target_].set(hadamard.previous_qubit_index_);
    mask_.set(hadamard.previous_qubit_index_);
}

int TparSynthesis::CheckDimension(int current_dimension)
{
    int new_dimension = 0;
    const int updated_dimension = util::ComputeRank(chr_.num_qubit(), chr_.num_data_qubit() + chr_.num_hadamard(), wires_);
    if (updated_dimension > current_dimension)
    {
        new_dimension = updated_dimension;
        oracle_.set_dim(new_dimension);
        tpar::repartition(floats_[0], chr_.phase_exponents(), oracle_);
        tpar::repartition(floats_[1], chr_.phase_exponents(), oracle_);
    }

    return new_dimension;
}

void TparSynthesis::ConstructFinalSubCircuit()
{
    circuit_.add_gate_list(builder_.Build(floats_[0], wires_, wires_));
    circuit_.add_gate_list(builder_.Build(floats_[1], wires_, chr_.outputs()));

    /*
     * Add the global phase
     */
    circuit_.add_gate_list(builder_.BuildGlobalPhase(chr_.num_qubit(), global_phase_, chr_.qubit_names()));
}


Circuit TparSynthesis::Execute()
{
    std::cout << "t-par running..." << std::endl;

    int dimension = chr_.num_data_qubit();

    /*
     * create an initial partition
     */
    CreatePartition();

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
        DetermineApplyPartition(hadamard);

        /*
         * Construct {CNOT, T} subcircuit for the frozen partitions
         */
        ConstructSubCircuit(hadamard);

        /*
         * Apply Hadamard gate
         */
        ApplyHadamard(hadamard);

        /*
         * Check for increases in dimension
         */
        dimension = CheckDimension(dimension);

        /*
         * Add new functions to the partition
         */
        CreatePartition();
    }

    /*
     * Construct the final {CNOT, T} subcircuit
     */
    ConstructFinalSubCircuit();

    return circuit_;
}

}