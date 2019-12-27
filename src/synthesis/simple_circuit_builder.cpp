#include <random>
#include <chrono>

#include "simple_circuit_builder.hpp"

#include "../matrix/matrix_reconstructor.hpp"

#include "../decomposer/gaussian_decomposer.hpp"

namespace tskd {

bool SimpleCircuitBuilder::Init(const std::vector<util::xor_func>& in,
                                const std::vector<util::xor_func>& out)
{
    bool is_io_different = true;

    bits_ = std::vector<util::xor_func>(qubit_num_);
    preparation_ = std::vector<util::xor_func>(qubit_num_);
    restoration_ = std::vector<util::xor_func>(qubit_num_);

    for (int i = 0; i < qubit_num_; i++)
    {
        is_io_different &= (in[i] == out[i]);
        preparation_[i] = util::xor_func(qubit_num_ + 1, 0);
        restoration_[i] = util::xor_func(qubit_num_ + 1, 0);
        preparation_[i].set(i);
        restoration_[i].set(i);
    }

    return is_io_different;
}

void SimpleCircuitBuilder::InitBits(const std::set<int>& phase_exponent_index_set,
                                    std::unordered_map<int, int>& target_phase_map)
{
    std::set<int>::iterator ti;
    int counter = 0;
    target_phase_map.clear();
    for (ti = phase_exponent_index_set.begin(), counter = 0; counter < qubit_num_; counter++)
    {
        if (counter < static_cast<int>(phase_exponent_index_set.size()))
        {
            bits_[counter] = phase_exponent_[*ti].second;
            target_phase_map.emplace(counter, *ti);
            ti++;
        }
        else
        {
            bits_[counter] = util::xor_func(dimension_ + 1, 0);
        }
    }
}

void SimpleCircuitBuilder::Prepare(std::list<Gate>& gate_list,
                                   const std::vector<util::xor_func>& in,
                                   const int num_partition,
                                   std::unordered_map<int, int>& target_phase_map)
{
    util::ToUpperEchelon(num_partition, dimension_, bits_, &restoration_, std::vector<std::string>());
    util::FixBasis(qubit_num_, dimension_, num_partition, in, bits_, &restoration_, std::vector<std::string>());
    util::Compose(qubit_num_, preparation_, restoration_);

    gate_list.splice(gate_list.end(), (*decomposer_)(qubit_num_, 0, preparation_, qubit_names_));
}

void SimpleCircuitBuilder::ApplyPhaseGates(std::list<Gate>& gate_list,
                                           const std::unordered_map<int, int>& target_phase_map)
{
    for (auto&& tp : target_phase_map)
    {
        const int target = tp.first;
        const int phase_exponent_index = tp.second;

        if (phase_exponent_[phase_exponent_index].first <= 4)
        {
            if (phase_exponent_[phase_exponent_index].first / 4 == 1)
            {
                gate_list.emplace_back("Z", qubit_names_[target]);
            }
            if (phase_exponent_[phase_exponent_index].first / 2 == 1)
            {
                gate_list.emplace_back("P", qubit_names_[target]);
            }
            if (phase_exponent_[phase_exponent_index].first % 2 == 1)
            {
                gate_list.emplace_back("T", qubit_names_[target]);
            }
        }
        else
        {
            if (phase_exponent_[phase_exponent_index].first == 5 || phase_exponent_[phase_exponent_index].first == 6)
            {
                gate_list.emplace_back("P*", qubit_names_[target]);
            }
            if (phase_exponent_[phase_exponent_index].first % 2 == 1)
            {
                gate_list.emplace_back("T*", qubit_names_[target]);
            }
        }
    }
}

void SimpleCircuitBuilder::UnPrepare()
{
    preparation_ = std::move(restoration_);
    restoration_ = std::vector<util::xor_func>(qubit_num_);
    // re-initialize
    for (int i = 0; i < qubit_num_; i++)
    {
        restoration_[i] = util::xor_func(qubit_num_ + 1, 0);
        restoration_[i].set(i);
    }
}

void SimpleCircuitBuilder::PrepareLastPart(std::list<Gate>& gate_list,
                                           const std::vector<util::xor_func>& in,
                                           const std::vector<util::xor_func>& out,
                                           MatrixReconstructor& sa)
{
    for (int i = 0; i < qubit_num_; i++)
    {
        bits_[i] = out[i];
    }
    std::unordered_map<int, int> none;
    if (option_.change_row_order())
    {
        bits_ = sa.Execute(static_cast<int>(bits_.size()), bits_, none);
    }
    util::ToUpperEchelon(qubit_num_, dimension_, bits_, &restoration_, std::vector<std::string>());
    util::FixBasis(qubit_num_, dimension_, qubit_num_, in, bits_, &restoration_, std::vector<std::string>());
    util::Compose(qubit_num_, preparation_, restoration_);

    gate_list.splice(gate_list.end(), (*decomposer_)(qubit_num_, 0, preparation_, qubit_names_));
}

std::list<Gate> SimpleCircuitBuilder::Build(const tpar::partitioning& partition,
                                            std::vector<util::xor_func>& in,
                                            const std::vector<util::xor_func>& out)
{
    std::list<Gate> ret;

    std::unordered_map<int, int> target_phase_map;

    if (Init(in, out) && partition.empty())
    {
        return ret;
    }

    /*
     * Reduce in to echelon form to decide on a basis
     */
    util::ToUpperEchelon(qubit_num_, dimension_, in, &preparation_, std::vector<std::string>());

    MatrixReconstructor sa(in, dimension_, qubit_num_);

    /*
     * For each partition... Compute *it, apply T gates, uncompute
     */
    for (auto&& it : partition)
    {
        /*
         * Initialize binary matrix
         */
        InitBits(it, target_phase_map);

        /*
         * Re-construct binary matrix
         */
        if (option_.change_row_order())
        {
            bits_ = sa.Execute(static_cast<int>(it.size()), bits_, target_phase_map);
        }

        /*
         * Prepare the bits
         */
        Prepare(ret, in, static_cast<int>(it.size()), target_phase_map);

        /*
         * Apply the phase gates
         */
        ApplyPhaseGates(ret, target_phase_map);

        /*
         * Unprepare the bits
         */
        UnPrepare();
    }

    /*
     * Reduce out to the basis of in
     */
    PrepareLastPart(ret, in, out, sa);

    return ret;
}

std::list<Gate> SimpleCircuitBuilder::BuildGlobalPhase(int qubit_num,
                                                       int phase,
                                                       const std::vector<std::string>& qubit_names)
{
    std::list<Gate> acc;
    int qubit = 0;

    if (phase % 2 == 1)
    {
        acc.splice(acc.end(), util::ComposeOM(qubit, qubit_names));
        qubit = (qubit + 1) % qubit_num;
    }
    for (int i = phase / 2; i > 0; i--)
    {
        acc.splice(acc.end(), util::ComposeImaginaryUnit(qubit, qubit_names));
        qubit = (qubit + 1) % qubit_num;
    }

    return acc;
}

}