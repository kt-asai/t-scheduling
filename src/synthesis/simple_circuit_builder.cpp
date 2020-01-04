#include <random>
#include <chrono>

#include "simple_circuit_builder.hpp"

#include "../matrix/matrix_reconstructor.hpp"

#include "../decomposer/gaussian_decomposer.hpp"

namespace tskd {

bool SimpleCircuitBuilder::init(const std::vector<util::xor_func>& in,
                                const std::vector<util::xor_func>& out)
{
    bool is_io_different = true;

    bits_ = std::vector<util::xor_func>(qubit_num_);
    preparation_ = std::vector<util::xor_func>(qubit_num_);
    restoration_ = std::vector<util::xor_func>(qubit_num_);
    identity_ = std::vector<util::xor_func>(qubit_num_);

    for (int i = 0; i < qubit_num_; i++)
    {
        is_io_different &= (in[i] == out[i]);
        preparation_[i] = util::xor_func(qubit_num_ + 1, 0);
        restoration_[i] = util::xor_func(qubit_num_ + 1, 0);
        identity_[i] = util::xor_func(qubit_num_ + 1, 0);
        preparation_[i].set(i);
        restoration_[i].set(i);
        identity_[i].set(i);
    }

    return is_io_different;
}

void SimpleCircuitBuilder::init_bits(const std::set<int>& phase_exponent_index_set,
                                     std::unordered_map<int, int>& target_phase_map,
                                     std::vector<util::xor_func>& in)
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

void SimpleCircuitBuilder::prepare(std::list<Gate>& gate_list,
                                   const std::vector<util::xor_func>& in,
                                   const int num_partition,
                                   std::unordered_map<int, int>& target_phase_map,
                                   MatrixReconstructor& sa)
{
    util::to_upper_echelon(num_partition, dimension_, bits_, &restoration_, std::vector<std::string>());
    util::fix_basis(qubit_num_, dimension_, num_partition, in, bits_, &restoration_, std::vector<std::string>());

    /*
     * Re-construct binary matrix
     */
    if (option_.change_row_order())
    {
        restoration_ = sa.execute(preparation_, restoration_, target_phase_map);
    }

    /*
     * set preparation before decompose matrix
     */
    std::vector<int> func_map;
    std::vector<util::xor_func> before_prep(identity_);
    util::compose(qubit_num_, before_prep, restoration_);

    util::compose(qubit_num_, preparation_, restoration_);

    gate_list.splice(gate_list.end(), (*decomposer_).execute(preparation_, func_map));
//    gate_list.splice(gate_list.end(), (*decomposer_)(layout_, qubit_num_, 0, preparation_, qubit_names_));

    /*
     * procedure after remove swap gate
     */
}

void SimpleCircuitBuilder::apply_phase_gates(std::list<Gate>& gate_list,
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

void SimpleCircuitBuilder::unprepare()
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

void SimpleCircuitBuilder::prepare_last_part(std::list<Gate>& gate_list,
                                             const std::vector<util::xor_func>& in,
                                             std::vector<util::xor_func>& out,
                                             MatrixReconstructor& sa)
{
    for (int i = 0; i < qubit_num_; i++)
    {
        bits_[i] = out[i];
    }

    std::unordered_map<int, int> bit_correspond_map;
    for (auto i = 0; i < out.size(); i++)
    {
        bit_correspond_map.emplace(i, i);
    }

    util::to_upper_echelon(qubit_num_, dimension_, bits_, &restoration_, std::vector<std::string>());
    util::fix_basis(qubit_num_, dimension_, qubit_num_, in, bits_, &restoration_, std::vector<std::string>());

    /*
     * Re-construct binary matrix
     */
    if (option_.change_row_order())
    {
        restoration_ = sa.execute(preparation_, restoration_, bit_correspond_map);
    }

    /*
     * set preparation before decompose matrix
     */
    std::vector<int> func_map;
    std::vector<util::xor_func> before_prep(identity_);
    util::compose(qubit_num_, before_prep, restoration_);

    util::compose(qubit_num_, preparation_, restoration_);


    // move bit place of out
    std::vector<util::xor_func> tmp_out(out.size());
    for (auto i = 0; i < out.size(); i++)
    {
        tmp_out[i] = out[bit_correspond_map[i]];
    }
    out = tmp_out;

    gate_list.splice(gate_list.end(), (*decomposer_).execute(preparation_, func_map));
//    gate_list.splice(gate_list.end(), (*decomposer_)(layout_, qubit_num_, 0, preparation_, qubit_names_));
}

std::list<Gate> SimpleCircuitBuilder::build(const tpar::partitioning& partition,
                                            std::vector<util::xor_func>& in,
                                            std::vector<util::xor_func>& out)
{
    std::list<Gate> ret;

    std::unordered_map<int, int> target_phase_map;

    if (init(in, out) && partition.empty())
    {
        return ret;
    }

    /*
     * Reduce in to echelon form to decide on a basis
     */
    util::to_upper_echelon(qubit_num_, dimension_, in, &preparation_, std::vector<std::string>());

    MatrixReconstructor sa(in, dimension_, qubit_num_);

    /*
     * For each partition... Compute *it, apply T gates, uncompute
     */
    for (auto&& it : partition)
    {
        /*
         * Initialize binary matrix
         */
        init_bits(it, target_phase_map, in);

        /*
         * Prepare the bits
         */
        prepare(ret, in, static_cast<int>(it.size()), target_phase_map, sa);

        /*
         * Apply the phase gates
         */
        apply_phase_gates(ret, target_phase_map);

        /*
         * Unprepare the bits
         */
        unprepare();
    }

    /*
     * Reduce out to the basis of in
     */
    prepare_last_part(ret, in, out, sa);

    return ret;
}

std::list<Gate> SimpleCircuitBuilder::build_global_phase(int qubit_num,
                                                       int phase,
                                                       const std::vector<std::string>& qubit_names)
{
    std::list<Gate> acc;
    int qubit = 0;

    if (phase % 2 == 1)
    {
        acc.splice(acc.end(), util::compose_om(qubit, qubit_names));
        qubit = (qubit + 1) % qubit_num;
    }
    for (int i = phase / 2; i > 0; i--)
    {
        acc.splice(acc.end(), util::compose_imaginary_unit(qubit, qubit_names));
        qubit = (qubit + 1) % qubit_num;
    }

    return acc;
}

}