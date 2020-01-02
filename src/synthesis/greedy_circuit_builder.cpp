#include <random>
#include <chrono>

#include "greedy_circuit_builder.hpp"

#include "../tpar/partition.hpp"

#include "../matrix/matrix_reconstructor.hpp"

namespace tskd {

bool GreedyCircuitBuilder::init(const std::vector<util::xor_func>& in,
                                const std::vector<util::xor_func>& out)
{
    bool is_io_different = true;

    bits_ = std::vector<util::xor_func>(qubit_num_);
    preparation_ = std::vector<util::xor_func>(qubit_num_);
    restoration_ = std::vector<util::xor_func>(qubit_num_);

    for (int i = 0; i < qubit_num_; i++)
    {
        is_io_different &= (in[i] == out[i]);
        for (int i = 0; i < qubit_num_; i++)
        {
            is_io_different &= (in[i] == out[i]);
            preparation_[i] = util::xor_func(qubit_num_ + 1, 0);
            restoration_[i] = util::xor_func(qubit_num_ + 1, 0);
            preparation_[i].set(i);
            restoration_[i].set(i);
        }
    }

    return is_io_different;
}

int GreedyCircuitBuilder::compute_time_step(const std::list<Gate>& gate_list)
{
    return static_cast<int>(gate_list.size()) * 2;
}

void GreedyCircuitBuilder::apply_phase_gates(std::list<Gate>& gate_list,
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

void GreedyCircuitBuilder::unprepare(const std::vector <util::xor_func>& restoration)
{
    preparation_ = std::move(restoration);
    restoration_ = std::vector<util::xor_func>(qubit_num_);
    // re-initialize
    for (int i = 0; i < qubit_num_; i++)
    {
        restoration_[i] = util::xor_func(qubit_num_ + 1, 0);
        restoration_[i].set(i);
    }
}

void GreedyCircuitBuilder::prepare_last_part(std::list<Gate>& gate_list,
                                             const std::vector<util::xor_func>& in,
                                             const std::vector<util::xor_func>& out,
                                             MatrixReconstructor& sa)
{
    for (int i = 0; i < qubit_num_; i++)
    {
        bits_[i] = out[i];
    }
    util::to_upper_echelon(qubit_num_, dimension_, bits_, &restoration_, std::vector<std::string>());
    util::fix_basis(qubit_num_, dimension_, qubit_num_, in, bits_, &restoration_, std::vector<std::string>());
    util::compose(qubit_num_, preparation_, restoration_);
//    std::cout << "-- preparation" << std::endl;
//    for (auto&& e : preparation_)
//    {
//        std::cout << e << std::endl;
//    }

    gate_list.splice(gate_list.end(), (*decomposer_)(layout_, qubit_num_, 0, preparation_, qubit_names_));
}

int GreedyCircuitBuilder::check_dimension(const Character& chr,
                                          std::vector <util::xor_func>& wires,
                                          int current_dimension)
{
    int new_dimension = 0;
    const int updated_dimension = util::compute_rank(chr.num_qubit(), chr.num_data_qubit() + chr.num_hadamard(), wires);
    if (updated_dimension > current_dimension)
    {
        new_dimension = updated_dimension;
        oracle_.set_dim(new_dimension);
    }

    return new_dimension;
}

std::list<Gate> GreedyCircuitBuilder::build(std::list<int>& index_list,
                                            std::list<int>& carry_index_list,
                                            std::vector<util::xor_func>& in,
                                            const std::vector<util::xor_func>& out)
{
    std::list<Gate> ret;

    std::vector<std::pair<int, int>> phase_target_list;

    if (init(in, out) && index_list.empty())
    {
        return ret;
    }

    /*
     * Reduce in to echelon form to decide on a basis
     */
    util::to_upper_echelon(qubit_num_, dimension_, in, &preparation_, std::vector<std::string>());
    init_prep_ = preparation_;

    MatrixReconstructor sa(in, dimension_, qubit_num_);

    while (!index_list.empty())
    {
        std::vector<util::xor_func> result_restoration(qubit_num_);
        for (int i = 0; i < qubit_num_; i++)
        {
            result_restoration[i] = util::xor_func(qubit_num_ + 1, 0);
            result_restoration[i].set(i);
        }
        std::list<Gate> result_gate_list;
        std::set<int> result_sub_part;
        std::list<int> delete_index_list;
        std::unordered_map<int, int> result_target_phase_map;

        /**
         * first build sub-circuit
         */
        for (auto it = index_list.begin(); it != index_list.end();)
        {
            std::list<Gate> tmp_gate_list;
            std::set<int> tmp_sub_part = result_sub_part;
            tmp_sub_part.insert(*it);

            std::vector<util::xor_func> tmp_bits(qubit_num_);
            std::vector<util::xor_func> tmp_preparation = preparation_;
            std::vector<util::xor_func> tmp_restoration = restoration_;

            std::unordered_map<int, int> tmp_target_phase_map;

            if (oracle_(phase_exponent_, tmp_sub_part))
            {
                /**
                 * create bits matrix
                 */
                std::set<int>::iterator ti;
                int counter = 0;
                for (ti = tmp_sub_part.begin(), counter = 0; counter < qubit_num_; counter++)
                {
                    if (counter < static_cast<int>(tmp_sub_part.size()))
                    {
                        tmp_bits[counter] = phase_exponent_[*ti].second;
                        tmp_target_phase_map.emplace(counter, *ti);
                        ti++;
                    }
                    else
                    {
                        tmp_bits[counter] = util::xor_func(dimension_ + 1, 0);
                    }
                }


                /**
                 * prepare preparation matrix
                 */
                const int num_partition = static_cast<int>(tmp_sub_part.size());
                util::to_upper_echelon(num_partition, dimension_, tmp_bits, &tmp_restoration, std::vector<std::string>());
                util::fix_basis(qubit_num_, dimension_, num_partition, in, tmp_bits, &tmp_restoration,
                               std::vector<std::string>());
                /**
                 * change row order in the matrix
                 */
                if (option_.change_row_order())
                {
                    tmp_restoration = sa.execute(tmp_preparation, tmp_restoration, tmp_target_phase_map);
                }

                util::compose(qubit_num_, tmp_preparation, tmp_restoration);
//                std::cout << "-- preparation" << std::endl;
//                for (auto&& e : tmp_preparation)
//                {
//                    std::cout << e << std::endl;
//                }

                /**
                 * create gate list
                 */
                tmp_gate_list.splice(tmp_gate_list.end(),
                                     (*decomposer_)(layout_, qubit_num_, 0, tmp_preparation, qubit_names_));

                /**
                 * check time step
                 */
                int upper_bound_time = option_.distillation_step() * static_cast<int>(tmp_sub_part.size());
                const int buffer_upper_bound = std::max(1, option_.num_buffer()) * option_.distillation_step();
                upper_bound_time = std::min(upper_bound_time, buffer_upper_bound);
                if (tmp_sub_part.size() == 1 || compute_time_step(tmp_gate_list) <= upper_bound_time)
                {
                    result_restoration = tmp_restoration;
                    result_sub_part = tmp_sub_part;
                    result_gate_list = tmp_gate_list;
                    result_target_phase_map = tmp_target_phase_map;
                    it = index_list.erase(it);
                }
                else
                {
                    it++;
                }
            }
            else
            {
                it++;
            }
        }

        /**
         * second build sub-circuit
         */
        for (auto it = carry_index_list.begin(); it != carry_index_list.end();)
        {
            std::list<Gate> tmp_gate_list;
            std::set<int> tmp_sub_part = result_sub_part;
            tmp_sub_part.insert(*it);

            std::vector<util::xor_func> tmp_bits(qubit_num_);
            std::vector<util::xor_func> tmp_preparation = preparation_;
            std::vector<util::xor_func> tmp_restoration = restoration_;

            std::unordered_map<int, int> tmp_target_phase_map;

            if (oracle_(phase_exponent_, tmp_sub_part))
            {
                /**
                 * create bits matrix
                 */
                std::set<int>::iterator ti;
                int counter = 0;
                for (ti = tmp_sub_part.begin(), counter = 0; counter < qubit_num_; counter++)
                {
                    if (counter < static_cast<int>(tmp_sub_part.size()))
                    {
                        tmp_bits[counter] = phase_exponent_[*ti].second;
                        tmp_target_phase_map.emplace(counter, *ti);
                        ti++;
                    }
                    else
                    {
                        tmp_bits[counter] = util::xor_func(dimension_ + 1, 0);
                    }
                }

                /**
                 * prepare preparation matrix
                 */
                const int num_partition = static_cast<int>(tmp_sub_part.size());
                util::to_upper_echelon(num_partition, dimension_, tmp_bits, &tmp_restoration, std::vector<std::string>());
                util::fix_basis(qubit_num_, dimension_, num_partition, in, tmp_bits, &tmp_restoration,
                               std::vector<std::string>());

                /**
                 * change row order in the matrix
                 */
                if (option_.change_row_order())
                {
                    tmp_restoration = sa.execute(tmp_preparation, tmp_restoration, tmp_target_phase_map);
                }

                util::compose(qubit_num_, tmp_preparation, tmp_restoration);

                /**
                 * create gate list
                 */
                tmp_gate_list.splice(tmp_gate_list.end(),
                                     (*decomposer_)(layout_, qubit_num_, 0, tmp_preparation, qubit_names_));

                /**
                 * check time step
                 */
                int upper_bound_time = option_.distillation_step() * static_cast<int>(tmp_sub_part.size());
                const int buffer_upper_bound = std::max(1, option_.num_buffer()) * option_.distillation_step();
                upper_bound_time = std::min(upper_bound_time, buffer_upper_bound);
                if (tmp_sub_part.size() == 1 || compute_time_step(tmp_gate_list) <= upper_bound_time)
                {
                    result_restoration = tmp_restoration;
                    result_sub_part = tmp_sub_part;
                    result_gate_list = tmp_gate_list;
                    result_target_phase_map = tmp_target_phase_map;
                    it = carry_index_list.erase(it);
                }
                else
                {
                    it++;
                }
            }
            else
            {
                it++;
            }
        }

        ret.splice(ret.end(), result_gate_list);

        /*
         * Apply the phase gates
         */
        apply_phase_gates(ret, result_target_phase_map);

        /*
         * Unprepare the bits
         */
        unprepare(result_restoration);
    }

    /*
     * Reduce out to the basis of in
     */
    prepare_last_part(ret, in, out, sa);

    return ret;
}

std::list<Gate> GreedyCircuitBuilder::build_global_phase(int qubit_num,
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
