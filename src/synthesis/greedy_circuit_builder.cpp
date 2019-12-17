#include "greedy_circuit_builder.hpp"

#include "../util/util.hpp"

#include "../circuit/gate.hpp"

#include "../tpar/partition.hpp"

namespace tskd {

bool GreedyCircuitBuilder::Init(const std::vector<util::xor_func>& in,
                                const std::vector<util::xor_func>& out)
{
    bool is_io_different = true;

//    bits_ = std::vector<util::xor_func>(qubit_num_);
//    preparation_ = std::vector<util::xor_func>(qubit_num_);
//    restoration_ = std::vector<util::xor_func>(qubit_num_);

    for (int i = 0; i < qubit_num_; i++)
    {
        is_io_different &= (in[i] == out[i]);
//        preparation_[i] = util::xor_func(qubit_num_ + 1, 0);
//        restoration_[i] = util::xor_func(qubit_num_ + 1, 0);
//        preparation_[i].set(i);
//        restoration_[i].set(i);
    }

    return is_io_different;
}

void GreedyCircuitBuilder::ChangeRowOrder(std::vector<std::pair<int, int>>& phase_target_list,
                                          std::vector<util::xor_func>& matrix)
{

}

int GreedyCircuitBuilder::ComputeTimeStep(const std::list<Gate>& gate_list)
{
    return static_cast<int>(gate_list.size());
}

void GreedyCircuitBuilder::ApplyPhaseGates(std::list<Gate>& gate_list,
                                           const std::set<int>& phase_exponent_index_set)
{
    int index = 0;
    for (auto&& phase_exponent_index : phase_exponent_index_set)
    {
        if (phase_exponent_[phase_exponent_index].first <= 4)
        {
            if (phase_exponent_[phase_exponent_index].first / 4 == 1)
            {
                gate_list.emplace_back("Z", qubit_names_[index]);
            }
            if (phase_exponent_[phase_exponent_index].first / 2 == 1)
            {
                gate_list.emplace_back("P", qubit_names_[index]);
            }
            if (phase_exponent_[phase_exponent_index].first % 2 == 1)
            {
                gate_list.emplace_back("T", qubit_names_[index]);
            }
        }
        else
        {
            if (phase_exponent_[phase_exponent_index].first == 5 || phase_exponent_[phase_exponent_index].first == 6)
            {
                gate_list.emplace_back("P*", qubit_names_[index]);
            }
            if (phase_exponent_[phase_exponent_index].first % 2 == 1)
            {
                gate_list.emplace_back("T*", qubit_names_[index]);
            }
        }
        index++;
    }
}

void GreedyCircuitBuilder::ApplyPhaseGates(std::list<Gate>& gate_list,
                                           const std::vector<std::pair<int, int>>& phase_target_list)
{
    for (auto&& p : phase_target_list)
    {
        const int phase_exponent_index = p.first;
        const int target = p.second;

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

int GreedyCircuitBuilder::CheckDimension(const Character& chr,
                                         std::vector <util::xor_func>& wires,
                                         int current_dimension)
{
    int new_dimension = 0;
    const int updated_dimension = util::ComputeRank(chr.num_qubit(), chr.num_data_qubit() + chr.num_hadamard(), wires);
    if (updated_dimension > current_dimension)
    {
        new_dimension = updated_dimension;
        oracle_.set_dim(new_dimension);
    }

    return new_dimension;
}

std::list<Gate> GreedyCircuitBuilder::Build(std::list<int>& index_list,
                                            std::list<int>& carry_index_list,
                                            std::vector<util::xor_func>& in,
                                            const std::vector<util::xor_func>& out)
{
//    std::cout << std::endl;
//    std::cout << "-->> greedy circuit builder" << std::endl;
//    std::cout << "# input" << std::endl;
//    std::cout << "-- index list" << std::endl;
//    for (auto&& e : index_list)
//    {
//        std::cout << e << ":" << phase_exponent_[e].second << std::endl;
//    }
//    std::cout << "-- carry index list" << std::endl;
//    for (auto&& e : carry_index_list)
//    {
//        std::cout << e << ":" << phase_exponent_[e].second << std::endl;
//    }

    std::list<Gate> ret;

    if (Init(in, out) && index_list.empty())
    {
        return ret;
    }

    std::vector<util::xor_func> bits(qubit_num_);
    std::vector<util::xor_func> preparation(qubit_num_);
    std::vector<util::xor_func> restoration(qubit_num_);
    for (int i = 0; i < qubit_num_; i++)
    {
        preparation[i] = util::xor_func(qubit_num_ + 1, 0);
        restoration[i] = util::xor_func(qubit_num_ + 1, 0);
        preparation[i].set(i);
        restoration[i].set(i);
    }
    std::vector<std::pair<int, int>> phase_target_list;

    /*
     * Reduce in to echelon form to decide on a basis
     */
    util::ToUpperEchelon(qubit_num_, dimension_, in, &preparation, std::vector<std::string>());

    while (!index_list.empty())
    {
//        std::cout << "-- index list size: " << index_list.size() << std::endl;
//        std::cout << "-- index list in loop" << std::endl;
//        for (auto&& e : index_list)
//        {
//            std::cout << e << ":" << phase_exponent_[e].second << std::endl;
//        }

        /**
         * first build sub-circuit
         */
        std::vector<util::xor_func> result_bits(qubit_num_);
        std::vector<util::xor_func> result_preparation(qubit_num_);
        std::vector<util::xor_func> result_restoration(qubit_num_);
        for (int i = 0; i < qubit_num_; i++)
        {
            result_preparation[i] = util::xor_func(qubit_num_ + 1, 0);
            result_restoration[i] = util::xor_func(qubit_num_ + 1, 0);
            result_preparation[i].set(i);
            result_restoration[i].set(i);
        }
        std::list<Gate> result_gate_list;
        std::set<int> result_sub_part;
        std::list<int> delete_index_list;
        std::vector<std::pair<int, int>> result_phase_target_list;

        for (auto it = index_list.begin(); it != index_list.end();)
        {
            std::list<Gate> tmp_gate_list;
            std::set<int> tmp_sub_part = result_sub_part;
            tmp_sub_part.insert(*it);

            std::vector<util::xor_func> tmp_bits(qubit_num_);
            std::vector<util::xor_func> tmp_preparation = preparation;
            std::vector<util::xor_func> tmp_restoration = restoration;

            std::vector<std::pair<int, int>> tmp_phase_target_list;

            if (oracle_(phase_exponent_, tmp_sub_part))
            {
//                std::cout << "oracle part num: " << tmp_sub_part.size() << std::endl;
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
//                std::cout << std::endl;
//                std::cout << "# prepare" << std::endl;
//                std::cout << "## num qubit: " << qubit_num_ << std::endl;
//                std::cout << "## dimension: " << dimension_ << std::endl;
//                std::cout << "## num partition: " << num_partition << std::endl;
//                std::cout << "## in" << std::endl;
//                for (auto&& e : in)
//                {
//                    std::cout << e << std::endl;
//                }
//                std::cout << "## tmp_bits" << std::endl;
//                for (auto&& e : tmp_bits)
//                {
//                    std::cout << e << std::endl;
//                }
//                std::cout << "## tmp_restoration" << std::endl;
//                for (auto&& e : tmp_restoration)
//                {
//                    std::cout << e << std::endl;
//                }
                util::ToUpperEchelon(num_partition, dimension_, tmp_bits, &tmp_restoration, std::vector<std::string>());
                util::FixBasis(qubit_num_, dimension_, num_partition, in, tmp_bits, &tmp_restoration,
                               std::vector<std::string>());
                util::Compose(qubit_num_, tmp_preparation, tmp_restoration);

                /**
                 * change row order in the matrix
                 */
                if (option_.change_row_order())
                {
                    // ChangeRowOrder();
                    ChangeRowOrder(tmp_phase_target_list, tmp_preparation);
                }

                /**
                 * create gate list
                 */
                if (option_.dec_type() == DecompositionType::kgauss)
                {
                    tmp_gate_list.splice(tmp_gate_list.end(),
                                         (*decomposer_)(qubit_num_, 0, tmp_preparation, qubit_names_));
                }

                /**
                 * check time step
                 */
                if (ComputeTimeStep(tmp_gate_list) < option_.distillation_step())
                {
                    result_bits = tmp_bits;
                    result_preparation = tmp_preparation;
                    result_restoration = tmp_restoration;
                    result_sub_part = tmp_sub_part;
                    result_gate_list = tmp_gate_list;
                    result_phase_target_list = tmp_phase_target_list;
//                    delete_index_list.push_back(*it);
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

//        index_list.splice(index_list.end(), std::move(delete_index_list));
//        index_list.unique();
//        std::cout << "after unique" << std::endl;
//        for (auto&& e : index_list)
//        {
//            std::cout << e << " ";
//        }
//        std::cout << std::endl;

        /**
         * second build sub-circuit
         */
        for (auto it = carry_index_list.begin(); it != carry_index_list.end();)
        {
            std::list<Gate> tmp_gate_list;
            std::set<int> tmp_sub_part = result_sub_part;
            tmp_sub_part.insert(*it);

            std::vector<util::xor_func> tmp_bits(qubit_num_);
            std::vector<util::xor_func> tmp_preparation = preparation;
            std::vector<util::xor_func> tmp_restoration = restoration;

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
                util::ToUpperEchelon(num_partition, dimension_, tmp_bits, &tmp_restoration, std::vector<std::string>());
                util::FixBasis(qubit_num_, dimension_, num_partition, in, tmp_bits, &tmp_restoration,
                               std::vector<std::string>());
                util::Compose(qubit_num_, tmp_preparation, tmp_restoration);

                /**
                 * change row order in the matrix
                 */
                if (option_.change_row_order())
                {
                    // ChangeRowOrder();
                }

                /**
                 * create gate list
                 */
                if (option_.dec_type() == DecompositionType::kgauss)
                {
                    tmp_gate_list.splice(tmp_gate_list.end(),
                                         (*decomposer_)(qubit_num_, 0, tmp_preparation, qubit_names_));
                }

                /**
                 * check time step
                 */
                if (ComputeTimeStep(tmp_gate_list) < option_.distillation_step())
                {
                    result_bits = tmp_bits;
                    result_preparation = tmp_preparation;
                    result_restoration = tmp_restoration;
                    result_sub_part = tmp_sub_part;
                    result_gate_list = tmp_gate_list;
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
        ApplyPhaseGates(ret, result_sub_part);

        /*
         * Unprepare the bits
         */
        preparation = std::move(result_restoration);
        restoration = std::vector<util::xor_func>(qubit_num_);
        // re-initialize
        for (int i = 0; i < qubit_num_; i++)
        {
            restoration[i] = util::xor_func(qubit_num_ + 1, 0);
            restoration[i].set(i);
        }
    }

    /*
     * Reduce out to the basis of in
     */
    for (int i = 0; i < qubit_num_; i++)
    {
        bits[i] = out[i];
    }
    util::ToUpperEchelon(qubit_num_, dimension_, bits, &restoration, std::vector<std::string>());
    util::FixBasis(qubit_num_, dimension_, qubit_num_, in, bits, &restoration, std::vector<std::string>());
    util::Compose(qubit_num_, preparation, restoration);

    if (option_.dec_type() == DecompositionType::kgauss)
    {
        ret.splice(ret.end(), (*decomposer_)(qubit_num_, 0, preparation, qubit_names_));
    }

    return ret;
}

std::list<Gate> GreedyCircuitBuilder::BuildGlobalPhase(int qubit_num,
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
