#include "tskd_synthesis.hpp"

namespace tskd {

void TskdSynthesis::Init(const Character& chr)
{
    global_phase_ = 0;
    remaining_.resize(2);
    part_.resize(2);

    /**
     * sort phase exponents
     */
    auto compare = [](auto& lhs, auto& rhs) -> bool
    {
        if (lhs.second.count() == rhs.second.count())
        {
            return lhs.second < rhs.second;
        }
        else
        {
            return lhs.second.count() < rhs.second.count();
        }
    };
    chr_.SortPhaseExponents(compare);

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

template<class T, typename oracle_type>
tpar::partitioning TskdSynthesis::AddPartition(std::list<int>& remaining_index_list_,
                                               const std::vector<T>& elts,
                                               const oracle_type& oracle)
{
    const int limited_step = option_.distillation_step();
    tpar::partitioning partition;

    bool update = true;
    while (update)
    {
        update = false;
        std::set<int> sub_part;
        for (auto it = remaining_index_list_.begin(); it != remaining_index_list_.end();)
        {
            util::xor_func tmp = (~mask_) & (chr_.phase_exponents()[*it].second);
            if (tmp.none())
            {
                sub_part.insert(*it);
                std::set<int> tmp_sub_part = sub_part;
                if (oracle_(elts, tmp_sub_part) && tmp_sub_part.size() < 3/* && build() > limited_step * tmmp_sub_part.size() */)
                {
                    sub_part = tmp_sub_part;
                    it = remaining_index_list_.erase(it);
                    update = true;
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
        partition.push_back(sub_part);
    }

    return partition;
}

void TskdSynthesis::CreatePartition()
{
    for (int i = 0; i < 2; ++i)
    {
        part_[i] = AddPartition(remaining_[i], chr_.phase_exponents(), oracle_);
    }
}

Circuit TskdSynthesis::Execute()
{
    std::cout << "t-scheduling running..." << std::endl;

    /**
     *
     */

    for (auto&& hadamard : chr_.hadamards())
    {

    }

    return tskd::Circuit();
}

}
