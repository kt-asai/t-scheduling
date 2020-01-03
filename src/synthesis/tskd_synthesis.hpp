#ifndef T_SCHEDULING_TSKD_SYNTHESIS_HPP
#define T_SCHEDULING_TSKD_SYNTHESIS_HPP

#include "synthesis.hpp"

#include "greedy_circuit_builder.hpp"

#include "../util/util.hpp"

#include "../layout/layout.hpp"

namespace tskd {

class TskdSynthesis : public Synthesis
{
private:
    util::Option option_;

    Character chr_;

    GreedyCircuitBuilder builder_;

    util::IndependentOracle oracle_;

    Circuit circuit_;

    int global_phase_;

    std::list<int> index_list_;
    std::list<int> carry_index_list_;

    util::xor_func mask_;

    std::vector<util::xor_func> wires_;
    std::list<int> remaining_;

    std::vector<int> bit_map_; // [from] = to

    void determine_apply_phase_set(Character::Hadamard& hadamard);

    void construct_subcircuit(Character::Hadamard& hadamard);

    void apply_hadamard(const Character::Hadamard& hadamard);

    void construct_final_subcircuit();

    inline void update_bit_map(const std::vector<util::xor_func>& original,
                               const std::vector<util::xor_func>& result)
    {
        std::set<int> used_bit_index_set;
        for (size_t i = 0; i < original.size(); i++)
        {
            for (size_t j = 0; j < result.size(); j++)
            {
                if (!used_bit_index_set.count(j) && original[i] == result[j])
                {
                    bit_map_[i] = j;
                    used_bit_index_set.insert(j);
                    break;
                }
            }
        }
    }

public:
    TskdSynthesis(const util::Option& option,
                  const Layout& layout,
                  const Character& chr)
            : option_(option),
              chr_(chr)
    {
        init(chr);

        oracle_ = util::IndependentOracle(chr.num_qubit(),
                                          chr.num_data_qubit(),
                                          chr.num_data_qubit() + chr.num_hadamard());

        builder_ = GreedyCircuitBuilder(option,
                                        layout,
                                        oracle_,
                                        chr.num_qubit(),
                                        chr.num_data_qubit() + chr.num_hadamard(),
                                        chr.qubit_names(),
                                        chr.phase_exponents());
    }

    void init(const Character& chr);

    Circuit execute() final;
};

}


#endif //T_SCHEDULING_TSKD_SYNTHESIS_HPP
