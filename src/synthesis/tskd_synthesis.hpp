#ifndef T_SCHEDULING_TSKD_SYNTHESIS_HPP
#define T_SCHEDULING_TSKD_SYNTHESIS_HPP

#include "synthesis.hpp"

#include "greedy_circuit_builder.hpp"

#include "../util/util.hpp"

namespace tskd {

class TskdSynthesis : public Synthesis
{
private:
    Character chr_;

    util::Option option_;

    GreedyCircuitBuilder builder_;

    util::IndependentOracle oracle_;

    Circuit circuit_;

    int global_phase_;

    std::vector<std::list<int>> index_list_;
    std::vector<std::list<int>> carry_index_list_;

    util::xor_func mask_;

    std::vector<util::xor_func> wires_;
    std::vector<std::list<int>> remaining_;

    void determine_apply_phase_set(Character::Hadamard& hadamard);

    void construct_subcircuit(const Character::Hadamard& hadamard);

    void apply_hadamard(const Character::Hadamard& hadamard);

    void construct_final_subcircuit();

public:
    TskdSynthesis(const Character& chr,
                  const util::Option& option)
            : chr_(chr),
              option_(option)
    {
        init(chr);

        oracle_ = util::IndependentOracle(chr.num_qubit(),
                                          chr.num_data_qubit(),
                                          chr.num_data_qubit() + chr.num_hadamard());

        builder_ = GreedyCircuitBuilder(option,
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
