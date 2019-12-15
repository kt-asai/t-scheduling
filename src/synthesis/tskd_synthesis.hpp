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

    util::xor_func mask_;

    std::vector<util::xor_func> wires_;
    std::vector<std::list<int>> remaining_;

public:
    TskdSynthesis(const Character& chr,
                  const util::Option& option)
            : chr_(chr),
              option_(option)
    {
        Init(chr);

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

    void Init(const Character& chr);

    Circuit Execute() final;
};

}


#endif //T_SCHEDULING_TSKD_SYNTHESIS_HPP
