#ifndef T_SCHEDULING_TSKD_SYNTHESIS_HPP
#define T_SCHEDULING_TSKD_SYNTHESIS_HPP

#include "synthesis.hpp"

#include "../circuit_builder.hpp"

#include "../util/util.hpp"

namespace tskd {

class TskdSynthesis : public Synthesis
{
private:
    Character chr_;

    util::Option option_;

    CircuitBuilder builder_;

    util::IndependentOracle oracle_;

    Circuit circuit_;

    int global_phase_;

    std::vector<tpar::partitioning> part_;

    util::xor_func mask_;

    std::vector<util::xor_func> wires_;
    std::vector<std::list<int>> remaining_;

    template<class T, typename oracle_type>
    tpar::partitioning AddPartition(std::list<int>& remaining_index_list_,
                                    const std::vector <T>& elts,
                                    const oracle_type& oracle);
    
    void CreatePartition();

public:
    TskdSynthesis(const Character& chr,
                  const util::Option& option)
            : chr_(chr),
              option_(option)
    {
        Init(chr);

        builder_ = CircuitBuilder(option, chr.num_qubit(), chr.num_data_qubit() + chr.num_hadamard(),
                                  chr.qubit_names(), chr.phase_exponents());

        oracle_ = util::IndependentOracle(chr.num_qubit(), chr.num_data_qubit(), chr.num_data_qubit() + chr.num_hadamard());
    }

    void Init(const Character& chr);

    Circuit Execute() final;
};

}


#endif //T_SCHEDULING_TSKD_SYNTHESIS_HPP
