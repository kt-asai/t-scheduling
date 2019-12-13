#ifndef T_SCHEDULING_TPAR_SYNTHESIS_HPP
#define T_SCHEDULING_TPAR_SYNTHESIS_HPP

#include "synthesis.hpp"

#include "../character.hpp"
#include "../circuit_builder.hpp"

#include "../util/util.hpp"
#include "../util/option.hpp"

#include "../tpar/partition.hpp"

namespace tskd {

class TparSynthesis : public Synthesis
{
private:
    Character chr_;

    util::Option option_;

    CircuitBuilder builder_;

    util::IndependentOracle oracle_;

    Circuit circuit_;

    int global_phase_;

    std::vector<tpar::partitioning> floats_;
    std::vector<tpar::partitioning> frozen_;

    util::xor_func mask_;

    std::vector<util::xor_func> wires_;
    std::vector<std::list<int>> remaining_;

    void CreatePartition();

    void DetermineApplyPartition(Character::Hadamard& hadamard);

    void ConstructSubCircuit(const Character::Hadamard& hadamard);

    void ApplyHadamard(const Character::Hadamard& hadamard);

    int CheckDimension(int dimension);

    void ConstructFinalSubCircuit();

public:
    TparSynthesis(const Character& chr,
                  const util::Option& option)
        : chr_(chr),
          option_(option)
    {
        init(chr);

        builder_ = CircuitBuilder(option, chr.num_qubit(), chr.num_data_qubit() + chr.num_hadamard(),
                                  chr.qubit_names(), chr.phase_exponents());

        oracle_ = util::IndependentOracle(chr.num_qubit(), chr.num_data_qubit(), chr.num_data_qubit() + chr.num_hadamard());
    }

    void init(const Character& chr);

    Circuit Execute() final;
};

}


#endif //T_SCHEDULING_TPAR_SYNTHESIS_HPP
