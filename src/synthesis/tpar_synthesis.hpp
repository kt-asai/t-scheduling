#ifndef T_SCHEDULING_TPAR_SYNTHESIS_HPP
#define T_SCHEDULING_TPAR_SYNTHESIS_HPP

#include "synthesis.hpp"

#include "simple_circuit_builder.hpp"

#include "../character/character.hpp"

#include "../util/util.hpp"
#include "../util/option.hpp"

#include "../tpar/partition.hpp"

namespace tskd {

class TparSynthesis : public Synthesis
{
private:
    Character chr_;

    util::Option option_;

    SimpleCircuitBuilder builder_;

    util::IndependentOracle oracle_;

    Circuit circuit_;

    int global_phase_;

    std::vector<tpar::partitioning> floats_;
    std::vector<tpar::partitioning> frozen_;

    util::xor_func mask_;

    std::vector<util::xor_func> wires_;
    std::vector<std::list<int>> remaining_;

    void create_partition();

    void determine_apply_partition(Character::Hadamard& hadamard);

    void construct_subcircuit(const Character::Hadamard& hadamard);

    void apply_hadamard(const Character::Hadamard& hadamard);

    int check_dimension(int dimension);

    void construct_final_subcircuit();

public:
    TparSynthesis(const Character& chr,
                  const util::Option& option)
        : chr_(chr),
          option_(option)
    {
        init(chr);

        builder_ = SimpleCircuitBuilder(option,
                                        chr.num_qubit(),
                                        chr.num_data_qubit() + chr.num_hadamard(),
                                        chr.qubit_names(),
                                        chr.phase_exponents());

        oracle_ = util::IndependentOracle(chr.num_qubit(),
                                          chr.num_data_qubit(),
                                          chr.num_data_qubit() + chr.num_hadamard());
    }

    void init(const Character& chr);

    Circuit execute() final;
};

}


#endif //T_SCHEDULING_TPAR_SYNTHESIS_HPP
