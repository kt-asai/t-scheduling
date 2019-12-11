#ifndef T_SCHEDULING_SYNTHESIS_HPP
#define T_SCHEDULING_SYNTHESIS_HPP

#include "character.hpp"
#include "util/util.hpp"
#include "util/option.hpp"

#include "tpar/partition.hpp"

namespace tskd {

class Synthesis
{
private:
    Character chr_;

    util::Option option_;

    Circuit circuit_;

    int global_phase_;

    std::vector<tpar::partitioning> floats_;
    std::vector<tpar::partitioning> frozen_;

    util::xor_func mask_;

    std::vector<util::xor_func> wires_;
    std::vector<std::list<int>> remaining_;

    template<typename oracle_type>
    void CreatePartition(const oracle_type& oracle);

    void ConstructCircuit();

    void ApplyHadamard(const Character::Hadamard& hadamard);

public:
    Synthesis(const Character& chr,
              const util::Option& option)
        : chr_(chr),
          option_(option)
    {
        init(chr);
    }

    void init(const Character& chr);

    Circuit Execute();
};

}


#endif //T_SCHEDULING_SYNTHESIS_HPP
