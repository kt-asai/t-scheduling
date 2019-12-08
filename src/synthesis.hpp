#ifndef T_SCHEDULING_SYNTHESIS_HPP
#define T_SCHEDULING_SYNTHESIS_HPP

#include "character.hpp"
#include "util.hpp"

#include "tpar/partition.hpp"

namespace tskd {

class Synthesis
{
private:
    Character chr_;

    Circuit circuit_;

    int global_phase_;

    std::vector<tpar::partitioning> floats_;
    std::vector<tpar::partitioning> frozen_;

    util::xor_func mask_;

    std::vector<util::xor_func> wires_;
    std::vector<std::list<int>> remaining_;

public:
    Synthesis(const Character& chr) : chr_(chr)
    {
        init(chr);
    }

    void init(const Character& chr);

    Circuit Execute();
};

}


#endif //T_SCHEDULING_SYNTHESIS_HPP
