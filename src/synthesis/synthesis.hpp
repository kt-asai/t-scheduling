#ifndef T_SCHEDULING_SYNTHESIS_HPP
#define T_SCHEDULING_SYNTHESIS_HPP

#include "../character.hpp"
#include "../circuit_builder.hpp"

#include "../util/option.hpp"
#include "../util/util.hpp"

namespace tskd {

class Synthesis
{
private:

public:
    Synthesis() = default;

    virtual ~Synthesis() = default;

    virtual Circuit Execute() = 0;
};

}

#endif //T_SCHEDULING_SYNTHESIS_HPP
