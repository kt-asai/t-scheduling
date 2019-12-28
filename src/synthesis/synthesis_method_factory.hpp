#ifndef T_SCHEDULING_SYNTHESIS_METHOD_FACTORY_HPP
#define T_SCHEDULING_SYNTHESIS_METHOD_FACTORY_HPP

#include "synthesis.hpp"
#include "tpar_synthesis.hpp"
#include "tskd_synthesis.hpp"

#include "../character/character.hpp"

#include "../util/option.hpp"

namespace tskd {

class SynthesisMethodFactory
{
private:
    Synthesis* synthesis_;

public:
    SynthesisMethodFactory() = default;

    Synthesis* create(const SynthesisMethod synthesis_method,
                      const Character& chr,
                      const util::Option& option)
    {
        switch (synthesis_method)
        {
            case SynthesisMethod::ktpar:
                synthesis_ = new TparSynthesis(chr, option);
                break;
            case SynthesisMethod::ktskd:
                synthesis_ = new TskdSynthesis(chr, option);
                break;
            default:
                synthesis_ = nullptr;
        }

        return synthesis_;
    }
};

}

#endif //T_SCHEDULING_SYNTHESIS_METHOD_FACTORY_HPP
