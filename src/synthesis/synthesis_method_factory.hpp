#ifndef T_SCHEDULING_SYNTHESIS_METHOD_FACTORY_HPP
#define T_SCHEDULING_SYNTHESIS_METHOD_FACTORY_HPP

#include "synthesis.hpp"
#include "tpar_synthesis.hpp"
#include "tskd_synthesis.hpp"

#include "../util/option.hpp"

#include "../layout/layout.hpp"

#include "../character/character.hpp"



namespace tskd {

class SynthesisMethodFactory
{
private:
    Synthesis* synthesis_;

public:
    SynthesisMethodFactory() = default;

    Synthesis* create(const SynthesisMethod synthesis_method,
                      const util::Option& option,
                      const Layout& layout,
                      const Character& chr)
    {
        switch (synthesis_method)
        {
            case SynthesisMethod::ktpar:
                synthesis_ = new TparSynthesis(option, layout, chr);
                break;
            case SynthesisMethod::ktskd:
                synthesis_ = new TskdSynthesis(option, layout, chr);
                break;
            default:
                synthesis_ = nullptr;
        }

        return synthesis_;
    }
};

}

#endif //T_SCHEDULING_SYNTHESIS_METHOD_FACTORY_HPP
