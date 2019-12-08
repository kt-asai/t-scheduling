#ifndef T_SCHEDULING_SYNTHESIS_HPP
#define T_SCHEDULING_SYNTHESIS_HPP

#include "character.hpp"

namespace tskd {

class Synthesis
{
private:
    Character character_;

    int num_qubit_;
    int num_ancilla_;
    int num_hadamard_;
    int dimension_;

    std::vector<std::string> qubit_names_;
    std::vector<bool> ancilla_list_;
    std::vector<util::phase_exponent> phase_exponents_;
    std::vector<Character::Hadamard> hadamards_;

public:
    Synthesis(Character character)
        : character_(character)
    {
        num_qubit_ = character.num_qubit();
        num_ancilla_ = character.num_ancilla_qubit();
        num_hadamard_ = character.num_hadamard();
        dimension_= character.num_qubit() - character.num_ancilla_qubit();

        qubit_names_ = character.qubit_names();
        ancilla_list_ = character.ancilla_list();
        phase_exponents_ = character.phase_exponents();
        hadamards_ = character.hadamards();
    }

    Circuit Execute();
};

}


#endif //T_SCHEDULING_SYNTHESIS_HPP
