#ifndef T_SCHEDULING_CHARACTER_HPP
#define T_SCHEDULING_CHARACTER_HPP

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <boost/dynamic_bitset.hpp>

#include "circuit.hpp"
#include "gate.hpp"


class Character
{
private:
    using xor_func = boost::dynamic_bitset<>;

    Circuit qc_;

    std::unordered_map<std::string, int> name_map_;
    std::unordered_map<std::string, int> result_val_map_;
    std::vector<std::string> phase_parity_list_;
    std::vector<std::string> hadamard_parity_list_;

    struct Hadamard
    {
        int target_;
        std::vector<xor_func> input_parity_;
    };

public:
    /**
     * constructor
     */
    Character() = default;

    /**
     * constructor
     * @param qc circuit info
     */
    Character(const Circuit& qc) : qc_(qc) { }

    /**
     * parse circuit info and calculate parity
     */
    void Parse();
};

#endif //T_SCHEDULING_CHARACTER_HPP
