#include "character.hpp"

void Character::Parse()
{
    // gate value map
    const std::unordered_map<std::string, int> gate_lookup{
            {"T", 1},
            {"T*", 7},
            {"P", 2},
            {"P*", 6},
            {"Z", 4},
            {"Y", 4}
    };

    std::vector<xor_func> parity_list;

    // initialization
    for (const std::string& str : qc_.qubit_names())
    {

    }

    for (const Gate& gate : qc_.gate_list())
    {

    }
}