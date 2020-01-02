#ifndef T_SCHEDULING_TPAR_SYNTHESIS_HPP
#define T_SCHEDULING_TPAR_SYNTHESIS_HPP

#include "synthesis.hpp"

#include "simple_circuit_builder.hpp"

#include "../character/character.hpp"

#include "../util/util.hpp"
#include "../util/option.hpp"

#include "../layout/layout.hpp"

#include "../tpar/partition.hpp"

namespace tskd {

class TparSynthesis : public Synthesis
{
private:
    util::Option option_;

    Character chr_;

    SimpleCircuitBuilder builder_;

    util::IndependentOracle oracle_;

    Circuit circuit_;

    int global_phase_;

    std::vector<tpar::partitioning> floats_;
    std::vector<tpar::partitioning> frozen_;

    util::xor_func mask_;

    std::vector<util::xor_func> wires_;
    std::vector<std::list<int>> remaining_;

    std::vector<int> bit_map_; // [from] = to

    void create_partition();

    void determine_apply_partition(Character::Hadamard& hadamard);

    void construct_subcircuit(Character::Hadamard& hadamard);

    void apply_hadamard(const Character::Hadamard& hadamard);

    int check_dimension(int dimension);

    void construct_final_subcircuit();

    inline void update_bit_map(const std::vector<util::xor_func>& original,
                               const std::vector<util::xor_func>& result)
    {
        std::set<int> used_bit_index_set;
        for (size_t i = 0; i < original.size(); i++)
        {
            for (size_t j = 0; j < result.size(); j++)
            {
                if (!used_bit_index_set.count(j) && original[i] == result[j])
                {
                    bit_map_[i] = j;
                    used_bit_index_set.insert(j);
                    break;
                }
            }
        }
    }

//    inline void update_bit_map(const std::vector<util::xor_func>& original,
//                               const std::vector<util::xor_func>& result)
//    {
////        std::cout << "before bit_map" << std::endl;
////        for (size_t i = 0; i < chr_.num_qubit(); i++)
////        {
////            std::cout << i << ":" << bit_map_[i] << std::endl;
////        }
//
//        std::vector<int> tmp_bit_map(chr_.num_qubit());
//        std::set<int> used_bit_index_set;
//        for (size_t i = 0; i < original.size(); i++)
//        {
//            for (size_t j = 0; j < result.size(); j++)
//            {
//                if (!used_bit_index_set.count(j) && original[i] == result[j])
//                {
//                    tmp_bit_map[i] = bit_map_[j];
//                    used_bit_index_set.insert(j);
//                    break;
//                }
//            }
//        }
//        bit_map_ = tmp_bit_map;
//
////        std::cout << "after bit_map" << std::endl;
////        for (size_t i = 0; i < chr_.num_qubit(); i++)
////        {
////            std::cout << i << ":" << bit_map_[i] << std::endl;
////        }
//    }

public:
    TparSynthesis(const util::Option& option,
                  const Layout& layout,
                  const Character& chr)
            : option_(option),
              chr_(chr)
    {
        init(chr);

        oracle_ = util::IndependentOracle(chr.num_qubit(),
                                          chr.num_data_qubit(),
                                          chr.num_data_qubit() + chr.num_hadamard());

        builder_ = SimpleCircuitBuilder(option,
                                        layout,
                                        chr.num_qubit(),
                                        chr.num_data_qubit() + chr.num_hadamard(),
                                        chr.qubit_names(),
                                        chr.phase_exponents());
    }

    void init(const Character& chr);

    Circuit execute() final;
};

}


#endif //T_SCHEDULING_TPAR_SYNTHESIS_HPP
