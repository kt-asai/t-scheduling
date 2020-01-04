#ifndef T_SCHEDULING_GREEDY_CIRCUIT_BUILDER_HPP
#define T_SCHEDULING_GREEDY_CIRCUIT_BUILDER_HPP

#include <vector>
#include <list>
#include <string>
#include <memory>

#include "../util/option.hpp"

#include "../layout/layout.hpp"

#include "../circuit/gate.hpp"
#include "../circuit/circuit.hpp"

#include "../character/character.hpp"

#include "../tpar/partition.hpp"

#include "../decomposer/matrix_decomposer.hpp"
#include "../decomposer/gaussian_decomposer.hpp"
#include "../decomposer/parallel_decomposer.hpp"

#include "../matrix/matrix_reconstructor.hpp"


namespace tskd {

/**
 * T-scheduling
 * this class build {CNOT, T} sub-circuit for given partitions greedy
 */
class GreedyCircuitBuilder
{
private:
    util::Option option_;

    Layout layout_;

    std::shared_ptr <MatrixDecomposer> decomposer_;

    util::IndependentOracle oracle_;

    int qubit_num_;
    int dimension_;

    std::vector<std::string> qubit_names_;
    std::vector<util::phase_exponent> phase_exponent_;

    std::vector<util::xor_func> bits_;
    std::vector<util::xor_func> preparation_;
    std::vector<util::xor_func> restoration_;

    std::vector<util::xor_func> identity_;

    bool init(const std::vector <util::xor_func>& in,
              const std::vector <util::xor_func>& out);

    int compute_time_step(const std::list<Gate>& gate_list);

    void apply_phase_gates(std::list<Gate>& gate_list,
                           const std::unordered_map<int, int>& target_phase_map);

    void unprepare(const std::vector <util::xor_func>& restoration);

    void prepare_last_part(std::list<Gate>& gate_list,
                           const std::vector<util::xor_func>& in,
                           std::vector<util::xor_func>& out,
                           MatrixReconstructor& sa);

public:
    GreedyCircuitBuilder() = default;

    template<typename oracle_type>
    GreedyCircuitBuilder(const util::Option& option,
                         const Layout& layout,
                         const oracle_type& oracle,
                         int qubit_num,
                         int dimension,
                         const std::vector <std::string>& qubit_names,
                         const std::vector <util::phase_exponent>& phase_exponent)
            : option_(option),
              layout_(layout),
              oracle_(oracle),
              qubit_num_(qubit_num),
              dimension_(dimension),
              qubit_names_(qubit_names),
              phase_exponent_(phase_exponent)
    {
        if (option.dec_type() == DecompositionType::kgauss)
        {
            decomposer_ = std::make_shared<GaussianDecomposer>(layout, qubit_num, 0, qubit_names);
        }
        else if (option.dec_type() == DecompositionType::kparallel)
        {
            decomposer_ = std::make_shared<ParallelDecomposer>(layout, qubit_num, 0, qubit_names);
        }
        else
        {
            std::cerr << "invalid type of decomposer" << std::endl;

            exit(1);
        }
    }

    std::list<Gate> build(std::list<int>& index_list,
                          std::list<int>& carry_index_list,
                          std::vector <util::xor_func>& in,
                          std::vector <util::xor_func>& out);

    static std::list<Gate> build_global_phase(int qubit_num,
                                              int phase,
                                              const std::vector <std::string>& qubit_names);

    int check_dimension(const Character& chr,
                        std::vector <util::xor_func>& wires,
                        int current_dimension);
};

}

#endif //T_SCHEDULING_GREEDY_CIRCUIT_BUILDER_HPP
