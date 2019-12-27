#ifndef T_SCHEDULING_SIMPLE_CIRCUIT_BUILDER_HPP
#define T_SCHEDULING_SIMPLE_CIRCUIT_BUILDER_HPP

#include <vector>
#include <list>
#include <string>
#include <memory>

#include "../util/option.hpp"

#include "../circuit/gate.hpp"
#include "../circuit/circuit.hpp"

#include "../tpar/partition.hpp"

#include "../decomposer/matrix_decomposer.hpp"
#include "../decomposer/gaussian_decomposer.hpp"
#include "../decomposer/parallel_decomposer.hpp"

#include "../matrix/matrix_reconstructor.hpp"

namespace tskd {

/**
 * this class build {CNOT, T} sub-circuit for given partitions
 */
class SimpleCircuitBuilder
{
private:
    util::Option option_;

    std::shared_ptr<MatrixDecomposer> decomposer_;

    int qubit_num_;
    int dimension_;

    std::vector<std::string> qubit_names_;
    std::vector<util::phase_exponent> phase_exponent_;

    std::vector<util::xor_func> bits_;
    std::vector<util::xor_func> preparation_;
    std::vector<util::xor_func> restoration_;

    bool Init(const std::vector<util::xor_func>& in,
              const std::vector<util::xor_func>& out);

    void InitBits(const std::set<int>& phase_exponent_index_set,
                  std::unordered_map<int, int>& target_phase_map);

    void Prepare(std::list<Gate>& gate_list,
                 const std::vector<util::xor_func>& in,
                 const int num_partition,
                 std::unordered_map<int, int>& target_phase_map);

    void ApplyPhaseGates(std::list<Gate>& gate_list,
                         const std::unordered_map<int, int>& target_phase_map);

    void UnPrepare();

    void PrepareLastPart(std::list<Gate>& gate_list,
                         const std::vector<util::xor_func>& in,
                         const std::vector<util::xor_func>& out,
                         MatrixReconstructor& sa);

public:
    SimpleCircuitBuilder() = default;

    SimpleCircuitBuilder(const util::Option& option,
                         int qubit_num,
                         int dimension,
                         const std::vector<std::string>& qubit_names,
                         const std::vector<util::phase_exponent>& phase_exponent)
        : option_(option),
          qubit_num_(qubit_num),
          dimension_(dimension),
          qubit_names_(qubit_names),
          phase_exponent_(phase_exponent)
    {
        if (option.dec_type() == DecompositionType::kgauss)
        {
            decomposer_ = std::make_shared<GaussianDecomposer>();
        }
        else if (option.dec_type() == DecompositionType::kparallel)
        {
            decomposer_ = std::make_shared<ParallelDecomposer>();
        }
        else
        {
            std::cerr << "invalid type of decomposer" << std::endl;

            exit(1);
        }
    }

    std::list<Gate> Build(const tpar::partitioning& partition,
                          std::vector<util::xor_func>& in,
                          const std::vector<util::xor_func>& out);

    static std::list<Gate> BuildGlobalPhase(int qubit_num,
                                            int phase,
                                            const std::vector<std::string>& qubit_names);
};

}

#endif //T_SCHEDULING_SIMPLE_CIRCUIT_BUILDER_HPP
