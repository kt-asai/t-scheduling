#ifndef T_SCHEDULING_PARALLELIZATION_ORACLE_HPP
#define T_SCHEDULING_PARALLELIZATION_ORACLE_HPP

#include <list>

#include "../circuit/gate.hpp"

namespace tskd {

class ParallelizationOracle
{
private:
    int num_qubit_;

    std::list<Gate> gate_list_;

public:
    ParallelizationOracle() = default;

    ParallelizationOracle(const int num_qubit,
                          const std::list<Gate>& mapped_gate_list,
                          const Gate& new_gate)
        : num_qubit_(num_qubit),
          gate_list_(mapped_gate_list)
    {
        gate_list_.push_back(new_gate);
    }

    bool check();
};

}

#endif //T_SCHEDULING_PARALLELIZATION_ORACLE_HPP
