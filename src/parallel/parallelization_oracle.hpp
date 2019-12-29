#ifndef T_SCHEDULING_PARALLELIZATION_ORACLE_HPP
#define T_SCHEDULING_PARALLELIZATION_ORACLE_HPP

#include <list>

#include "../circuit/gate.hpp"

namespace tskd {

class ParallelizationOracle
{
private:
    int num_qubit_;
    int num_distillation_;

    void init();

public:
    ParallelizationOracle() = default;

    ParallelizationOracle(const int num_qubit,
                          const int num_distillation)
        : num_qubit_(num_qubit),
          num_distillation_(num_distillation)
      {
        init();
      }

    bool check(const std::vector<Gate>& gate_list);
};

}

#endif //T_SCHEDULING_PARALLELIZATION_ORACLE_HPP
