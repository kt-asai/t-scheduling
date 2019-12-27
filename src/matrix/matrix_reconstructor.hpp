#ifndef T_SCHEDULING_MATRIX_RECONSTRUCTOR_HPP
#define T_SCHEDULING_MATRIX_RECONSTRUCTOR_HPP

#include <random>
#include <unordered_map>

#include "../util/util.hpp"

namespace tskd {

class MatrixReconstructor
{
private:
    std::vector<util::xor_func> input_;

    int dimension_;
    int num_qubit_;

    // SA parameters
    std::random_device seed_generator_;
    std::mt19937 engine_;

    int rate_;
    std::chrono::milliseconds req_time_;

public:
    MatrixReconstructor(std::vector<util::xor_func>& input,
                       const int dimension,
                       const int num_qubit)
        : input_(input),
          dimension_(dimension),
          num_qubit_(num_qubit)
    {
        engine_ = std::mt19937(seed_generator_());
        rate_ = 10000;
        req_time_ = std::chrono::milliseconds(500);
    }

    std::vector<util::xor_func> Execute(const int num_partition,
                                        const std::vector<util::xor_func>& bits,
                                        std::unordered_map<int, int>& target_phase_map);
};

}


#endif //T_SCHEDULING_MATRIX_RECONSTRUCTOR_HPP
