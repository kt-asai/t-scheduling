#ifndef T_SCHEDULING_SIMULATOR_HPP
#define T_SCHEDULING_SIMULATOR_HPP

#include <set>

#include "../util/option.hpp"

#include "../circuit/circuit.hpp"

namespace tskd {

class Simulator
{
private:
    util::Option option_;

    Circuit circuit_;

    bool use_buffer_;
    int buffer_size_;
    int buffer_capacity_;

    std::vector<int> distillations_;

    int result_time_step_;

    enum GateType
    {
        kcnot,
        khadamard,
        kphase,
    };

    inline bool capacity_is_full()
    {
        return buffer_capacity_ == buffer_size_;
    }

    void update_buffer_capacity(int inc_step,
                                std::set<int>& used_index_set);

    void update_distillations(int inc_step,
                              std::set<int>& used_index_set);

    int get_magic_state(int required_magic_states,
                        std::set<int>& used_index_set);

    void calculate();

public:
    Simulator() = default;

    Simulator(const util::Option& option,
              const Circuit& circuit)
        : option_(option),
          circuit_(circuit),
          use_buffer_(false),
          buffer_capacity_(0)
    {
        distillations_ = std::vector<int>(option.num_distillation(), 0);

        if (option.num_buffer() > 0)
        {
            buffer_size_ = option.num_buffer();
            use_buffer_ = true;
        }
        else
        {
            buffer_size_ = option.num_distillation();
        }


        calculate();
    }

    int result_time_step() const
    {
        return result_time_step_;
    }

    void print()
    {
        circuit_.print();

        std::cout << "# total time step: " << result_time_step_ << std::endl;
    }

};

}

#endif //T_SCHEDULING_SIMULATOR_HPP
