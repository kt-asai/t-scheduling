#include "simulator.hpp"

namespace tskd {

void Simulator::update_distillations(int inc_step,
                                     std::set<int>& used_index_set)
{
    if (use_buffer_)
    {
        if (buffer_capacity_ < buffer_size_)
        {
            for (size_t i = 0; i < distillations_.size(); i++)
            {
                distillations_[i] = std::min(distillations_[i] + inc_step, option_.distillation_step());
            }
        }
    }
    else
    {
        for (size_t i = 0; i < distillations_.size(); i++)
        {
            if (used_index_set.count(i) == 0)
            {
                distillations_[i] = std::min(distillations_[i] + inc_step, option_.distillation_step());
            }
        }
    }
}

void Simulator::update_buffer_capacity(int inc_step,
                                       std::set<int>& used_index_set)
{
    update_distillations(inc_step, used_index_set);

    if (use_buffer_)
    {
        for (size_t i = 0; i < distillations_.size(); i++)
        {
            if (buffer_capacity_ == buffer_size_)
            {
                break;
            }

            if (distillations_[i] == option_.distillation_step())
            {
                buffer_capacity_++;
                distillations_[i] = 0;
            }
        }

        // if capacity of buffer is full, stop distillation procedure
        if (buffer_capacity_ == buffer_size_)
        {
            for (size_t i = 0; i < distillations_.size(); i++)
            {
                distillations_[i] = 0;
            }
        }
    }
    else
    {
        buffer_capacity_ = 0;
        for (size_t i = 0; i < distillations_.size(); i++)
        {
            if (distillations_[i] == option_.distillation_step())
            {
                buffer_capacity_++;
            }
        }
    }
}

int Simulator::get_magic_state(int required_magic_states,
                               std::set<int>& used_index_set)
{
    int ret_magic_state = 0;
    if (use_buffer_)
    {
        ret_magic_state = std::min(buffer_capacity_, std::min(required_magic_states, option_.num_distillation() * 2));
        buffer_capacity_ -= ret_magic_state;
    }
    else
    {
        for (size_t i = 0; i < distillations_.size(); i++)
        {
            if (ret_magic_state == required_magic_states)
            {
                break;
            }

            if (distillations_[i] == option_.distillation_step())
            {
                ret_magic_state++;
                used_index_set.insert(i);
                distillations_[i] = 0;
                buffer_capacity_--;
            }
        }
    }

    return ret_magic_state;
}

void Simulator::calculate()
{
    constexpr int inc_time_step = 2;
    int current_time_step = 0;
    GateType previous_gate_type = GateType::kphase;

    int num_t_gate = 0;
    int num_p_gate = 0;
    int num_z_gate = 0;

    std::set<int> used_index_set;
    std::set<int> none;

    for (auto&& gate : circuit_.gate_list())
    {
        if (previous_gate_type == GateType::kphase
            && (gate.type() == "H" ||gate.type() == "cnot" || gate.type() == "tof"))
        {
            if (num_p_gate > 0 && num_t_gate == 0)
            {
                update_buffer_capacity(inc_time_step, none);
                current_time_step += inc_time_step;
            }
            if (num_t_gate > 0)
            {
                // calculate buffer capacity
                while (num_t_gate > 0)
                {
                    num_t_gate -= get_magic_state(num_t_gate, used_index_set);
                    update_buffer_capacity(3, used_index_set);
                    used_index_set.clear();
                    current_time_step += 3;
                }
            }

            num_t_gate = 0;
            num_p_gate = 0;
            num_z_gate = 0;
        }

        if (gate.type() == "H")
        {
            if (previous_gate_type != GateType::khadamard)
            {
                previous_gate_type = GateType::khadamard;
                update_buffer_capacity(inc_time_step, none);
                current_time_step += inc_time_step;
            }

            continue;
        }

        if (gate.type() == "cnot" || gate.type() == "tof")
        {
            previous_gate_type = GateType::kcnot;
            update_buffer_capacity(inc_time_step, none);
            current_time_step += inc_time_step;

            continue;
        }

        if (gate.type() == "T" || gate.type() == "T*")
        {
            previous_gate_type = GateType::kphase;
            num_t_gate++;

            continue;
        }

        if (gate.type() == "P" || gate.type() == "P*")
        {
            previous_gate_type = GateType::kphase;
            num_p_gate++;

            continue;
        }

        if (gate.type() == "Z")
        {
            previous_gate_type = GateType::kphase;
            num_z_gate++;
        }
    }

    result_time_step_ = current_time_step;
}

}