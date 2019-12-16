#include "tskd_synthesis.hpp"

namespace tskd {

void TskdSynthesis::Init(const Character& chr)
{
    global_phase_ = 0;
    remaining_.resize(2);

    /**
     * sort phase exponents
     */
//    auto compare = [](auto& lhs, auto& rhs) -> bool
//    {
//        if (lhs.second.count() == rhs.second.count())
//        {
//            return lhs.second < rhs.second;
//        }
//        else
//        {
//            return lhs.second.count() < rhs.second.count();
//        }
//    };
//    chr_.SortPhaseExponents(compare);

    /*
     * initialize some stuff
     */
    mask_ = util::xor_func(chr.num_data_qubit() + chr.num_hadamard() + 1, 0);
    mask_.set(chr.num_data_qubit() + chr.num_hadamard());
    for (int i = 0, j = 0; i < chr.num_qubit(); i++)
    {
        circuit_.add_qubit(chr.qubit_names()[i]);
        circuit_.set_ancilla(chr.qubit_names()[i], chr.ancilla_list()[i]);
        wires_.emplace_back(chr.num_data_qubit() + chr.num_hadamard()+ 1, 0);
        if (!chr.ancilla_list()[i])
        {
            wires_[i].set(j);
            mask_.set(j);
            j++;
        }
    }

    /*
     * initialize the remaining list
     */
    int index = 0;
    for (auto&& phase_exponent : chr_.phase_exponents())
    {
        if (phase_exponent.second == util::xor_func(chr.num_data_qubit() + chr.num_hadamard() + 1, 0))
        {
            global_phase_ = phase_exponent.first;
        }
        else if (phase_exponent.first % 2 == 1)
        {
            remaining_[0].push_back(index);
        }
        else if (phase_exponent.first != 0)
        {
            remaining_[1].push_back(index);
        }
        index++;
    }
}

Circuit TskdSynthesis::Execute()
{
    std::cout << "t-scheduling running..." << std::endl;

    int dimension = chr_.num_data_qubit();

    /**
     *
     */

//    for (int i = 0; i < 2; ++i)
//    {
//        for (auto it = remaining_[i].begin(); it != remaining_[i].end();)
//        {
//            util::xor_func tmp = (~mask_) & (chr_.phase_exponents()[*it].second);
//            if (tmp.none())
//            {
//                index_list[i].push_back(*it);
//            }
//        }
//    }

    std::vector<std::list<int>> index_list(2);
    std::vector<std::list<int>> carry_index_list(2);

    for (auto&& hadamard : chr_.hadamards())
    {
//        std::cout << std::endl;
//        std::cout << "------- hadamard -------" << std::endl;
//        std::cout << "--- remining[0]" << std::endl;
//        for (auto&& r : remaining_[0])
//        {
//            std::cout << r << std::endl;
//        }
//        std::cout << "--- remining[1]" << std::endl;
//        for (auto&& r : remaining_[1])
//        {
//            std::cout << r << std::endl;
//        }

        /*
         * determine apply (carry) index list
         */
//        std::cout << std::endl;
//        std::cout << "----- before hadamard -----" << std::endl;
//        std::cout << "--- index list[0]" << std::endl;
//        for (auto&& e : index_list[0])
//        {
//            std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//        }
//        std::cout << "--- carry index list[0]" << std::endl;
//        for (auto&& e : carry_index_list[0])
//        {
//            std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//        }
//
//        std::cout << "--- index list[1]" << std::endl;
//        for (auto&& e : index_list[1])
//        {
//            std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//        }
//        std::cout << "--- carry index list[1]" << std::endl;
//        for (auto&& e : carry_index_list[1])
//        {
//            std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//        }

        std::vector<std::list<int>> tmp_index_list(2);
        std::vector<std::list<int>> tmp_carry_index_list(2);
        for (int i = 0; i < 2; ++i)
        {
            for (auto it = remaining_[i].begin(); it != remaining_[i].end();)
            {
                util::xor_func tmp = (~mask_) & (chr_.phase_exponents()[*it].second);
                if (tmp.none())
                {
                    auto ti = hadamard.in_.find(*it);
                    if (ti != hadamard.in_.end())
                    {
                        tmp_index_list[i].push_back(*it);
                    }
                    else
                    {
                        tmp_carry_index_list[i].push_back(*it);
                    }
                    it = remaining_[i].erase(it);
                }
                else
                {
                    it++;
                }
            }
        }

        index_list[0] = tmp_index_list[0];
        index_list[1] = tmp_index_list[1];
        carry_index_list[0] = tmp_carry_index_list[0];
        carry_index_list[1] = tmp_carry_index_list[1];

        /**
         * Construct sub-circuit
         */
//        std::cout << std::endl;
//        std::cout << "----- after hadamard -----" << std::endl;
//        std::cout << "--- index list[0]" << std::endl;
//        for (auto&& e : index_list[0])
//        {
//            std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//        }
//        std::cout << "--- carry index list[0]" << std::endl;
//        for (auto&& e : carry_index_list[0])
//        {
//            std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//        }

        circuit_.add_gate_list(builder_.Build(index_list[0], carry_index_list[0], wires_, wires_));

//        std::cout << "--- index list[1]" << std::endl;
//        for (auto&& e : index_list[1])
//        {
//            std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//        }
//        std::cout << "--- carry index list[1]" << std::endl;
//        for (auto&& e : carry_index_list[1])
//        {
//            std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//        }

        circuit_.add_gate_list(builder_.Build(index_list[1], carry_index_list[1], wires_, hadamard.input_wires_parity_));
        remaining_[0].splice(remaining_[0].begin(), carry_index_list[0]);
        remaining_[1].splice(remaining_[1].begin(), carry_index_list[1]);
        for (int i = 0; i < chr_.num_qubit(); i++)
        {
            wires_[i] = hadamard.input_wires_parity_[i];
        }

//        std::cout << std::endl;
//        std::cout << "----- after construct -----" << std::endl;
//        std::cout << "--- index list[0]" << std::endl;
//        for (auto&& e : index_list[0])
//        {
//            std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//        }
//        std::cout << "--- carry index list[0]" << std::endl;
//        for (auto&& e : carry_index_list[0])
//        {
//            std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//        }
//        circuit_.add_gate_list(builder_.Build(index_list[0], carry_index_list[0], wires_, wires_));

//        std::cout << "--- index list[1]" << std::endl;
//        for (auto&& e : index_list[1])
//        {
//            std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//        }
//        std::cout << "--- carry index list[1]" << std::endl;
//        for (auto&& e : carry_index_list[1])
//        {
//            std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//        }

        /*
         * Apply Hadamard gate
        */
        circuit_.add_gate("H", chr_.qubit_names()[hadamard.target_]);
        wires_[hadamard.target_].reset();
        wires_[hadamard.target_].set(hadamard.previous_qubit_index_);
        mask_.set(hadamard.previous_qubit_index_);

        /*
         * Check for increases in dimension
         */
        dimension = builder_.CheckDimension(chr_, wires_, dimension);
    }

    /*
     * Construct the final {CNOT, T} subcircuit
     */
    std::vector<std::list<int>> final_index_list(2);
    for (int i = 0; i < 2; ++i)
    {
        for (auto it = remaining_[i].begin(); it != remaining_[i].end(); it++)
        {
            final_index_list[i].push_back(*it);
        }
    }

//    std::cout << std::endl;
//    std::cout << "----- final part -----" << std::endl;
//    std::cout << "-- final index list[0]" << std::endl;
//    for (auto&& e : final_index_list[0])
//    {
//        std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//    }
//    std::cout << "-- carry index list[1]" << std::endl;
//    for (auto&& e : final_index_list[1])
//    {
//        std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//    }

    index_list[0] = std::move(carry_index_list[0]);
    index_list[1] = std::move(carry_index_list[1]);
    std::list<int> none_list;

    circuit_.add_gate_list(builder_.Build(final_index_list[0], none_list, wires_, wires_));
    circuit_.add_gate_list(builder_.Build(final_index_list[1], none_list, wires_, chr_.outputs()));

//    std::cout << std::endl;
//    std::cout << "----- final part -----" << std::endl;
//    std::cout << "-- carry index list[0]" << std::endl;
//    for (auto&& e : carry_index_list[0])
//    {
//        std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//    }
//    std::cout << "-- carry index list[1]" << std::endl;
//    for (auto&& e : carry_index_list[1])
//    {
//        std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//    }
//    index_list[0] = std::move(carry_index_list[0]);
//    index_list[1] = std::move(carry_index_list[1]);
//    std::list<int> none_list;
//
//    std::cout << std::endl;
//    std::cout << "----- final part -----" << std::endl;
//    std::cout << "-- index list[0]" << std::endl;
//    for (auto&& e : index_list[0])
//    {
//        std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//    }
//    std::cout << "-- index list[1]" << std::endl;
//    for (auto&& e : index_list[1])
//    {
//        std::cout << e << ":" << chr_.phase_exponents()[e].second << std::endl;
//    }
//    std::cout << std::endl;
//
//    circuit_.add_gate_list(builder_.Build(index_list[0], none_list, wires_, wires_));
//    circuit_.add_gate_list(builder_.Build(index_list[1], none_list, wires_, chr_.outputs()));

    /*
     * Add the global phase
     */
//    circuit_.add_gate_list(builder_.BuildGlobalPhase(chr_.num_qubit(), global_phase_, chr_.qubit_names()));

    return circuit_;
}

}
