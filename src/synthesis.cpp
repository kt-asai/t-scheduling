#include "tpar/partition.hpp"
#include "tpar/matroid.hpp"
#include "synthesis.hpp"
#include "character.hpp"
#include "circuit_builder.hpp"

namespace tskd {

void Synthesis::init(const Character& chr)
{
    floats_.resize(2);
    frozen_.resize(2);
    remaining_.resize(2);

    /*
     * initialize some stuff
     */
    global_phase_ = 0;
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
    for (int i = 0; i < chr.phase_exponents().size(); i++)
    {
        if (chr.phase_exponents()[i].second == util::xor_func(chr.num_data_qubit() + chr.num_hadamard() + 1, 0))
        {
            global_phase_ = chr.phase_exponents()[i].first;
        }
        else if (chr.phase_exponents()[i].first % 2 == 1)
        {
            remaining_[0].push_back(i);
        }
        else if (chr.phase_exponents()[i].first != 0)
        {
            remaining_[1].push_back(i);
        }
    }
}

Circuit Synthesis::Execute()
{
    std::cout << "running..." << std::endl;

    int dimension = chr_.num_data_qubit();
    int applied = 0;
    util::IndependentOracle oracle(chr_.num_qubit(), dimension, chr_.num_data_qubit() + chr_.num_hadamard());
    CircuitBuilder builder(chr_.num_qubit(), dimension, chr_.qubit_names(), chr_.phase_exponents());

    /*
    std::cout << "- wires" << std::endl;
    for (auto&& w : wires_)
    {
        std::cout << w << std::endl;
    }
    std::cout << "- mask" << std::endl;
    std::cout << mask_ << std::endl;
    std::cout << "- remaining[0]" << std::endl;
    for (auto&& r : remaining_[0])
    {
        std::cout << r << ":" << chr_.phase_exponents()[r].second << std::endl;
    }
    std::cout << "- remaining[1]" << std::endl;
    for (auto&& r : remaining_[1])
    {
        std::cout << r << ":" << chr_.phase_exponents()[r].second << std::endl;
    }
    */

    /*
     * create an initial partition
     */
    for (int j = 0; j < 2; j++)
    {
        for (auto it = remaining_[j].begin(); it != remaining_[j].end();)
        {
            util::xor_func tmp = (~mask_) & (chr_.phase_exponents()[*it].second);
            if (tmp.none())
            {
                tpar::add_to_partition(floats_[j], *it, chr_.phase_exponents(), oracle);
                it = remaining_[j].erase(it);
            }
            else
            {
                it++;
            }
        }
    }

    /*
    std::cout << "- floats[0]" << std::endl;
    for (auto f : floats[0])
    {
        std::cout << "--" << std::endl;
        for (auto ff : f)
        {
            std::cout << ff << ":" << phase_exponents_[ff].second << std::endl;
        }
    }
    std::cout << "- floats[1]" << std::endl;
    for (auto f : floats[1])
    {
        std::cout << "--" << std::endl;
        for (auto ff : f)
        {
            std::cout << ff << ":" << phase_exponents_[ff].second << std::endl;
        }
    }

    if (disp_log) cerr << "  " << phase_expts.size() - (remaining[0].size() + remaining[1].size())
                       << "/" << phase_expts.size() << " phase rotations partitioned\n" << flush;
   */

    /**
     * Synthesize the circuit by applying H-gate to greed
     * 1. freeze partitions that are not disjoint from the hadamard input
     * 2. construct CNOT+T circuit
     * 3. apply the hadamard gate
     * 4. add new functions to the partition
     */
    for (auto&& hadamard : chr_.hadamards())
    {
        /*
        std::cout << h_count << ".番目--" << std::endl;
        std::cout << hadamard.target_ << std::endl;
        std::cout << hadamard.previous_qubit_index_ << std::endl;
        std::cout << "- input parity" << std::endl;
        for (auto&& e : hadamard.input_wires_parity_)
        {
            std::cout << e << std::endl;
        }
        std::cout << "- in" << std::endl;
        for (auto&& e : hadamard.in_)
        {
            std::cout << e << std::endl;
        }
        */

        /*
         * determine frozen partitions
         */
        for (int i = 0; i < 2; i++)
        {
            frozen_[i] = tpar::freeze_partitions(floats_[i], hadamard.in_);
            applied += tpar::num_elts(frozen_[i]);
        }

        /*
         * Construct {CNOT, T} subcircuit for the frozen partitions
         */
//        ret.circ.splice(ret.circ.end(),
//                        construct_circuit(phase_expts, frozen[0], wires, wires, n + m, n + h, names));
//        ret.circ.splice(ret.circ.end(),
//                        construct_circuit(phase_expts, frozen[1], wires, it->wires, n + m, n + h, names));
        for (int i = 0; i < chr_.num_qubit(); i++)
        {
            wires_[i] = hadamard.input_wires_parity_[i];
        }

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
        const int updated_dimension = util::ComputeRank(chr_.num_qubit(), chr_.num_data_qubit() + chr_.num_hadamard(), wires_);
        if (updated_dimension > dimension)
        {
            dimension = updated_dimension;
            oracle.set_dim(dimension);
            tpar::repartition(floats_[0], chr_.phase_exponents(), oracle);
            tpar::repartition(floats_[1], chr_.phase_exponents(), oracle);
        }

        /*
         * Add new functions to the partition
         */
        for (int i = 0; i < 2; i++)
        {
            for (auto it = remaining_[i].begin(); it != remaining_[i].end();)
            {
                util::xor_func tmp = (~mask_) & (chr_.phase_exponents()[*it].second);
                if (tmp.none())
                {
                    tpar::add_to_partition(floats_[i], *it, chr_.phase_exponents(), oracle);
                    it = remaining_[i].erase(it);
                }
                else
                {
                    it++;
                }
            }
        }
//        if (disp_log) cerr << "    " << phase_expts.size() - (remaining[0].size() + remaining[1].size())
//                           << "/" << phase_expts.size() << " phase rotations partitioned\n" << flush;
    }

//    applied += num_elts(floats[0]) + num_elts(floats[1]);
//    // Construct the final {CNOT, T} subcircuit
//    ret.circ.splice(ret.circ.end(),
//                    construct_circuit(phase_expts, floats[0], wires, wires, n + m, n + h, names));
//    ret.circ.splice(ret.circ.end(),
//                    construct_circuit(phase_expts, floats[1], wires, outputs, n + m, n + h, names));
//    if (disp_log) cerr << "  " << applied << "/" << phase_expts.size() << " phase rotations applied\n" << flush;

//    // Add the global phase
//    ret.circ.splice(ret.circ.end(), global_phase_synth(n + m, global_phase_, names));

    return circuit_;
}

}