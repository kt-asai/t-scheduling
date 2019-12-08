#include "tpar/partition.hpp"
#include "tpar/matroid.hpp"
#include "synthesis.hpp"
#include "character.hpp"

namespace tskd {

Circuit Synthesis::Execute()
{
    std::cout << "running..." << std::endl;

    Circuit circuit;
    util::IndependentOracle oracle(num_qubit_, dimension_, (num_qubit_ - num_ancilla_) + num_hadamard_);

    auto floats = std::vector<tpar::partitioning>(2);
    auto frozen = std::vector<tpar::partitioning>(2);

    util::xor_func mask((num_qubit_ - num_ancilla_) + num_hadamard_ + 1, 0);
    std::vector<util::xor_func> wires;
    std::vector<std::list<int>> remaining(2);
    int tmp, h_count = 1, applied = 0, j;
    std::list<std::pair<std::string, std::list<std::string> > > circ;
    int global_phase = 0;

    /*
     * initialize some stuff
     */
    mask.set((num_qubit_ - num_ancilla_) + num_hadamard_);
    for (int i = 0, j = 0; i < num_qubit_; i++)
    {
        circuit.add_qubit(qubit_names_[i]);
        circuit.set_ancilla(qubit_names_[i], ancilla_list_[i]);
        wires.emplace_back((num_qubit_ - num_ancilla_) + num_hadamard_ + 1, 0);
        if (!ancilla_list_[i])
        {
            wires[i].set(j);
            mask.set(j);
            j++;
        }
    }

    /*
     * initialize the remaining list
     */
    for (int i = 0; i < phase_exponents_.size(); i++)
    {
        if (phase_exponents_[i].second == util::xor_func((num_qubit_ - num_ancilla_) + num_hadamard_ + 1, 0))
        {
            global_phase = phase_exponents_[i].first;
        }
        else if (phase_exponents_[i].first % 2 == 1)
        {
            remaining[0].push_back(i);
        }
        else if (phase_exponents_[i].first != 0)
        {
            remaining[1].push_back(i);
        }
    }

    /*
    std::cout << "- wires" << std::endl;
    for (auto&& w : wires)
    {
        std::cout << w << std::endl;
    }
    std::cout << "- mask" << std::endl;
    std::cout << mask << std::endl;
    std::cout << "- remaining[0]" << std::endl;
    for (auto&& r : remaining[0])
    {
        std::cout << r << ":" << phase_exponents_[r].second << std::endl;
    }
    std::cout << "- remaining[1]" << std::endl;
    for (auto&& r : remaining[1])
    {
        std::cout << r << ":" << phase_exponents_[r].second << std::endl;
    }
     */

    /*
     * create an initial partition
     */
    for (j = 0; j < 2; j++)
    {
        for (auto it = remaining[j].begin(); it != remaining[j].end();)
        {
            util::xor_func tmp = (~mask) & (phase_exponents_[*it].second);
            if (tmp.none())
            {
                tpar::add_to_partition(floats[j], *it, phase_exponents_, oracle);
                it = remaining[j].erase(it);
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
     * 1. freeze partitions that are not disjoint from the hadamard input
     * 2. construct CNOT+T circuit
     * 3. apply the hadamard gate
     * 4. add new functions to the partition
     */
    for (auto it = hadamards_.begin(); it != hadamards_.end(); it++, h_count++)
    {
        /*
        std::cout << h_count << ".番目--" << std::endl;
        std::cout << it->target_ << std::endl;
        std::cout << it->previous_qubit_index_ << std::endl;
        std::cout << "- input parity" << std::endl;
        for (auto&& e : it->input_wires_parity_)
        {
            std::cout << e << std::endl;
        }
        std::cout << "- in" << std::endl;
        for (auto&& e : it->in_)
        {
            std::cout << e << std::endl;
        }
        if (disp_log) cerr << "  Hadamard " << h_count << "/" << hadamards.size() << "\n" << flush;
         */

        /*
         * determine frozen partitions
         */
        for (j = 0; j < 2; j++)
        {
            frozen[j] = tpar::freeze_partitions(floats[j], it->in_);
            applied += tpar::num_elts(frozen[j]);
        }

        /*
         * Construct {CNOT, T} subcircuit for the frozen partitions
         */
//        ret.circ.splice(ret.circ.end(),
//                        construct_circuit(phase_expts, frozen[0], wires, wires, n + m, n + h, names));
//        ret.circ.splice(ret.circ.end(),
//                        construct_circuit(phase_expts, frozen[1], wires, it->wires, n + m, n + h, names));
        for (int i = 0; i < num_qubit_; i++)
        {
            wires[i] = it->input_wires_parity_[i];
        }

        /*
         * Apply Hadamard gate
         */
        circuit.add_gate("H", qubit_names_[it->target_]);
        wires[it->target_].reset();
        wires[it->target_].set(it->previous_qubit_index_);
        mask.set(it->previous_qubit_index_);

        /*
         * Check for increases in dimension
         */
        tmp = util::ComputeRank(num_qubit_, (num_qubit_ - num_ancilla_) + num_hadamard_, wires);
        if (tmp > dimension_)
        {
            dimension_ = tmp;
            oracle.set_dim(dimension_);
            tpar::repartition(floats[0], phase_exponents_, oracle);
            tpar::repartition(floats[1], phase_exponents_, oracle);
        }

        /*
         * Add new functions to the partition
         */
        for (j = 0; j < 2; j++)
        {
            for (auto it = remaining[j].begin(); it != remaining[j].end();)
            {
                util::xor_func tmp = (~mask) & (phase_exponents_[*it].second);
                if (tmp.none())
                {
                    tpar::add_to_partition(floats[j], *it, phase_exponents_, oracle);
                    it = remaining[j].erase(it);
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
//    ret.circ.splice(ret.circ.end(), global_phase_synth(n + m, global_phase, names));

    return circuit;
}

}