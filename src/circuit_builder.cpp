#include "circuit_builder.hpp"
#include "gate.hpp"

namespace tskd {

std::list<Gate> CircuitBuilder::build(const tpar::partitioning& partition,
                                      std::vector<util::xor_func>& in,
                                      const std::vector<util::xor_func>& out)
{
    std::list<Gate> ret;

    auto bits = std::vector<util::xor_func>(qubit_num_);
    auto pre = std::vector<util::xor_func>(qubit_num_);
    auto post = std::vector<util::xor_func>(qubit_num_);
    bool is_io_different = true;
    int i = 0;

    for (int i = 0; i < qubit_num_; i++)
    {
        is_io_different &= (in[i] == out[i]);
        pre[i] = util::xor_func(qubit_num_ + 1, 0);
        post[i] = util::xor_func(qubit_num_ + 1, 0);
        pre[i].set(i);
        post[i].set(i);
    }
    if (is_io_different && (partition.empty()))
    {
        return ret;
    }

    /*
     * Reduce in to echelon form to decide on a basis
     */
    util::ToUpperEchelon(qubit_num_, dimension_, in, &pre, std::vector<std::string>());

    /*
     * For each partition... Compute *it, apply T gates, uncompute
     */
    for (tpar::partitioning::const_iterator it = partition.begin(); it != partition.end(); it++)
    {
        std::set<int>::iterator ti;
        int counter = 0;
        for (ti = it->begin(), counter = 0; counter < qubit_num_; counter++)
        {
            if (counter < it->size())
            {
                bits[counter] = phase_exponent_[*ti].second;
                ti++;
            }
            else
            {
                bits[counter] = util::xor_func(dimension_ + 1, 0);
            }
        }

        /*
         * Prepare the bits
         */
        util::ToUpperEchelon(it->size(), dimension_, bits, &post, std::vector<std::string>());
        util::FixBasis(qubit_num_, dimension_, it->size(), in, bits, &post, std::vector<std::string>());
        util::Compose(qubit_num_, pre, post);

        std::cout << "---- pre" << std::endl;
        for (auto&& e : pre)
        {
            std::cout << e << std::endl;
        }
        std::cout << "---- post" << std::endl;
        for (auto&& e : post)
        {
            std::cout << e << std::endl;
        }

//        ret.splice(ret.end(), gauss_CNOT_synth(num, 0, pre, names));

//        if (synth_method == GAUSS)
//        {
//            ret.splice(ret.end(), gauss_CNOT_synth(num, 0, pre, names));
//        }
//        else if (synth_method == PMH)
//        {
//            ret.splice(ret.end(), CNOT_synth(num, pre, names));
//        }

        /*
         * Apply the T gates
         */
        for (ti = it->begin(), i = 0; ti != it->end(); ti++, i++)
        {
            if (phase_exponent_[*ti].first <= 4)
            {
                if (phase_exponent_[*ti].first / 4 == 1)
                {
                    ret.emplace_back("Z", qubit_names_[i]);
                }
                if (phase_exponent_[*ti].first / 2 == 1)
                {
                    ret.emplace_back("P", qubit_names_[i]);
                }
                if (phase_exponent_[*ti].first % 2 == 1)
                {
                    ret.emplace_back("T", qubit_names_[i]);
                }
            }
            else
            {
                if (phase_exponent_[*ti].first == 5 || phase_exponent_[*ti].first == 6)
                {
                    ret.emplace_back("P*", qubit_names_[i]);
                }
                if (phase_exponent_[*ti].first % 2 == 1)
                {
                    ret.emplace_back("T*", qubit_names_[i]);
                }
            }
        }

        // unprepare the bits
        pre = std::move(post);
        post = std::vector<util::xor_func>(qubit_num_);
        // re-initialize
        for (i = 0; i < qubit_num_; i++)
        {
            post[i] = util::xor_func(qubit_num_ + 1, 0);
            post[i].set(i);
        }
    }

    // Reduce out to the basis of in
    for (i = 0; i < qubit_num_; i++)
    {
        bits[i] = out[i];
    }
    util::ToUpperEchelon(qubit_num_, dimension_, bits, &post, std::vector<std::string>());
    util::FixBasis(qubit_num_, dimension_, qubit_num_, in, bits, &post, std::vector<std::string>());
    util::Compose(qubit_num_, pre, post);
//    ret.splice(ret.end(), gauss_CNOT_synth(num, 0, pre, names));
;
    return ret;
}

}