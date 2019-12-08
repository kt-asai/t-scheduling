#include "circuit_builder.hpp"

namespace tskd {

std::vector<Gate> CircuitBuilder::build(const tpar::partitioning& partition,
                                        std::vector<util::xor_func>& in,
                                        const std::vector<util::xor_func>& out)
{
    std::vector<Gate> circuit;

//    gatelist ret, tmp, rev;
    auto bits = std::vector<util::xor_func>(qubit_num_);
    auto pre = std::vector<util::xor_func>(qubit_num_);
    auto post = std::vector<util::xor_func>(qubit_num_);
    bool flg = true;

    for (int i = 0; i < qubit_num_; i++)
    {
        flg &= (in[i] == out[i]);
    }
    for (int i = 0; i < qubit_num_; i++)
    {
        flg &= (in[i] == out[i]);
        pre[i] = util::xor_func(qubit_num_ + 1, 0);
        post[i] = util::xor_func(qubit_num_ + 1, 0);
        pre[i].set(i);
        post[i].set(i);
    }
    if (flg && (partition.empty()))
    {
        return circuit;
    }

    /*
     * Reduce in to echelon form to decide on a basis
     */
    util::ToUpperEchelon(qubit_num_, dimension_, in, &pre, std::vector<std::string>());

    /*
     * For each partition... Compute *it, apply T gates, uncompute
     */
//    for (partitioning::const_iterator it = part.begin(); it != part.end(); it++)
//    {
//        for (ti = it->begin(), i = 0; i < num; i++)
//        {
//            if (i < it->size())
//            {
//                bits[i] = phase[*ti].second;
//                ti++;
//            }
//            else
//            {
//                bits[i] = xor_func(dim + 1, 0);
//            }
//        }
//
//        // prepare the bits
//        to_upper_echelon(it->size(), dim, bits, &post, vector<string>());
//        fix_basis(num, dim, it->size(), in, bits, &post, vector<string>());
//
//        if (synth_method == GAUSS)
//        {
//            ret.splice(ret.end(), gauss_CNOT_synth(num, 0, pre, names));
//        }
//        else if (synth_method == PMH)
//        {
//            ret.splice(ret.end(), CNOT_synth(num, pre, names));
//        }
//
//        // apply the T gates
//        list <string> tmp_lst;
//        for (ti = it->begin(), i = 0; ti != it->end(); ti++, i++)
//        {
//            tmp_lst.clear();
//            tmp_lst.push_back(names[i]);
//            if (phase[*ti].first <= 4)
//            {
//                if (phase[*ti].first / 4 == 1) ret.push_back(make_pair("Z", tmp_lst));
//                if (phase[*ti].first / 2 == 1) ret.push_back(make_pair("P", tmp_lst));
//                if (phase[*ti].first % 2 == 1) ret.push_back(make_pair("T", tmp_lst));
//            }
//            else
//            {
//                if (phase[*ti].first == 5 || phase[*ti].first == 6) ret.push_back(make_pair("P*", tmp_lst));
//                if (phase[*ti].first % 2 == 1) ret.push_back(make_pair("T*", tmp_lst));
//            }
//        }
//
//        // unprepare the bits
//        pre = std::move(post);
//        post = vector<xor_func>(num);
//        // re-initialize
//        for (i = 0; i < num; i++)
//        {
//            post[i] = xor_func(num + 1, 0);
//            post[i].set(i);
//        }
//    }
//
//    // Reduce out to the basis of in
//    for (i = 0; i < num; i++)
//    {
//        bits[i] = out[i];
//    }
//    to_upper_echelon(num, dim, bits, &post, vector<string>());
//    fix_basis(num, dim, num, in, bits, &post, vector<string>());
//    compose(num, pre, post);
//    if (synth_method == GAUSS) ret.splice(ret.end(), gauss_CNOT_synth(num, 0, pre, names));
//    else if (synth_method == PMH) ret.splice(ret.end(), CNOT_synth(num, pre, names));

    return circuit;
}

}