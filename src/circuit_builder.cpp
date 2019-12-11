#include "circuit_builder.hpp"
#include "gate.hpp"
#include "gaussian_decomposer.hpp"

namespace tskd {

bool CircuitBuilder::Init(const std::vector<util::xor_func>& in,
                          const std::vector<util::xor_func>& out)
{
    bool is_io_different = true;

    bits_ = std::vector<util::xor_func>(qubit_num_);
    preparation_ = std::vector<util::xor_func>(qubit_num_);
    restoration_ = std::vector<util::xor_func>(qubit_num_);

    for (int i = 0; i < qubit_num_; i++)
    {
        is_io_different &= (in[i] == out[i]);
        preparation_[i] = util::xor_func(qubit_num_ + 1, 0);
        restoration_[i] = util::xor_func(qubit_num_ + 1, 0);
        preparation_[i].set(i);
        restoration_[i].set(i);
    }

    return is_io_different;
}

void CircuitBuilder::InitBits(const std::set<int>& phase_exponent_index_set)
{
    std::set<int>::iterator ti;
    int counter = 0;
    for (ti = phase_exponent_index_set.begin(), counter = 0; counter < qubit_num_; counter++)
    {
        if (counter < static_cast<int>(phase_exponent_index_set.size()))
        {
            bits_[counter] = phase_exponent_[*ti].second;
            ti++;
        }
        else
        {
            bits_[counter] = util::xor_func(dimension_ + 1, 0);
        }
    }
}

void CircuitBuilder::Prepare(std::list<Gate>& gate_list,
                             const std::vector<util::xor_func>& in,
                             const int num_partition)
{
    util::ToUpperEchelon(num_partition, dimension_, bits_, &restoration_, std::vector<std::string>());
    util::FixBasis(qubit_num_, dimension_, num_partition, in, bits_, &restoration_, std::vector<std::string>());
    util::Compose(qubit_num_, preparation_, restoration_);

    /*
    std::cout << "---- preparation_" << std::endl;
    for (auto&& e : preparation_)
    {
        std::cout << e << std::endl;
    }
    std::cout << "---- restoration_" << std::endl;
    for (auto&& e : restoration_)
    {
        std::cout << e << std::endl;
    }
     */

    if (option_.dec_type() == DecompositionType::kgauss)
    {
        gate_list.splice(gate_list.end(), GaussianDecomposer()(qubit_num_, 0, preparation_, qubit_names_));
    }
}

void CircuitBuilder::ApplyPhaseGates(std::list<Gate>& gate_list,
                                     const std::set<int>& phase_exponent_index_set)
{
    int index = 0;
    for (auto&& phase_exponent_index : phase_exponent_index_set)
    {
        if (phase_exponent_[phase_exponent_index].first <= 4)
        {
            if (phase_exponent_[phase_exponent_index].first / 4 == 1)
            {
                gate_list.emplace_back("Z", qubit_names_[index]);
            }
            if (phase_exponent_[phase_exponent_index].first / 2 == 1)
            {
                gate_list.emplace_back("P", qubit_names_[index]);
            }
            if (phase_exponent_[phase_exponent_index].first % 2 == 1)
            {
                gate_list.emplace_back("T", qubit_names_[index]);
            }
        }
        else
        {
            if (phase_exponent_[phase_exponent_index].first == 5 || phase_exponent_[phase_exponent_index].first == 6)
            {
                gate_list.emplace_back("P*", qubit_names_[index]);
            }
            if (phase_exponent_[phase_exponent_index].first % 2 == 1)
            {
                gate_list.emplace_back("T*", qubit_names_[index]);
            }
        }
        index++;
    }
}

void CircuitBuilder::UnPrepare()
{
    preparation_ = std::move(restoration_);
    restoration_ = std::vector<util::xor_func>(qubit_num_);
    // re-initialize
    for (int i = 0; i < qubit_num_; i++)
    {
        restoration_[i] = util::xor_func(qubit_num_ + 1, 0);
        restoration_[i].set(i);
    }
}

void CircuitBuilder::PrepareLastPart(std::list<Gate>& gate_list,
                                     const std::vector<util::xor_func>& in,
                                     const std::vector<util::xor_func>& out)
{
    for (int i = 0; i < qubit_num_; i++)
    {
        bits_[i] = out[i];
    }
    util::ToUpperEchelon(qubit_num_, dimension_, bits_, &restoration_, std::vector<std::string>());
    util::FixBasis(qubit_num_, dimension_, qubit_num_, in, bits_, &restoration_, std::vector<std::string>());
    util::Compose(qubit_num_, preparation_, restoration_);

    if (option_.dec_type() == DecompositionType::kgauss)
    {
        gate_list.splice(gate_list.end(), GaussianDecomposer()(qubit_num_, 0, preparation_, qubit_names_));
    }
}

std::list<Gate> CircuitBuilder::Build(const tpar::partitioning& partition,
                                      std::vector<util::xor_func>& in,
                                      const std::vector<util::xor_func>& out)
{
    std::list<Gate> ret;

    if (Init(in, out) && partition.empty())
    {
        return ret;
    }

    /*
    std::cout << "-->> build circuit" << std::endl;
    std::cout << "# partition size: " << static_cast<int>(partition.size()) << std::endl;
    for (auto&& e : partition)
    {
        std::cout << "# --" << std::endl;
        for (auto&& ee : e)
        {
            std::cout << "# " << phase_exponent_[ee].second << std::endl;
        }
    }
    std::cout << "# input " << std::endl;
    for (auto&& e : in)
    {
        std::cout << "# " << e << std::endl;
    }
    std::cout << "# output " << std::endl;
    for (auto&& e : out)
    {
        std::cout << "# " << e << std::endl;
    }
     */

    /*
     * Reduce in to echelon form to decide on a basis
     */
    util::ToUpperEchelon(qubit_num_, dimension_, in, &preparation_, std::vector<std::string>());

    /*
     * For each partition... Compute *it, apply T gates, uncompute
     */
    for (auto&& it : partition)
    {
        InitBits(it);
//        std::set<int>::iterator ti;
//        int counter = 0;
//        for (ti = it.begin(), counter = 0; counter < qubit_num_; counter++)
//        {
//            if (counter < static_cast<int>(it.size()))
//            {
//                bits_[counter] = phase_exponent_[*ti].second;
//                ti++;
//            }
//            else
//            {
//                bits_[counter] = util::xor_func(dimension_ + 1, 0);
//            }
//        }
        /*
         * Prepare the bits
         */
        Prepare(ret, in, static_cast<int>(it.size()));

        /*
         * Apply the T gates
         */
        ApplyPhaseGates(ret, it);

        /*
         * Unprepare the bits
         */
        UnPrepare();
    }

    /*
     * Reduce out to the basis of in
     */
    PrepareLastPart(ret, in, out);

    return ret;
}

std::list<Gate> CircuitBuilder::BuildGlobalPhase(int qubit_num,
                                                 int phase,
                                                 const std::vector<std::string>& qubit_names)
{
    std::list<Gate> acc;
    int qubit = 0;

    if (phase % 2 == 1)
    {
        acc.splice(acc.end(), util::ComposeOM(qubit, qubit_names));
        qubit = (qubit + 1) % qubit_num;
    }
    for (int i = phase / 2; i > 0; i--)
    {
        acc.splice(acc.end(), util::ComposeImaginaryUnit(qubit, qubit_names));
        qubit = (qubit + 1) % qubit_num;
    }

    return acc;
}

}