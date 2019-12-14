#ifndef T_SCHEDULING_OPTION_HPP
#define T_SCHEDULING_OPTION_HPP

#include <iostream>
#include <string>

enum SynthesisMethod
{
    ktpar,
    ktskd
};


enum PartitionType
{
    kmatroid,
    kgreed
};

enum DecompositionType
{
    kgauss,
    kparallel
};

namespace tskd {
namespace util {

class Option {
private:
    std::string input_path_;
    std::string output_path_;

    int num_distillation_;
    int distillation_step_;

    bool change_row_order_;

    SynthesisMethod syn_method_;
    PartitionType part_type_;
    DecompositionType dec_type_;

public:
    Option() = default;

    Option(std::string  input_path_,
           int num_distillation,
           int distillation_step,
           bool change_row_order,
           SynthesisMethod syn_method,
           PartitionType part_type,
           DecompositionType dec_type)
            : input_path_(std::move(input_path_)),
              num_distillation_(num_distillation),
              distillation_step_(distillation_step),
              change_row_order_(change_row_order),
              syn_method_(syn_method),
              part_type_(part_type),
              dec_type_(dec_type)
    {
        output_path_ = input_path_ + "-opt";
    }

    std::string input_path() const
    {
        return input_path_;
    }

    std::string output_path() const
    {
        return output_path_;
    }

    int num_distillation() const
    {
        return num_distillation_;
    }

    int distillation_step() const
    {
        return distillation_step_;
    }

    bool change_row_order() const
    {
        return change_row_order_;
    }

    SynthesisMethod syn_method() const
    {
        return syn_method_;
    }

    PartitionType part_type() const
    {
        return part_type_;
    }

    DecompositionType dec_type() const
    {
        return dec_type_;
    }

    void set_input_path(const std::string& input_path)
    {
        input_path_ = input_path;
        output_path_ = input_path_ + "-opt";
    }

    void set_num_distillation(int num_distillation)
    {
        num_distillation_ = num_distillation;
    }

    void set_distillation_step(int distillation_step)
    {
        distillation_step_ = distillation_step;
    }

    void set_change_row_order(bool change_row_order)
    {
        change_row_order_ = change_row_order;
    }

    void set_syn_method(const SynthesisMethod& syn_method)
    {
        syn_method_ = syn_method;
    }

    void set_part_type(const PartitionType& part_type)
    {
        part_type_ = part_type;
    }

    void set_dec_type(const DecompositionType& dec_type)
    {
        dec_type_ = dec_type;
    }

    void show()
    {
        std::cout << "# Option list" << std::endl;
        std::cout << "## input file: " << input_path_ << std::endl;
        std::cout << "## output file: " << output_path_ << std::endl;
        std::cout << "## number of distillatoion: " << num_distillation_ << std::endl;
        std::cout << "## required distillation step: " << distillation_step_ << std::endl;
        std::cout << "## synthesis method: ";
        switch (syn_method_)
        {
            case SynthesisMethod::ktpar:
                std::cout << "t-par" << std::endl;
                break;
            case SynthesisMethod::ktskd:
                std::cout << "t-scheduling" << std::endl;
                break;
        }
        std::cout << "## change row order: " << std::boolalpha << change_row_order_ << std::endl;
        std::cout << "## partition type: ";
        switch (part_type_)
        {
            case PartitionType::kmatroid:
                std::cout << "matroid" << std::endl;
                break;
            case PartitionType::kgreed:
                std::cout << "greed" << std::endl;
                break;
        }
        std::cout << "## decomposition type: ";
        switch (dec_type_)
        {
            case DecompositionType::kgauss:
                std::cout << "gauss" << std::endl;
                break;
            case DecompositionType::kparallel:
                std::cout << "parallel" << std::endl;
                break;
        }
    }
};

}
}




#endif //T_SCHEDULING_OPTION_HPP