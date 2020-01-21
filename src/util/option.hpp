#ifndef T_SCHEDULING_OPTION_HPP
#define T_SCHEDULING_OPTION_HPP

#include <iostream>
#include <string>

enum SynthesisMethod
{
    ktpar,
    ktskd
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

    int num_buffer_;
    int num_buffer_row_;

    bool change_row_order_;

    SynthesisMethod syn_method_;
    DecompositionType dec_type_;

public:
    Option() = default;

    Option(std::string  input_path_,
           int num_distillation,
           int distillation_step,
           int num_buffer,
           bool change_row_order,
           SynthesisMethod syn_method,
           DecompositionType dec_type)
            : input_path_(std::move(input_path_)),
              num_distillation_(num_distillation),
              distillation_step_(distillation_step),
              num_buffer_(num_buffer),
              change_row_order_(change_row_order),
              syn_method_(syn_method),
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

    int num_buffer() const
    {
        return num_buffer_;
    }

    int num_buffer_row() const
    {
        return num_buffer_row_;
    }

    bool change_row_order() const
    {
        return change_row_order_;
    }

    SynthesisMethod syn_method() const
    {
        return syn_method_;
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

    void set_num_buffer(int num_buffer)
    {
        num_buffer_ = num_buffer;
    }

    void set_num_buffer_row(int num_row)
    {
        num_buffer_row_ = num_row;
    }

    void set_change_row_order(bool change_row_order)
    {
        change_row_order_ = change_row_order;
    }

    void set_syn_method(const SynthesisMethod& syn_method)
    {
        syn_method_ = syn_method;
    }

    void set_dec_type(const DecompositionType& dec_type)
    {
        dec_type_ = dec_type;
    }

    void show()
    {
        std::cout << "# Option list" << std::endl;
        std::cout << "# input file: " << input_path_ << std::endl;
        std::cout << "# output file: " << output_path_ << std::endl;
        std::cout << "# number of distillatoion: " << num_distillation_ << std::endl;
        std::cout << "# required distillation step: " << distillation_step_ << std::endl;
        std::cout << "# number of buffer: " << num_buffer_ << std::endl;
        std::cout << "# synthesis method: ";
        switch (syn_method_)
        {
            case SynthesisMethod::ktpar:
                std::cout << "t-par" << std::endl;
                break;
            case SynthesisMethod::ktskd:
                std::cout << "t-scheduling" << std::endl;
                break;
        }
        std::cout << "# change row order: " << std::boolalpha << change_row_order_ << std::endl;
        std::cout << "# decomposition type: ";
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
