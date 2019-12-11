#ifndef T_SCHEDULING_OPTION_HPP
#define T_SCHEDULING_OPTION_HPP

#include <iostream>
#include <string>

enum PartitionType
{
    kmatroid,
    kgreed
};

enum DecompositionType
{
    kgauss,
    ktskd
};

namespace tskd {
namespace util {

class Option {
private:
    std::string input_path_;
    std::string output_path_;

    bool change_row_order_;

    PartitionType part_type_;
    DecompositionType dec_type_;

public:
    Option() = default;

    Option(std::string  input_path_,
           bool change_row_order,
           PartitionType part_type,
           DecompositionType dec_type)
            : input_path_(std::move(input_path_)),
              change_row_order_(change_row_order),
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

    bool change_row_order() const
    {
        return change_row_order_;
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

    void set_change_row_order(bool change_row_order)
    {
        change_row_order_ = change_row_order;
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
        std::cout << "## change row order: " << change_row_order_ << std::endl;
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
            case DecompositionType::ktskd:
                std::cout << "tskd" << std::endl;
                break;
        }
    }
};

}
}




#endif //T_SCHEDULING_OPTION_HPP
