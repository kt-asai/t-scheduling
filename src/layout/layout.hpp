#ifndef T_SCHEDULING_LAYOUT_HPP
#define T_SCHEDULING_LAYOUT_HPP

#include <vector>

#include "../util/option.hpp"

#include "../circuit/circuit.hpp"

namespace tskd {

enum BitType
{
    kdata,
    kancilla
};

enum LayoutType
{
    kuv,  // up and vertical
    kuh,  // up and horizontal
    kudv, // up, down and vertical
    kudh, // up, down and horizontal
};

class Layout
{
public:
    struct LogicalBit
    {
        std::string name_;
        int x_;
        int y_;
        BitType type_;

        LogicalBit() { }

        LogicalBit(const std::string& name,
                   int x,
                   int y,
                   BitType type)
            : name_(name),
              x_(x),
              y_(y),
              type_(type) { }
    };

private:
    util::Option option_;

    Circuit circuit_;

    LayoutType type_;

    int width_;
    int height_;

    std::vector<std::vector<LogicalBit>> grid_;

    void init();

public:
    Layout() = default;

    Layout(const util::Option& option,
           const Circuit& circuit)
        : option_(option),
          circuit_(circuit),
          type_(LayoutType::kuv),
          width_(0),
          height_(0)
    {
        init();
    }

    LayoutType type() const
    {
        return type_;
    }

    int width() const
    {
        return width_;
    }

    int height() const
    {
        return height_;
    }

    std::vector<std::vector<LogicalBit>> grid() const
    {
        return grid_;
    }

    void print() const
    {
        std::cout << "# Layout info" << std::endl;
        std::cout << "# type:  ";
        switch (type_)
        {
            case LayoutType::kuv:
                std::cout << "up vertical" << std::endl;
                break;
            case LayoutType::kuh:
                std::cout << "up horizontal" << std::endl;
                break;
            case LayoutType::kudv:
                std::cout << "up, down vertical" << std::endl;
                break;
            case LayoutType::kudh:
                std::cout << "up, down horizontal" << std::endl;
                break;
        }
        std::cout << "# width: " << width_ << std::endl;
        std::cout << "# height:" << height_ << std::endl;
        for (int y = 0; y < height_; y++)
        {
            for (int x = 0; x < width_; x++)
            {
                std::cout << grid_[y][x].name_;
            }
            std::cout << std::endl;
        }
    }
};

}

#endif //T_SCHEDULING_LAYOUT_HPP
