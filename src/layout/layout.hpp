#ifndef T_SCHEDULING_LAYOUT_HPP
#define T_SCHEDULING_LAYOUT_HPP

#include <vector>

#include "../util/option.hpp"

namespace tskd {

enum BitType
{
    kdata,
    kancilla
};

class Layout
{
public:
    struct LogicalBit
    {
        int id_;
        int x_;
        int y_;
        BitType type_;

        LogicalBit() { }

        LogicalBit(int id,
                   int x,
                   int y,
                   BitType type)
            : id_(id),
              x_(x),
              y_(y),
              type_(type) { }
    };

private:
    util::Option option_;

    int width_;
    int height_;

    std::vector<std::vector<LogicalBit>> grid_;

    void init();

public:
    Layout(const util::Option& option) : option_(option) { }

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
};

}

#endif //T_SCHEDULING_LAYOUT_HPP
