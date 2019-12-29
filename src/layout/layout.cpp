#include "layout.hpp"

namespace tskd {

void Layout::init()
{
    constexpr int patch_size = 2;
    int id = 1;

    // TODO: calculate width and height of gird
    width_ = 10;
    height_ = 10;

    grid_ = std::vector<std::vector<LogicalBit>>(height_, std::vector<LogicalBit>(width_));

    for (int y = 0; y < height_; y++)
    {
        for (int x = 0; x < width_; x++)
        {
            if (x % patch_size == 0 && y % patch_size == 0)
            {
                grid_[y][x] = LogicalBit(id, x, y, BitType::kdata);
            }
            else
            {
                grid_[y][x] = LogicalBit(id, x, y, BitType::kancilla);
            }
        }
    }
}

}