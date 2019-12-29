#include "layout.hpp"

namespace tskd {

void Layout::init()
{
    // patch size 2x2
    constexpr int patch_size = 2;

    // distillation circuit size
    constexpr int ds_width = 4;
    constexpr int ds_height = 8;
    const std::vector<int> ds_length = {ds_width, ds_height};
    const std::vector<LayoutType> types = {LayoutType::kuv, LayoutType::kuh, LayoutType::kudv, LayoutType::kudh};

    /*
     * 0.up and vertical
     * 1.up and horizontal
     * 2.up, down vertical
     * 3.up, down and horizontal
     */
    double ratio = 100000;
    for (int i = 0; i < 4; i++)
    {
        const int rotation = i % 2;
        const int width = i < 2 ? ds_length[rotation] * option_.num_distillation()
                                : ds_length[rotation] * ((option_.num_distillation() + 1) / 2);
        const int horizontal_qubit_num = width / patch_size;
        const int height = (((circuit_.qubit_num() - 1) / horizontal_qubit_num) + 1) * patch_size;
        const std::pair<int, int> size = std::minmax(width, height);
        const double new_ration = size.second / size.first;
        if (ratio > new_ration)
        {
            ratio = new_ration;
            type_ = types[i];
            width_ = width;
            height_ = height;
        }
    }

    /*
     * construct grid
     */
    grid_ = std::vector<std::vector<LogicalBit>>(height_, std::vector<LogicalBit>(width_));
    int index = 0;
    for (int y = 0; y < height_; y++)
    {
        for (int x = 0; x < width_; x++)
        {
            if (x % patch_size == 0 && y % patch_size == 0)
            {
                if (index < circuit_.qubit_num())
                {
                    grid_[y][x] = LogicalBit(circuit_.qubit_names()[index], x, y, BitType::kdata);
                    index++;
                }
                else
                {
                    grid_[y][x] = LogicalBit("A", x, y, BitType::kancilla);
                }
            }
            else
            {
                grid_[y][x] = LogicalBit("A", x, y, BitType::kancilla);
            }
        }
    }
}

}