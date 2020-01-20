#include "layout.hpp"

namespace tskd {

void Layout::init()
{
    /*
     * determine data block size
     */
    calculate_width_and_height();

    /*
     * construct node list and grid
     */
    make_node_list_and_grid();

    /*
     * construct edge list
     */
    make_edge_list();

}

void Layout::calculate_width_and_height()
{
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
    for (int i = 0; i < 1; i++)
    {
        const int rotation = i % 2;
        const int width = i < 2 ? ds_length[rotation] * option_.num_distillation()
                                : ds_length[rotation] * ((option_.num_distillation() + 1) / 2);
        const int horizontal_qubit_num = width / patch_size_;
        const int height = (((circuit_.qubit_num() - 1) / horizontal_qubit_num) + 1) * patch_size_;
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
}

void Layout::make_node_list_and_grid()
{
    grid_ = std::vector<std::vector<std::shared_ptr<Node>>>(height_, std::vector<std::shared_ptr<Node>>(width_));
    int index = 0;
    for (int y = 0; y < height_; y++)
    {
        for (int x = 0; x < width_; x++)
        {
            if (x % patch_size_ == 0 && y % patch_size_ == 0)
            {
                if (index < circuit_.qubit_num())
                {
                    // grid_[y][x] =
                    new_node(circuit_.qubit_names()[index], x, y, NodeType::kdata);
                    index++;
                }
                else
                {
                    // grid_[y][x] =
                    new_node("A", x, y, NodeType::kancilla);
                }
            }
            else
            {
                // grid_[y][x] =
                new_node("A", x, y, NodeType::kancilla);
            }
        }
    }
}

void Layout::make_edge_list()
{
    /*
     * make horizontal edge
     */
    for (int y = 0; y < height_; y++)
    {
        for (int x = 0; x < width_ - 1; x++)
        {
            new_edge(grid_[y][x], grid_[y][x + 1], true);
        }
    }

    /*
     * make vertical edge
     */
    for (int y = 0; y < height_ - 1; y++)
    {
        for (int x = 0; x < width_; x++)
        {
            new_edge(grid_[y][x], grid_[y + 1][x], false);
        }
    }
}

}