#ifndef T_SCHEDULING_LAYOUT_HPP
#define T_SCHEDULING_LAYOUT_HPP

#include <vector>
#include <memory>

#include "../util/option.hpp"

#include "../circuit/circuit.hpp"

namespace tskd {

enum NodeType
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

class Edge;

class Node
{
private:
    int id_;
    std::string name_;
    int x_;
    int y_;
    NodeType type_;

    std::shared_ptr<Edge> right_edge_;
    std::shared_ptr<Edge> left_edge_;
    std::shared_ptr<Edge> upper_edge_;
    std::shared_ptr<Edge> lower_edge_;
    std::vector<std::shared_ptr<Edge>> edge_list_;

public:
    Node() = default;

    Node(int id,
         const std::string& name,
         int x,
         int y,
         NodeType type)
        : id_(id),
          name_(name),
          x_(x),
          y_(y),
          type_(type),
          right_edge_(nullptr),
          left_edge_(nullptr),
          upper_edge_(nullptr),
          lower_edge_(nullptr) { }

    int id() const
    {
        return id_;
    }

    std::string name() const
    {
        return name_;
    }

    int x() const
    {
        return x_;
    }

    int y() const
    {
        return y_;
    }

    std::vector<std::shared_ptr<Edge>> edge_list() const
    {
        return edge_list_;
    };

    /*
     * dir means that following
     * - 0: right
     * - 1: left
     * - 2: upper
     * - 3: lower
     */
    inline void add_edge(const std::shared_ptr<Edge>& edge,
                         int dir)
    {
        edge_list_.push_back(edge);
        if (dir == 0)
        {
            right_edge_ = edge;
        }
        else if (dir == 1)
        {
            left_edge_ = edge;
        }
        else if (dir == 2)
        {
            upper_edge_ = edge;
        }
        else if (dir == 3)
        {
            lower_edge_ = edge;
        }
    }

    void print() const
    {
        std::cout << "(" << id_ << "," << name_ << "," << x_ << "," << y_ << ")" << std::endl;
    }
};

class Edge
{
private:
    int id_;
    std::shared_ptr<Node> node_a_;
    std::shared_ptr<Node> node_b_;
public:
    Edge() = default;

    Edge(int id,
         const std::shared_ptr<Node>& node_a,
         const std::shared_ptr<Node>& node_b)
        : id_(id),
          node_a_(node_a),
          node_b_(node_b) { }

    int id() const
    {
        return id_;
    }

    std::shared_ptr<Node> node_a() const
    {
        return node_a_;
    }

    std::shared_ptr<Node> node_b() const
    {
        return node_b_;
    }

    void print() const
    {
        std::cout << id_ << ": ";
        std::cout << "(" << node_a_->id() << "," << node_a_->name() << "," << node_a_->x() << "," << node_a_->y() << ")";
        std::cout << " <--> ";
        std::cout << "(" << node_b_->id() << "," << node_b_->name() << "," << node_b_->x() << "," << node_b_->y() << ")";
        std::cout << std::endl;
    }
};

class Layout
{
private:
    util::Option option_;

    Circuit circuit_;

    LayoutType type_;

    int patch_size_;
    int width_;
    int height_;

//    std::vector<Node> node_list_;
    std::vector<std::shared_ptr<Node>> node_list_;
    std::vector<std::shared_ptr<Edge>> horizontal_edge_list_;
    std::vector<std::shared_ptr<Edge>> vertical_edge_list_;
    std::vector<std::shared_ptr<Edge>> edge_list_;
    std::vector<std::vector<std::shared_ptr<Node>>> grid_;

    void init();

    void calculate_width_and_height();

    void make_node_list_and_grid();

    void make_edge_list();

    inline void new_node(const std::string& name,
                         int x,
                         int y,
                         const NodeType& type)
    {
        const int id = static_cast<int>(node_list_.size());
        auto node = std::make_shared<Node>(id, name, x, y, type);
        node_list_.push_back(node);
        grid_[y][x] = node;
    }

    inline void new_edge(std::shared_ptr<Node>& node_a,
                         std::shared_ptr<Node>& node_b,
                         bool horizontal)
    {
        const int id = static_cast<int>(edge_list_.size());
        auto edge = std::make_shared<Edge>(id, node_a, node_b);
        edge_list_.push_back(edge);

        int dir_a = 0;
        int dir_b = 0;
        if (horizontal)
        {
            dir_a = 0;
            dir_b = 1;
            horizontal_edge_list_.push_back(edge);
        }
        else
        {
            dir_a = 3;
            dir_b = 2;
            vertical_edge_list_.push_back(edge);
        }
        node_a->add_edge(edge, dir_a);
        node_b->add_edge(edge, dir_b);
    }

public:
    Layout() = default;

    Layout(const util::Option& option,
           const Circuit& circuit)
        : option_(option),
          circuit_(circuit),
          type_(LayoutType::kuv),
          patch_size_(2),
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

    std::vector<std::shared_ptr<Node>> node_list() const
    {
        return node_list_;
    }

    std::vector<std::shared_ptr<Edge>> horizontal_edge_list() const
    {
        return horizontal_edge_list_;
    }

    std::vector<std::shared_ptr<Edge>> vertical_edge_list() const
    {
        return vertical_edge_list_;
    }

    std::vector<std::shared_ptr<Edge>> edge_list() const
    {
        return edge_list_;
    }

    std::vector<std::vector<std::shared_ptr<Node>>> grid() const
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
//        for (int y = 0; y < height_; y++)
//        {
//            for (int x = 0; x < width_; x++)
//            {
//                std::cout << grid_[y][x]->name();
//            }
//            std::cout << std::endl;
//        }
    }

    void print_edge() const
    {
        for (auto&& edge : edge_list_)
        {
            edge->print();
        }
    }
};

}

#endif //T_SCHEDULING_LAYOUT_HPP
