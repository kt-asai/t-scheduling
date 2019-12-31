#ifndef T_SCHEDULING_PARALLELIZATION_ORACLE_HPP
#define T_SCHEDULING_PARALLELIZATION_ORACLE_HPP

#include <list>

#include "z3++.h"

#include "../circuit/gate.hpp"

#include "../layout/layout.hpp"

namespace tskd {

class ParallelizationOracle
{
private:
    Layout layout_;

    int num_gate_;

    std::unordered_map<std::string, int> control_qubit_map_;    // <control qubit name, id of gate>
    std::unordered_map<std::string, int> target_qubit_map_;     // <target  qubit name, id of gate>
    std::unordered_map<std::string, int> related_qubit_id_map_; // <related qubit name, id of gate>

    int var_count_;

    void init(const std::list<Gate>& gate_list);

    void make_vars(z3::context& context,
                   z3::expr_vector& edge_expr_array,
                   z3::expr_vector& node_expr_array);

    void make_base_constraint(z3::solver& solver,
                              z3::expr_vector& edge_expr_array,
                              z3::expr_vector& node_expr_array);

    void make_range_constraint(z3::solver& solver,
                               z3::expr_vector& edge_expr_array,
                               z3::expr_vector& node_expr_array);

    void make_edge_constraints(z3::solver& solver,
                               z3::expr_vector& edge_expr_array);

    void make_adj_nodes_constraint(z3::solver& solver,
                                   z3::expr_vector& edge_expr_array,
                                   z3::expr_vector& node_expr_array);

    void make_remove_edge_constraint(z3::solver& solver,
                                     z3::expr_vector& edge_expr_array);

    inline int new_variable()
    {
        var_count_ += 1;

        return var_count_;
    }

    void verification(z3::solver& s);

public:
    ParallelizationOracle() = default;

    ParallelizationOracle(const Layout& layout)
        : layout_(layout),
          var_count_(0) { }

    bool check(std::list<Gate>& gate_list,
               const Gate& new_gate);

    bool check(const std::list<Gate>& gate_list);
};

}

#endif //T_SCHEDULING_PARALLELIZATION_ORACLE_HPP
