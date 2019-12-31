#include "parallelization_oracle.hpp"

namespace tskd {

void ParallelizationOracle::init(const std::list<Gate>& gate_list)
{
    int id = 1;
    for (auto&& gate : gate_list)
    {
        // store control qubit name
        for (auto&& name : gate.control_list())
        {
            control_qubit_map_.emplace(name, id);
            related_qubit_id_map_.emplace(name, id);
        }

        // store target qubit name
        for (auto&& name : gate.target_list())
        {
            target_qubit_map_.emplace(name, id);
            related_qubit_id_map_.emplace(name, id);
        }
        id++;
    }
}

void ParallelizationOracle::make_vars(z3::context& context,
                                      z3::expr_vector& edge_expr_array,
                                      z3::expr_vector& node_expr_array)
{
    // make edge variables
    for (size_t i = 0; i < layout_.edge_list().size(); i++)
    {
        edge_expr_array.push_back(context.int_const(std::to_string(new_variable()).c_str()));
    }

    // make node id variables
    for (size_t i = 0; i < layout_.node_list().size(); i++)
    {
        node_expr_array.push_back(context.int_const(std::to_string(new_variable()).c_str()));
    }
}

void ParallelizationOracle::make_base_constraint(z3::solver& solver,
                                                 z3::expr_vector& edge_expr_array,
                                                 z3::expr_vector& node_expr_array)
{
    make_range_constraint(solver, edge_expr_array, node_expr_array);

    make_edge_constraints(solver, edge_expr_array);

    make_adj_nodes_constraint(solver, edge_expr_array, node_expr_array);

    make_remove_edge_constraint(solver, edge_expr_array);
}

void ParallelizationOracle::make_range_constraint(z3::solver& solver,
                                                    z3::expr_vector& edge_expr_array,
                                                    z3::expr_vector& node_expr_array)
{
    for (auto&& node : layout_.node_list())
    {
        if (node->name() == "A")
        {
            z3::expr e =  node_expr_array[node->id()];
            solver.add(1 <= e && e <= num_gate_);
        }
        else
        {
            if (related_qubit_id_map_.count(node->name()))
            {
                z3::expr e =  node_expr_array[node->id()];
                solver.add(e == related_qubit_id_map_[node->name()]);
            }
            else
            {
                z3::expr e =  node_expr_array[node->id()];
                solver.add(e == 0);
            }
        }
    }

    for (auto&& edge : layout_.edge_list())
    {
        z3::expr e =  edge_expr_array[edge->id()];
        solver.add(0 <= e && e <= 1);
    }
}

void ParallelizationOracle::make_edge_constraints(z3::solver& solver,
                                                  z3::expr_vector& edge_expr_array)
{
    for (auto&& node : layout_.node_list())
    {
        const int num_connected_edge = static_cast<int>(node->edge_list().size());
        if (node->name() != "A")
        {
            if (related_qubit_id_map_.count(node->name()))
            {
                if (num_connected_edge == 4)
                {
                    z3::expr e1 =  edge_expr_array[node->edge_list()[0]->id()];
                    z3::expr e2 =  edge_expr_array[node->edge_list()[1]->id()];
                    z3::expr e3 =  edge_expr_array[node->edge_list()[2]->id()];
                    z3::expr e4 =  edge_expr_array[node->edge_list()[3]->id()];

                    solver.add((e1 + e2 + e3 + e4) == 1);
                }
                if (num_connected_edge == 3)
                {
                    z3::expr e1 =  edge_expr_array[node->edge_list()[0]->id()];
                    z3::expr e2 =  edge_expr_array[node->edge_list()[1]->id()];
                    z3::expr e3 =  edge_expr_array[node->edge_list()[2]->id()];

                    solver.add((e1 + e2 + e3) == 1);
                }
                if (num_connected_edge == 2)
                {
                    z3::expr e1 =  edge_expr_array[node->edge_list()[0]->id()];
                    z3::expr e2 =  edge_expr_array[node->edge_list()[1]->id()];

                    solver.add((e1 + e2) == 1);
                }
            }
            else
            {
                if (num_connected_edge == 4)
                {
                    z3::expr e1 =  edge_expr_array[node->edge_list()[0]->id()];
                    z3::expr e2 =  edge_expr_array[node->edge_list()[1]->id()];
                    z3::expr e3 =  edge_expr_array[node->edge_list()[2]->id()];
                    z3::expr e4 =  edge_expr_array[node->edge_list()[3]->id()];

                    solver.add((e1 + e2 + e3 + e4) == 0);
                }
                if (num_connected_edge == 3)
                {
                    z3::expr e1 =  edge_expr_array[node->edge_list()[0]->id()];
                    z3::expr e2 =  edge_expr_array[node->edge_list()[1]->id()];
                    z3::expr e3 =  edge_expr_array[node->edge_list()[2]->id()];

                    solver.add((e1 + e2 + e3) == 0);
                }
                if (num_connected_edge == 2)
                {
                    z3::expr e1 =  edge_expr_array[node->edge_list()[0]->id()];
                    z3::expr e2 =  edge_expr_array[node->edge_list()[1]->id()];

                    solver.add((e1 + e2) == 0);
                }
            }
        }
        else
        {
            if (num_connected_edge == 4)
            {
                z3::expr e1 =  edge_expr_array[node->edge_list()[0]->id()];
                z3::expr e2 =  edge_expr_array[node->edge_list()[1]->id()];
                z3::expr e3 =  edge_expr_array[node->edge_list()[2]->id()];
                z3::expr e4 =  edge_expr_array[node->edge_list()[3]->id()];

                solver.add((e1 + e2 + e3 + e4) == 0
                ||  (e1 + e2 + e3 + e4) == 2
                || (e1 + e2 + e3 + e4) == 3
                || (e1 + e2 + e3 + e4) == 4);
            }
            if (num_connected_edge == 3)
            {
                z3::expr e1 =  edge_expr_array[node->edge_list()[0]->id()];
                z3::expr e2 =  edge_expr_array[node->edge_list()[1]->id()];
                z3::expr e3 =  edge_expr_array[node->edge_list()[2]->id()];

                solver.add((e1 + e2 + e3) == 0 ||  (e1 + e2 + e3) == 2 || (e1 + e2 + e3) == 3);
            }
            if (num_connected_edge == 2)
            {
                z3::expr e1 =  edge_expr_array[node->edge_list()[0]->id()];
                z3::expr e2 =  edge_expr_array[node->edge_list()[1]->id()];

                solver.add((e1 + e2) == 0 ||  (e1 + e2) == 2);
            }
        }
    }
}

void ParallelizationOracle::make_adj_nodes_constraint(z3::solver& solver,
                                                      z3::expr_vector& edge_expr_array,
                                                      z3::expr_vector& node_expr_array)
{
    for (auto&& edge : layout_.edge_list())
    {
        z3::expr edge_expr = edge_expr_array[edge->id()];
        z3::expr node_expr1 = node_expr_array[edge->node_a()->id()];
        z3::expr node_expr2 = node_expr_array[edge->node_b()->id()];

        solver.add(z3::implies(edge_expr == 1, node_expr1 == node_expr2));
    }
}

void ParallelizationOracle::make_remove_edge_constraint(z3::solver& solver,
                                                        z3::expr_vector& edge_expr_array)
{
    for (auto&& edge : layout_.horizontal_edge_list())
    {
        const int count = control_qubit_map_.count(edge->node_a()->name()) + control_qubit_map_.count(edge->node_b()->name());
        if ((edge->node_a()->name() != "A" || edge->node_b()->name() != "A") && count > 0)
        {
            z3::expr e = edge_expr_array[edge->id()];

            solver.add(e == 0);
        }
    }

    for (auto&& edge : layout_.vertical_edge_list())
    {
        const int count = target_qubit_map_.count(edge->node_a()->name()) + target_qubit_map_.count(edge->node_b()->name());
        if ((edge->node_a()->name() != "A" || edge->node_b()->name() != "A") && count > 0)
        {
            z3::expr e = edge_expr_array[edge->id()];

            solver.add(e == 0);
        }
    }
}

void ParallelizationOracle::verification(z3::solver& s)
{
    z3::model m = s.get_model();

    std::vector<int> node_var(m.size() + 1);
    std::vector<int> edge_var(m.size() + 1);

    // traversing the model
    for (unsigned i = 0; i < m.size(); i++) {
        z3::func_decl v = m[i];
        assert(v.arity() == 0);
        node_var[std::stoi(v.name().str())] = m.get_const_interp(v).get_numeral_int();
        edge_var[std::stoi(v.name().str())] = m.get_const_interp(v).get_numeral_int();
    }

    for (int y = 0; y < layout_.height(); y++)
    {
        for (int x = 0; x < layout_.width(); x++)
        {
            int index = layout_.edge_list().size() + 1 + (y * layout_.width()) + x;
            std::cout << node_var[index];
        }
        std::cout << std::endl;
    }
}

bool ParallelizationOracle::check(std::list<Gate>& gate_list,
                                  const Gate& new_gate)
{
    gate_list.push_back(new_gate);

    const bool result = check(gate_list);

    gate_list.pop_back();

    return result;
}

bool ParallelizationOracle::check(const std::list<Gate>& gate_list)
{
    num_gate_ = static_cast<int>(gate_list.size());

    if (num_gate_ == 1)
    {
        return true;
    }

    init(gate_list);

    z3::context c;
    z3::expr_vector edge_expr_array(c);
    z3::expr_vector node_expr_array(c);

    make_vars(c, edge_expr_array, node_expr_array);

    z3::solver s(c);

    make_base_constraint(s, edge_expr_array, node_expr_array);

    return s.check() == z3::sat;
}

}