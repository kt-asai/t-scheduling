#ifndef T_SCHEDULING_GATE_H
#define T_SCHEDULING_GATE_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace tskd {

class Gate
{
private:
    std::string type_;
    std::vector<std::string> control_list_;
    std::vector<std::string> target_list_;

public:
    /**
     * constructor
     */
    Gate() = default;

    Gate(const std::string& type)
            : type_(type) { }

    /**
     * constructor for single qubit gate
     * @param type gate name
     * @param target targat qubit name
     */
    Gate(const std::string& type,
         const std::string& target)
            : type_(type)
    {
        target_list_.push_back(target);
    }

    /**
     * constructor for a 2 qubit gate (e.g.cnot)
     * @param type gate name
     * @param control control qubit name
     * @param target target qubit name
     */
    Gate(const std::string& type,
         const std::string& control,
         const std::string& target)
            : type_(type)
    {
        control_list_.push_back(control);
        target_list_.push_back(target);
    }

    /**
     * constructor for a multi-target gate (e.g.multi-target cnot)
     * @param type gate name
     * @param control control qubit name
     * @param target_list target qubit names
     */
    Gate(const std::string& type,
         const std::string& control,
         const std::vector<std::string>& target_list)
            : type_(type),
              target_list_(target_list)
    {
        control_list_.push_back(control);
    }

    /**
     * constructor of a multi-control gate (e.g.toffoli)
     * @param type gate name
     * @param control_list control qubit names
     * @param target target qubit name
     */
    Gate(const std::string& type,
         const std::vector<std::string>& control_list,
         const std::string& target)
            : type_(type),
              control_list_(control_list)
    {
        target_list_.push_back(target);
    }

    /**
     * constructor for a multi-control and multi-target gate
     * @param type gate name
     * @param control_list control qubit names
     * @param target_list target qubit names
     */
    Gate(const std::string& type,
         const std::vector<std::string>& control_list,
         const std::vector<std::string>& target_list)
            : type_(type),
              control_list_(control_list),
              target_list_(target_list)
    { }

    /**
     * retunr gate name
     * @return gate name
     */
    std::string type() const
    {
        return type_;
    }

    /**
     * return control qubit list
     * @return control list
     */
    std::vector<std::string> control_list() const
    {
        return control_list_;
    }

    /**
     * return target qubit list
     * @return target list
     */
    std::vector<std::string> target_list() const
    {
        return target_list_;
    }

    /**
     * print gate status
     */
    void print() const
    {
        std::cout << type_ << " ";
        for (auto& e : control_list_)
        {
            std::cout << e << " ";
        }
        for (auto& e : target_list_)
        {
            std::cout << e << " ";
        }
        std::cout << std::endl;
    }
};

}

#endif //T_SCHEDULING_GATE_H
