#include "partition.hpp"

namespace tskd {
namespace tpar {

template<typename T>
bool is_disjoint(const std::set <T>& A, const std::set <T>& B)
{
    typename std::set<T>::iterator itA = A.begin(), itB = B.begin();

    while (itA != A.end() && itB != B.end())
    {
        if (*itA == *itB)
        {
            return false;
        }
        else if (*itA < *itB)
        {
            itA++;
        }
        else
        {
            itB++;
        }
    }

    return true;
}

std::ostream& operator<<(std::ostream& output, const partitioning& part)
{
    partitioning::const_iterator Si;
    std::set<int>::const_iterator yi;

    for (Si = part.begin(); Si != part.end(); Si++)
    {
        output << "{";
        for (yi = Si->begin(); yi != Si->end(); yi++)
        {
            output << *yi << ",";
        }
        output << "}";
    }

    return output;
}

// Take a partition and a set of ints, and return all partitions that are not
//   disjoint with the set, also removing them from the partition
partitioning freeze_partitions(partitioning& part, std::set<int>& st)
{
    partitioning ret;
    partitioning::iterator it, tmp;

    for (it = part.begin(); it != part.end();)
    {
        if (!is_disjoint(*it, st))
        {
            tmp = it;
            it++;
            ret.splice(ret.begin(), part, tmp);
        }
        else
        {
            it++;
        }
    }
    return ret;
}

int num_elts(partitioning& part)
{
    int tot = 0;
    partitioning::iterator it;
    for (it = part.begin(); it != part.end(); it++)
    {
        tot += it->size();
    }
    return tot;
}

partitioning create(std::set<int>& st)
{
    return partitioning(1, st);
}

}
}
