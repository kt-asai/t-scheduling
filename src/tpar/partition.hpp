#ifndef T_SCHEDULING_PARTITION_HPP
#define T_SCHEDULING_PARTITION_HPP

#include <list>
#include <set>
#include <iostream>

namespace tskd {
namespace tpar {

using partitioning = std::list<std::set<int>>;
using path_iterator = std::list<std::pair<int, partitioning::iterator>>::iterator;

std::ostream& operator<<(std::ostream& output,
                         const partitioning& part);

partitioning freeze_partitions(partitioning& part,
                               std::set<int>& st);

int num_elts(partitioning& part);

partitioning create(std::set<int>& st);

}
}

#endif //T_SCHEDULING_PARTITION_HPP