#ifndef T_SCHEDULING_UTIL_HPP
#define T_SCHEDULING_UTIL_HPP

#include <boost/dynamic_bitset.hpp>

namespace tskd {
namespace util {

using xor_func = boost::dynamic_bitset<>;
using phase_exponent = std::pair<int, xor_func>;

/**
* make triangular to determine the rank (destructive)
* @param num_qubit number of qubit in the circuit
* @param num_qubit_and_hadamard number of qubit and hadamard in the circuit
* @param bits parity matrix
* @return rank of matrix
*/
int ComputeRankDestructive(int num_qubit,
                           int num_qubit_and_hadamard,
                           std::vector<xor_func>& bits);

/**
* make triangular to determine the rank
* @param num_qubit number of qubit in the circuit
* @param num_qubit_and_hadamard number of qubit and hadamard in the circuit
* @param bits parity matrix
* @return rank of matrix
*/
int ComputeRank(int num_qubit,
                int num_qubit_and_hadamard,
                std::vector<xor_func>& bits);

/**
* check linear independence of one vector wrt a matrix (destructive)
* @param num_qubit number of qubit in the circuit
* @param bits  parity matrix
* @param parity parity to check
* @return whether parity is linearly independent of bits
*/
bool IsIndependentDestructive(int num_qubit,
                              const std::vector<xor_func>& bits,
                              xor_func& parity);

/**
* check linear independence of one vector wrt a matrix
* @param num_qubit number of qubit in the circuit
* @param bits  parity matrix
* @param parity parity to check
* @return whether parity is linearly independent of bits
*/
bool IsIndependent(int num_qubit,
                   const std::vector<xor_func>& bits,
                   const xor_func& parity);

}
}

#endif //T_SCHEDULING_UTIL_HPP
