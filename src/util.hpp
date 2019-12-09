#ifndef T_SCHEDULING_UTIL_HPP
#define T_SCHEDULING_UTIL_HPP

#include <set>
#include <list>
#include <boost/dynamic_bitset.hpp>

#include "gate.hpp"

namespace tskd {
namespace util {

using xor_func = boost::dynamic_bitset<>;
using phase_exponent = std::pair<int, xor_func>;

class IndependentOracle
{
private:
    int num_;
    int dim_;
    int length_;

public:
    /**
     * constructor
     */
    IndependentOracle()
            : num_(0),
              dim_(0),
              length_(0) { }

    /**
     * contructor
     * @param numin
     * @param dimin
     * @param lengthin
     */
    IndependentOracle(int numin, int dimin, int lengthin)
            : num_(numin),
              dim_(dimin),
              length_(lengthin) { }

    /**
     * matroid oracle
     * @param expnts
     * @param lst
     * @return
     */
    bool operator()(const std::vector<phase_exponent>& expnts,
                    const std::set<int>& lst) const;

    /**
     * set new dimension
     * @param newdim new dimension
     */
    void set_dim(int newdim) { dim_ = newdim; }

    /**
     * Shortcut to find a linearly dependent element faster
     * @param expnts
     * @param lst
     * @return
     */
    int retrieve_lin_dep(const std::vector<phase_exponent>& expnts,
                         const std::set<int>& lst) const;
};

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

std::list<Gate> ToUpperEchelon(int m,
                               int n,
                               std::vector<xor_func>& bits,
                               std::vector<xor_func> *mat,
                               const std::vector<std::string>& qubit_names);

std::list<Gate> ToLowerEchelon(int m,
                               int n,
                               std::vector<xor_func>& bits,
                               std::vector<xor_func>* mat,
                               const std::vector<std::string>& qubit_names);

std::list<Gate> FixBasis(int m,
                         int n,
                         int k,
                         const std::vector<xor_func>& fst,
                         std::vector<xor_func>& snd,
                         std::vector<xor_func>* mat,
                         const std::vector<std::string>& qubit_names);

/*
 * A := B^{-1} A
 */
void Compose(int num,
             std::vector<xor_func>& A,
             const std::vector<xor_func>& B);

std::list<Gate> ComposeX(int target,
                         const std::vector<std::string>& qubit_names);

std::list<Gate> ComposeSwap(int a,
                            int b,
                            const std::vector<std::string>& qubit_names);

std::list<Gate> ComposeCNOT(int target,
                            int control,
                            const std::vector<std::string>& qubit_names);

}
}

#endif //T_SCHEDULING_UTIL_HPP
