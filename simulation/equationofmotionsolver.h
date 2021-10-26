#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing tbb headers
#ifndef EQUATIONOFMOTIONSOLVER_H
#define EQUATIONOFMOTIONSOLVER_H
#include <tbb/concurrent_vector.h>
#include <Eigen/Core>
#include <map>
#include <unordered_map>
#include <initializer_list>
#include <vector>
#include <memory>
#include "mosek.h"

class EquationOfMotionSolver
{
public:
    EquationOfMotionSolver();
    ~EquationOfMotionSolver();
    EquationOfMotionSolver& operator=(EquationOfMotionSolver&) = delete;

    void ClearAndResize(std::size_t N);     // size N must be set; return execution time

    void AddEntriesToStructure(const int* idx_begin, const int* idx_end); // insert nxn matrix of indices of non-zero entries
    void CreateStructure();

    // add values to non-zero elements
    void AddToEquation(const double *linearEntries, const double *quadraticEntries, const std::initializer_list<int> ids);

    constexpr static unsigned dofs = 2; // number of degrees of freedom per node (size of per-node blocks)
    constexpr static unsigned dofssq = dofs*dofs;

    // MOSEK
    MSKrealt objective_value;    // value of the optimized expression (should be near zero)
    double solution_norm;
    bool Solve();   // true if successful
    void GetTentativeResult(int idx, Eigen::Vector2d &vec);  // solution => convenient vector form

private:
    // nonzero values of the Q-matrix in the term 1/2 xt.Q.x; size is nnz*DOFs
    std::vector<int> qosubi, qosubj;
    std::vector<double> qoval;

    std::vector<int> csr_rows, csr_cols;

    // linear term
    std::vector<int> csubj;
    std::vector<double> cval, sln;

    unsigned N;      // number of variables (divided by DOF)
    unsigned nnz;    // number of non-zero entries in Q (lower triangle)

    std::vector<std::unique_ptr<tbb::concurrent_vector<unsigned>>> rows_neighbors;  // list of indices of nz-columns per row

    void AddNNZEntry(int row, int column);    // reserve non-zero positions one-by-one (thread-safe)

    void ResizeRows();

    void AddToQ(const int row, const int column, const double v11, const double v12, const double v21, const double v22);
    void AddToC(const int idx, const double v1, const double v2);
    unsigned get_offset(const int row, const int column) const;

    // MOSEK
    MSKenv_t env  = NULL;
    static void MSKAPI printstr(void *, const char str[]);
};

#endif // EQUATIONOFMOTIONSOLVER_H
#endif
