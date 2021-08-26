#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing tbb headers
#ifndef EQUATIONOFMOTIONSOLVER_H
#define EQUATIONOFMOTIONSOLVER_H
#include <tbb/concurrent_vector.h>
#include <Eigen/Core>
#include <map>
#include <unordered_map>
#include <vector>
#include "mosek.h"
#define DOFS 2

class EquationOfMotionSolver
{
public:
    EquationOfMotionSolver();
    ~EquationOfMotionSolver();

    void ClearAndResize(std::size_t N);     // size N must be set; return execution time
    void AddEntriesToStructure(int idx1, int idx2, int idx3);   // mark entries as non-zero for combinations of these indices
    void CreateStructure();

    // creating the values array
    void AddToEquation(const double &constTerm, const Eigen::Matrix<double,6,1> &linearTerm, const Eigen::Matrix<double,6,6> &quadraticTerm, const int (&ids)[3]);
    void AddToEquation(const double &constTerm, const Eigen::Vector2d &linearTerm, const Eigen::Matrix2d &quadraticTerm, int id);

    bool Solve();   // true if successful
    void TestSolve(); // test solve using sample data

    void GetTentativeResult(int idx, Eigen::Vector2d &vec);  // solution => convenient vector form

    MSKrealt objective_value;    // value of the optimized expression (should be near zero)
    double solution_norm;

private:
    MSKenv_t     env  = NULL;

    // nonzero values of the Q-matrix in the term 1/2 xt.Q.x; size is nnz*DOFs
    std::vector<MSKint32t> qosubi, qosubj;
    std::vector<MSKrealt> qoval;

    // linear term
    std::vector<MSKint32t> csubj;
    std::vector<MSKrealt> cval, sln;

    MSKrealt cfix;    // constant term

    unsigned N;      // number of variables (divided by DOF)
    unsigned nnz;    // number of non-zero entries in Q (lower triangle)

    std::vector<tbb::concurrent_vector<unsigned>*> rows_Neighbors;
    std::vector<std::vector<unsigned>*> rows_pcsr;   // per row mappings between columns and offset in "values"

    static void MSKAPI printstr(void *, const char str[]);
    void AddNNZEntry(int row, int column);    // reserve non-zero positions one-by-one (thread-safe)

    void ResizeRows();

    void AddToQ(const int row, const int column, const Eigen::Matrix2d &mat);
    void AddToC(const int idx, const Eigen::Vector2d &vec);
    void AddToConstTerm(const double &c);
};

#endif // EQUATIONOFMOTIONSOLVER_H
#endif
