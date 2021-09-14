#include "equationofmotionsolver.h"
#include <iostream>
#include <algorithm>
#include "spdlog/spdlog.h"

EquationOfMotionSolver::EquationOfMotionSolver()
{
    MSKrescodee  r;

    r = MSK_makeenv(&env, NULL);
    if (r != MSK_RES_OK) throw std::runtime_error("makeenv");

}

EquationOfMotionSolver::~EquationOfMotionSolver()
{
    MSKrescodee  r;
    r = MSK_deleteenv(&env);
    if (r != MSK_RES_OK) spdlog::critical("~EquationOfMotionSolver: MSK_deleteenv returned {}",r);
    else spdlog::info("~EquationOfMotionSolver() done");
}


void MSKAPI EquationOfMotionSolver::printstr(void *, const char str[])
{
    spdlog::info("{}",str);
}


void EquationOfMotionSolver::ClearAndResize(std::size_t N_)
{
    this->N=N_;
    cfix = 0;
    if(csubj.size() < N*DOFS) csubj.resize(N*DOFS*1.5);
    if(cval.size() < N*DOFS) cval.resize(N*DOFS*1.5);
    if(sln.size() < N*DOFS) sln.resize(N*DOFS*1.5);

    std::fill(cval.begin(), cval.begin()+N*DOFS, 0);

    while(rows_Neighbors.size()<N)
        rows_Neighbors.push_back(std::make_unique<tbb::concurrent_vector<unsigned>>(10));

    while(rows_pcsr.size()<N)
        rows_pcsr.push_back(std::make_unique<std::vector<unsigned>>(10));

#pragma omp parallel for
    for(unsigned i=0;i<N;i++)
    {
        rows_Neighbors[i]->clear();
        rows_Neighbors[i]->push_back(i);    // diagonal elements must be non-zero
        rows_pcsr[i]->clear(); // clear the mapping of (i,j)->offset
        csubj[i*DOFS+0]=i*DOFS+0;
        csubj[i*DOFS+1]=i*DOFS+1;
    }
}

void EquationOfMotionSolver::AddNNZEntry(int row, int column)
{
    if(row < 0 || column < 0) return; // the element does not belong in the matrix
    if(row < column) std::swap(row,column);    // enforce lower-triangular matrix
    if((unsigned)row >= N) throw std::runtime_error("trying to insert an element beyond the matrix size");
    rows_Neighbors[row]->push_back(column);
}

void EquationOfMotionSolver::AddEntriesToStructure(int idx1, int idx2, int idx3)
{
    AddNNZEntry(idx1,idx2);
    AddNNZEntry(idx1,idx3);
    AddNNZEntry(idx3,idx2);
}

void EquationOfMotionSolver::CreateStructure()
{
    // CREATE STRUCTURE ARRAYS

    // sort the neighbor list of each row
#pragma omp parallel for
    for(unsigned i=0;i<N;i++)
    {
        tbb::concurrent_vector<unsigned> &rn = *rows_Neighbors[i];
        std::sort(rn.begin(),rn.end());
        auto unique_res = std::unique(rn.begin(), rn.end());
        unsigned newSize = std::distance(rn.begin(),unique_res);
        rn.resize(newSize);
    }

    // count non-zero entries
    nnz = 0;
#pragma omp parallel for reduction(+:nnz)
    for(unsigned i=0;i<N;i++) nnz+=(rows_Neighbors[i]->size()*DOFS*DOFS-1);

    // ensure that the arrays are of sufficient sizes
    if(qosubi.size() < nnz) qosubi.resize(nnz*1.5);
    if(qosubj.size() < nnz) qosubj.resize(nnz*1.5);
    if(qoval.size() < nnz) qoval.resize(nnz*1.5);
    std::fill(qoval.begin(), qoval.begin()+nnz, 0);

    // enumerate entries
    unsigned count=0;
    for(unsigned row=0;row<N;row++)
    {
        tbb::concurrent_vector<unsigned> &sorted_vec = *rows_Neighbors[row];
        if(sorted_vec.size() == 0) throw std::runtime_error("matrix row contains no entries");

        int previous_column = -1;
        for(unsigned int const &local_column : sorted_vec)
        {
            rows_pcsr[row]->push_back(count);
            if(row > local_column) {
                qosubi[count+0]=row*DOFS+0;
                qosubj[count+0]=local_column*DOFS+0;
                qosubi[count+1]=row*DOFS+0;
                qosubj[count+1]=local_column*DOFS+1;
                qosubi[count+2]=row*DOFS+1;
                qosubj[count+2]=local_column*DOFS+0;
                qosubi[count+3]=row*DOFS+1;
                qosubj[count+3]=local_column*DOFS+1;

                count+=DOFS*DOFS;
            }
            else if(row == local_column)
            {
                qosubi[count+0]=row*DOFS+0;
                qosubj[count+0]=local_column*DOFS+0;
                qosubi[count+1]=row*DOFS+1;
                qosubj[count+1]=local_column*DOFS+0;
                qosubi[count+2]=row*DOFS+1;
                qosubj[count+2]=local_column*DOFS+1;

                count+=DOFS*DOFS-1;
            }
            else throw std::runtime_error("matrix is not lower-triangular");
            if((int)local_column <= previous_column) throw std::runtime_error("column entries are not sorted");
            previous_column = local_column;
        }
    }

    if(nnz!=count)
    {
        spdlog::critical("csr_rows[{}]=={}, whereas count=={}",N,nnz,count);
        throw std::runtime_error("nnz != count");
    }
}

// creating the values array
void EquationOfMotionSolver::AddToQ(const int row, const int column, const Eigen::Matrix2d &mat)
{
    if (row < 0 || column < 0 || row < column) return;
    else if((unsigned)row >= N || (unsigned)column >= N) throw std::runtime_error("AddToQ: out of range");

    // find the value array offset corresponding to the "row/column" entry
    int offset = -1;
    tbb::concurrent_vector<unsigned>&vec = *rows_Neighbors[row];

    for(unsigned count = 0;count<vec.size();count++)
    {
        if(vec.at(count) == (unsigned)column)
        {
            offset = rows_pcsr[row]->at(count);
            break;
        }
    }

    if(offset<0) throw std::runtime_error("AddToQ: column index not found");
    else if((unsigned)offset >= nnz) throw std::runtime_error("AddToQ: offset >= nnz");

    if(row > column)
    {
#pragma omp atomic
            qoval[offset+0] += mat.coeff(0,0);
#pragma omp atomic
            qoval[offset+1] += mat.coeff(0,1);
#pragma omp atomic
            qoval[offset+2] += mat.coeff(1,0);
#pragma omp atomic
            qoval[offset+3] += mat.coeff(1,1);
    }
    else if(row == column)
    {
#pragma omp atomic
            qoval[offset+0] += mat.coeff(0,0);
#pragma omp atomic
            qoval[offset+1] += mat.coeff(1,0);
#pragma omp atomic
            qoval[offset+2] += mat.coeff(1,1);
    }
}

void EquationOfMotionSolver::AddToC(const int idx, const Eigen::Vector2d &vec)
{
    if(idx < 0) return;
    if((unsigned)idx >= N) throw std::runtime_error("AddToC: index out of range");

#pragma omp atomic
    cval[idx*DOFS+0]+=vec[0];
#pragma omp atomic
    cval[idx*DOFS+1]+=vec[1];
}

void EquationOfMotionSolver::AddToConstTerm(const double &c)
{
#pragma omp atomic
    cfix+=c;
}

void EquationOfMotionSolver::AddToEquation(const double &constTerm,
                                           const Eigen::Matrix<double,6,1> &linearTerm,
                                           const Eigen::Matrix<double,6,6> &quadraticTerm,
                                           const int (&ids)[3])
{
    AddToConstTerm(constTerm);

    for(int i=0;i<3;i++)
    {
        int row = ids[i];
        if(row < 0) continue;
        AddToC(row, linearTerm.block(i*2,0,2,1));
        for(int j=0;j<3;j++)
        {
            int col = ids[j];
            if(col < 0) continue;
            AddToQ(row, col, quadraticTerm.block(i*2,j*2,2,2));
        }
    }
}

void EquationOfMotionSolver::AddToEquation(const double &constTerm,
                                           const Eigen::Vector2d &linearTerm,
                                           const Eigen::Matrix2d &quadraticTerm, int id)
{
    AddToConstTerm(constTerm);
    AddToC(id, linearTerm);
    AddToQ(id, id, quadraticTerm);
}


bool EquationOfMotionSolver::Solve()
{
    MSKrescodee  r;
    MSKtask_t    task = NULL;
    r = MSK_maketask(env, 0, N*DOFS, &task);
    if (r != MSK_RES_OK) throw std::runtime_error("maketask");

    //    r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
    //    if (r != MSK_RES_OK) throw std::runtime_error("linkfunctotaskstream");

    int numvar = N*DOFS;
    r = MSK_appendvars(task, numvar);
    if (r != MSK_RES_OK) throw std::runtime_error("appendvars");

    for (int j = 0; j < numvar; j++)
    {
        r = MSK_putvarbound(task, j, MSK_BK_FR, -MSK_DPAR_DATA_TOL_BOUND_INF, MSK_DPAR_DATA_TOL_BOUND_INF);
        if (r != MSK_RES_OK) throw std::runtime_error("MSK_putvarbound");
    }

    r = MSK_putclist(task, numvar, csubj.data(), cval.data());
    if (r != MSK_RES_OK)
    {
        spdlog::critical("MSK_putclist returns {}",r);
        throw std::runtime_error("MSK_putclist");
    }

    r = MSK_putqobj(task, nnz, qosubi.data(), qosubj.data(), qoval.data());
    if (r != MSK_RES_OK) throw std::runtime_error("MSK_putqobj");

    r = MSK_putcfix(task, cfix);
    if (r != MSK_RES_OK) throw std::runtime_error("MSK_putcfix");

    r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
    if (r != MSK_RES_OK) throw std::runtime_error("MSK_putobjsense");


    MSKrescodee trmcode;

    r = MSK_optimizetrm(task, &trmcode);
    if(r == MSK_RES_ERR_OBJ_Q_NOT_PSD)
    {
//        spdlog::info("EquationOfMotionSolver: The quadratic coefficient matrix in the objective is not positive semidefinite");
        r = MSK_deletetask(&task);
        if (r != MSK_RES_OK) spdlog::warn("MSK_deletetask error");
        return false;
    }
    if (r != MSK_RES_OK)
    {
        spdlog::critical("EquationOfMotionSolver: MSK_optimizetrm returns", r);
        throw std::runtime_error("MSK_optimizetrm");
    }

    //     MSK_solutionsummary(task, MSK_STREAM_LOG);

    MSKsolstae solsta;
    r = MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
    if (r != MSK_RES_OK) throw std::runtime_error("MSK_getsolsta result");
    if (solsta != MSK_SOL_STA_OPTIMAL) throw std::runtime_error("solsta != MSK_SOL_STA_OPTIMAL");

    MSK_getxx(task, MSK_SOL_ITR, sln.data());

    MSK_getdualobj(task, MSK_SOL_ITR, &objective_value);
    //     std::cout << "\nMSK_SOL_ITR sol = " << objective_value << std::endl;

    r = MSK_deletetask(&task);
    if (r != MSK_RES_OK) spdlog::warn("MSK_deletetask error");

    solution_norm = 0;
#pragma omp parallel for reduction(+:solution_norm)
    for(unsigned i=0;i<N*DOFS;i++) solution_norm+=(double)(sln[i]*sln[i]);

    solution_norm = sqrt(solution_norm);

    return true;
}


void EquationOfMotionSolver::GetTentativeResult(int idx, Eigen::Vector2d &vec)
{
    vec[0] = sln[idx*DOFS+0];
    vec[1] = sln[idx*DOFS+1];
}
