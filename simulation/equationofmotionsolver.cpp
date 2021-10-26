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
    if(csubj.size() < N*dofs)
    {
        csubj.resize(N*dofs*2);
        cval.resize(N*dofs*2);
        sln.resize(N*dofs*2);
    }

    if(csr_rows.size() < N+1) csr_rows.resize(N*2+1);

    std::fill(cval.begin(), cval.begin()+N*dofs, 0);

    if(rows_neighbors.size()<N)
    {
        rows_neighbors.reserve(N*2);
        while(rows_neighbors.size()<N*2)
            rows_neighbors.push_back(std::make_unique<tbb::concurrent_vector<unsigned>>(20));
    }

#pragma omp parallel for
    for(unsigned i=0;i<N;i++)
    {
        rows_neighbors[i]->clear();
        csubj[i*dofs+0]=i*dofs+0;
        csubj[i*dofs+1]=i*dofs+1;
    }
}

void EquationOfMotionSolver::AddNNZEntry(int row, int column)
{
    if(row < 0 || column < 0) return; // the element does not belong in the matrix
    if(row < column) std::swap(row,column);    // enforce lower-triangular matrix
    if((unsigned)row >= N) throw std::runtime_error("trying to insert an element beyond the matrix size");
    rows_neighbors[row]->push_back(column);
}

void EquationOfMotionSolver::AddEntriesToStructure(const int* idx_begin, const int* idx_end)
{
    for(auto iter=(idx_begin+1); iter!=idx_end; ++iter)
        for(auto j=idx_begin; j!=iter; ++j)
            AddNNZEntry(*iter,*j);
}


void EquationOfMotionSolver::CreateStructure()
{
    // CREATE STRUCTURE ARRAYS

    // sort the neighbor list of each row
    nnz = 0;
#pragma omp parallel for reduction(+:nnz)
    for(unsigned i=0;i<N;i++)
    {
        tbb::concurrent_vector<unsigned> &rn = *rows_neighbors[i];
        rn.push_back(i);    // add diagonal entry
        std::sort(rn.begin(),rn.end());
        rn.resize(std::distance(rn.begin(),std::unique(rn.begin(), rn.end())));     // remove duplicates
        nnz+=(rn.size());
    }

    csr_rows[N] = nnz;
    if(csr_cols.size() < nnz) csr_cols.resize(nnz*1.5);
    nnz*=dofssq;

    // ensure that the arrays are of sufficient sizes
    if(qosubi.size() < nnz) qosubi.resize(nnz*1.5);
    if(qosubj.size() < nnz) qosubj.resize(nnz*1.5);
    if(qoval.size() < nnz) qoval.resize(nnz*1.5);
    std::fill(qoval.begin(), qoval.begin()+nnz, 0);

    // enumerate entries
    unsigned count=0;
    for(unsigned row=0;row<N;row++)
    {
        csr_rows[row] = count;

        tbb::concurrent_vector<unsigned> &rn = *rows_neighbors[row];

        for(unsigned int const &local_column : rn)
        {
            if(row < local_column) throw std::runtime_error("matrix is not lower-triangular");

            csr_cols[count] = local_column;

            qosubi[dofssq*count+0]=row*dofs+0;
            qosubj[dofssq*count+0]=local_column*dofs+0;
            qosubi[dofssq*count+1]=row*dofs+0;
            qosubj[dofssq*count+1]=local_column*dofs+1;
            qosubi[dofssq*count+2]=row*dofs+1;
            qosubj[dofssq*count+2]=local_column*dofs+0;
            qosubi[dofssq*count+3]=row*dofs+1;
            qosubj[dofssq*count+3]=local_column*dofs+1;

            if(row==local_column)
            {
                qosubi[dofssq*count+1]=row*dofs;
                qosubj[dofssq*count+1]=local_column*dofs;
            }

            count++;
        }
    }

    if(nnz != dofssq*count)
    {
        spdlog::critical("csr_rows[{}]=={}, whereas count=={}",N,nnz,count);
        throw std::runtime_error("nnz != count");
    }
}

unsigned EquationOfMotionSolver::get_offset(const int row, const int column) const
{

    int col_offset_begin = csr_rows[row];
    int col_offset_end = csr_rows[row+1];

    const int *start_pt = &csr_cols[col_offset_begin];
    const int *end_pt = &csr_cols[col_offset_end];

    auto it = std::lower_bound(start_pt,end_pt,column);
    if(it == end_pt || *it != column)
        throw std::runtime_error("EquationOfMotionSolver::get_offset(): (i,j) index not found");
    unsigned offset = std::distance(start_pt,it)+col_offset_begin;
    offset*=dofssq;

/*

    // find the value array offset corresponding to the "row/column" entry
    int offset2 = -1;
    tbb::concurrent_vector<unsigned>&vec = *rows_neighbors[row];

    for(unsigned count = 0;count<vec.size();count++)
    {
        if(vec.at(count) == (unsigned)column)
        {
            offset2 = rows_pcsr[row]->at(count);
            break;
        }
    }

    if(offset2<0) throw std::runtime_error("AddToQ: column index not found");
    else if((unsigned)offset2 >= nnz) throw std::runtime_error("AddToQ: offset >= nnz");

    if(offset2 != offset)
    {
        spdlog::info("offset {}; offset2 {}",offset,offset2);
        throw std::runtime_error("offsets don't match");
    }
*/
    return offset;
}

void EquationOfMotionSolver::AddToQ(const int row, const int column, const double v11, const double v12, const double v21, const double v22)
{
    if (row < 0 || column < 0 || row < column) return;
    else if((unsigned)row >= N || (unsigned)column >= N) throw std::runtime_error("AddToQ: out of range");
    int offset = get_offset(row,column);



    if(row > column)
    {
#pragma omp atomic
        qoval[offset+0] += v11;
#pragma omp atomic
        qoval[offset+1] += v12;
#pragma omp atomic
        qoval[offset+2] += v21;
#pragma omp atomic
        qoval[offset+3] += v22;
    }
    else if(row == column)
    {
#pragma omp atomic
        qoval[offset+0] += v11;
#pragma omp atomic
        qoval[offset+2] += v21;
#pragma omp atomic
        qoval[offset+3] += v22;
    }
}

void EquationOfMotionSolver::AddToC(const int idx, const double v1, const double v2)
{
    if(idx < 0) return;
    if((unsigned)idx >= N) throw std::runtime_error("AddToC: index out of range");

#pragma omp atomic
    cval[idx*dofs+0] += v1;
#pragma omp atomic
    cval[idx*dofs+1] += v2;
}



void EquationOfMotionSolver::AddToEquation(const double *lE, const double *qE, const std::initializer_list<int> ids)
{
    unsigned n = ids.size();
    unsigned n2 = n*dofs;
    unsigned i_ct=0;
    for(auto i=ids.begin();i!=ids.end();i++,i_ct+=dofs)
    {
        int row = *i;
        if(row < 0) continue;
        AddToC(row, *(lE+i_ct), *(lE+i_ct+1));
        unsigned j_ct=0;
        for(auto j=ids.begin();j!=ids.end();j++,j_ct+=dofs)
        {
            int col = *j;
            if(col < 0) continue;
            double m11 = *(qE+j_ct*n2+i_ct);
            double m12 = *(qE+(j_ct+1)*n2+i_ct);
            double m21 = *(qE+j_ct*n2+i_ct+1);
            double m22 = *(qE+(j_ct+1)*n2+i_ct+1);
            AddToQ(row, col, m11, m12, m21, m22);
        }
    }
}






// MOSEK

bool EquationOfMotionSolver::Solve()
{
    MSKrescodee  r;
    MSKtask_t    task = NULL;
    r = MSK_maketask(env, 0, N*dofs, &task);
    if (r != MSK_RES_OK) throw std::runtime_error("maketask");

    //    r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
    //    if (r != MSK_RES_OK) throw std::runtime_error("linkfunctotaskstream");

    int numvar = N*dofs;
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

    //r = MSK_putcfix(task, cfix);
    //if (r != MSK_RES_OK) throw std::runtime_error("MSK_putcfix");

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
    for(unsigned i=0;i<N*dofs;i++) solution_norm+=(double)(sln[i]*sln[i]);

    solution_norm = sqrt(solution_norm);

    return true;
}


void EquationOfMotionSolver::GetTentativeResult(int idx, Eigen::Vector2d &vec)
{
    vec[0] = sln[idx*dofs+0];
    vec[1] = sln[idx*dofs+1];
}
