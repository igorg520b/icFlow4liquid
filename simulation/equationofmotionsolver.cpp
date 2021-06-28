#include "equationofmotionsolver.h"
#include <iostream>
#include <algorithm>

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
    if (r != MSK_RES_OK) std::cout << "~EquationOfMotionSolver MSK_deleteenv" << std::endl;

    // free the elements of rows_Neighbors and rows_pcsr
    for(std::size_t i=0;i<rows_Neighbors.size();i++) delete rows_Neighbors[i];
    for(std::size_t i=0;i<rows_pcsr.size();i++) delete rows_pcsr[i];
    std::cout << "~EquationOfMotionSolver() done" << std::endl;
}


void MSKAPI EquationOfMotionSolver::printstr(void *, const char str[])
{
    std::cout << str << std::endl;
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
        rows_Neighbors.push_back(new tbb::concurrent_vector<unsigned>(10));

    while(rows_pcsr.size()<N)
        rows_pcsr.push_back(new std::vector<unsigned>(10));

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

void EquationOfMotionSolver::AddElementToStructure(int row, int column)
{
    if(row < 0 || column < 0) return; // the element does not belong in the matrix
    else if((unsigned)row >= N || (unsigned)column >= N) throw std::runtime_error("trying to insert element beyond the matrix size");

    // lower-triangular matrix, so enforce row>=column
    if(row>=column) rows_Neighbors[row]->push_back(column);
    else rows_Neighbors[column]->push_back(row);
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
        if(rows_Neighbors[row]->size() == 0) throw std::runtime_error("matrix row contains no entries");
        tbb::concurrent_vector<unsigned> &sorted_vec = *rows_Neighbors[row];

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



            if((int)local_column <= previous_column) {
                std::cout << "sorted_vec: ";
                for(unsigned int const &idx : sorted_vec) std::cout << idx << ", ";
                throw std::runtime_error("column entries are not sorted");
            }
            previous_column = local_column;
        }
    }

    if(nnz!=count)
    {
        std::cout << "csr_rows["<<N<<"]="<<nnz<< std::endl;
        std::cout << "nnz="<<nnz<<"; count="<<count<< std::endl;
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
    tbb::concurrent_vector<unsigned>*vec = rows_Neighbors[row];
    for(unsigned count = 0;count<vec->size();count++)
    {
        if(vec->at(count) == (unsigned)column)
        {
            offset = rows_pcsr[row]->at(count);
            break;
        }
    }
    if(offset<0) throw std::runtime_error("AddToQ: column index not found");
    else if((unsigned)offset >= nnz) throw std::runtime_error("offset >= nnz");

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
    cval[idx*DOFS+0]+=vec.x();
#pragma omp atomic
    cval[idx*DOFS+1]+=vec.y();
}

void EquationOfMotionSolver::AddToConstTerm(const double c)
{
#pragma omp atomic
    cfix+=c;
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
        std::cout << "MSK_putclist returns " << r << std::endl;
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
        std::cout << "The quadratic coefficient matrix in the objective is not positive semidefinite" << std::endl;
        return false;
    }
    if (r != MSK_RES_OK)
    {
        std::cout << "MSK_optimizetrm returns " << r << std::endl;
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
    if (r != MSK_RES_OK) std::cout << "MSK_deletetask error" << std::endl;

    // calculate solution norm
    solution_norm_prev = solution_norm;
    solution_norm = 0;
#pragma omp parallel for reduction(+:solution_norm)
    for(int i=0;i<N*DOFS;i++) solution_norm+=(double)(sln[i]*sln[i]);

    solution_norm = sqrt(solution_norm);

    return true;
}


void EquationOfMotionSolver::GetTentativeResult(int idx, Eigen::Vector2d &vec)
{
    vec.x() = sln[idx*DOFS+0];
    vec.y() = sln[idx*DOFS+1];
}

void EquationOfMotionSolver::TestSolve()
{
    MSKrescodee  r;

    MSKtask_t    task = NULL;
    r = MSK_maketask(env, 0, 10000, &task);
    if (r != MSK_RES_OK) throw std::runtime_error("maketask");

//    r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
//    if (r != MSK_RES_OK) throw std::runtime_error("linkfunctotaskstream");

    // TASK1

    int numvar = 2;
    r = MSK_appendvars(task, numvar);
     if (r != MSK_RES_OK) throw std::runtime_error("appendvars");

     for (int j = 0; j < numvar; j++)
     {
         r = MSK_putvarbound(task, j, MSK_BK_FR, -MSK_DPAR_DATA_TOL_BOUND_INF, MSK_DPAR_DATA_TOL_BOUND_INF);
         if (r != MSK_RES_OK) throw std::runtime_error("MSK_putvarbound");
     }

     const MSKint32t csubj[] = {0,1};
     const MSKrealt cval[] = {5.0, -4.0};
     r = MSK_putclist(task, 2, csubj, cval);
     if (r != MSK_RES_OK) throw std::runtime_error("MSK_putclist");

     const MSKint32t qosubi[] = {0,1,1};
     const MSKint32t qosubj[] = {0,0,1};
     const MSKrealt qoval[] = {2.0, 1.0, 8.0};

     r = MSK_putqobj(task, 3, qosubi, qosubj, qoval);
     if (r != MSK_RES_OK) throw std::runtime_error("MSK_putqobj");

     r = MSK_putcfix(task, 10);
     if (r != MSK_RES_OK) throw std::runtime_error("MSK_putcfix");

     r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
     if (r != MSK_RES_OK) throw std::runtime_error("MSK_putobjsense");


     MSKrescodee trmcode;

     r = MSK_optimizetrm(task, &trmcode);
     if (r != MSK_RES_OK) throw std::runtime_error("MSK_optimizetrm");

//     MSK_solutionsummary(task, MSK_STREAM_LOG);

     MSKsolstae solsta;
     r = MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
     if (r != MSK_RES_OK) throw std::runtime_error("MSK_getsolsta result");
     if (solsta != MSK_SOL_STA_OPTIMAL) throw std::runtime_error("solsta != MSK_SOL_STA_OPTIMAL");

     MSKrealt result[numvar];
     MSK_getxx(task, MSK_SOL_ITR, result);

     std::cout << "\nx0 = " << result[0];
     std::cout << "\nx1 = " << result[1] << std::endl;

     MSKrealt dualobj;
     MSK_getdualobj(task, MSK_SOL_ITR, &dualobj);
     std::cout << "\nMSK_SOL_ITR sol = " << dualobj << std::endl;

     r = MSK_deletetask(&task);
     if (r != MSK_RES_OK) std::cout << "MSK_deletetask error" << std::endl;





     // TASK 2
     r = MSK_maketask(env, 0, 10000, &task);
     if (r != MSK_RES_OK) throw std::runtime_error("maketask");

//     r = MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
//     if (r != MSK_RES_OK) throw std::runtime_error("linkfunctotaskstream");

     numvar = 3;
     r = MSK_appendvars(task, numvar);
      if (r != MSK_RES_OK) throw std::runtime_error("appendvars");

     for (int j = 0; j < numvar; j++)
     {
         r = MSK_putvarbound(task, j, MSK_BK_FR, -MSK_DPAR_DATA_TOL_BOUND_INF, MSK_DPAR_DATA_TOL_BOUND_INF);
         if (r != MSK_RES_OK) throw std::runtime_error("MSK_putvarbound");
     }

     const MSKint32t csubj2[] = {0,1,2};
     const MSKrealt cval2[] = {5.0, -4.0, 7.0};
     r = MSK_putclist(task, numvar, csubj2, cval2);
     if (r != MSK_RES_OK) throw std::runtime_error("MSK_putclist");

     const MSKint32t qosubi2[] = {0,   1,   2,   2,   2};
     const MSKint32t qosubj2[] = {0,   1,   0,   1,   2};
     const MSKrealt qoval2[] =   {2.0, 3.0, 1.0, 1.0, 4.0};

     r = MSK_putqobj(task, 5, qosubi2, qosubj2, qoval2);
     if (r != MSK_RES_OK) throw std::runtime_error("MSK_putqobj");

     r = MSK_putcfix(task, 10);
     if (r != MSK_RES_OK) throw std::runtime_error("MSK_putcfix");

     r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
     if (r != MSK_RES_OK) throw std::runtime_error("MSK_putobjsense");


     r = MSK_optimizetrm(task, &trmcode);
     if (r != MSK_RES_OK) throw std::runtime_error("MSK_optimizetrm");

//     MSK_solutionsummary(task, MSK_STREAM_LOG);

     r = MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
     if (r != MSK_RES_OK) throw std::runtime_error("MSK_getsolsta result");
     if (solsta != MSK_SOL_STA_OPTIMAL) throw std::runtime_error("solsta != MSK_SOL_STA_OPTIMAL");

     MSKrealt result2[3];
     MSK_getxx(task, MSK_SOL_ITR, result2);

     std::cout << "\nx0 = " << result2[0];
     std::cout << "\nx1 = " << result2[1];
     std::cout << "\nx2 = " << result2[2] << std::endl;

     MSK_getdualobj(task, MSK_SOL_ITR, &dualobj);
     std::cout << "\nMSK_SOL_ITR sol = " << dualobj << std::endl;

     r = MSK_deletetask(&task);
     if (r != MSK_RES_OK) std::cout << "MSK_deletetask error" << std::endl;


}
