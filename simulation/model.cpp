#include <vtkPointData.h>
#include <QtGlobal>
#include "model.h"
#include "mesh.h"

icy::Model::Model()
{
    mesh = new icy::Mesh();
    Reset(prms);
}

icy::Model::~Model()
{
    delete mesh;
}

void icy::Model::Reset(SimParams &prms)
{
    mesh->Reset(prms.CharacteristicLength, prms.InteractionDistance);
    UnsafeUpdateGeometry();
}

void icy::Model::Prepare(void)
{
    abortRequested = false;
//    timeStepFactor = 1;
}

bool icy::Model::Step(void)
{
    int iter, attempt = 0;
    bool converges=false;
    bool sln_res, ccd_res; // false if matrix is not PSD
    double h;
    do
    {
        iter = 0;
        h = prms.InitialTimeStep*timeStepFactor; // time step
        InitialGuess(prms, h, timeStepFactor);
        std::pair<bool, double> ccd_result = mesh->EnsureNoIntersectionViaCCD();
        ccd_res = ccd_result.first;
        std::cout << std::scientific << std::setprecision(1);
        std::cout << "\nSTEP: " << currentStep << "-" << attempt << " TCF " << timeStepFactor << std::endl;
        sln_res=true;

        while(ccd_res && sln_res && iter < prms.MaxIter && (iter < prms.MinIter || !converges))
        {
            if(abortRequested) {Aborting(); return false;}
            sln_res = AssembleAndSolve(prms, h);
            ccd_result = mesh->EnsureNoIntersectionViaCCD();
            ccd_res = ccd_result.first;

            double ratio = iter == 0 ? 0 : eqOfMotion.solution_norm/eqOfMotion.solution_norm_prev;
            converges = (eqOfMotion.solution_norm < prms.ConvergenceCutoff || ratio < prms.ConvergenceEpsilon);

            std::cout << "IT "<< std::left << std::setw(2) << iter;
            std::cout << " obj " << std::setw(10) << eqOfMotion.objective_value;
            std::cout << " sln " << std::setw(10) << eqOfMotion.solution_norm;
            if(iter) std::cout << " ra " << std::setw(10) << ratio;
            else std::cout << "tsf " << std::setw(20) << timeStepFactor;
            std::cout << std::endl;
            iter++;
        }

        if(!ccd_res)
        {
            qDebug() << "intersection detected";
            attempt++;
            timeStepFactor*=(ccd_result.second*0.8);
        }
        else if(!sln_res)
        {
            qDebug() << "could not solve";
            attempt++;
            timeStepFactor*=0.5;
        }
        else if(!converges)
        {
            qDebug() << "sln did not converge";
            attempt++;
            timeStepFactor*=0.5;
        }
        if(attempt > 20) throw std::runtime_error("could not solve");
    } while (!ccd_res || !sln_res || !converges);

    if(timeStepFactor < 1) timeStepFactor *= 1.2;
    if(timeStepFactor > 1) timeStepFactor=1;
    // accept step
    AcceptTentativeValues(h);
    currentStep++;

    emit stepCompleted();
    return(currentStep < prms.MaxSteps);
}

void icy::Model::RequestAbort(void)
{
    abortRequested = true;
}

void icy::Model::Aborting()
{
    //perform any cleanup if step was aborted
    qDebug() << "icy::ModelController::Aborting()";
    abortRequested = false;
    emit stepAborted();
}

void icy::Model::UnsafeUpdateGeometry()
{
    vtk_update_mutex.lock();
    mesh->UnsafeUpdateGeometry();
    vtk_update_mutex.unlock();
}

void icy::Model::ChangeVisualizationOption(icy::Model::VisOpt option)
{
    vtk_update_mutex.lock();
    mesh->ChangeVisualizationOption((int)option);
    vtk_update_mutex.unlock();
}

void icy::Model::PositionIndenter(double offset)
{
    vtk_update_mutex.lock();
    unsigned n = mesh->indenter.nodes.size();
    Eigen::Vector2d y_direction = Eigen::Vector2d(0,-1.0);
    for(unsigned i=0;i<n;i++)
    {
        icy::Node &nd = mesh->indenter.nodes[i];
        nd.intended_position = nd.x_initial + offset*y_direction;
    }
    vtk_update_mutex.unlock();
}


void icy::Model::InitialGuess(SimParams &prms, double timeStep, double timeStepFactor)
{
    std::size_t nNodes = mesh->allNodes.size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh->allNodes[i];
        if(nd->pinned)
        {
            nd->xt = (timeStepFactor)*nd->intended_position + (1-timeStepFactor)*nd->xn;
            nd->vn=Eigen::Vector2d::Zero();
            nd->x_hat = nd->xn;
        }
        else
        {
            nd->x_hat = nd->xn + timeStep*nd->vn;
            nd->x_hat.y() -= prms.Gravity*timeStep*timeStep;
            nd->xt = nd->xn + timeStep*nd->vn;
        }
    }
}


bool icy::Model::AssembleAndSolve(SimParams &prms, double timeStep)
{
    eqOfMotion.ClearAndResize(mesh->freeNodeCount);

    unsigned nElems = mesh->allElems.size();
    unsigned nNodes = mesh->allNodes.size();

#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++) mesh->allElems[i]->AddToSparsityStructure(eqOfMotion);

    mesh->UpdateTree(prms.InteractionDistance);
    vtk_update_mutex.lock();
    mesh->DetectContactPairs(prms.InteractionDistance);
    vtk_update_mutex.unlock();

    unsigned nInteractions = mesh->collision_interactions.size();

#pragma omp parallel for
    for(unsigned i=0;i<nInteractions;i++)
        mesh->collision_interactions[i].AddToSparsityStructure(eqOfMotion);

    eqOfMotion.CreateStructure();


    // assemble
    bool mesh_iversion_detected = false;
#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++)
    {
        bool elem_entry_ok = mesh->allElems[i]->ComputeEquationEntries(eqOfMotion, prms, timeStep);
        if(!elem_entry_ok) mesh_iversion_detected = true;
    }

    if(mesh_iversion_detected) return false; // mesh inversion

#pragma omp parallel for
    for(unsigned i=0;i<nNodes;i++) mesh->allNodes[i]->ComputeEquationEntries(eqOfMotion, prms, timeStep);

#pragma omp parallel for
    for(unsigned i=0;i<nInteractions;i++) mesh->collision_interactions[i].Evaluate(eqOfMotion, prms, timeStep);

    // solve
    bool solve_result = eqOfMotion.Solve();
    if(!solve_result) return false;

    // pull
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh->allNodes[i];
        Eigen::Vector2d delta_x;
        if(!nd->pinned)
        {
            eqOfMotion.GetTentativeResult(nd->eqId, delta_x);
            nd->xt+=delta_x;
        }
    }

    // for visualization
    vtk_update_mutex.lock();
#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++) mesh->allElems[i]->EvaluateVelocityDivergence();
    vtk_update_mutex.unlock();

    return true;
}

void icy::Model::AcceptTentativeValues(double timeStep)
{
    vtk_update_mutex.lock();
    unsigned nNodes = mesh->allNodes.size();
#pragma omp parallel for
    for(unsigned i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh->allNodes[i];
        Eigen::Vector2d dx = nd->xt-nd->xn;
        nd->vn = dx/timeStep;
        nd->xn = nd->xt;
    }
    vtk_update_mutex.unlock();
}

