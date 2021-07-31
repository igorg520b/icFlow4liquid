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
    currentStep = 0;
    timeStepFactor = 1;
    simulationTime = 0;
}

void icy::Model::Prepare()
{
    abortRequested = false;
}

bool icy::Model::Step()
{
    mesh->UpdateTree(prms.InteractionDistance);

    int attempt = 0;
    bool converges=false;
    bool sln_res, ccd_res; // false if matrix is not PSD
    double h;
    do
    {
        int iter = 0;
        h = prms.InitialTimeStep*timeStepFactor; // time step
        InitialGuess(prms, h, timeStepFactor);
        std::pair<bool, double> ccd_result = mesh->EnsureNoIntersectionViaCCD();
        ccd_res = ccd_result.first;
        std::cout << std::scientific << std::setprecision(1);
        std::cout << "STEP: " << currentStep << "-" << attempt << " TCF " << timeStepFactor << std::endl;
        sln_res=true;
        double first_solution_norm = 0;

        while(ccd_res && sln_res && iter < prms.MaxIter && (iter < prms.MinIter || !converges))
        {
            if(abortRequested) {Aborting(); return false;}
            sln_res = AssembleAndSolve(prms, h);
            ccd_result = mesh->EnsureNoIntersectionViaCCD();
            ccd_res = ccd_result.first;

            double ratio = 0;
            if(iter == 0) { first_solution_norm = eqOfMotion.solution_norm; }
            else if(first_solution_norm > prms.ConvergenceCutoff) ratio = eqOfMotion.solution_norm/first_solution_norm;
            converges = (eqOfMotion.solution_norm < prms.ConvergenceCutoff ||
                         (ratio > 0 && ratio < prms.ConvergenceEpsilon));

            std::cout << "IT "<< std::left << std::setw(2) << iter;
            std::cout << " obj " << std::setw(10) << eqOfMotion.objective_value;
            std::cout << " sln " << std::setw(10) << eqOfMotion.solution_norm;
            if(iter) std::cout << " ra " << std::setw(10) << ratio;
            else std::cout << "tsf " << std::setw(20) << timeStepFactor;
            std::cout << '\n';
            iter++;
        }

        if(!ccd_res)
        {
            std::cout << "intersection detected - discarding this attempt\n";
            attempt++;
            timeStepFactor*=(ccd_result.second*0.8);
        }
        else if(!sln_res)
        {
            std::cout << "could not solve - discarding this attempt\n";
            attempt++;
            timeStepFactor*=0.5;
        }
        else if(!converges)
        {
            std::cout << "did not converge - discarding this attempt\n";
            attempt++;
            timeStepFactor*=0.5;
        }
        if(attempt > 20) throw std::runtime_error("could not solve");
        std::cout << std::endl;
    } while (!ccd_res || !sln_res || !converges);

    // accept step
    bool plasticDeformation = AcceptTentativeValues(h);
    if(plasticDeformation) GetNewMaterialPosition();
    currentStep++;

    // gradually increase the time step
    if(timeStepFactor < 1) timeStepFactor *= 1.2;
    if(timeStepFactor > 1) timeStepFactor=1;

    emit stepCompleted();
    return(currentStep < prms.MaxSteps);
}

void icy::Model::RequestAbort(void) { abortRequested = true; }

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
        icy::Node* nd = mesh->indenter.nodes[i];
        nd->intended_position = nd->x_initial + offset*y_direction;
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


bool icy::Model::AssembleAndSolve(SimParams &prms, double timeStep, bool restShape)
{
    eqOfMotion.ClearAndResize(mesh->freeNodeCount);

    unsigned nElems = mesh->allElems.size();
    unsigned nNodes = mesh->allNodes.size();

#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++) mesh->allElems[i]->AddToSparsityStructure(eqOfMotion);

    unsigned nInteractions;
    if(!restShape)
    {
        mesh->UpdateTree(prms.InteractionDistance);
        vtk_update_mutex.lock();
        mesh->DetectContactPairs(prms.InteractionDistance);
        vtk_update_mutex.unlock();
        nInteractions = mesh->collision_interactions.size();

#pragma omp parallel for
        for(unsigned i=0;i<nInteractions;i++) mesh->collision_interactions[i].AddToSparsityStructure(eqOfMotion);
    }

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
/*
#pragma omp parallel for
    for(unsigned i=0;i<nNodes;i++)
    {
        Node *nd = mesh->allNodes[i];
        nd->ComputeEquationEntries(eqOfMotion, prms, timeStep);
    }
*/
    if(!restShape)
    {
#pragma omp parallel for
        for(unsigned i=0;i<nNodes;i++)
        {
            Node *nd = mesh->allNodes[i];
            nd->AddSpringEntries(eqOfMotion, prms, timeStep, spring);
        }

#pragma omp parallel for
        for(unsigned i=0;i<nInteractions;i++) mesh->collision_interactions[i].Evaluate(eqOfMotion, prms, timeStep);
    }



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

    return true;
}

bool icy::Model::AcceptTentativeValues(double timeStep)
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

    // plastic behavior
    unsigned nElems = mesh->allElems.size();
    bool plasticDeformation = false;
#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++)
    {
        icy::Element *elem = mesh->allElems[i];
        bool result = elem->PlasticDeformation(prms, timeStep);
        if(result) plasticDeformation = true;
        elem->ComputeVisualizedVariables();
    }

    simulationTime+=timeStep;
    mesh->area_current = std::accumulate(mesh->allElems.begin(), mesh->allElems.end(),0.0,
                                         [](double a, Element* m){return a+m->area_current;});
    return plasticDeformation;
}

void icy::Model::GetNewMaterialPosition()
{
    // relax each mesh fragment separately
    for(MeshFragment *mf : mesh->allFragments)
    {
        if(!mf->deformable) continue;

        unsigned freeNodes = 0;
        for(Node *nd : mf->nodes)
        {
            if(nd->pinned) nd->eqId = -1;
            else nd->eqId = freeNodes++;
        }

        unsigned nNodes = mf->nodes.size();
        unsigned nElems = mf->elems.size();


        // solve the equation iteratively (as usual, but without collisions)

        int attempt = 0;
        bool converges=false;
        bool sln_res; // false if matrix is not PSD
        double h = prms.InitialTimeStep*timeStepFactor;;
        do
        {
            int iter = 0;
            h = prms.InitialTimeStep*timeStepFactor; // time step

            // initial coordinates -> tentative
#pragma omp parallel for
            for(std::size_t i=0;i<nNodes;i++)
            {
                icy::Node *nd = mf->nodes[i];
                nd->xt = nd->x_hat = nd->x_initial;
            }

            std::cout << std::scientific << std::setprecision(1);
            std::cout << "\nRELAXING STEP: " << attempt << " h " << h << std::endl;
            sln_res=true;
            double first_solution_norm = 0;

            while(sln_res && iter < prms.MaxIter*5 && (iter < prms.MinIter || !converges))
            {
                sln_res = AssembleAndSolve(prms, h, true);

                double ratio = 0;
                if(iter == 0) { first_solution_norm = eqOfMotion.solution_norm; }
                else if(first_solution_norm > prms.ConvergenceCutoff) ratio = eqOfMotion.solution_norm/first_solution_norm;
                converges = (eqOfMotion.solution_norm < prms.ConvergenceCutoff ||
                             (ratio > 0 && ratio < prms.ConvergenceEpsilon));

                std::cout << "IT "<< std::left << std::setw(2) << iter;
                std::cout << " obj " << std::setw(10) << eqOfMotion.objective_value;
                std::cout << " sln " << std::setw(10) << eqOfMotion.solution_norm;
                if(iter) std::cout << " ra " << std::setw(10) << ratio;
                else std::cout << "h " << std::setw(20) << h;
                std::cout << '\n';
                iter++;
            }

            if(!sln_res)
            {
                std::cout << "could not solve - discarding this attempt\n";
                attempt++;
                h*=0.5;
            }
            else if(!converges)
            {
                std::cout << "did not converge - discarding this attempt\n";
                attempt++;
                h*=0.5;
            }
            if(attempt > 20) throw std::runtime_error("relaxation step: could not solve");
            std::cout << std::endl;
        } while (!sln_res || !converges);


        // pull
#pragma omp parallel for
        for(std::size_t i=0;i<nNodes;i++)
        {
            icy::Node *nd = mf->nodes[i];
            Eigen::Vector2d delta_x;
            if(!nd->pinned)
            {
                eqOfMotion.GetTentativeResult(nd->eqId, delta_x);
                nd->xt+=delta_x;
            }
        }

        // infer new PiMultipliers
#pragma omp parallel for
        for(unsigned i=0;i<nElems;i++)
        {
            icy::Element *elem = mf->elems[i];
            Eigen::Matrix2d DmPrime, Dm;
            DmPrime << elem->nds[0]->xt-elem->nds[2]->xt, elem->nds[1]->xt-elem->nds[2]->xt;
            Dm << elem->nds[0]->x_initial-elem->nds[2]->x_initial, elem->nds[1]->x_initial-elem->nds[2]->x_initial;
            elem->PiMultiplier = DmPrime*Dm.inverse()*elem->PiMultiplier;
        }
        // tentative -> initial
#pragma omp parallel for
        for(std::size_t i=0;i<nNodes;i++)
        {
            icy::Node *nd = mf->nodes[i];
            if(!nd->pinned)
            {
                nd->x_initial = nd->xt;
            }
            nd->area = 0;
        }

        for(unsigned i=0;i<nElems;i++)
        {
            icy::Element *elem = mf->elems[i];
            elem->PrecomputeInitialArea();
            for(int j=0;j<3;j++) elem->nds[j]->area += elem->area_initial/3;
        }
    }

    mesh->freeNodeCount=0;
    for(unsigned i=0;i<mesh->allNodes.size();i++)
    {
        Node *nd = mesh->allNodes[i];
        if(nd->pinned) nd->eqId=-1;
        else nd->eqId=mesh->freeNodeCount++;
    }
}



void icy::Model::AttachSpring(double X, double Y, double radius)
{
    std::cout << "icy::Model::AttachSpring " << X << "; " << Y << "; " << radius << std::endl;
    Eigen::Vector2d attachmentPos(X,Y);
    vtk_update_mutex.lock();
#pragma omp parallel for
    for(unsigned i=0;i<mesh->allNodes.size();i++)
    {
        icy::Node *nd = mesh->allNodes[i];
        if((nd->xn-attachmentPos).norm()<radius)
        {
            nd->spring_attached = 1;
            nd->spring_attachment_position = nd->xn;
        }
        else
        {
            nd->spring_attached = 0;
        }
    }
    vtk_update_mutex.unlock();
}

void icy::Model::ReleaseSpring()
{
    std::cout << "icy::Model::ReleaseSpring " << std::endl;
    vtk_update_mutex.lock();
#pragma omp parallel for
    for(unsigned i=0;i<mesh->allNodes.size();i++)
    {
        icy::Node *nd = mesh->allNodes[i];
        nd->spring_attached=0;
    }
    vtk_update_mutex.unlock();
}

void icy::Model::AdjustSpring(double dX, double dY)
{
    std::cout << "icy::Model::AdjustSpring " << dX << "; " << dY << std::endl;
    spring << dX,dY;
}

