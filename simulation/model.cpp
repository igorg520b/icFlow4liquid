#include <vtkPointData.h>
#include <QtGlobal>
#include "model.h"
#include "mesh.h"
#include "spdlog/spdlog.h"


void icy::Model::Reset(unsigned setup)
{
    currentStep = 0;
    timeStepFactor = 1;
    simulationTime = 0;

    vtk_update_mutex.lock();
    mesh.Reset(prms.CharacteristicLength, prms.InteractionDistance, setup);
    vtk_update_mutex.unlock();
    std::clog << "icy::Model::Reset() done\n";
}

void icy::Model::Prepare()
{
    abortRequested = false;
}

bool icy::Model::Step()
{
    mesh.UpdateTree(prms.InteractionDistance);

    int attempt = 0;
    bool converges=false;
    bool sln_res, ccd_res; // false if matrix is not PSD
    double h;
    do
    {
        int iter = 0;
        h = prms.InitialTimeStep*timeStepFactor; // time step
        InitialGuess(h, timeStepFactor);
        std::pair<bool, double> ccd_result = mesh.EnsureNoIntersectionViaCCD();
        ccd_res = ccd_result.first;
        spdlog::info("\n┌{2:─^{1}}┬{3:─^{1}}┬{4:─^{1}}┬{5:─^{1}}┬{6:─^{1}}┐ ST {7:>}-{8:<2}"
                     ,"",colWidth, " it "," obj "," sln "," tsf "," ra ", currentStep,attempt);

        sln_res=true;
        double first_solution_norm = 0;

        while(ccd_res && sln_res && iter < prms.MaxIter && (iter < prms.MinIter || !converges))
        {
            if(abortRequested) {Aborting(); return false;}
            sln_res = AssembleAndSolve(h);
            ccd_result = mesh.EnsureNoIntersectionViaCCD();
            ccd_res = ccd_result.first;

            double ratio = 0;
            if(iter == 0) { first_solution_norm = eqOfMotion.solution_norm; }
            else if(first_solution_norm > prms.ConvergenceCutoff) ratio = eqOfMotion.solution_norm/first_solution_norm;
            converges = (eqOfMotion.solution_norm < prms.ConvergenceCutoff ||
                         (ratio > 0 && ratio < prms.ConvergenceEpsilon));

            spdlog::info("│{2: ^{1}d}│{3: ^{1}.3e}│{4: ^{1}.3e}│{5: ^{1}.3e}│{6: ^{1}.3e}│",
    "",colWidth, iter, eqOfMotion.objective_value,eqOfMotion.solution_norm, timeStepFactor, ratio);
            iter++;
        }
        spdlog::info("└{0:─^{1}}┴{0:─^{1}}┴{0:─^{1}}┴{0:─^{1}}┴{0:─^{1}}┘","",colWidth);

        if(!ccd_res)
        {
            spdlog::info("intersection detected - discarding attempt {}",attempt);
            attempt++;
            timeStepFactor*=(std::max((1-ccd_result.second)*0.8,0.1));
        }
        else if(!sln_res)
        {
            spdlog::info("matrix not PSD - discarding attempt {}",attempt);
            attempt++;
            timeStepFactor*=0.5;
        }
        else if(!converges)
        {
            spdlog::info("did not converge - discarding attempt {}",attempt);
            attempt++;
            timeStepFactor*=0.5;
        }
        if(attempt > 20) throw std::runtime_error("could not solve");
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
    spdlog::info("icy::ModelController::Aborting()");
    abortRequested = false;
    emit stepAborted();
}

void icy::Model::UnsafeUpdateGeometry()
{
    vtk_update_mutex.lock();
    mesh.UnsafeUpdateGeometry();
    vtk_update_mutex.unlock();
}

void icy::Model::ChangeVisualizationOption(icy::Model::VisOpt option)
{
    vtk_update_mutex.lock();
    mesh.ChangeVisualizationOption((int)option);
    vtk_update_mutex.unlock();
}

void icy::Model::SetIndenterPosition(double position)
{
    vtk_update_mutex.lock();
    mesh.SetIndenterPosition(position);
    vtk_update_mutex.unlock();
}


void icy::Model::InitialGuess(double timeStep, double timeStepFactor)
{
    std::size_t nNodes = mesh.allNodes.size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh.allNodes[i];
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


bool icy::Model::AssembleAndSolve(double timeStep, bool restShape)
{
    eqOfMotion.ClearAndResize(mesh.freeNodeCount);

    unsigned nElems = mesh.allElems.size();
    unsigned nNodes = mesh.allNodes.size();

#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++) mesh.allElems[i]->AddToSparsityStructure(eqOfMotion);

    unsigned nInteractions;
    if(!restShape)
    {
        mesh.UpdateTree(prms.InteractionDistance);
        vtk_update_mutex.lock();
        mesh.DetectContactPairs(prms.InteractionDistance);
        vtk_update_mutex.unlock();
        nInteractions = mesh.collision_interactions.size();

#pragma omp parallel for
        for(unsigned i=0;i<nInteractions;i++) mesh.collision_interactions[i].AddToSparsityStructure(eqOfMotion);
    }

    eqOfMotion.CreateStructure();


    // assemble
    bool mesh_iversion_detected = false;
#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++)
    {
        bool elem_entry_ok = mesh.allElems[i]->ComputeEquationEntries(eqOfMotion, prms, timeStep);
        if(!elem_entry_ok) mesh_iversion_detected = true;
    }

    if(mesh_iversion_detected) return false; // mesh inversion

    if(!restShape)
    {
#pragma omp parallel for
        for(unsigned i=0;i<nNodes;i++)
        {
            Node *nd = mesh.allNodes[i];
            nd->AddSpringEntries(eqOfMotion, prms, timeStep, spring);
        }

#pragma omp parallel for
        for(unsigned i=0;i<nInteractions;i++) mesh.collision_interactions[i].Evaluate(eqOfMotion, prms, timeStep);
    }

    // solve
    bool solve_result = eqOfMotion.Solve();
    if(!solve_result) return false;

    // pull
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh.allNodes[i];
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
    unsigned nNodes = mesh.allNodes.size();
#pragma omp parallel for
    for(unsigned i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh.allNodes[i];
        Eigen::Vector2d dx = nd->xt-nd->xn;
        nd->vn = dx/timeStep;
        nd->xn = nd->xt;
    }
    vtk_update_mutex.unlock();

    // plastic behavior
    unsigned nElems = mesh.allElems.size();
    bool plasticDeformation = false;
#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++)
    {
        icy::Element *elem = mesh.allElems[i];
        bool result = elem->PlasticDeformation(prms, timeStep);
        if(result) plasticDeformation = true;
        elem->ComputeVisualizedVariables();
    }

    simulationTime+=timeStep;
    mesh.area_current = std::accumulate(mesh.allElems.begin(), mesh.allElems.end(),0.0,
                                         [](double a, Element* m){return a+m->area_current;});
    return plasticDeformation;
}

void icy::Model::GetNewMaterialPosition()
{
    // relax each mesh fragment separately
    for(MeshFragment &mf : mesh.fragments)
    {
        if(!mf.deformable) continue;

        unsigned freeNodes = 0;
        for(Node *nd : mf.nodes)
        {
            if(nd->pinned) nd->eqId = -1;
            else nd->eqId = freeNodes++;
        }

        unsigned nNodes = mf.nodes.size();
        unsigned nElems = mf.elems.size();


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
                icy::Node *nd = mf.nodes[i];
                nd->xt = nd->x_hat = nd->x_initial;
            }

            spdlog::info("\n╔{0:─^{1}}┬{0:─^{1}}┬{0:─^{1}}┬{0:─^{1}}┬{0:─^{1}}┐  R-STEP: {7:<3}\n"
                         "│{2: ^{1}}│{3: ^{1}}│{4: ^{1}}│{5: ^{1}}│{6: ^{1}}│","",colWidth, "it","obj","sln","h","ra", attempt);

            sln_res=true;
            double first_solution_norm = 0;

            while(sln_res && iter < prms.MaxIter*5 && (iter < prms.MinIter || !converges))
            {
                sln_res = AssembleAndSolve(h, true);

                double ratio = 0;
                if(iter == 0) { first_solution_norm = eqOfMotion.solution_norm; }
                else if(first_solution_norm > prms.ConvergenceCutoff) ratio = eqOfMotion.solution_norm/first_solution_norm;
                converges = (eqOfMotion.solution_norm < prms.ConvergenceCutoff ||
                             (ratio > 0 && ratio < prms.ConvergenceEpsilon));

                spdlog::info("│{2: ^{1}d}│{3: ^{1}.3e}│{4: ^{1}.3e}│{5: ^{1}.3e}│{6: ^{1}.3e}│",
        "",colWidth, iter, eqOfMotion.objective_value,eqOfMotion.solution_norm, h, ratio);

                iter++;
            }
            spdlog::info("└{0:─^{1}}┴{0:─^{1}}┴{0:─^{1}}┴{0:─^{1}}┴{0:─^{1}}┘","",colWidth);

            if(!sln_res)
            {
                spdlog::info("Matrix not PSD - discarding this attempt");
                attempt++;
                h*=0.5;
            }
            else if(!converges)
            {
                spdlog::info("did not converge - discarding this attempt");
                attempt++;
                h*=0.5;
            }
            if(attempt > 20) throw std::runtime_error("relaxation step: could not solve");
        } while (!sln_res || !converges);


        // pull
#pragma omp parallel for
        for(std::size_t i=0;i<nNodes;i++)
        {
            icy::Node *nd = mf.nodes[i];
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
            icy::Element *elem = mf.elems[i];
            Eigen::Matrix2d DmPrime;
            DmPrime << elem->nds[0]->xt-elem->nds[2]->xt, elem->nds[1]->xt-elem->nds[2]->xt;
            elem->PiMultiplier = DmPrime*elem->Dm.inverse()*elem->PiMultiplier;
        }
        // tentative -> initial
#pragma omp parallel for
        for(std::size_t i=0;i<nNodes;i++)
        {
            icy::Node *nd = mf.nodes[i];
            if(!nd->pinned) nd->x_initial = nd->xt;
        }

        mf.PostMeshingEvaluations();
    }

    mesh.freeNodeCount=0;
    for(unsigned i=0;i<mesh.allNodes.size();i++)
    {
        Node *nd = mesh.allNodes[i];
        if(nd->pinned) nd->eqId=-1;
        else nd->eqId=mesh.freeNodeCount++;
    }
}



void icy::Model::AttachSpring(double X, double Y, double radius)
{
    spdlog::debug("icy::Model::AttachSpring ({},{}); radius {}",X,Y,radius);
    Eigen::Vector2d attachmentPos(X,Y);
    vtk_update_mutex.lock();
#pragma omp parallel for
    for(unsigned i=0;i<mesh.allNodes.size();i++)
    {
        icy::Node *nd = mesh.allNodes[i];
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
    spdlog::debug("icy::Model::ReleaseSpring");
    vtk_update_mutex.lock();
#pragma omp parallel for
    for(unsigned i=0;i<mesh.allNodes.size();i++)
    {
        icy::Node *nd = mesh.allNodes[i];
        nd->spring_attached=0;
    }
    vtk_update_mutex.unlock();
}

void icy::Model::AdjustSpring(double dX, double dY)
{
    spdlog::debug("icy::Model::AdjustSpring {}-{}",dX,dY);
    spring << dX,dY;
}

