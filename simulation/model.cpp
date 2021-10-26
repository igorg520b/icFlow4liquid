#include <vtkPointData.h>
#include <QtGlobal>
#include "model.h"
#include "mesh.h"
#include "spdlog/spdlog.h"

// CONTROLLER - RESET, STEP, ABORT

void icy::Model::Reset(unsigned setup)
{
    currentStep = 0;
    timeStepFactor = 1;
    simulationTime = 0;

    mesh.Reset(prms.CharacteristicLength, prms.InteractionDistance, setup);
    topologyInvalid = true;
    spdlog::info("icy::Model::Reset() done");
}

void icy::Model::Prepare()
{
    abortRequested = false;
}

bool icy::Model::Step()
{
    if(prms.EnableCollisions)
    {
        mesh.UpdateTree(prms.InteractionDistance);
    }

    int attempt = 0;
    bool converges=false;
    bool sln_res, ccd_res=true; // false if matrix is not PSD
    double h;
    do
    {
        int iter = 0;
        h = prms.InitialTimeStep*timeStepFactor; // time step
        InitialGuess(h, timeStepFactor);
        if(prms.EnableCollisions) ccd_res = mesh.EnsureNoIntersectionViaCCD();
        spdlog::info("\n┌{2:─^{1}}┬{3:─^{1}}┬{4:─^{1}}┬{5:─^{1}}┬{6:─^{1}}┐ ST {7:>}-{8:<2}"
                     ,"",colWidth, " it "," obj "," sln "," tsf "," ra ", currentStep,attempt);

        sln_res=true;
        double first_solution_norm = 0;

        while((!prms.EnableCollisions || ccd_res) && sln_res && iter < prms.MaxIter && (iter < prms.MinIter || !converges))
        {
            if(abortRequested) {Aborting(); return false;}
            sln_res = AssembleAndSolve(h, prms.EnableCollisions, true, mesh.allNodes, mesh.allElems, mesh.allCZs);
            if(prms.EnableCollisions) ccd_res = mesh.EnsureNoIntersectionViaCCD();

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

        if(!ccd_res || !sln_res || !converges)
        {
            spdlog::info("discarding attempt {}; intersection {}; PSD {}; converges {}",attempt, ccd_res,sln_res,converges);
            attempt++;
            timeStepFactor*=0.5;
        }
        if(attempt > 20) throw std::runtime_error("Model::Step() could not solve");
    } while (!ccd_res || !sln_res || !converges);

    // accept step
    bool plasticDeformation = AcceptTentativeValues(h);

    if(plasticDeformation) GetNewMaterialPosition();
    Fracture(h);
    currentStep++;

    // gradually increase the time step
    if(timeStepFactor < 1) timeStepFactor *= 1.2;
    if(timeStepFactor > 1) timeStepFactor = 1;

    emit stepCompleted();
    return(currentStep < prms.MaxSteps);
}

void icy::Model::RequestAbort(void) { abortRequested = true; }

void icy::Model::Aborting() { abortRequested = false; emit stepAborted(); }



// CONTROLLER - VTK VISUALIZATION

void icy::Model::UnsafeSynchronizeVTK()
{
    vtk_update_mutex.lock();
    {
        if(topologyInvalid) mesh.RegenerateVisualizedGeometry();
        else if(displacementsInvalid) mesh.UnsafeUpdateGeometry();
    }
    vtk_update_mutex.unlock();
    topologyInvalid = false;
    displacementsInvalid = false;
}

void icy::Model::ChangeVisualizationOption(icy::Model::VisOpt option)
{
    UnsafeSynchronizeVTK();
    vtk_update_mutex.lock();
    mesh.ChangeVisualizationOption((int)option);
    vtk_update_mutex.unlock();
}



// CONTROLLER - GUI SPRING AND INDENTER

void icy::Model::SetIndenterPosition(double position)
{
    vtk_update_mutex.lock();
    mesh.SetIndenterPosition(position);
    vtk_update_mutex.unlock();
    displacementsInvalid = true;
}

void icy::Model::AttachSpring(double X, double Y, double radius)
{
    spdlog::info("icy::Model::AttachSpring ({},{}); radius {}",X,Y,radius);
    spring = Eigen::Vector2d::Zero();
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
    spdlog::info("icy::Model::ReleaseSpring");
    vtk_update_mutex.lock();
    for(unsigned i=0;i<mesh.allNodes.size();i++) mesh.allNodes[i]->spring_attached=0;
    vtk_update_mutex.unlock();
}

void icy::Model::AdjustSpring(double dX, double dY)
{
    spdlog::info("icy::Model::AdjustSpring {}-{}",dX,dY);
    vtk_update_mutex.lock();
    spring << dX,dY;
    vtk_update_mutex.unlock();
}



// SOLVING - ELASTICITY AND FRACTURE MODEL

void icy::Model::InitialGuess(double timeStep, double timeStepFactor)
{
    std::size_t nNodes = mesh.allNodes.size();
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = mesh.allNodes[i];
        if(nd==nullptr) { spdlog::critical("InitialGuess: nd==nullptr for i=={}",i); throw std::runtime_error("InitialGuess: nd==nullptr");}
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

bool icy::Model::AssembleAndSolve(double timeStep, bool enable_collisions, bool enable_spring,
                                  std::vector<icy::Node*> &nodes, std::vector<icy::Element*> &elems,
                                  std::vector<icy::CohesiveZone*> &czs)
{
    for(Node *nd : mesh.allNodes) nd->eqId = -1;

    unsigned nElems = elems.size();
    unsigned nNodes = nodes.size();
    unsigned nCzs = czs.size();

    // assign sequential indices to free nodes
    unsigned freeNodeCount = 0;
    for(unsigned i=0;i<nNodes;i++)
        if(!nodes[i]->pinned) nodes[i]->eqId = freeNodeCount++;

    eqOfMotion.ClearAndResize(freeNodeCount);

#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++) elems[i]->AddToSparsityStructure(eqOfMotion);

    if(prms.EnableCZs)
    {
#pragma omp parallel for
        for(unsigned i=0;i<nCzs;i++) czs[i]->AddToSparsityStructure(eqOfMotion);
    }

    if(enable_collisions)
    {
        mesh.UpdateTree(prms.InteractionDistance);
        vtk_update_mutex.lock();
        mesh.DetectContactPairs(prms.InteractionDistance);
        vtk_update_mutex.unlock();
#pragma omp parallel for
        for(unsigned i=0;i<mesh.contacts_final_list.size();i++) mesh.contacts_final_list[i].AddToSparsityStructure(eqOfMotion);
    }

    eqOfMotion.CreateStructure();

    // assemble
    bool mesh_inversion_detected = false;
#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++)
    {
        bool elem_entry_ok = elems[i]->ComputeEquationEntries(eqOfMotion, prms, timeStep);
        if(!elem_entry_ok) mesh_inversion_detected = true;
    }
    if(mesh_inversion_detected)
    {
        spdlog::info("mesh inversion detected while assembling elements");
        return false; // mesh inversion
    }

    if(enable_spring)
    {
#pragma omp parallel for
        for(unsigned i=0;i<nNodes;i++) nodes[i]->AddSpringEntries(eqOfMotion, prms, timeStep, spring);
    }

    if(enable_collisions)
    {
#pragma omp parallel for
        for(unsigned i=0;i<mesh.contacts_final_list.size();i++) mesh.contacts_final_list[i].Evaluate(eqOfMotion, prms, timeStep);
    }

    if(prms.EnableCZs)
    {
        bool non_smooth_loading_detected = false;
#pragma omp parallel for
        for(unsigned i=0;i<nCzs;i++)
        {
            bool result = czs[i]->ComputeEquationEntries(eqOfMotion, prms, timeStep);
            if(!result) non_smooth_loading_detected = true;
        }
        if(non_smooth_loading_detected)
        {
            spdlog::info("non-smooth cz loading");
            return false; // mesh inversion
        }
    }


    // solve
    bool solve_result = eqOfMotion.Solve();
    if(!solve_result) return false;

    // pull into Node::xt
#pragma omp parallel for
    for(std::size_t i=0;i<nNodes;i++)
    {
        icy::Node *nd = nodes[i];
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
#pragma omp parallel for
    for(unsigned i=0;i<mesh.allCZs.size();i++)
    {
        CohesiveZone *cz = mesh.allCZs[i];
        cz->AcceptValues();
    }

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

    // remove failed czs from the list
    auto result = std::remove_if(mesh.allCZs.begin(),mesh.allCZs.end(),[](CohesiveZone *cz){return !cz->isActive;});

    if(result != mesh.allCZs.end())
    {
//        for(auto iter = result; iter!=mesh.allCZs.end(); iter++) (*iter)->Disconnect();     // "inform" the nodes that CZ disappeared
        mesh.allCZs.erase(result,mesh.allCZs.end());
        topologyInvalid = true;
    }
    vtk_update_mutex.unlock();



    // plastic behavior
    unsigned nElems = mesh.allElems.size();
    bool plasticDeformation = false;

    if(prms.EnablePlasticity)
    {
#pragma omp parallel for
        for(unsigned i=0;i<nElems;i++)
        {
            icy::Element *elem = mesh.allElems[i];
            bool result = elem->PlasticDeformation(prms, timeStep);
            if(result) plasticDeformation = true;
            elem->ComputeVisualizedVariables();
        }
    }
    else
    {
#pragma omp parallel for
        for(unsigned i=0;i<nElems;i++) mesh.allElems[i]->ComputeVisualizedVariables();
    }

    simulationTime+=timeStep;
    mesh.area_current = std::accumulate(mesh.allElems.begin(), mesh.allElems.end(),0.0,
                                         [](double a, Element* m){return a+m->area_current;});
    displacementsInvalid = true;
    return plasticDeformation;
}

void icy::Model::GetNewMaterialPosition()
{
    // relax each mesh fragment separately
    for(auto &mf : mesh.fragments)
    {
        if(!mf->isDeformable) continue;

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
        double h = prms.InitialTimeStep*timeStepFactor;
        do
        {
            int iter = 0;

            // initial coordinates -> tentative
#pragma omp parallel for
            for(std::size_t i=0;i<nNodes;i++)
            {
                icy::Node *nd = mf->nodes[i];
                nd->xt = nd->x_hat = nd->x_initial;
            }

            spdlog::info("\n╔{0:─^{1}}┬{0:─^{1}}┬{0:─^{1}}┬{0:─^{1}}┬{0:─^{1}}┐  R-STEP: {7:<3}\n"
                         "│{2: ^{1}}│{3: ^{1}}│{4: ^{1}}│{5: ^{1}}│{6: ^{1}}│","",colWidth, "it","obj","sln","h","ra", attempt);

            sln_res=true;
            double first_solution_norm = 0;

            while(sln_res && iter < prms.MaxIter*5 && (iter < prms.MinIter || !converges))
            {
                sln_res = AssembleAndSolve(h, false, false, mesh.allNodes, mesh.allElems, mesh.allCZs);

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
            Eigen::Matrix2d DmPrime;
            DmPrime << elem->nds[0]->xt-elem->nds[2]->xt, elem->nds[1]->xt-elem->nds[2]->xt;
            elem->PiMultiplier = DmPrime*elem->DmInv*elem->PiMultiplier;
        }
        // tentative -> initial
#pragma omp parallel for
        for(std::size_t i=0;i<nNodes;i++)
        {
            icy::Node *nd = mf->nodes[i];
            if(!nd->pinned) nd->x_initial = nd->xt;
        }

        mf->PostMeshingEvaluations();
    }
}

void icy::Model::Fracture(double timeStep)
{
    vtk_update_mutex.lock();
    mesh.ComputeFractureDirections(prms, timeStep, true);   // even if fracture is disabled, compute for visualization
    vtk_update_mutex.unlock();

    if(!prms.EnableFracture) return;

    mesh.updateMinMax = false;  // temporarily disable the adaptive scale when visualizing variables (to avoid flickering)
    int fracture_step_count = 0;

    while(mesh.maxNode != nullptr && fracture_step_count < prms.FractureMaxSubsteps && !abortRequested)
    {
        spdlog::info("Step {}; FR-Step {}; maxNode {}", this->currentStep, fracture_step_count, mesh.maxNode->globId);

        // perform the FractureStep - identify the fracturing node, change topology and do local relaxation
        vtk_update_mutex.lock();
        mesh.PropagateCrack(prms);
        mesh.CreateLeaves();
        topologyInvalid = true;
        vtk_update_mutex.unlock();

        Fracture_LocalSubstep();

        vtk_update_mutex.lock();
        mesh.ComputeFractureDirections(prms, 0, false);
        vtk_update_mutex.unlock();

        fracture_step_count++;
        emit fractureProgress();    // if needed, update the VTK representation and render from the main thread
    }
    mesh.updateMinMax = true;


    if(fracture_step_count > 0)
    {
        // any work that occurs after fracture, e.g. identify separated fragments
    }
}

void icy::Model::Fracture_LocalSubstep()
{
    mesh.ResetFractureTimer(prms);
    mesh.InferLocalSupport(prms);

    if(mesh.local_support.size() == 0) throw std::runtime_error("Fracture_LocalSubstep: local node range is zero");

    double substepping_timestep_factor = 0.1;
    int attempt = 0;
    bool converges=false;
    bool sln_res, ccd_res=true; // false if matrix is not PSD
    double h;
    if(prms.EnableCollisions) mesh.UpdateTree(prms.InteractionDistance);

    do
    {
        int iter = 0;
        h = prms.InitialTimeStep*timeStepFactor*substepping_timestep_factor; // time step
        InitialGuess(h, timeStepFactor*substepping_timestep_factor);
        if(prms.EnableCollisions) ccd_res = mesh.EnsureNoIntersectionViaCCD();
        spdlog::info("\n┏{2:━^{1}}┬{3:━^{1}}┬{4:━^{1}}┬{5:━^{1}}┬{6:━^{1}}┐ FR-ST {7:>}-{8:<2}"
                     ,"",colWidth, " it "," obj "," sln "," tsf "," ra ", currentStep,attempt);

        sln_res=true;
        double first_solution_norm = 0;

        while((!prms.EnableCollisions || ccd_res) && sln_res && iter < prms.MaxIter && (iter < prms.MinIter || !converges))
        {
            if(abortRequested) {Aborting(); return;}
            sln_res = AssembleAndSolve(h, prms.EnableCollisions, true, mesh.local_support, mesh.local_elems, mesh.local_czs);
            if(prms.EnableCollisions) ccd_res = mesh.EnsureNoIntersectionViaCCD();

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

        if(!ccd_res || !sln_res || !converges)
        {
            substepping_timestep_factor*=0.2;
            spdlog::info("discarding attempt {}; ccd_res {}; sln_res {}; converges {}; stf {:.3e}",
                         attempt, ccd_res, sln_res, converges, substepping_timestep_factor);
            attempt++;
        }
        if(attempt > 20) throw std::runtime_error("Fracture_LocalSubstep: could not solve");
    } while (!ccd_res || !sln_res || !converges);

#pragma omp parallel for
        for(unsigned i=0;i<mesh.local_elems.size();i++) mesh.local_elems[i]->ComputeVisualizedVariables();
}

