#include "mesh.h"
#include "model.h"
#include <numeric>
#include <algorithm>
#include <iterator>
#include <unordered_set>
#include <queue>
#include "spdlog/spdlog.h"

void icy::Mesh::Reset(double MeshSizeMax, double offset, unsigned typeOfSetup_)
{
    maxNode = nullptr;
    typeOfSetup = typeOfSetup_;
    fragments.clear();
    movableBoundary.clear();
    movableNodes.clear();
    collision_interactions.clear();

    switch(typeOfSetup)
    {
    // indentation
    case 0:
    {
        fragments.resize(3);
        fragments[0].GenerateIndenter(MeshSizeMax/2, 0, 1+0.15*1.1, 0.15, 2);
        fragments[1].GenerateBrick(MeshSizeMax,2,1);
        fragments[2].GenerateContainer(MeshSizeMax,offset);

        movableBoundary.resize(fragments[0].boundary_edges.size());
        std::copy(fragments[0].boundary_edges.begin(),fragments[0].boundary_edges.end(),movableBoundary.begin());
    }
        break;

        // shear
    case 1:
    {
        const double radius = 0.1;
        const double height = 0.8;
        const double width = 2.0;
        fragments.resize(5);
        fragments[0].GenerateBrick(MeshSizeMax, width, height);
        for(Node *nd : fragments[0].nodes) if(nd->group.test(2)) nd->pinned=true;
        fragments[1].GenerateIndenter(MeshSizeMax/2, width*0.1/2, height+radius+offset*1.5, radius, 1);
        fragments[2].GenerateIndenter(MeshSizeMax/2, -width*0.9/2, height+radius+offset*1.5, radius, 1);
        fragments[3].GenerateIndenter(MeshSizeMax/2, -width*0.1/2, -(radius+offset*1.5), radius, 1);
        fragments[4].GenerateIndenter(MeshSizeMax/2, width*0.9/2, -(radius+offset*1.5), radius, 1);
        std::copy(fragments[1].boundary_edges.begin(),fragments[1].boundary_edges.end(),std::back_inserter(movableBoundary));
    }
        break;

        // stretch
    case 2:
        fragments.resize(1);
        fragments[0].GenerateBrick(MeshSizeMax,2,1);
        for(Node *nd : fragments[0].nodes) if(nd->group.test(0) || nd->group.test(1)) nd->pinned=true;
        movableBoundary.resize(fragments[0].special_boundary.size());
        std::copy(fragments[0].special_boundary.begin(),fragments[0].special_boundary.end(),movableBoundary.begin());
        break;

        // self-collision test
    case 3:
        fragments.resize(1);
        fragments[0].GenerateSelfCollisionBrick(MeshSizeMax,2,1);
        movableBoundary.resize(fragments[0].special_boundary.size());
        std::copy(fragments[0].special_boundary.begin(),fragments[0].special_boundary.end(),movableBoundary.begin());
        break;
    }

    // make a list of nodes that only belong to the movable boudnary
    std::unordered_set<Node*> s;
    for(auto const &p : movableBoundary) { s.insert(p.first); s.insert(p.second); }
    movableNodes.resize(s.size());
    std::copy(s.begin(),s.end(),movableNodes.begin());
    for(unsigned i=0;i<movableNodes.size();i++) movableNodes[i]->indId=i;


    RegenerateVisualizedGeometry();
    ChangeVisualizationOption(VisualizingVariable);

    tree_update_counter=0;
    area_initial = area_current = std::accumulate(allElems.begin(),
                                                  allElems.end(),0.0,
                                                  [](double a, Element* m){return a+m->area_initial;});
    updateMinMax = true;
    std::clog << "icy::Mesh::Reset() done\n";
}

icy::Mesh::Mesh()
{
    // initialize LUT
    const int nLut = 51;
    hueLut->SetNumberOfTableValues(nLut);
    for ( int i=0; i<nLut; i++)
            hueLut->SetTableValue(i, (double)lutArrayTemperatureAdj[i][0],
                    (double)lutArrayTemperatureAdj[i][1],
                    (double)lutArrayTemperatureAdj[i][2], 1.0);

    visualized_values->SetName("visualized_values");

    ugrid_deformable->SetPoints(points_deformable);
    dataSetMapper_deformable->SetInputData(ugrid_deformable);
    dataSetMapper_deformable->UseLookupTableScalarRangeOn();
    dataSetMapper_deformable->SetLookupTable(hueLut);

    actor_mesh_deformable->SetMapper(dataSetMapper_deformable);
    actor_mesh_deformable->GetProperty()->VertexVisibilityOff();
    actor_mesh_deformable->GetProperty()->EdgeVisibilityOn();
    actor_mesh_deformable->GetProperty()->SetColor(0.84938, 0.872213, 0.848103);
    actor_mesh_deformable->GetProperty()->SetEdgeColor(90.0/255.0, 90.0/255.0, 97.0/255.0);
    actor_mesh_deformable->GetProperty()->LightingOff();
    actor_mesh_deformable->GetProperty()->ShadingOff();
    actor_mesh_deformable->GetProperty()->SetInterpolationToFlat();
    actor_mesh_deformable->PickableOff();
    actor_mesh_deformable->GetProperty()->SetLineWidth(1);

    // boundary that encompasses all objects
    ugrid_boundary_all->SetPoints(points_deformable);
    dataSetMapper_boundary_all->SetInputData(ugrid_boundary_all);
    actor_boundary_all->SetMapper(dataSetMapper_boundary_all);
    actor_boundary_all->GetProperty()->EdgeVisibilityOn();
    actor_boundary_all->GetProperty()->VertexVisibilityOn();
    actor_boundary_all->GetProperty()->SetColor(0,0,0.4);
    actor_boundary_all->GetProperty()->SetEdgeColor(0,0,0.4);
    actor_boundary_all->GetProperty()->SetVertexColor(0.4,0,0);
    actor_boundary_all->GetProperty()->SetPointSize(5);
    actor_boundary_all->GetProperty()->SetLineWidth(3);
    actor_boundary_all->PickableOff();

    // indenter's intended position
    ugrid_indenter_intended->SetPoints(points_indenter_intended);
    dataSetMapper_indenter_intended->SetInputData(ugrid_indenter_intended);
    actor_boundary_intended_indenter->SetMapper(dataSetMapper_indenter_intended);
    actor_boundary_intended_indenter->GetProperty()->EdgeVisibilityOn();
    actor_boundary_intended_indenter->GetProperty()->VertexVisibilityOn();
    actor_boundary_intended_indenter->GetProperty()->SetColor(0.3,0,0.4);
    actor_boundary_intended_indenter->GetProperty()->SetEdgeColor(0.6,0,0.04);
    actor_boundary_intended_indenter->GetProperty()->SetVertexColor(0.04,0,0);
    actor_boundary_intended_indenter->GetProperty()->SetPointSize(2);
    actor_boundary_intended_indenter->GetProperty()->SetLineWidth(1);

    // collisions
    ugrid_collisions->SetPoints(points_collisions);
    mapper_collisions->SetInputData(ugrid_collisions);
    actor_collisions->SetMapper(mapper_collisions);
    actor_collisions->GetProperty()->EdgeVisibilityOn();
    actor_collisions->GetProperty()->VertexVisibilityOn();
    actor_collisions->GetProperty()->SetColor(0.3,0,0.4);
    actor_collisions->GetProperty()->SetEdgeColor(0.6,0,0.04);
    actor_collisions->GetProperty()->SetVertexColor(0.04,0,0);
    actor_collisions->GetProperty()->SetPointSize(2);
    actor_collisions->GetProperty()->SetLineWidth(1);

    collision_interactions.reserve(10000);
    broadlist_ccd.reserve(100000);
    broadlist_contact.reserve(100000);

    mesh_root_contact.isLeaf = mesh_root_ccd.isLeaf = false;
    mesh_root_contact.test_self_collision = mesh_root_ccd.test_self_collision = true;
}

void icy::Mesh::SetIndenterPosition(double position)
{
    switch(typeOfSetup)
    {
    case 0:
    case 1:
    {
        Eigen::Vector2d y_direction = Eigen::Vector2d(0,-1.1);
        for(unsigned i=0;i<movableNodes.size();i++)
        {
            icy::Node* nd = movableNodes[i];
            nd->intended_position = nd->x_initial + position*y_direction;
        }
    }
        break;
    case 2:
    case 3:
        Eigen::Vector2d x_direction = Eigen::Vector2d(0.5,0);
        for(unsigned i=0;i<movableNodes.size();i++)
        {
            icy::Node* nd = movableNodes[i];
            nd->intended_position = nd->x_initial + (position-0.3)*x_direction*(1+nd->x_initial.y()/3);
        }
        break;
    }
}



// VTK VISUALIZATION

void icy::Mesh::RegenerateVisualizedGeometry()
{
    tree_update_counter = 0;
    allNodes.clear();
    allElems.clear();
    allBoundaryEdges.clear();
    unsigned count = 0;

    global_leaves_ccd.clear();
    global_leaves_contact.clear();
    fragmentRoots_ccd.clear();
    fragmentRoots_contact.clear();

    for(MeshFragment &mf : fragments)
    {
        for(unsigned i=0;i<mf.nodes.size();i++)
        {
            Node *nd = mf.nodes[i];
            nd->globId = count++;
            allNodes.push_back(nd);
        }
        mf.GenerateLeaves(allBoundaryEdges.size());

        for(unsigned i=0;i<mf.elems.size();i++) allElems.push_back(mf.elems[i]);
        allBoundaryEdges.insert(allBoundaryEdges.end(), mf.boundary_edges.begin(), mf.boundary_edges.end());
        global_leaves_ccd.insert(global_leaves_ccd.end(), mf.leaves_for_ccd.begin(), mf.leaves_for_ccd.end());
        global_leaves_contact.insert(global_leaves_contact.end(),mf.leaves_for_contact.begin(), mf.leaves_for_contact.end());
        fragmentRoots_ccd.push_back(&mf.root_ccd);
        fragmentRoots_contact.push_back(&mf.root_contact);
    }

    points_deformable->SetNumberOfPoints(allNodes.size());
    cellArray_deformable->Reset();

    // deformable material - elements
    for(icy::Element *tr : allElems)
    {
        vtkIdType pts[3] = {tr->nds[0]->globId, tr->nds[1]->globId, tr->nds[2]->globId};
        cellArray_deformable->InsertNextCell(3, pts);
    }
    ugrid_deformable->SetCells(VTK_TRIANGLE, cellArray_deformable);

    // all boundaries
    cellArray_boundary_all->Reset();
    for(auto edge : allBoundaryEdges)
    {
        vtkIdType pts[2] = {edge.first->globId, edge.second->globId};
        cellArray_boundary_all->InsertNextCell(2, pts);
    }
    ugrid_boundary_all->SetCells(VTK_LINE, cellArray_boundary_all);


    // intended position of the indenter
    points_indenter_intended->SetNumberOfPoints(movableNodes.size());
    cellArray_indenter_intended->Reset();
    for(const auto &edge : movableBoundary)
    {
        vtkIdType pts[2] = {edge.first->indId, edge.second->indId};
        cellArray_indenter_intended->InsertNextCell(2, pts);
    }
    ugrid_indenter_intended->SetCells(VTK_LINE, cellArray_indenter_intended);

    UnsafeUpdateGeometry();
}

void icy::Mesh::UnsafeUpdateGeometry()
{
    double x[3]={};

    if(showDeformation == ShowDeformationOption::initial)
    {
        for(icy::Node *nd : allNodes)
        {
            x[0] = nd->x_initial[0];
            x[1] = nd->x_initial[1];
            points_deformable->SetPoint((vtkIdType)nd->globId, x);
        }
    }
    else if(showDeformation == ShowDeformationOption::current)
    {
        for(icy::Node *nd : allNodes)
        {
            x[0] = nd->xn[0];
            x[1] = nd->xn[1];
            points_deformable->SetPoint((vtkIdType)nd->globId, x);
        }
    }
    points_deformable->Modified();

    for(const icy::Node* nd : movableNodes)
    {
        x[0]=nd->intended_position[0];
        x[1]=nd->intended_position[1];
        points_indenter_intended->SetPoint((vtkIdType)nd->indId, x);
    }
    points_indenter_intended->Modified();

    if(VisualizingVariable != icy::Model::VisOpt::none) UpdateValues();


    // collisions
    points_collisions->SetNumberOfPoints(collision_interactions.size()*2);
    cellArray_collisions->Reset();
    if(showDeformation == ShowDeformationOption::current)
    {
        for(unsigned i=0;i<collision_interactions.size();i++)
        {
            vtkIdType pts[2] = {2*i, 2*i+1};
            cellArray_collisions->InsertNextCell(2, pts);
            Interaction &intr = collision_interactions[i];
            x[0] = intr.ndP->xn[0];
            x[1] = intr.ndP->xn[1];
            points_collisions->SetPoint((vtkIdType)2*i, x);
            x[0] = intr.D[0];
            x[1] = intr.D[1];
            points_collisions->SetPoint((vtkIdType)2*i+1, x);
        }
    }
    points_indenter_intended->Modified();
    ugrid_collisions->SetCells(VTK_LINE, cellArray_collisions);
    actor_collisions->Modified();
}

void icy::Mesh::ChangeVisualizationOption(int option)
{
    spdlog::info("icy::Model::ChangeVisualizationOption {}", option);
    VisualizingVariable = option;

    switch(VisualizingVariable)
    {
    case icy::Model::VisOpt::none:
        dataSetMapper_deformable->ScalarVisibilityOff();
        ugrid_deformable->GetPointData()->RemoveArray("visualized_values");
        ugrid_deformable->GetCellData()->RemoveArray("visualized_values");
        return;

    case icy::Model::VisOpt::node_group:
    case icy::Model::VisOpt::vel_mag:
    case icy::Model::VisOpt::adj_elems_count_nd:
    case icy::Model::VisOpt::nd_max_normal_traction:
    case icy::Model::VisOpt::nd_isBoundary:
        ugrid_deformable->GetCellData()->RemoveArray("visualized_values");
        ugrid_deformable->GetPointData()->AddArray(visualized_values);
        ugrid_deformable->GetPointData()->SetActiveScalars("visualized_values");
        dataSetMapper_deformable->SetScalarModeToUsePointData();
        dataSetMapper_deformable->ScalarVisibilityOn();
        break;

    default:
        ugrid_deformable->GetPointData()->RemoveArray("visualized_values");
        ugrid_deformable->GetCellData()->AddArray(visualized_values);
        ugrid_deformable->GetCellData()->SetActiveScalars("visualized_values");
        dataSetMapper_deformable->SetScalarModeToUseCellData();
        dataSetMapper_deformable->ScalarVisibilityOn();
        break;
    }
    UpdateValues();
}

void icy::Mesh::UpdateValues()
{
    if(allNodes.size()==0)
    {
        dataSetMapper_deformable->ScalarVisibilityOff();
        ugrid_deformable->GetPointData()->RemoveArray("visualized_values");
        ugrid_deformable->GetCellData()->RemoveArray("visualized_values");
        return;
    }

    visualized_values->SetNumberOfValues(allElems.size());

    switch((icy::Model::VisOpt)VisualizingVariable)
    {
    // per node
    case icy::Model::VisOpt::node_group:
        visualized_values->SetNumberOfValues(allNodes.size());
        for(std::size_t i=0;i<allNodes.size();i++) visualized_values->SetValue(i, allNodes[i]->group.to_ulong());
        break;

    case icy::Model::VisOpt::vel_mag:
        visualized_values->SetNumberOfValues(allNodes.size());
        for(std::size_t i=0;i<allNodes.size();i++) visualized_values->SetValue(i, allNodes[i]->vn.norm());
        break;

    case icy::Model::VisOpt::adj_elems_count_nd:
        visualized_values->SetNumberOfValues(allNodes.size());
        for(std::size_t i=0;i<allNodes.size();i++) visualized_values->SetValue(i, allNodes[i]->adj_elems.size());
        break;

    case icy::Model::VisOpt::nd_max_normal_traction:
        visualized_values->SetNumberOfValues(allNodes.size());
        for(std::size_t i=0;i<allNodes.size();i++) visualized_values->SetValue(i, allNodes[i]->max_normal_traction);
        break;

    case icy::Model::VisOpt::nd_isBoundary:
        visualized_values->SetNumberOfValues(allNodes.size());
        for(std::size_t i=0;i<allNodes.size();i++) visualized_values->SetValue(i, allNodes[i]->isBoundary ? 1 : 0);
        break;

        // plasticity
    case icy::Model::VisOpt::plasticity_norm:
        for(std::size_t i=0;i<allElems.size();i++)
            visualized_values->SetValue(i, (allElems[i]->PiMultiplier-Eigen::Matrix2d::Identity()).norm());
        break;

    case icy::Model::VisOpt::plasticity_gamma:
        for(std::size_t i=0;i<allElems.size();i++)
            visualized_values->SetValue(i, allElems[i]->plasticity_gamma);
        break;

    case icy::Model::VisOpt::plasticity_tau_ratio:
        for(std::size_t i=0;i<allElems.size();i++)
            visualized_values->SetValue(i, std::clamp(allElems[i]->plasticity_tau_ratio,0.0,1.0));
        break;


        // stress
    case icy::Model::VisOpt::stress_xx:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->CauchyStress(0,0));
        break;

    case icy::Model::VisOpt::stress_yy:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->CauchyStress(1,1));
        break;

    case icy::Model::VisOpt::stress_hydrostatic:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->hydrostatic_stress);
        break;

    case icy::Model::VisOpt::ps1:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->principal_stress1);
        break;

    case icy::Model::VisOpt::ps2:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->principal_stress2);
        break;

    case icy::Model::VisOpt::shear_stress:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->max_shear_stress);
        break;


        // strain
    case icy::Model::VisOpt::Green_strain_xx:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->GreenStrain(0,0));
        break;

    case icy::Model::VisOpt::Green_strain_yy:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->GreenStrain(1,1));
        break;

    case icy::Model::VisOpt::Green_strain_xy:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->GreenStrain(0,1));
        break;


        // mesh quality measures
    case icy::Model::VisOpt::QM1:
        for(std::size_t i=0;i<allElems.size();i++)
            visualized_values->SetValue(i, allElems[i]->quality_measure_Wicke);
        break;

    case icy::Model::VisOpt::avg_edge_len:
        for(std::size_t i=0;i<allElems.size();i++)
        {
            Element *elem = allElems[i];
            double e0 = (elem->nds[1]->x_initial - elem->nds[2]->x_initial).norm();
            double e1 = (elem->nds[0]->x_initial - elem->nds[2]->x_initial).norm();
            double e2 = (elem->nds[1]->x_initial - elem->nds[0]->x_initial).norm();
            visualized_values->SetValue(i, (e0+e1+e2)/3.0);
        }
        break;


    // other per-element values
    case icy::Model::VisOpt::elem_area:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->area_initial);
        break;

    case icy::Model::VisOpt::energy_density:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->strain_energy_density);
        break;


    case icy::Model::VisOpt::volume_change:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->volume_change);
        break;

    case icy::Model::VisOpt::velocity_div:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->velocity_divergence);
        break;

    case icy::Model::VisOpt::elem_group:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->group);
        break;


    default:
        break;
    }

    visualized_values->Modified();

    if(updateMinMax)
    {
        double minmax[2];
        visualized_values->GetValueRange(minmax);
        if(VisualizingVariable == icy::Model::VisOpt::volume_change)
        {
            double range = minmax[1]-minmax[0];
            hueLut->SetTableRange(1-range*0.75,1+range*0.75);
        }
        else if(VisualizingVariable == icy::Model::VisOpt::velocity_div)
        {
            double upper_boundary = std::max(10.0, minmax[1]);
            double lower_boundary = std::min(-10.0, minmax[0]);
            double boundary = std::max(upper_boundary, std::abs(lower_boundary));
            hueLut->SetTableRange(-boundary, boundary);
        }
        else if(VisualizingVariable == icy::Model::VisOpt::stress_hydrostatic)
        {
            double absVal = std::max(std::abs(minmax[0]), std::abs(minmax[1]));
            hueLut->SetTableRange(-absVal, absVal);
        }
        else if(VisualizingVariable == icy::Model::VisOpt::plasticity_norm)
        {
            hueLut->SetTableRange(0, 0.5);
        }
        else if(VisualizingVariable == icy::Model::VisOpt::QM1)
        {
            hueLut->SetTableRange(0, 1);
        }
        else
        {
            hueLut->SetTableRange(minmax[0], minmax[1]);
        }
    }
}



// COLLISION DETECTION

void icy::Mesh::UpdateTree(float distance_threshold)
{
    // update leafs
    unsigned nLeafs = global_leaves_ccd.size();
#pragma omp parallel for
    for(unsigned i=0;i<nLeafs;i++)
    {
        BVHN *leaf_ccd = global_leaves_ccd[i];
        Node *nd1, *nd2;
        std::tie(nd1,nd2) = allBoundaryEdges[leaf_ccd->feature_idx];
        kDOP8 &box_ccd = leaf_ccd->box;
        box_ccd.Reset();
        box_ccd.Expand(nd1->xn[0], nd1->xn[1]);
        box_ccd.Expand(nd2->xn[0], nd2->xn[1]);
        box_ccd.Expand(nd1->xt[0], nd1->xt[1]);
        box_ccd.Expand(nd2->xt[0], nd2->xt[1]);

        BVHN *leaf_contact = global_leaves_contact[i];
        std::tie(nd1,nd2) = allBoundaryEdges[leaf_contact->feature_idx];

        kDOP8 &box_contact = leaf_contact->box;
        box_contact.Reset();
        box_contact.Expand(nd1->xt[0], nd1->xt[1]);
        box_contact.Expand(nd2->xt[0], nd2->xt[1]);
        box_contact.ExpandBy(distance_threshold);
    }

    // update or build the rest of the tree
    // TODO: parallel
    if(tree_update_counter%10 != 0)
    {
        mesh_root_ccd.Update();
        mesh_root_contact.Update();
    }
    else
    {
        BVHN::BVHNFactory.releaseAll(); // does not release the leaves and roots

        if(fragmentRoots_ccd.size()>1)
        {
            for(MeshFragment &mf : fragments)
            {
                mf.root_ccd.Build(&mf.leaves_for_ccd,0);
                mf.root_contact.Build(&mf.leaves_for_contact,0);
            }
            mesh_root_ccd.Build(&fragmentRoots_ccd,0);
            mesh_root_contact.Build(&fragmentRoots_contact,0);
        }
        else
        {
            mesh_root_ccd.Build(&fragments.front().leaves_for_ccd,0);
            mesh_root_contact.Build(&fragments.front().leaves_for_contact,0);
        }
    }
    tree_update_counter++;
}

void icy::Mesh::AddToNarrowListIfNeeded(unsigned edge_idx, unsigned node_idx, double distance_threshold)
{
    Node *ndA, *ndB, *ndP;
    std::tie(ndA,ndB) = allBoundaryEdges[edge_idx];
    ndP = allNodes[node_idx];
    if(ndA==ndP || ndB==ndP) return;    // a node can't collide with itself
    Eigen::Vector2d D;
    double dist = icy::Interaction::SegmentPointDistance(ndA->xt, ndB->xt, ndP->xt, D);
    if(dist < distance_threshold)
    {
        auto result = narrow_list_contact.insert((long long) edge_idx << 32 | node_idx);
        bool inserted = result.second;
        if(inserted)
        {
            Interaction i;
            i.ndA = ndA;
            i.ndB = ndB;
            i.ndP = ndP;
            i.D = D;
            collision_interactions.push_back(i);
        }
    }
    //        unsigned edge_idx = (unsigned)(contact_entry >> 32);
    //        int node_idx = (unsigned)(contact_entry & 0xffffffff);
}

void icy::Mesh::DetectContactPairs(double distance_threshold)
{
    // BROAD PHASE
    broadlist_contact.clear();
//    if(fragments.size()>1) //TODO: allow self-collisions
    mesh_root_contact.SelfCollide(broadlist_contact);

    unsigned nBroadListContact = broadlist_contact.size();

    narrow_list_contact.clear();
    collision_interactions.clear();

    // NARROW PHASE
#pragma omp parallel for
    for(unsigned i=0;i<nBroadListContact/2;i++)
    {
        unsigned edge1_idx = broadlist_contact[i*2];
        unsigned edge2_idx = broadlist_contact[i*2+1];
        Node *nd1, *nd2, *nd3, *nd4;
        std::tie(nd1,nd2) = allBoundaryEdges[edge1_idx];
        std::tie(nd3,nd4) = allBoundaryEdges[edge2_idx];

        AddToNarrowListIfNeeded(edge1_idx, nd3->globId, distance_threshold);
        AddToNarrowListIfNeeded(edge1_idx, nd4->globId, distance_threshold);
        AddToNarrowListIfNeeded(edge2_idx, nd1->globId, distance_threshold);
        AddToNarrowListIfNeeded(edge2_idx, nd2->globId, distance_threshold);
    }
}

std::pair<bool, double> icy::Mesh::EnsureNoIntersectionViaCCD()
{
    broadlist_ccd.clear();

    mesh_root_ccd.SelfCollide(broadlist_ccd);

    // CCD
    ccd_results.clear();

    unsigned nBroadListCCD = broadlist_ccd.size();
//    qDebug() << "nBroadListCCD size " << nBroadListCCD;

    bool final_state_contains_edge_intersection = false;

#pragma omp parallel for
    for(unsigned i=0;i<nBroadListCCD/2;i++)
    {
        unsigned edge1_idx = broadlist_ccd[i*2];
        unsigned edge2_idx = broadlist_ccd[i*2+1];

        Node *nd1, *nd2, *nd3, *nd4;
        std::tie(nd1,nd2) = allBoundaryEdges[edge1_idx];
        std::tie(nd3,nd4) = allBoundaryEdges[edge2_idx];

        if(EdgeIntersection(edge1_idx, edge2_idx)==true) final_state_contains_edge_intersection = true;

        bool result;
        double t;
        std::tie(result, t) = CCD(edge1_idx, nd3->globId);
        if(result) ccd_results.push_back(t);
        std::tie(result, t) = CCD(edge1_idx, nd4->globId);
        if(result) ccd_results.push_back(t);
        std::tie(result, t) = CCD(edge2_idx, nd1->globId);
        if(result) ccd_results.push_back(t);
        std::tie(result, t) = CCD(edge2_idx, nd2->globId);
        if(result) ccd_results.push_back(t);
    }

    if(ccd_results.size() > 0)
    {
        // find the smallest t; Model will discard the step
        auto iter = std::min_element(ccd_results.begin(), ccd_results.end());
        // std::cout << "CCD algorithm detected intersection;" << ccd_results.size() << "; min " << *iter << std::endl;
        return std::make_pair(false, *iter);
    }
    else if(final_state_contains_edge_intersection)
    {
        return std::make_pair(false, 0.5);
    }
    else
    {
        // proceed without intersections
        return std::make_pair(true, 0);
    }
}

bool icy::Mesh::EdgeIntersection(unsigned edgeIdx1, unsigned edgeIdx2)
{
    Node *ndA, *ndB, *ndC, *ndD;
    std::tie(ndA,ndB) = allBoundaryEdges[edgeIdx1];
    std::tie(ndC,ndD) = allBoundaryEdges[edgeIdx2];

    // if edges are adjacent, consider them non-intersecting
    if(ndA==ndC || ndA == ndD || ndB==ndC || ndB == ndD) return false;

    gte::Segment2<double> seg1;
    seg1.p[0] = {ndA->xt.x(), ndA->xt.y()};
    seg1.p[1] = {ndB->xt.x(), ndB->xt.y()};

    gte::Segment2<double> seg2;
    seg2.p[0] = {ndC->xt.x(), ndC->xt.y()};
    seg2.p[1] = {ndD->xt.x(), ndD->xt.y()};

    bool result = mTIQuery(seg1, seg2).intersect;
    return result;
}

std::pair<bool, double> icy::Mesh::CCD(unsigned edge_idx, unsigned node_idx)
{
    Node *ndA, *ndB, *ndP;
    std::tie(ndA,ndB) = allBoundaryEdges[edge_idx];
    ndP = allNodes[node_idx];

    if(ndA == ndP || ndB == ndP) return std::make_pair(false,0);

    double t;
    double px,py,ptx,pty;
    double v1x,v1y,v2x,v2y,v1tx,v1ty,v2tx,v2ty;

    px=ndP->xn[0];
    ptx=ndP->xt[0];
    py=ndP->xn[1];
    pty=ndP->xt[1];

    v1x=ndA->xn[0];
    v1tx=ndA->xt[0];
    v1y=ndA->xn[1];
    v1ty=ndA->xt[1];

    v2x=ndB->xn[0];
    v2tx=ndB->xt[0];
    v2y=ndB->xn[1];
    v2ty=ndB->xt[1];

    double t1, t2;

    double expr0t = pty*v1x - ptx*v1y + v1y*v2tx - v1x*v2ty - pty*v2x + v1ty*v2x - 2*v1y*v2x + py*(v1tx - 2*v1x - v2tx + 2*v2x) +
            px*(-v1ty + 2*v1y + v2ty - 2*v2y) + ptx*v2y - v1tx*v2y + 2*v1x*v2y;

    double expr1t= expr0t*expr0t -
            4*(-(ptx*v1ty) + px*v1ty + ptx*v1y - px*v1y + v1ty*v2tx - v1y*v2tx + ptx*v2ty - px*v2ty - v1tx*v2ty + v1x*v2ty +
               py*(-v1tx + v1x + v2tx - v2x) - v1ty*v2x + v1y*v2x + pty*(v1tx - v1x - v2tx + v2x) - ptx*v2y + px*v2y + v1tx*v2y -
               v1x*v2y)*(py*(v1x - v2x) + v1y*v2x - v1x*v2y + px*(-v1y + v2y));

    if(expr1t>=0)
    {
        double sqrt1t = sqrt(expr1t);

        t1=(py*v1tx - px*v1ty + pty*v1x - 2*py*v1x - ptx*v1y + 2*px*v1y - py*v2tx + v1y*v2tx + px*v2ty - v1x*v2ty - pty*v2x + 2*py*v2x +
                   v1ty*v2x - 2*v1y*v2x + ptx*v2y - 2*px*v2y - v1tx*v2y + 2*v1x*v2y +
                   sqrt1t)/
                (2.*(ptx*v1ty - px*v1ty - ptx*v1y + px*v1y - v1ty*v2tx + v1y*v2tx - ptx*v2ty + px*v2ty + v1tx*v2ty - v1x*v2ty +
                     pty*(-v1tx + v1x + v2tx - v2x) + v1ty*v2x - v1y*v2x +
                     py*(v1tx - v1x - v2tx + v2x) + ptx*v2y - px*v2y - v1tx*v2y + v1x*v2y));

        t2 = (-(py*v1tx) + px*v1ty - pty*v1x + 2*py*v1x + ptx*v1y - 2*px*v1y +
              py*v2tx - v1y*v2tx - px*v2ty + v1x*v2ty + pty*v2x - 2*py*v2x -
              v1ty*v2x + 2*v1y*v2x - ptx*v2y + 2*px*v2y + v1tx*v2y - 2*v1x*v2y +
              sqrt1t)/
           (2.*(-(ptx*v1ty) + px*v1ty + ptx*v1y - px*v1y + v1ty*v2tx - v1y*v2tx + ptx*v2ty - px*v2ty - v1tx*v2ty + v1x*v2ty +
                py*(-v1tx + v1x + v2tx - v2x) - v1ty*v2x + v1y*v2x + pty*(v1tx - v1x - v2tx + v2x) -
                ptx*v2y + px*v2y + v1tx*v2y - v1x*v2y));

        bool result = false;
        double result_s;
        if(t1>0 && t1<1)
        {
            double denom_x = -v1x + t1*(-v1tx + v1x + v2tx - v2x) + v2x;
            if(denom_x != 0)
            {
                double s = (px*(-1 + t1) - ptx*t1 + t1*v2tx + v2x - t1*v2x)/denom_x;
                if(s>=0 && s<=1) { result = true; t=t1; result_s = s; }
            }
            else
            {
                double denom_y = -v1y + t1*(-v1ty + v1y + v2ty - v2y) + v2y;
                if(denom_y != 0)
                {
                    double s = (py*(-1 + t1) - pty*t1 + t1*v2ty + v2y - t1*v2y)/denom_y;
                    if(s>=0 && s<=1) { result = true; t=t1; result_s=s;}
                }
            }
        }

        if(!result && t2>0 && t2<1)
        {
            double denom_x = -v1x + t2*(-v1tx + v1x + v2tx - v2x) + v2x;
            if(denom_x != 0)
            {
                double s = (px*(-1 + t2) - ptx*t2 + t2*v2tx + v2x - t2*v2x)/denom_x;
                if(s>=0 && s<=1) {result = true; t=t2;result_s=s;}
            }
            else
            {
                double denom_y = -v1y + t2*(-v1ty + v1y + v2ty - v2y) + v2y;
                if(denom_y != 0)
                {
                    double s = (py*(-1 + t2) - pty*t2 + t2*v2ty + v2y - t2*v2y)/denom_y;
                    if(s>=0 && s<=1) {result = true; t=t2;result_s=s;}
                }
            }
        }
        if(result)
        {
            if(t<1e-5) spdlog::warn("A({},{})-({},{}); B({},{})-({},{}); P({},{})-({},{}); t={}, s={}",
                                 v1x, v1y, v1tx, v1ty, v2x, v2y, v2tx, v2ty, px, py, ptx, pty, t, result_s);
            return std::make_pair(true,t);
        }
    }
    return std::make_pair(false,0);
}



// FRACTURE MODEL

void icy::Mesh::ComputeFractureDirections(const icy::SimParams &prms, double timeStep, bool startingFracture)
{
    const double threashold = prms.FractureTractionThreshold;
    const double temporal_attenuation = prms.FractureTemporalAttenuation;

#pragma omp parallel for
    for(unsigned i=0;i<allNodes.size();i++)
    {
        icy::Node *nd = allNodes[i];
        if(nd->pinned) continue;
        nd->ComputeFanVariables(prms);
    }

    maxNode = nullptr;
    if(!prms.EnableFracture) return;

    if(startingFracture)
    {
        breakable_range.clear();
        for(unsigned i=0;i<allNodes.size();i++)
        {
            icy::Node *nd = allNodes[i];
            if(nd->pinned) continue;
            if(nd->time_loaded_above_threshold >= temporal_attenuation)
            {
                if(nd->max_normal_traction > threashold) breakable_range.push_back(nd);
                else nd->time_loaded_above_threshold = 0;
            }
            else nd->time_loaded_above_threshold+=timeStep;
        }
        std::sort(breakable_range.begin(), breakable_range.end(), [](Node *nd1, Node *nd2)
            {return nd1->max_normal_traction > nd2->max_normal_traction;});

        // set some reasonable limit on the breakable_range
        const unsigned breakable_range_limit = std::max(prms.FractureMaxSubsteps/10,10);
        if(breakable_range.size()>breakable_range_limit) breakable_range.resize(breakable_range_limit);
    }
    else
    {
        // insert the recently created crack tips into the breakable range
        for(Node *nct : new_crack_tips)
        {
            nct->ComputeFanVariables(prms);
            nct->time_loaded_above_threshold = temporal_attenuation;
            auto find_result = std::find(breakable_range.begin(), breakable_range.end(),nct);
            bool already_contains = find_result!=breakable_range.end();
            if(!already_contains) breakable_range.push_back(nct);
        }
        // new_crack_tips.clear();

        // remove the nodes that were affected by the crack on the previous step
        breakable_range.erase(std::remove_if(breakable_range.begin(), breakable_range.end(),
                                          [temporal_attenuation](Node *nd)
                {return nd->max_normal_traction==0 || (nd->time_loaded_above_threshold < temporal_attenuation && !nd->isCrackTip);}),
                breakable_range.end());

        for(Node *nd : breakable_range) nd->ComputeFanVariables(prms);


    }

    if(breakable_range.size() > 0)
    {
        // take out maximal node from breakable_range
        auto it_nd = std::max_element(breakable_range.begin(), breakable_range.end(),
                                      [](Node *nd1, Node *nd2) {
                if(nd2->isCrackTip && nd2->max_normal_traction>0 && !nd1->isCrackTip) return true;
                return nd1->max_normal_traction < nd2->max_normal_traction; });

        if((*it_nd)->max_normal_traction > 0)
        {
            maxNode = *it_nd;
            maxNode->time_loaded_above_threshold = 0;
            breakable_range.erase(it_nd);
        }
    }
}

void icy::Mesh::EstablishSplittingEdge(Edge &splitEdge, Node* nd, const double phi, const double theta,
                            const Edge e0, const Edge e1, Element *elem)
{
    icy::Node *nd0, *nd1;
    std::tie(nd0,nd1) = elem->CW_CCW_Node(nd);

    const Eigen::Vector2d &nd_vec = nd->xn;
    const Eigen::Vector2d &nd0_vec = nd0->xn;
    const Eigen::Vector2d &nd1_vec = nd1->xn;

    double factor0 = sin(phi)*(nd0_vec-nd_vec).norm();
    double factor1 = sin(theta)*(nd1_vec-nd_vec).norm();
    double whereToSplit = factor1/(factor0+factor1);  // ~1 means the split is near nd0, ~0 means it is near nd1

    constexpr double epsilon = 1e-7;
    if((whereToSplit < epsilon && elem->isCCWBoundary(nd)) || (whereToSplit > 1-epsilon && elem->isCWBoundary(nd)))
    {
        // this should not happen
        spdlog::critical("icy::Mesh::EstablishSplittingEdge: degenerate element; nd {}, param {}", nd->locId, whereToSplit);
        spdlog::critical("phi {}; theta {}; factor0 {}; factor1 {}", phi, theta, factor0, factor1);
        spdlog::critical("elem {}-{}-{}", elem->nds[0]->locId, elem->nds[1]->locId, elem->nds[2]->locId);
        nd->PrintoutFan();
        throw std::runtime_error("degenerate element 1");
    }

    if(whereToSplit < fracture_epsilon && !elem->isCCWBoundary(nd))
    {
        splitEdge = e1;
    }
    else if(whereToSplit > 1.0-fracture_epsilon && !elem->isCWBoundary(nd))
    {
        splitEdge = e0;
    }
    else
    {
        icy::Element *elem_adj = elem->getAdjacentElementOppositeToNode(nd);
        if(elem_adj==nullptr)
            SplitBoundaryElem(elem, nd, nd0, nd1, whereToSplit, splitEdge);
        else
            SplitNonBoundaryElem(elem, elem_adj, nd, nd0, nd1, whereToSplit, splitEdge);
    }
    splitEdge.toSplit = true;
    splitEdge.elems[0]->edges[splitEdge.edge_in_elem_idx[0]]=splitEdge;
    splitEdge.elems[1]->edges[splitEdge.edge_in_elem_idx[1]]=splitEdge;
}

void icy::Mesh::SplitNode(const SimParams &prms)
{
    spdlog::info("SplitNode {}",maxNode->globId);
    if(maxNode == nullptr) throw std::runtime_error("SplitNode: trying to split nullptr");

    new_crack_tips.clear();
    affected_elements_during_split.clear();

    icy::Node* nd = maxNode;
    nd->time_loaded_above_threshold = 0;
    for(Element *e : nd->adj_elems) affected_elements_during_split.insert(e);

    // subsequent calculations are based on the fracture direction where the traction is maximal
    icy::Node::SepStressResult &ssr = nd->result_with_max_traction;

    // ensure that the interior node has two split faces
    bool isBoundary = (ssr.faces[1] == nullptr);
    if(isBoundary != nd->isBoundary) std::runtime_error("SplitNode: isBoundary != nd->isBoundary");

    // assert that the splitted node belongs to the element
    if(!ssr.faces[0]->containsNode(nd)) throw std::runtime_error("SplitNode: mesh toplogy error 0");
    icy::Edge splitEdge_fw;
    EstablishSplittingEdge(splitEdge_fw, nd, ssr.phi[0], ssr.theta[0], ssr.e[0], ssr.e[1], ssr.faces[0]);

    Node *split0 = splitEdge_fw.getOtherNode(nd);
    Node *split1 = nullptr;

    icy::Edge splitEdge_bw;
    if(!isBoundary)
    {
        // determine splitEdge_bw (create if necessary)
        if(!ssr.faces[1]->containsNode(nd)) throw std::runtime_error("SplitNode: mesh toplogy error 1");
        EstablishSplittingEdge(splitEdge_bw, nd, ssr.phi[1], ssr.theta[1], ssr.e[2], ssr.e[3], ssr.faces[1]);
        split1=splitEdge_bw.getOtherNode(nd);
    }

    Fix_X_Topology(nd); // split as fan.front().e[0] --- splitEdge_fw --- fan.back().e[1]

    if(!isBoundary)
    {
        if(split1==nullptr) throw std::runtime_error("SplitNode: split1==nullptr");
        if(split1->isBoundary)
        {
            Fix_X_Topology(split1);
        }
        else
        {
            split1->isCrackTip = true;
            new_crack_tips.push_back(split1);
            split1->weakening_direction = (split1->xn - nd->xn).normalized();
            for(Element *e : split1->adj_elems) affected_elements_during_split.insert(e);
        }
    }

    if(split0->isBoundary)
    {
        Fix_X_Topology(split0);
    }
    else
    {
        split0->isCrackTip = true;
        new_crack_tips.push_back(split0);
        split0->weakening_direction = (split0->xn - nd->xn).normalized();
        for(Element *e : split0->adj_elems) affected_elements_during_split.insert(e);
    }

    UpdateEdges();

    nd->weakening_direction = Eigen::Vector2d::Zero();
    nd->isCrackTip = false;
}

void icy::Mesh::Fix_X_Topology(Node *nd)
{
    nd->fan.clear();
    for(unsigned k=0;k<nd->adj_elems.size();k++)
    {
        icy::Element *elem = nd->adj_elems[k];

        Node::Sector s;
        s.face = elem;
        Eigen::Vector2d tcv = elem->getCenter() - nd->x_initial;
        s.centerAngle = atan2(tcv.y(), tcv.x());

        short thisIdx, CWIdx, CCWIdx;
        elem->getIdxs(nd, thisIdx, CWIdx, CCWIdx);

        s.nd[0] = elem->nds[CWIdx];
        s.nd[1] = elem->nds[CCWIdx];

        nd->fan.push_back(s);
    }

    std::sort(nd->fan.begin(), nd->fan.end(),
              [](const Node::Sector &f0, const Node::Sector &f1)
    {return f0.centerAngle < f1.centerAngle; });


    auto cw_boundary = std::find_if(nd->fan.begin(), nd->fan.end(),
                                    [nd](const Node::Sector &f){return f.face->CWEdge(nd).isBoundary;});
    if(cw_boundary != nd->fan.end()) std::rotate(nd->fan.begin(), cw_boundary, nd->fan.end());

    Node *split = nullptr;
    bool replacing = false;
    for(unsigned i=1;i<nd->fan.size();i++)
    {
        Node::Sector &s = nd->fan[i];
        if(s.face->CWEdge(nd).toSplit || s.face->CWEdge(nd).isBoundary) replacing=!replacing;
        if(replacing) {
            if(split==nullptr) {
                split = nd->fragment->AddNode();
                allNodes.push_back(split);
                split->Initialize(nd);
            }
            s.face->ReplaceNode(nd, split);
        }
    }

    if(split==nullptr) spdlog::warn("Fix_X_Topology: nothing to split; locId {}",nd->locId);

    for(Node::Sector &s : nd->fan) affected_elements_during_split.insert(s.face);
}

void icy::Mesh::UpdateEdges()
{
    std::unordered_set<Node*> affected_nodes_during_split; // their neighbors are also affected
    std::unordered_set<Element *> expanded_set_elems1;
    std::unordered_set<Element *> expanded_set_elems2;

    for(Element *elem : affected_elements_during_split)
    {
        for(int k=0;k<3;k++)
        {
            affected_nodes_during_split.insert(elem->nds[k]);
            if(elem->incident_elems[k]!=nullptr) expanded_set_elems1.insert(elem->incident_elems[k]);
            for(Element *elem2 : elem->nds[k]->adj_elems)
            {
                expanded_set_elems1.insert(elem2);
                for(int m=0;m<3;m++) affected_nodes_during_split.insert(elem2->nds[m]);
            }
        }
        expanded_set_elems1.insert(elem);
    }

    for(Node *nd : affected_nodes_during_split)
    {
        for(Element *elem : nd->adj_elems) expanded_set_elems2.insert(elem);
        nd->adj_elems.clear();
        nd->area=0;
        nd->isBoundary=false;
    }

    for(Element *elem : expanded_set_elems1)
    {
        for(int k=0;k<3;k++)
        {
            if(elem->incident_elems[k]!=nullptr) expanded_set_elems2.insert(elem->incident_elems[k]);
            elem->incident_elems[k]=nullptr;
        }
        expanded_set_elems2.insert(elem);
    }

    std::unordered_map<uint64_t, Edge> edges_map;

    for(Element *elem : expanded_set_elems2)
    {
        for(int i=0;i<3;i++)
        {
            Node *nd = elem->nds[i];
            // only work with the nodes in the affected_nodes set
            if(affected_nodes_during_split.find(nd)!=affected_nodes_during_split.end())
            {
                nd->adj_elems.push_back(elem);
            }
            // process edges
            icy::Node *nd0 = elem->nds[i];
            icy::Node *nd1 = elem->nds[(i+1)%3];
            uint64_t key = Node::make_key(nd0,nd1);
            auto result = edges_map.try_emplace(key, nd0, nd1);
            Edge &e = result.first->second;
            e.AddElement(elem, (i+2)%3);
        }
    }

    std::unordered_map<uint64_t, Edge> correctly_inferred_edges;
    // only take edges in the affected_elements set
    for(Element *elem : expanded_set_elems1)
    {
        for(int i=0;i<3;i++)
        {
            icy::Node *nd0 = elem->nds[i];
            icy::Node *nd1 = elem->nds[(i+1)%3];
            uint64_t key = Node::make_key(nd0,nd1);
            icy::Edge &existing_edge = edges_map.at(key);
            correctly_inferred_edges.insert({key,existing_edge});
        }
    }
/*
    allBoundaryEdges.erase(std::remove_if(allBoundaryEdges.begin(), allBoundaryEdges.end(),
                   [affected_nodes_during_split](std::pair<Node*,Node*> e)
    {return (affected_nodes_during_split.find(e.first)!=affected_nodes_during_split.end() &&
                affected_nodes_during_split.find(e.second)!=affected_nodes_during_split.end());}),
            allBoundaryEdges.end());
*/
    for(auto kvp : correctly_inferred_edges)
    {
        Edge &existing_edge = kvp.second;
        icy::Element *elem_of_edge0 = existing_edge.elems[0];
        icy::Element *elem_of_edge1 = existing_edge.elems[1];
        short idx0 = existing_edge.edge_in_elem_idx[0];
        short idx1 = existing_edge.edge_in_elem_idx[1];

        if(elem_of_edge0 == nullptr && elem_of_edge1 == nullptr) throw std::runtime_error("icy::Mesh::UpdateEdges: disconnected edge");
        existing_edge.isBoundary = (existing_edge.elems[0] == nullptr || existing_edge.elems[1] == nullptr);

        if(elem_of_edge0 != nullptr) elem_of_edge0->edges[idx0] = existing_edge;
        if(elem_of_edge1 != nullptr) elem_of_edge1->edges[idx1] = existing_edge;

        if(!existing_edge.isBoundary)
        {
            elem_of_edge0->incident_elems[idx0] = elem_of_edge1;
            elem_of_edge1->incident_elems[idx1] = elem_of_edge0;
        }

//        if(existing_edge.isBoundary) allBoundaryEdges.emplace_back(existing_edge.nds[0],existing_edge.nds[1]);
    }

    for(Node *nd : affected_nodes_during_split) nd->PrepareFan();
}

void icy::Mesh::SplitBoundaryElem(Element *originalElem, Node *nd, Node *nd0, Node *nd1, double where, Edge &insertedEdge)
{
    if(!originalElem->containsNode(nd)) throw std::runtime_error("SplitBoundaryElem: originalElem does not contain nd");
    if(!originalElem->containsNode(nd0)) throw std::runtime_error("SplitBoundaryElem: originalElem does not contain nd0");
    if(!originalElem->containsNode(nd1)) throw std::runtime_error("SplitBoundaryElem: originalElem does not contain nd1");
    short ndIdx = originalElem->getNodeIdx(nd);
    short nd0Idx = originalElem->getNodeIdx(nd0);
    short nd1Idx = originalElem->getNodeIdx(nd1);

    MeshFragment *fragment = nd->fragment;

    Element *insertedFace = fragment->AddElement();
    insertedFace->PiMultiplier = originalElem->PiMultiplier;
    allElems.push_back(insertedFace);
    nd->adj_elems.push_back(insertedFace);

    Node *split = fragment->AddNode();
    allNodes.push_back(split);
    split->InitializeLERP(nd0, nd1, where);
    split->isBoundary=true;
    split->adj_elems.push_back(originalElem);
    split->adj_elems.push_back(insertedFace);

    originalElem->nds[nd1Idx] = split;
    insertedFace->nds[ndIdx] = nd;
    insertedFace->nds[nd1Idx] = nd1;
    insertedFace->nds[nd0Idx] = split;
    insertedFace->edges[nd0Idx] = originalElem->edges[nd0Idx];
    insertedFace->PiMultiplier = originalElem->PiMultiplier; // TODO: test if this works

    insertedEdge = Edge(nd, split);
    insertedEdge.AddElement(insertedFace, nd1Idx);
    insertedEdge.AddElement(originalElem, nd0Idx);
    insertedEdge.isBoundary = false;
    insertedEdge.toSplit = true;

    Edge exteriorEdge1 = Edge(split, nd1);
    exteriorEdge1.isBoundary = true;
    exteriorEdge1.AddElement(insertedFace, ndIdx);
    insertedFace->edges[ndIdx] = exteriorEdge1;

    Edge exteriorEdge2 = Edge(split, nd0);
    exteriorEdge2.isBoundary = true;
    exteriorEdge2.AddElement(originalElem, ndIdx);
    originalElem->edges[ndIdx] = exteriorEdge2;

    originalElem->PrecomputeInitialArea();
    insertedFace->PrecomputeInitialArea();

    affected_elements_during_split.insert({originalElem,insertedFace});
}

void icy::Mesh::SplitNonBoundaryElem(Element *originalElem, Element *adjElem, Node *nd,
                                 Node *nd0, Node *nd1, double where, Edge &insertedEdge)
{
    originalElem->AssertEdges();
    if(!originalElem->containsNode(nd)) throw std::runtime_error("originalElem does not contain nd");
    if(!originalElem->containsNode(nd0)) throw std::runtime_error("originalElem does not contain nd0");
    if(!originalElem->containsNode(nd1)) throw std::runtime_error("originalElem does not contain nd1");

    short ndIdx_orig = originalElem->getNodeIdx(nd);
    short nd0Idx_orig = originalElem->getNodeIdx(nd0);
    short nd1Idx_orig = originalElem->getNodeIdx(nd1);

    Node *oppositeNode = adjElem->getOppositeNode(nd0, nd1);
    if(!adjElem->containsNode(oppositeNode)) throw std::runtime_error("adjElem does not contain opposite node");
    if(!adjElem->containsNode(nd0)) throw std::runtime_error("adjElem does not contain nd0");
    if(!adjElem->containsNode(nd1)) throw std::runtime_error("adjElem does not contain nd1");
    short nd0Idx_adj = adjElem->getNodeIdx(nd0);
    short nd1Idx_adj = adjElem->getNodeIdx(nd1);
    short oppIdx_adj = adjElem->getNodeIdx(oppositeNode);

    MeshFragment *fragment = nd->fragment;

    Element *insertedFace = fragment->AddElement();
    allElems.push_back(insertedFace);
    Element *insertedFace_adj = fragment->AddElement();
    allElems.push_back(insertedFace_adj);
    insertedFace->PiMultiplier = insertedFace_adj->PiMultiplier = originalElem->PiMultiplier;
    Node *split = fragment->AddNode();
    allNodes.push_back(split);

    split->InitializeLERP(nd0, nd1, where);

    split->adj_elems.push_back(originalElem);
    split->adj_elems.push_back(insertedFace);
    split->adj_elems.push_back(adjElem);
    split->adj_elems.push_back(insertedFace_adj);

    originalElem->nds[nd1Idx_orig] = split;
    insertedFace->nds[ndIdx_orig] = nd;
    insertedFace->nds[nd1Idx_orig] = nd1;
    insertedFace->nds[nd0Idx_orig] = split;
    insertedFace->edges[nd0Idx_orig] = originalElem->edges[nd0Idx_orig];
    insertedFace->PiMultiplier = originalElem->PiMultiplier; // TODO: test if this works

    nd->adj_elems.push_back(insertedFace);
    nd1->adj_elems.push_back(insertedFace);
    split->adj_elems.push_back(insertedFace);

    adjElem->nds[nd1Idx_adj] = split;
    insertedFace_adj->nds[oppIdx_adj] = oppositeNode;
    insertedFace_adj->nds[nd1Idx_adj] = nd1;
    insertedFace_adj->nds[nd0Idx_adj] = split;
    insertedFace_adj->edges[nd0Idx_adj] = adjElem->edges[nd0Idx_adj];
    insertedFace_adj->PiMultiplier = adjElem->PiMultiplier; // TODO: test if this works

    oppositeNode->adj_elems.push_back(insertedFace_adj);
    nd1->adj_elems.push_back(insertedFace_adj);
    split->adj_elems.push_back(insertedFace_adj);

    insertedEdge = Edge(nd, split);
    insertedEdge.AddElement(insertedFace, nd1Idx_orig);
    insertedEdge.AddElement(originalElem, nd0Idx_orig);
    insertedEdge.isBoundary = false;
    insertedEdge.toSplit = true;

    Edge insertedEdge_adj = Edge(oppositeNode, split);
    insertedEdge_adj.AddElement(insertedFace_adj, nd1Idx_adj);
    insertedEdge_adj.AddElement(adjElem, nd0Idx_adj);
    insertedEdge_adj.isBoundary = false;
    insertedEdge_adj.toSplit = false;
    insertedFace_adj->edges[nd1Idx_adj] = insertedEdge_adj;
    adjElem->edges[nd0Idx_adj] = insertedEdge_adj;

    Edge exteriorEdge1 = Edge(split, nd1);
    exteriorEdge1.isBoundary = false;
    exteriorEdge1.AddElement(insertedFace, ndIdx_orig);
    exteriorEdge1.AddElement(insertedFace_adj, oppIdx_adj);
    insertedFace->edges[ndIdx_orig] = exteriorEdge1;
    insertedFace_adj->edges[oppIdx_adj] = exteriorEdge1;

    Edge exteriorEdge2 = Edge(split, nd0);
    exteriorEdge2.isBoundary = false;
    exteriorEdge2.AddElement(originalElem, ndIdx_orig);
    exteriorEdge2.AddElement(adjElem, oppIdx_adj);
    originalElem->edges[ndIdx_orig] = exteriorEdge2;
    adjElem->edges[oppIdx_adj] = exteriorEdge2;

    originalElem->PrecomputeInitialArea();
    insertedFace->PrecomputeInitialArea();
    adjElem->PrecomputeInitialArea();
    insertedFace_adj->PrecomputeInitialArea();

    affected_elements_during_split.insert({originalElem,insertedFace,adjElem,insertedFace_adj});
}

void icy::Mesh::InferLocalSupport(SimParams &prms)
{
    if(maxNode==nullptr) throw std::runtime_error("InferLocalSupport: maxNode==nullptr");

    local_elems.clear();
    std::copy(maxNode->adj_elems.begin(),maxNode->adj_elems.end(),std::back_inserter(local_elems));
    CreateSupportRange(prms.FractureSubstepLevels, local_elems);

    std::unordered_set<Node*> local_support_set;
    for(Element *elem : local_elems) for(int k=0;k<3;k++) local_support_set.insert(elem->nds[k]);
    local_support.clear();
    std::copy(local_support_set.begin(), local_support_set.end(),std::back_inserter(local_support));

    // reset the loading timer in the vicinity of the crack
    local_elems2.clear();
    std::copy(maxNode->adj_elems.begin(),maxNode->adj_elems.end(),std::back_inserter(local_elems2));
    CreateSupportRange(prms.FractureTimerLevels, local_elems2);
    for(icy::Node *nd : allNodes) nd->reset_timing = nd->isSupportNode = false;
    for(Element *e : local_elems2)
        for(int k=0;k<3;k++)
        {
            e->nds[k]->time_loaded_above_threshold = 0;
            e->nds[k]->reset_timing=true;
        }

    // for visualization - mark support range (stored in breakable_range)
    for(icy::Node *nd : local_support) nd->isSupportNode = true; // for visualization
}

void icy::Mesh::CreateSupportRange(int neighborLevel, std::vector<Element*> &initial_set)
{
    for(unsigned i=0;i<allElems.size();i++) allElems[i]->traversal=0;

    std::queue<Element*> q_wave;
    for(Element *e : initial_set)
    {
        e->traversal=1;
        q_wave.push(e);
    }
    initial_set.clear();

    while(q_wave.size() > 0)
    {
        icy::Element *elem = q_wave.front();
        q_wave.pop();
        initial_set.push_back(elem);

        unsigned short level = elem->traversal;
        if(level < neighborLevel)
        {
            for(int i=0;i<3;i++)
            {
                icy::Element *adj_e = elem->incident_elems[i];
                if(adj_e!= nullptr && adj_e->traversal==0)
                {
                    adj_e->traversal=level+1;
                    q_wave.push(adj_e);
                }
            }
        }
    }
}
