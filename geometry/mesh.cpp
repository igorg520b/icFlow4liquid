#include "mesh.h"
#include "model.h"
#include <numeric>
#include <algorithm>
#include <iterator>
#include <unordered_set>
#include <queue>
#include "spdlog/spdlog.h"
#include <boost/container/flat_set.hpp>
#include <array>
#include <set>

void icy::Mesh::Reset(double MeshSizeMax, double offset, unsigned typeOfSetup_)
{
    maxNode = nullptr;
    typeOfSetup = typeOfSetup_;
    fragments.clear();
    movableBoundary.clear();
    movableNodes.clear();
    contacts_narrow_set.clear();
    contacts_final_list.clear();
    allNodes.clear();
    allElems.clear();
    allCZs.clear();
    allNodes.reserve(reserve_param);
    allElems.reserve(reserve_param);
    allCZs.reserve(reserve_param);

    switch(typeOfSetup)
    {
    // indentation
    case 0:
        for(int k=0;k<3;k++)
            fragments.push_back(std::make_unique<MeshFragment>(this));

        fragments[0]->GenerateIndenter(MeshSizeMax/2, 0, 1+0.15*1.1, 0.15, 2);
        fragments[1]->GenerateBrick(MeshSizeMax,2,1);
        fragments[2]->GenerateContainer(MeshSizeMax,offset);
        for(const auto &b : fragments[0]->boundaryEdges) movableBoundary.emplace_back(b->nds[0],b->nds[1]);
        break;

        // shear
    case 1:
    {
        const double radius = 0.1;
        const double height = 0.8;
        const double width = 2.0;
        for(int k=0;k<5;k++)
            fragments.push_back(std::make_unique<MeshFragment>(this));

        fragments[0]->GenerateBrick(MeshSizeMax, width, height);
        for(Node *nd : fragments[0]->nodes) if(nd->group.test(3)) nd->pinned=true;
        fragments[1]->GenerateIndenter(MeshSizeMax/2, width*0.1/2, height+radius+offset*1.5, radius, 1);
        fragments[2]->GenerateIndenter(MeshSizeMax/2, -width*0.9/2, height+radius+offset*1.5, radius, 1);
        fragments[3]->GenerateIndenter(MeshSizeMax/2, -width*0.1/2, -(radius+offset*1.5), radius, 1);
        fragments[4]->GenerateIndenter(MeshSizeMax/2, width*0.9/2, -(radius+offset*1.5), radius, 1);
        for(const auto &b : fragments[1]->boundaryEdges) movableBoundary.emplace_back(b->nds[0],b->nds[1]);
    }
        break;

        // stretch
    case 2:
        fragments.push_back(std::make_unique<MeshFragment>(this));
        fragments[0]->GenerateBrick(MeshSizeMax,2,1);
        for(Node *nd : fragments[0]->nodes) if(nd->group.test(1) || nd->group.test(2)) nd->pinned=true;
        movableBoundary.reserve(fragments[0]->special_boundary.size());
        for(const auto &b : fragments[0]->special_boundary) movableBoundary.push_back(b);
        break;

        // self-collision test
    case 3:
        fragments.push_back(std::make_unique<MeshFragment>(this));
        fragments[0]->GenerateSelfCollisionBrick(MeshSizeMax,2,1);
        movableBoundary.reserve(fragments[0]->special_boundary.size());
        for(const auto &b : fragments[0]->special_boundary) movableBoundary.push_back(b);
        break;

        //CZs
    case 4:
        fragments.push_back(std::make_unique<MeshFragment>(this));
        fragments[0]->GenerateCZBrick(MeshSizeMax,2,1);
        movableBoundary.reserve(fragments[0]->special_boundary.size());
        for(const auto &b : fragments[0]->special_boundary) movableBoundary.push_back(b);
        break;

        //CZs - fracture test
    case 5:
        fragments.push_back(std::make_unique<MeshFragment>(this));
        fragments[0]->GenerateCZBrickFractureTest(MeshSizeMax,2,1);
        movableBoundary.reserve(fragments[0]->special_boundary.size());
        for(const auto &b : fragments[0]->special_boundary) movableBoundary.push_back(b);
        break;
    }

    // make a list of nodes that only belong to the movable boudnary
    std::unordered_set<Node*> s;
    for(auto const &p : movableBoundary) { s.insert(p.first); s.insert(p.second); }
    movableNodes.resize(s.size());
    std::copy(s.begin(),s.end(),movableNodes.begin());
    for(unsigned i=0;i<movableNodes.size();i++) movableNodes[i]->indId=i;


    tree_update_counter=0;
    updateMinMax = true;

    area_initial = area_current = std::accumulate(allElems.begin(), allElems.end(),0.0,
                                                  [](double a, Element* m){return a+m->area_initial;});
    CreateLeaves();
    RegenerateVisualizedGeometry();
    ChangeVisualizationOption(VisualizingVariable);
}

icy::Mesh::Mesh()
{
    // initialize LUT
    const int nLut = 51;
    hueLut->SetNumberOfTableValues(nLut);
    for ( int i=0; i<nLut; i++)
            hueLut->SetTableValue(i, lutArrayTemperatureAdj[i][0],
                    lutArrayTemperatureAdj[i][1],
                    lutArrayTemperatureAdj[i][2], 1.0);

    hueLutBoundary->SetNumberOfTableValues(4);
    hueLutBoundary->SetTableValue(0,0,     0,     0, 1.0);
    hueLutBoundary->SetTableValue(1,0,     0,   0.2, 1.0);
    hueLutBoundary->SetTableValue(2,0,   0.5,   0.5, 1.0);
    hueLutBoundary->SetTableValue(3,0.6,   0,     0, 1.0);
    hueLutBoundary->SetTableRange(0,3);

    hueLut_czs->SetNumberOfTableValues(2);
    hueLut_czs->SetTableValue(1,0.95, 0.97, 0.91, 1.0);
    hueLut_czs->SetTableValue(0,0.45, 0.47, 0.41, 1.0);
    hueLut_czs->SetTableRange(0,1);

    // initialize various VTK objects
    visualized_values->SetName("visualized_values");
    visualized_values_edges->SetName("visualized_values_edges");
    visualized_values_czs->SetName("visualized_values_czs");

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
    actor_mesh_deformable->GetProperty()->SetLineWidth(0.6);

    // boundary
    ugrid_boundary_all->SetPoints(points_deformable);
    dataSetMapper_boundary_all->SetInputData(ugrid_boundary_all);

    dataSetMapper_boundary_all->UseLookupTableScalarRangeOn();
    dataSetMapper_boundary_all->SetLookupTable(hueLutBoundary);
    ugrid_boundary_all->GetCellData()->AddArray(visualized_values_edges);
    ugrid_boundary_all->GetCellData()->SetActiveScalars("visualized_values_edges");
    dataSetMapper_boundary_all->SetScalarModeToUseCellData();

    actor_boundary_all->SetMapper(dataSetMapper_boundary_all);
    actor_boundary_all->GetProperty()->EdgeVisibilityOn();
    actor_boundary_all->GetProperty()->VertexVisibilityOn();
    actor_boundary_all->GetProperty()->SetColor(0,0,0.4);
    actor_boundary_all->GetProperty()->SetEdgeColor(0,0,0.4);
    actor_boundary_all->GetProperty()->SetVertexColor(0.4,0,0);
    actor_boundary_all->GetProperty()->SetPointSize(3);
    actor_boundary_all->GetProperty()->SetLineWidth(1.4);
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

    //czs
    ugrid_czs->SetPoints(points_deformable);
    mapper_czs->SetInputData(ugrid_czs);
    mapper_czs->SetLookupTable(hueLut_czs);

    ugrid_czs->GetCellData()->AddArray(visualized_values_czs);
    ugrid_czs->GetCellData()->SetActiveScalars("visualized_values_czs");
    mapper_czs->SetScalarModeToUseCellData();

    actor_czs->SetMapper(mapper_czs);
    actor_czs->GetProperty()->VertexVisibilityOff();
    actor_czs->GetProperty()->EdgeVisibilityOn();
    actor_czs->GetProperty()->SetColor(0.95, 0.97, 0.91);
    actor_czs->GetProperty()->SetEdgeColor(110.0/255.0, 120.0/255.0, 12.0/255.0);
    actor_czs->GetProperty()->LightingOff();
    actor_czs->GetProperty()->ShadingOff();
    actor_czs->GetProperty()->SetInterpolationToFlat();
    actor_czs->PickableOff();
    actor_czs->GetProperty()->SetLineWidth(3);

    constexpr int initial_size =1000;
    contacts_final_list.reserve(initial_size*10);
    broadlist_ccd.reserve(initial_size*100);
    local_czs.reserve(initial_size);
    local_elems.reserve(initial_size);
    local_support.reserve(initial_size);

    mesh_root_ccd.test_self_collision = true;
    mesh_root_ccd.boundaryEdge = nullptr;
    mesh_root_ccd.isLeaf = false;
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
    case 4:
    case 5:
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
    points_deformable->SetNumberOfPoints(allNodes.size());

    // deformable material - elements
    cellArray_deformable->Reset();
    for(icy::Element *tr : allElems)
    {
        vtkIdType pts[3] = {tr->nds[0]->globId, tr->nds[1]->globId, tr->nds[2]->globId};
        cellArray_deformable->InsertNextCell(3, pts);
    }
    ugrid_deformable->SetCells(VTK_TRIANGLE, cellArray_deformable);

    // czs
    cellArray_czs->Reset();
    for(const icy::CohesiveZone *cz : allCZs)
    {
        vtkIdType pts[4] = {cz->nds[0]->globId, cz->nds[1]->globId, cz->nds[3]->globId, cz->nds[2]->globId};
        cellArray_czs->InsertNextCell(4, pts);
    }
    ugrid_czs->SetCells(VTK_QUAD,cellArray_czs);
    visualized_values_czs->SetNumberOfValues(allCZs.size());

    // VTK-displayed (thick) boundaries
    cellArray_boundary_all->Reset();
    visualized_values_edges->SetNumberOfValues(globalBoundaryEdges.size());

    unsigned count = 0;
    for(const auto b : globalBoundaryEdges)
    {
        vtkIdType pts[2] = {b->nds[0]->globId, b->nds[1]->globId};
        cellArray_boundary_all->InsertNextCell(2, pts);
        visualized_values_edges->SetValue(count++, b->status);
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
    points_collisions->SetNumberOfPoints(contacts_final_list.size()*2);
    cellArray_collisions->Reset();
    if(showDeformation == ShowDeformationOption::current)
    {
        for(unsigned i=0;i<contacts_final_list.size();i++)
        {
            vtkIdType pts[2] = {2*i, 2*i+1};
            cellArray_collisions->InsertNextCell(2, pts);
            Interaction &intr = contacts_final_list[i];
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
    case icy::Model::VisOpt::node_traversal:
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

    case icy::Model::VisOpt::node_traversal:
        visualized_values->SetNumberOfValues(allNodes.size());
        for(std::size_t i=0;i<allNodes.size();i++) visualized_values->SetValue(i, allNodes[i]->traversal);
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
    case icy::Model::VisOpt::quality_measure:
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
    case icy::Model::VisOpt::elem_group:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->group.to_ulong());
        break;

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

    case icy::Model::VisOpt::elem_isBoundary:
        for(std::size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->isBoundary());
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
        else if(VisualizingVariable == icy::Model::VisOpt::quality_measure)
        {
            hueLut->SetTableRange(0, 1);
        }
        else
        {
            hueLut->SetTableRange(minmax[0], minmax[1]);
        }
    }

    for(std::size_t i=0;i<allCZs.size();i++) visualized_values_czs->SetValue(i, allCZs[i]->isDamaged ? 1 : 0);
    visualized_values_czs->Modified();

}



// COLLISION DETECTION

void icy::Mesh::CreateLeaves()
{
    tree_update_counter = 0;
    global_leaves_ccd.clear();
    fragmentRoots_ccd.clear();
    globalBoundaryEdges.clear();

    for(auto &mf : fragments)
    {
        mf->CreateLeaves();
        globalBoundaryEdges.insert(globalBoundaryEdges.end(), mf->boundaryEdges.begin(),mf->boundaryEdges.end());
        global_leaves_ccd.insert(global_leaves_ccd.end(), mf->leaves_for_ccd.begin(), mf->leaves_for_ccd.end());
        fragmentRoots_ccd.push_back(&mf->root_ccd);
    }
}

void icy::Mesh::UpdateTree(float distance_threshold)
{
    // update leafs
    unsigned nLeafs = global_leaves_ccd.size();
#pragma omp parallel for
    for(unsigned i=0;i<nLeafs;i++)
    {
        BVHN *leaf_ccd = global_leaves_ccd[i];
        leaf_ccd->Expand_CCD(distance_threshold);
    }

    // update or build the rest of the tree
    if(tree_update_counter%10 != 0)
    {
        mesh_root_ccd.Update();
    }
    else
    {
        BVHN::BVHNFactory.releaseAll(); // does not release the leaves and roots

        if(fragmentRoots_ccd.size()>1)
        {
            for(auto &mf : fragments)
            {
                mf->root_ccd.Build(&mf->leaves_for_ccd,0);
            }
            mesh_root_ccd.Build(&fragmentRoots_ccd,0);
        }
        else
        {
            mesh_root_ccd.Build(&fragments.front()->leaves_for_ccd,0);
        }
    }
    tree_update_counter++;
}

void icy::Mesh::DetectContactPairs(const double distance_threshold)
{
    // BROAD PHASE of contact detection
    broadlist_ccd.clear();
    mesh_root_ccd.SelfCollide(broadlist_ccd);

    // NARROW PHASE
    contacts_narrow_set.clear();
    unsigned nBroadList = broadlist_ccd.size();
#pragma omp parallel for
    for(unsigned i=0;i<nBroadList;i++)
    {
        BVHN *bvhn1, *bvhn2;
        std::tie(bvhn1,bvhn2) = broadlist_ccd[i];

        if(bvhn1->boundaryEdge->isDeformable() && !bvhn2->boundaryEdge->isDeformable()) std::swap(bvhn1,bvhn2);

        if(!bvhn1->boundaryEdge->isDeformable() && bvhn2->boundaryEdge->isDeformable())
        {
            // indenter-deformable interaction
            Node *nd1 = bvhn1->boundaryEdge->nds[0];
            Node *nd2 = bvhn1->boundaryEdge->nds[1];
            Node *nd3 = bvhn2->boundaryEdge->nds[0];
            Node *nd4 = bvhn2->boundaryEdge->nds[1];
            AddToNarrowSet_NodeVsEdge(nd1, nd2, nd3, distance_threshold);
            AddToNarrowSet_NodeVsEdge(nd1, nd2, nd4, distance_threshold);
            AddToNarrowSet_NodeVsEdge(nd3, nd4, nd1, distance_threshold);
            AddToNarrowSet_NodeVsEdge(nd3, nd4, nd2, distance_threshold);
        }
        else if(bvhn1->boundaryEdge->isDeformable() && bvhn2->boundaryEdge->isDeformable())
        {
            // deformable-deformable interaction
        }
    }

    contacts_final_list.resize(contacts_narrow_set.size());
    std::copy(contacts_narrow_set.begin(),contacts_narrow_set.end(),contacts_final_list.begin());
}


bool icy::Mesh::PointInTriangle(Eigen::Vector2d pt, Eigen::Vector2d v1, Eigen::Vector2d v2, Eigen::Vector2d v3)
{
    auto sign = [](Eigen::Vector2d p1, Eigen::Vector2d p2, Eigen::Vector2d p3)
    {return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]);};

    double d1 = sign(pt, v1, v2);
    double d2 = sign(pt, v2, v3);
    double d3 = sign(pt, v3, v1);
    bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);
    return !(has_neg && has_pos);
}

void icy::Mesh::AddToNarrowSet_NodeVsEdge(Node* ndA, Node* ndB, Node *ndP, const double distance_threshold)
{
    if(ndA==ndP || ndB==ndP) return;    // a node can't collide with its incident edge
    Eigen::Vector2d D;
    double dist = icy::Interaction::SegmentPointDistance(ndA->xt, ndB->xt, ndP->xt, D);
    if(dist < distance_threshold)
    {
        Interaction i(ndA, ndB, ndP, D);
        contacts_narrow_set.insert(i);
    }
}

bool icy::Mesh::EnsureNoIntersectionViaCCD()
{
    broadlist_ccd.clear();
    mesh_root_ccd.SelfCollide(broadlist_ccd);

    unsigned nBroadList = broadlist_ccd.size();

    volatile bool intersection_detected = false;

#pragma omp parallel for shared(intersection_detected)
    for(unsigned i=0;i<nBroadList;i++)
    {
        if(intersection_detected) continue;
        BVHN *bvhn1, *bvhn2;
        std::tie(bvhn1,bvhn2) = broadlist_ccd[i];

        if(bvhn1->boundaryEdge->isDeformable() && !bvhn2->boundaryEdge->isDeformable()) std::swap(bvhn1,bvhn2);

        if(!bvhn1->boundaryEdge->isDeformable() && bvhn2->boundaryEdge->isDeformable())
        {
            Node *nd1 = bvhn1->boundaryEdge->nds[0];
            Node *nd2 = bvhn1->boundaryEdge->nds[1];
            Node *nd3 = bvhn2->boundaryEdge->nds[0];
            Node *nd4 = bvhn2->boundaryEdge->nds[1];
            if(EdgeIntersection(nd1,nd2,nd3,nd4)) intersection_detected = true;
            if(CCD(nd1, nd2, nd3)) intersection_detected = true;
            if(CCD(nd1, nd2, nd4)) intersection_detected = true;
            if(CCD(nd3, nd4, nd1)) intersection_detected = true;
            if(CCD(nd3, nd4, nd2)) intersection_detected = true;
        }
        else if(bvhn1->boundaryEdge->isDeformable() && bvhn2->boundaryEdge->isDeformable())
        {
            // check for possible deformable-deformable intersection
        }
    }

    return !intersection_detected;
}

bool icy::Mesh::EdgeIntersection(const Node* ndA, const Node* ndB, const Node* ndC, const Node* ndD)
{
    // if edges are adjacent, consider them non-intersecting
    if(ndA==ndC || ndA == ndD || ndB==ndC || ndB == ndD) return false;

    gte::Segment2<double> seg1;
    seg1.p[0] = {ndA->xt[0], ndA->xt[1]};
    seg1.p[1] = {ndB->xt[0], ndB->xt[1]};

    gte::Segment2<double> seg2;
    seg2.p[0] = {ndC->xt[0], ndC->xt[1]};
    seg2.p[1] = {ndD->xt[0], ndD->xt[1]};

    gte::TIQuery<double, gte::Segment2<double>, gte::Segment2<double>> mTIQuery;
    return mTIQuery(seg1, seg2).intersect;
}

bool icy::Mesh::CCD(const Node* ndA, const Node* ndB, const Node* ndP)
{
    if(ndA == ndP || ndB == ndP) return false;

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
            return true; // t contains the collision time
        }
    }
    return false;
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

void icy::Mesh::PropagateCrack(const SimParams &prms)
{
    icy::Node* nd = maxNode;
    if(nd == nullptr) throw std::runtime_error("PropagateCrack: trying to split nullptr");

    MeshFragment *fr = nd->fragment;
    nd->time_loaded_above_threshold = 0;
    new_crack_tips.clear();

    // subsequent calculations are based on the fracture direction where the traction is maximal
    icy::Node::SepStressResult &ssr = nd->result_with_max_traction;

    // ensure that the interior node has two split faces
    bool isBoundary = (ssr.faces[1] == nullptr);
    if(isBoundary != nd->isBoundary) std::runtime_error("PropagateCrack: isBoundary != nd->isBoundary");
    if(!ssr.faces[0]->containsNode(nd)) throw std::runtime_error("PropagateCrack: mesh toplogy error 0");

    Node *edge_split0, *edge_split1 = nullptr;

    EstablishSplittingEdge(nd, ssr.phi[0], ssr.theta[0], ssr.faces[0], edge_split0);

    // add boundaries to ssr.faces[0] and its incident elem along the edge nd--edge_split0
    uint8_t edge_idx = ssr.faces[0]->getEdgeIdx(nd,edge_split0);
    Element* incidentElem1 = dynamic_cast<Element*>(ssr.faces[0]->incident_elems[edge_idx]);
    if(incidentElem1 == nullptr) throw std::runtime_error("PropagateCrack: insertedElem1 == nullptr");
    if(prms.EnableInsertCZs)
    {
        CohesiveZone *cz = AddCZ();
        cz->Initialize(ssr.faces[0],edge_idx,incidentElem1,incidentElem1->getEdgeIdx(nd,edge_split0));
    }
    else
    {
        fr->AddBoundary(ssr.faces[0],edge_idx,2);
        fr->AddBoundary(incidentElem1,incidentElem1->getEdgeIdx(nd,edge_split0),2);
    }

    if(!isBoundary)
    {
        if(!ssr.faces[1]->containsNode(nd)) throw std::runtime_error("PropagateCrack: mesh toplogy error 1");
        EstablishSplittingEdge(nd, ssr.phi[1], ssr.theta[1], ssr.faces[1], edge_split1);

        // add new boundaries
        uint8_t edge_idx2 = ssr.faces[1]->getEdgeIdx(nd,edge_split1);
        Element* insertedElem2 = dynamic_cast<Element*>(ssr.faces[1]->incident_elems[edge_idx2]);
        if(insertedElem2 == nullptr) throw std::runtime_error("PropagateCrack: insertedElem2 == nullptr");
        if(prms.EnableInsertCZs)
        {
            CohesiveZone *cz = AddCZ();
            cz->Initialize(ssr.faces[1],edge_idx2,insertedElem2,insertedElem2->getEdgeIdx(nd,edge_split1));
        }
        else
        {
            fr->AddBoundary(ssr.faces[1],edge_idx2,2);
            fr->AddBoundary(insertedElem2,insertedElem2->getEdgeIdx(nd,edge_split1),2);
        }
    }

    // SPLIT nd
    Node *nd_split = fr->AddNode();
    {
        nd->CreateUnrotatedFan();

        // find the sector whose CW boundary is nd--edge_split0
        auto cw_boundary_sector = std::find_if(nd->fan.begin(), nd->fan.end(),
                                               [edge_split0](const Node::Sector &f){return f.nd[0]==edge_split0;});
        if(cw_boundary_sector == nd->fan.end()) throw std::runtime_error("PropagateCrack: cw boundary not found 1");
        std::rotate(nd->fan.begin(), cw_boundary_sector, nd->fan.end());


        nd->adj_elems.clear();
        nd_split->Initialize(nd);
        bool other_side = false;

        for(Node::Sector &s : nd->fan)
        {
            if(!other_side &&
                    ((isBoundary && s.face->isCWBoundary(nd) && s.nd[0] != edge_split0) ||
                     (!isBoundary && s.nd[0] == edge_split1)))
                other_side = true;

            if(other_side)
            {
                nd_split->adj_elems.push_back(s.face);
                s.face->ReplaceNode(nd, nd_split);
            }
            else
            {
                nd->adj_elems.push_back(s.face);
            }
        }

        nd->weakening_direction.setZero();
        nd->isCrackTip = false;
        nd->PrepareFan();
        nd_split->PrepareFan();
    }


    if(edge_split0->isBoundary)
    {
        Node *edge_split0_split = Fix_X_Topology(edge_split0, nd_split);
    }
    else
    {
        edge_split0->isCrackTip = true;
        edge_split0->weakening_direction = (edge_split0->xn - nd->xn).normalized();
        new_crack_tips.push_back(edge_split0);
        edge_split0->PrepareFan();
    }

    if(!isBoundary)
    {
        if(edge_split1->isBoundary)
        {
            Node *edge_split1_split = Fix_X_Topology(edge_split1, nd);
        }
        else
        {
            edge_split1->isCrackTip = true;
            edge_split1->weakening_direction = (edge_split1->xn - nd->xn).normalized();
            new_crack_tips.push_back(edge_split1);
            edge_split1->PrepareFan();
        }
    }

}

void icy::Mesh::EstablishSplittingEdge(Node* nd, const double phi, const double theta, Element *elem, Node* &adjacentNode)
{
    icy::Node *nd0, *nd1;
    std::tie(nd0,nd1) = elem->CW_CCW_Node(nd);

    const Eigen::Vector2d &nd_vec = nd->xn;
    const Eigen::Vector2d &nd0_vec = nd0->xn;
    const Eigen::Vector2d &nd1_vec = nd1->xn;

    double factor0 = sin(phi)*(nd0_vec-nd_vec).norm();
    double factor1 = sin(theta)*(nd1_vec-nd_vec).norm();
    double whereToSplit = factor1/(factor0+factor1);  // ~1 means the split is near nd0, ~0 means it is near nd1

    if((theta == 0 && elem->isCCWBoundary(nd)) || (phi == 0 && elem->isCWBoundary(nd)))
        throw std::runtime_error("trying to split the boundary");

    if(theta == 0 && !elem->isCCWBoundary(nd))
        adjacentNode = nd1;
    else if(phi == 0 && !elem->isCWBoundary(nd))
        adjacentNode = nd0;
    else
    {
        Node* additionalNode;
        std::tie(adjacentNode,additionalNode) = elem->SplitElem(nd,nd0,nd1,whereToSplit);
        if(additionalNode!=nullptr) new_crack_tips.push_back(additionalNode);
    }
}




icy::Node* icy::Mesh::Fix_X_Topology(Node *nd, Node *alignment_node)
{
    // create a new "fan" with X-topology allowed
    nd->CreateUnrotatedFan();

    // in the fan, find the entry with clock-wise boundary
    auto cw_bdry = std::find_if(nd->fan.begin(), nd->fan.end(),
                                [nd,alignment_node](const Node::Sector &f){return f.face->isCWBoundary(nd) && f.nd[0]==alignment_node;});
    if(cw_bdry == nd->fan.end()) throw std::runtime_error("Fix_X_Topology: could not rotate");
    std::rotate(nd->fan.begin(), cw_bdry, nd->fan.end());

    Node *split = nd->fragment->AddNode();
    split->Initialize(nd);

    bool other_side = false;
    nd->adj_elems.clear();

    for(Node::Sector &s : nd->fan)
    {
        if(s.nd[0] != alignment_node && s.face->isCWBoundary(nd)) other_side = true;
        if(other_side)
        {
            split->adj_elems.push_back(s.face);
            s.face->ReplaceNode(nd, split);
        }
        else
        {
            nd->adj_elems.push_back(s.face);
        }
    }

    if(!other_side) { nd->PrintoutFan(); throw std::runtime_error("Fix_X_Topology: could not split"); }

    split->PrepareFan();
    nd->PrepareFan();
    return split;
}







void icy::Mesh::InferLocalSupport(SimParams &prms)
{
    if(maxNode==nullptr) throw std::runtime_error("InferLocalSupport: maxNode==nullptr");

    local_elems.clear();
    local_czs.clear();

    std::copy(maxNode->adj_elems.begin(),maxNode->adj_elems.end(),std::back_inserter(local_elems));
    CreateSupportRange(prms.FractureSubstepLevels);
}

void icy::Mesh::ResetFractureTimer(SimParams &prms)
{
    // reset the loading timer in the vicinity of the crack
    CreateSupportRange(prms.FractureTimerLevels);
    for(Node *nd : local_support) nd->time_loaded_above_threshold = 0;
}



void icy::Mesh::CreateSupportRange(const int neighborLevel)
{
    local_czs.clear();
    local_support.clear();
    local_elems.clear();

    for(auto &elem : allElems) elem->traversed=false;
    for(auto &cz : allCZs) cz->traversed = false;
    for(auto &nd : allNodes) nd->traversal = 0;


    std::queue<Node*> q_wave;
    q_wave.push(maxNode);
    maxNode->traversal = 1;

    while(q_wave.size() > 0)
    {
        Node *nd = q_wave.front();
        q_wave.pop();
        local_support.push_back(nd);
        unsigned short level = nd->traversal;

        for(Element *e : nd->adj_elems)
        {
            if(!e->traversed)
            {
                local_elems.push_back(e);
                e->traversed = true;
                if(level<neighborLevel)
                {
                    for(Node *nd2 : e->nds)
                    {
                        if(nd2->traversal == 0)
                        {
                            nd2->traversal = level+1;
                            q_wave.push(nd2);
                        }
                    }
                }

                for(BaseElement *be : e->incident_elems)
                    if(be->type == BaseElement::ElementType::CZ)
                    {
                        CohesiveZone *cz = static_cast<CohesiveZone*>(be);
                        cz->traversed = true;
                        local_czs.push_back(cz);
                        for(Node *nd2 : cz->nds)
                        {
                            if(nd2->traversal == 0)
                            {
                                nd2->traversal = level+1;
                                q_wave.push(nd2);
                            }
                        }
                    }
            }
        }
    }
}

icy::CohesiveZone* icy::Mesh::AddCZ()
{
    CohesiveZone *cz = MeshFragment::CZFactory.take();
    cz->parentMesh = this;
    allCZs.push_back(cz);
    return cz;
}
