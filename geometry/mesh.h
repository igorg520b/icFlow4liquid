#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing tbb headers
#ifndef FL333_H
#define FL333_H

#include <gmsh.h>

#include <vector>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>

#include "meshfragment.h"
#include "element.h"
#include "interaction.h"
#include "bvh/bvhn.h"

#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkNamedColors.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>

#include "Mathematics/IntrSegment2Segment2.h"

namespace icy { class Mesh; class Model; }

class icy::Mesh
{
public:
    MeshFragment indenter;
    std::vector<MeshFragment*> allMeshes;   // including the indenter
    std::vector<icy::Node*> allNodes;
    std::vector<icy::Element*> allElems;
    std::vector<std::pair<Node*,Node*>> allBoundaryEdges; // for visualization
    unsigned freeNodeCount;

    Mesh();
    void Reset(double CharacteristicLengthMax, double offset);


private:
    void UpdateValues();
    void UnsafeUpdateGeometry();
    void RegenerateVisualizedGeometry();    // from the collection of individual meshes, build allNodes, allElems, etc.

    // Collision detection
public:
    tbb::concurrent_vector<Interaction> collision_interactions;
    void DetectContactPairs(double distance_threshold);
    std::pair<bool, double> EnsureNoIntersectionViaCCD();

private:
    BVHN root_ccd, root_contact;
    std::vector<BVHN*> global_leafs_ccd, global_leafs_contact, fragmentRoots_ccd, fragmentRoots_contact;
    std::vector<unsigned> broadlist_ccd, broadlist_contact; // indices of potentially colliding edges
    tbb::concurrent_unordered_set<long long> narrow_list_contact;
    tbb::concurrent_vector<double> ccd_results; // if not empty, time step is multiplied by the minimal value on the list

    void AddToNarrowListIfNeeded(unsigned edge_idx, unsigned node_idx, double distance_threshold);
    std::pair<bool, double> CCD(unsigned edge_idx, unsigned node_idx);  // if intersects, return [true, time]
    bool EdgeIntersection(unsigned edgeIdx1, unsigned edgeIdx2); // true if edges intersect
    void UpdateTree(float distance_threshold);
    void BuildTree(float distance_threshold);
    unsigned tree_update_counter = 0;

    gte::TIQuery<double, gte::Segment2<double>, gte::Segment2<double>> mTIQuery;

    // VTK
public:
    void ChangeVisualizationOption(int option);  // called from the main thread
    vtkNew<vtkLookupTable> hueLut;
    vtkNew<vtkActor> actor_collisions;
    vtkNew<vtkActor> actor_mesh_deformable;
    vtkNew<vtkActor> actor_boundary_all;
    vtkNew<vtkActor> actor_boundary_intended_indenter;

private:

    int VisualizingVariable = 0;

    vtkNew<vtkPoints> points_deformable;
    vtkNew<vtkPoints> points_indenter_intended;   // prescribed indenter location
    vtkNew<vtkDoubleArray> visualized_values;

    // elements
    vtkNew<vtkUnstructuredGrid> ugrid_deformable;
    vtkNew<vtkCellArray> cellArray_deformable;
    vtkNew<vtkDataSetMapper> dataSetMapper_deformable;

    // boundary
    vtkNew<vtkUnstructuredGrid> ugrid_boundary_all;
    vtkNew<vtkCellArray> cellArray_boundary_all;
    vtkNew<vtkDataSetMapper> dataSetMapper_boundary_all;

    // boundary-intended
    vtkNew<vtkUnstructuredGrid> ugrid_indenter_intended;
    vtkNew<vtkCellArray> cellArray_indenter_intended;
    vtkNew<vtkDataSetMapper> dataSetMapper_indenter_intended;

    // collisions
    vtkNew<vtkPoints> points_collisions;
    vtkNew<vtkUnstructuredGrid> ugrid_collisions;
    vtkNew<vtkDataSetMapper> mapper_collisions;
    vtkNew<vtkCellArray> cellArray_collisions;

    static constexpr float lutArrayTemperatureAdj[51][3] =
    {{0.770938, 0.951263, 0.985716}, {0.788065, 0.959241, 0.986878},
     {0.805191, 0.96722, 0.98804}, {0.822318, 0.975199, 0.989202},
     {0.839445, 0.983178, 0.990364}, {0.856572, 0.991157, 0.991526},
     {0.872644, 0.995552, 0.98386}, {0.887397, 0.995466, 0.965157},
     {0.902149, 0.99538, 0.946454}, {0.916902, 0.995294, 0.927751},
     {0.931655, 0.995208, 0.909049}, {0.946408, 0.995123, 0.890346},
     {0.961161, 0.995037, 0.871643}, {0.975913, 0.994951, 0.85294},
     {0.990666, 0.994865, 0.834237}, {0.996257, 0.991758, 0.815237},
     {0.994518, 0.986234, 0.795999}, {0.992779, 0.98071, 0.77676},
     {0.99104, 0.975186, 0.757522}, {0.989301, 0.969662, 0.738283},
     {0.987562, 0.964138, 0.719045}, {0.985823, 0.958614, 0.699807},
     {0.984084, 0.953089, 0.680568}, {0.982345, 0.947565, 0.66133},
     {0.97888, 0.936201, 0.641773}, {0.974552, 0.921917, 0.622058},
     {0.970225, 0.907633, 0.602342}, {0.965897, 0.893348, 0.582626},
     {0.961569, 0.879064, 0.562911}, {0.957242, 0.86478, 0.543195},
     {0.952914, 0.850496, 0.52348}, {0.948586, 0.836212, 0.503764},
     {0.944259, 0.821927, 0.484048}, {0.939066, 0.801586, 0.464871},
     {0.933626, 0.779513, 0.445847}, {0.928186, 0.757441, 0.426823},
     {0.922746, 0.735368, 0.4078}, {0.917306, 0.713296, 0.388776},
     {0.911866, 0.691223, 0.369752}, {0.906426, 0.669151, 0.350728},
     {0.900986, 0.647078, 0.331704}, {0.895546, 0.625006, 0.312681},
     {0.889975, 0.597251, 0.298625}, {0.884388, 0.568785, 0.285191},
     {0.8788, 0.54032, 0.271756}, {0.873212, 0.511855, 0.258322},
     {0.867625, 0.483389, 0.244888}, {0.862037, 0.454924, 0.231453},
     {0.856449, 0.426459, 0.218019}, {0.850862, 0.397993, 0.204584},
     {0.845274, 0.369528, 0.19115}};

    friend class icy::Model;
};
#endif
#endif
