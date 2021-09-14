#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing tbb headers
#ifndef FL333_H
#define FL333_H

#include <gmsh.h>

#include <vector>
#include <unordered_set>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>

#include "meshfragment.h"
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

#include <Eigen/Core>

namespace icy { class Mesh; class Model; }

class icy::Mesh
{
public:
    unsigned typeOfSetup; // 0:intentation; 1:shear; 2:stretching

    std::vector<MeshFragment> fragments;   // including the indenter
    std::vector<icy::Node*> allNodes;
    std::vector<icy::Element*> allElems;
    std::unordered_map<uint64_t,BoundaryEdge> globalBoundaryEdges;

    std::vector<std::pair<Node*,Node*>> movableBoundary;    // controlled via GUI
    std::vector<icy::Node*> movableNodes;

    Mesh();
    void Reset(double MeshSizeMax, double offset, unsigned typeOfSetup_);
    void RegenerateVisualizedGeometry();    // from the collection of individual meshes, build allNodes, allElems, etc.
    void SetIndenterPosition(double position);
    double area_initial, area_current;


    // COLLISION DETECTION
public:
    tbb::concurrent_unordered_set<Interaction,Interaction::MyHash> contacts_narrow_set;
    std::vector<Interaction> contacts_final_list;
    void DetectContactPairs(const double distance_threshold);
    std::pair<bool, double> EnsureNoIntersectionViaCCD();

private:
    BVHN mesh_root_ccd, mesh_root_contact;  // only used when the mesh contains more than one fragment

    std::vector<BVHN*> global_leaves_ccd, global_leaves_contact, fragmentRoots_ccd, fragmentRoots_contact;
    std::vector<uint64_t> broadlist_ccd, broadlist_contact; // keys of potentially colliding edges
    tbb::concurrent_vector<double> ccd_results; // if not empty, time step is multiplied by the minimal value on the list

    void AddToNarrowListIfNeeded(Node* ndA, Node* ndB, bool edgeIsActive, Node *ndP, const double distance_threshold);
    static std::pair<bool, double> CCD(const Node* ndA, const Node* ndB, const bool isActive, const Node* ndP);  // if intersects, return [true, time]
    static bool EdgeIntersection(const Node* e1n1, const Node* e1n2,const Node* e2n1, const Node* e2n2); // true if edges intersect
    void CreateLeaves();
    void UpdateTree(float distance_threshold);
    unsigned tree_update_counter = 0;



    // FRACTURE
private:
    std::vector<Node*> breakable_range;     // populated in ComputeFractureDirections() when startingFracture==true
    std::vector<Node*> new_crack_tips;      // populated in SplitNode(), then used when startingFracture==false
    icy::Node *maxNode;
    constexpr static double fracture_epsilon = 0.1;   // if an edge splits too close to its vertex, then just go through the vertex
    void ComputeFractureDirections(const SimParams &prms, double timeStep, bool startingFracture);
    void SplitNode(const SimParams &prms);
    void EstablishSplittingEdge(Node* nd, const double phi, const double theta, Element *elem, Node* &adjacentNode);
    void SplitBoundaryElem(Element *originalElem, Node *nd, Node *nd0, Node *nd1, double where, Node*& insertedNode);
    void SplitNonBoundaryElem(Element *originalElem, Element *adjElem, Node *nd,
                                     Node *nd0, Node *nd1, double where, Node*& insertedNode);
    Node* Fix_X_Topology(Node *nd);

    void InferLocalSupport(SimParams &prms);
    void CreateSupportRange(int neighborLevel, std::vector<Element*> &initial_set);

    static double get_angle(const Eigen::Vector2d u, const Eigen::Vector2d v)
    { return (180.0/M_PI)*abs(acos(std::clamp((double)u.normalized().dot(v.normalized()),-1.0,1.0))); };

    std::vector<Element*> local_elems, local_elems2; // elems corresponding to breakable_range;
    std::vector<Node*> local_support;
    std::vector<BoundaryEdge> boundaries_created;
    void RemoveAdjBoundaries(Node *nd);
    void InsertAdjBoundaries(Node *nd);
public:

    // VTK
public:
    void ChangeVisualizationOption(int option);  // called from the main thread
    vtkNew<vtkLookupTable> hueLut;
    vtkNew<vtkActor> actor_collisions;
    vtkNew<vtkActor> actor_mesh_deformable;
    vtkNew<vtkActor> actor_boundary_all;
    vtkNew<vtkActor> actor_boundary_intended_indenter;

    enum ShowDeformationOption { initial, current };
    ShowDeformationOption showDeformation = ShowDeformationOption::current;
    bool updateMinMax;

private:
    vtkNew<vtkLookupTable> hueLutBlackRed;

    void UpdateValues();
    void UnsafeUpdateGeometry();

    int VisualizingVariable = 0;

    vtkNew<vtkPoints> points_deformable;
    vtkNew<vtkPoints> points_indenter_intended;   // prescribed indenter location
    vtkNew<vtkDoubleArray> visualized_values;
    vtkNew<vtkIntArray> visualized_values_edges;

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

    static constexpr double lutArrayTemperatureAdj[51][3] =
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
