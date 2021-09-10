#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing tbb headers
#ifndef MESHFRAGMENT_H
#define MESHFRAGMENT_H

#include <vector>
#include <unordered_map>
#include <string>
#include "element.h"
#include "bvh/ConcurrentPool.h"
#include "bvh/bvhn.h"
#include "spdlog/spdlog.h"

namespace icy { class MeshFragment; class Mesh; }

class icy::MeshFragment
{    
public:
    ~MeshFragment();

    bool deformable;
    std::vector<icy::Node*> nodes;
    std::vector<icy::Element*> elems;

//    std::vector<std::pair<icy::Node*,icy::Node*>> boundary_edges;
    std::unordered_map<uint64_t,std::pair<Node*,Node*>> boundaryEdgesMap;

    std::vector<std::pair<icy::Node*,icy::Node*>> special_boundary, fixed_boundary; // these are used when some nodes are pinned
    std::vector<std::pair<unsigned,unsigned>> inner_boundary_edges; // for the experiments with fluid material

    void SaveFragment(std::string fileName);

    void GenerateBrick(double ElementSize, double width, double height);
    void GenerateBrick_Simple(double ElementSize, double width, double height);
    void GenerateSelfCollisionBrick(double ElementSize, double width, double height);
    void GenerateSpecialBrick(double ElementSize);
    void GenerateIndenter(double ElementSize, double cx, double cy, double radius, double aspect);
    void GenerateBall(double x, double y, double r1, double r2, double ElementSize);
    void GenerateContainer(double ElementSize, double offset);

    void RemeshSpecialBrick(double ElementSize);
    void RemeshWithBackgroundMesh(double ElementSize);
    void PostMeshingEvaluations();  // element/node area and connectivity information
    void Swap();

private:
    unsigned nFirstGroupElems, nFirstGroupNodes, nInnerBoundaryNodes;
    std::vector<icy::Node*> nodes_tmp;  // for remeshing
    std::vector<icy::Element*> elems_tmp;// for remeshing
    std::vector<std::pair<icy::Node*,icy::Node*>> boundary_edges_tmp;
    void GetFromGmsh();     // populate "nodes" and "elems"
    void KeepGmshResult();  // save gmsh representation into gmshResult

    static ConcurrentPool<Node> NodeFactory;
    static ConcurrentPool<Element> ElementFactory;
    icy::Node* AddNode();   // automatically assign locId and fragment, then add to nodes vector; also used by Mesh class
    icy::Element* AddElement();

    struct GmshEntity
    {
        std::pair<int, int> dimTags;
        std::vector<std::pair<int,int>> outDimTags_boundary;
        std::vector<int> boundary;  // same as outDimTags_boundary, but only the second element of the pair
        std::vector<std::size_t> nodeTags_nodes;
        std::vector<double> coord_nodes;
        std::vector<double> parametricCoord_nodes;
        std::vector<int> elementTypes_elements;
        std::vector<std::vector<std::size_t>> elementTags_elements;
        std::vector<std::vector<std::size_t>> nodeTags_elements;
    };
    std::vector<GmshEntity> gmshResult; // store the result returned by the gmsh, to reuse later as a background mesh
    std::unordered_map<std::size_t, std::size_t> mtags; // gmsh nodeTag -> sequential position in nodes[]

// BROAD PHASE
public:
    BVHN root_ccd, root_contact;
    std::vector<BVHN*> leaves_for_ccd, leaves_for_contact;
    void CreateLeaves();
private:
    static ConcurrentPool<BVHN> BVHNLeafFactory;

// FRACTURE MODEL
private:
    void ConnectIncidentElements();     // infer adjacency information and create the "fan", i.e. sorted vector of sectors per node

    friend icy::Mesh;
};

#endif // MESHFRAGMENT_H
#endif
