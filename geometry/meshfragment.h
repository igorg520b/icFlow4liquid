#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing tbb headers
#ifndef MESHFRAGMENT_H
#define MESHFRAGMENT_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
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

    bool isDeformable;
    std::vector<icy::Node*> nodes;
    std::vector<icy::Element*> elems;
    std::vector<std::pair<icy::Node*,icy::Node*>> special_boundary; // used when some nodes are pinned (movable via GUI)

    void SaveFragment(std::string fileName);

    void GenerateBrick(double ElementSize, double width, double height);
    void GenerateSelfCollisionBrick(double ElementSize, double width, double height);
    void GenerateIndenter(double ElementSize, double cx, double cy, double radius, double aspect);
    void GenerateContainer(double ElementSize, double offset);
    void PostMeshingEvaluations();  // element/node area and connectivity information

private:
    unsigned nFirstGroupElems, nFirstGroupNodes, nInnerBoundaryNodes;
    void GetFromGmsh();     // populate "nodes" and "elems"

    static ConcurrentPool<Node> NodeFactory;
    static ConcurrentPool<Element> ElementFactory;
    icy::Node* AddNode();       // add to nodes vector from NodeFactory; used in the fracture algorithm
    icy::Element* AddElement(); // add to elems vector from ElementFactory

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


// COLLISION DETECTION AND CONTACT MODEL
public:
    BVHN root_ccd;
    std::vector<BVHN*> leaves_for_ccd;
    void CreateLeaves();

    std::vector<std::pair<Node*, Node*>> boundaryEdges; // for collision with rigid shapes

private:
    static ConcurrentPool<BVHN> BVHNLeafFactory;


// FRACTURE MODEL
private:
    void ConnectIncidentElements();     // infer adjacency information and create the "fan", i.e. sorted vector of sectors per node

    friend icy::Mesh;
};

#endif // MESHFRAGMENT_H
#endif
