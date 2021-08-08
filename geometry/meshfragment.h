#ifndef MESHFRAGMENT_H
#define MESHFRAGMENT_H

#include <vector>
#include <string>
#include "element.h"
#include "bvh/ConcurrentPool.h"
#include "bvh/bvhn.h"

namespace icy { class MeshFragment; }

class icy::MeshFragment
{    
public:
    ~MeshFragment();

    bool deformable;
    std::vector<icy::Node*> nodes;
    std::vector<icy::Element*> elems;

    std::vector<std::pair<icy::Node*,icy::Node*>> boundary_edges;
    std::vector<std::pair<unsigned,unsigned>> inner_boundary_edges; // for the experiments with fluid material

    void SaveFragment(std::string fileName);

    void GenerateBrick(double ElementSize);
    void GenerateSpecialBrick(double ElementSize);
    void GenerateIndenter(double ElementSize);
    void GenerateBall(double x, double y, double r1, double r2, double ElementSize);
    void GenerateContainer(double ElementSize, double offset);

    void RemeshSpecialBrick(double ElementSize);
    void RemeshWithBackgroundMesh(double ElementSize);
    void PostMeshingEvaluations(bool inferConnectivity=true);  // element/node area and connectivity information
    void Swap();

private:
    unsigned nFirstGroupElems, nFirstGroupNodes, nInnerBoundaryNodes;
    std::vector<icy::Node*> nodes_tmp;  // for remeshing
    std::vector<icy::Element*> elems_tmp;// for remeshing
    std::vector<std::pair<icy::Node*,icy::Node*>> boundary_edges_tmp;
    void GetFromGmsh();     // populate "nodes" and "elems" vectors
    void KeepGmshResult();  // save into gmshResult
    static ConcurrentPool<Node> NodeFactory;
    static ConcurrentPool<Element> ElementFactory;

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
    void GenerateLeafs(unsigned edge_idx);
private:
    static ConcurrentPool<BVHN> BVHNLeafFactory;

};

#endif // MESHFRAGMENT_H
