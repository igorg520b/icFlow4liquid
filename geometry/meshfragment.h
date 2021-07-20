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
    // TODO: destructor to allow removal of the fragments (release BVHNs)

public:
    bool deformable;
    std::vector<icy::Node*> nodes;
    std::vector<icy::Element*> elems;

    std::vector<std::pair<icy::Node*,icy::Node*>> boundary_edges;
    std::vector<std::pair<unsigned,unsigned>> inner_boundary_edges; // for the experiments with fluid material

    unsigned freeNodeCount;
    void SaveFragment(std::string fileName);

    void GenerateBrick(double ElementSize);
    void GenerateSpecialBrick(double ElementSize);
    void GenerateIndenter(double ElementSize);
    void GenerateBall(double x, double y, double r1, double r2, double ElementSize);
    void GenerateContainer(double ElementSize, double offset);

private:
    std::vector<icy::Node*> nodes_tmp;  // for remeshing
    std::vector<icy::Element*> elems_tmp;// for remeshing
    void GetFromGmsh();
    static ConcurrentPool<Node> NodeFactory;
    static ConcurrentPool<Element> ElementFactory;

// BROAD PHASE

public:
    BVHN root_ccd, root_contact;
    std::vector<BVHN*> leafs_for_ccd, leafs_for_contact;
    void GenerateLeafs(unsigned edge_idx);
private:
    static ConcurrentPool<BVHN> BVHNLeafFactory;

};

#endif // MESHFRAGMENT_H
