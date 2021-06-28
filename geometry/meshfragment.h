#ifndef MESHFRAGMENT_H
#define MESHFRAGMENT_H

#include "element.h"
#include "bvh/ConcurrentPool.h"
#include "bvh/bvhn.h"
#include <vector>

namespace icy { class MeshFragment; }

class icy::MeshFragment
{    
    // TODO: destructor to allow removal of the fragments (release BVHNs)

public:
    bool deformable;
    std::vector<icy::Node> nodes;
    std::vector<icy::Element> elems;
    std::vector<int> boundary_nodes;
    std::vector<std::pair<icy::Node*,icy::Node*>> boundary_edges;
    unsigned freeNodeCount;

    void GenerateBrick(double ElementSize);
    void GenerateIndenter(double ElementSize);
    void GenerateCup(double ElementSize);
    void GenerateBall(double x, double y, double r1, double r2, double ElementSize);
    void GenerateContainer(double ElementSize, double offset);
//    void GenerateSelfCollisionTest(double ElementSize);

private:
    void GetFromGmsh();

// BROAD PHASE

public:
    BVHN root_ccd, root_contact;
    std::vector<BVHN*> leafs_for_ccd, leafs_for_contact;
    void GenerateLeafs(unsigned edge_idx);
private:
    static ConcurrentPool<BVHN> BVHNLeafFactory;


};

#endif // MESHFRAGMENT_H
