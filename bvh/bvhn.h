#ifndef BVHN_H
#define BVHN_H

#include <vector>
#include <utility>
#include "kdop8.h"
#include "ConcurrentPool.h"

#include "element.h"
#include "boundaryedge.h"

namespace icy { class BVHN; }

class icy::BVHN
{
public:
    static ConcurrentPool<std::vector<BVHN*>> VectorFactory;
    static ConcurrentPool<BVHN> BVHNFactory;

    kDOP8 box;
    bool isLeaf() {return boundaryEdge!=nullptr;};
    bool test_self_collision;   // can disable self-collision tests on fragments
    int level;

    // reference to the geometrical feature (edge) that this BVHN envelopes (if leaf)
    const BoundaryEdge *boundaryEdge;

    BVHN();
    void Build(std::vector<BVHN*> *bvs, int level_);
    void Update();
    void SelfCollide(std::vector<std::pair<BVHN*,BVHN*>> &broad_list);
    void Collide(BVHN *b, std::vector<std::pair<BVHN*,BVHN*>> &broad_list);

    void Expand_CCD(float distance_threshold);

private:
    BVHN *child1, *child2;
};

#endif // BVHN_H
