#ifndef BVHN_H
#define BVHN_H

#include <vector>
#include <utility>
#include "kdop8.h"
#include "ConcurrentPool.h"

#include "element.h"

namespace icy { class BVHN; }

class icy::BVHN
{
public:
    static ConcurrentPool<std::vector<BVHN*>> VectorFactory;
    static ConcurrentPool<BVHN> BVHNFactory;

    kDOP8 box;
    bool isLeaf;
    bool test_self_collision;   // can disable self-collision tests on fragments
    int level;

    // reference to the geometrical feature (edge) that this BVHN envelopes (if leaf)
    // only one of these pointers will be active
    const std::pair<Node*, Node*> *boundaryEdge;
    const Element *elem;
    bool isEdge() const {return boundaryEdge!=nullptr;}
    bool isElem() const {return elem!=nullptr;}

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
