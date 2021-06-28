#ifndef BVHN_H
#define BVHN_H

#include <vector>
#include "kdop8.h"
#include "ConcurrentPool.h"

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

    unsigned feature_idx;

    BVHN();
    void Build(std::vector<BVHN*> *bvs, int level_);
    void Update();
    void SelfCollide(std::vector<unsigned> &broad_list);
    void Collide(BVHN *b, std::vector<unsigned> &broad_list);

private:
    BVHN *child1, *child2;
};

#endif // BVHN_H
