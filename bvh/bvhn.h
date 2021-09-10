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

//    unsigned feature_idx;
    uint64_t feature_key;   // reference to the geometrical feature (edge) that this BVHN envelopes (if leaf)

    BVHN();
    void Build(std::vector<BVHN*> *bvs, int level_);
    void Update();
    void SelfCollide(std::vector<uint64_t> &broad_list);
    void Collide(BVHN *b, std::vector<uint64_t> &broad_list);

private:
    BVHN *child1, *child2;
};

#endif // BVHN_H
