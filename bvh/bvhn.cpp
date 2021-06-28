#include <stdexcept>
#include <cfloat>
#include <algorithm>
#include "bvhn.h"

icy::ConcurrentPool<std::vector<icy::BVHN*>> icy::BVHN::VectorFactory(50);
icy::ConcurrentPool<icy::BVHN> icy::BVHN::BVHNFactory(50000);

icy::BVHN::BVHN() {}

void icy::BVHN::Build(std::vector<BVHN*> *bvs, int level_)
{
    // TODO: ensure that the algorithm will also work with just one element in bvs

    if(level_ > 100) throw std::runtime_error("BVH level is over 100");
    level = level_;
    auto count = bvs->size();
    if(count == 0) throw new std::runtime_error("bvs->size==0 in BVHN::Initialize");
    else if(count == 1) throw new std::runtime_error("bvs->size==1 in BVHN::Initialize");

    isLeaf = false;

    box.Reset();
    for(auto const &bv : *bvs) box.Expand(bv->box); // expand box to the size of bvs collection

    std::vector<icy::BVHN*> *left = VectorFactory.take();
    std::vector<icy::BVHN*> *right = VectorFactory.take();
    left->clear();
    right->clear();
    left->reserve(bvs->size());
    right->reserve(bvs->size());

    if (box.dX >= box.dY)
    {
        for(auto const &bv : *bvs)
        {
            if(bv->box.ctrX < box.ctrX) left->push_back(bv);
            else right->push_back(bv);
        }

        // make sure that there is at least one element on each side
        if(left->size() == 0)
        {
            auto iter = std::min_element(right->begin(), right->end(),
                                         [](BVHN* b1, BVHN* b2) {return b1->box.ctrX < b2->box.ctrX;});

            // move "selected" from left to right
            left->push_back(*iter);
            right->erase(iter);
        }
        else if(right->size() == 0)
        {
            auto iter = std::max_element(left->begin(), left->end(),
                                         [](BVHN* b1, BVHN* b2) {return b1->box.ctrX < b2->box.ctrX;});

            // move selected from left to right
            right->push_back(*iter);
            left->erase(iter);
        }
    }
    else
    {
        for(auto const &bv : *bvs)
        {
            if(bv->box.ctrY < box.ctrY) left->push_back(bv);
            else right->push_back(bv);
        }

        // make sure that there is at least one element on each side
        if(left->size() == 0)
        {
            auto iter = std::min_element(right->begin(), right->end(),
                                         [](BVHN* b1, BVHN* b2) {return b1->box.ctrY < b2->box.ctrY;});

            // move "selected" from left to right
            left->push_back(*iter);
            right->erase(iter);
        }
        else if(right->size() == 0)
        {
            auto iter = std::max_element(left->begin(), left->end(),
                                         [](BVHN* b1, BVHN* b2) {return b1->box.ctrY < b2->box.ctrY;});

            // move selected from left to right
            right->push_back(*iter);
            left->erase(iter);
        }
    }


    if(left->size() == 1)
    {
        child1 = left->front();
        child1->level=level+1;
//        if(!child1->isLeaf) throw std::runtime_error("lone child is not leaf");
    }
    else if(left->size() > 1)
    {
        child1 = BVHNFactory.take();
        child1->test_self_collision = this->test_self_collision;
        child1->Build(left, level+1);
    }
    else throw std::runtime_error("left.size < 1");
    VectorFactory.release(left);

    if(right->size() == 1)
    {
        child2 = right->front();
        child2->level=level+1;
//        if(!child2->isLeaf) throw std::runtime_error("lone child is not leaf");
    }
    else if(right->size() > 1)
    {
        child2 = BVHNFactory.take();
        child2->test_self_collision = this->test_self_collision;
        child2->Build(right, level+1);
    }
    else throw std::runtime_error("right.size < 1");
    VectorFactory.release(right);
}

void icy::BVHN::Update()
{
    // assume that leaves are already updated
    if(isLeaf) throw std::runtime_error("Updated called on BVHN leaf");
    if(!child1->isLeaf) child1->Update();
    if(!child2->isLeaf) child2->Update();
    box.Reset();
    box.Expand(child1->box);
    box.Expand(child2->box);
}

void icy::BVHN::SelfCollide(std::vector<unsigned> &broad_list)
{
    if (isLeaf) return;
    if(child1->test_self_collision) child1->SelfCollide(broad_list);
    if(child2->test_self_collision) child2->SelfCollide(broad_list);
    child1->Collide(child2, broad_list);
}

void icy::BVHN::Collide(BVHN *b, std::vector<unsigned> &broad_list)
{
    if(!box.Overlaps(b->box)) return;
    if (this->isLeaf && b->isLeaf)
    {
        broad_list.push_back(feature_idx);
        broad_list.push_back(b->feature_idx);
    }
    else if (this->isLeaf)
    {
        Collide(b->child1, broad_list);
        Collide(b->child2, broad_list);
    }
    else
    {
        b->Collide(child1, broad_list);
        b->Collide(child2, broad_list);
    }
}

