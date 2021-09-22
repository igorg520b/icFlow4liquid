#include <stdexcept>
#include <cfloat>
#include <algorithm>
#include <iostream>
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
    if(count == 0) throw std::runtime_error("bvs->size==0 in icy::BVHN::Build");
    else if(count == 1) throw std::runtime_error("bvs->size==1 in icy::BVHN::Build");

    isLeaf = false;

    box.Reset();
    for(auto const &bv : *bvs) box.Expand(bv->box); // expand box to the size of bvs collection

    std::vector<icy::BVHN*> *left = VectorFactory.take();
    std::vector<icy::BVHN*> *right = VectorFactory.take();
    left->clear();
    right->clear();

    std::vector<icy::BVHN*>::iterator iter;
    if (box.dX >= box.dY)
    {
        float ctrX = box.ctrX;
        iter = std::partition(bvs->begin(),bvs->end(),[ctrX](const BVHN *bv){return bv->box.ctrX < ctrX;});
    }
    else
    {
        float ctr = box.ctrY;
        iter = std::partition(bvs->begin(),bvs->end(),[ctr](const BVHN *bv){return bv->box.ctrY < ctr;});
    }
    if(iter==bvs->begin()) iter++;
    else if(iter==bvs->end()) iter--;
    left->resize(std::distance(bvs->begin(),iter));
    right->resize(std::distance(iter,bvs->end()));
    std::copy(bvs->begin(),iter,left->begin());
    std::copy(iter,bvs->end(),right->begin());

    if(left->size() == 1)
    {
        child1 = left->front();
        child1->level=level+1;
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
    if(isLeaf) throw std::runtime_error("Update called on BVHN leaf");
    if(!child1->isLeaf) child1->Update();
    if(!child2->isLeaf) child2->Update();
    box.Reset();
    box.Expand(child1->box);
    box.Expand(child2->box);
}

void icy::BVHN::SelfCollide(std::vector<std::pair<BVHN*,BVHN*>> &broad_list)
{
    if (isLeaf || !test_self_collision) return;
    child1->SelfCollide(broad_list);
    child2->SelfCollide(broad_list);
    child1->Collide(child2, broad_list);
}

void icy::BVHN::Collide(BVHN *b, std::vector<std::pair<BVHN*,BVHN*>> &broad_list)
{
    if(!box.Overlaps(b->box)) return;
    if (this->isLeaf && b->isLeaf)
    {
        broad_list.emplace_back(this,b);
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

void icy::BVHN::Expand_CCD(float distance_threshold)
{
    box.Reset();
    if(!isLeaf) throw std::runtime_error("Expand_CCD: not leaf");
    if(boundaryEdge != nullptr)
    {
        Node *nd1 = boundaryEdge->first;
        Node *nd2 = boundaryEdge->second;
        box.Expand(nd1->xt[0], nd1->xt[1]);
        box.Expand(nd2->xt[0], nd2->xt[1]);
        box.Expand(nd1->xn[0], nd1->xn[1]);
        box.Expand(nd2->xn[0], nd2->xn[1]);
    }
    else if(elem != nullptr)
    {
        for(int i=0;i<3;i++)
        {
            Node *nd = elem->nds[i];
            box.Expand(nd->xt[0], nd->xt[1]);
            box.Expand(nd->xn[0], nd->xn[1]);
        }
    }
    else throw std::runtime_error("Expand_CCD: feature not set");
    box.ExpandBy(distance_threshold);
}
