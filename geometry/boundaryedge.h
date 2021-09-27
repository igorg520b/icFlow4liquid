#ifndef BOUNDARYEDGE_H
#define BOUNDARYEDGE_H

#include <utility>
#include "node.h"

namespace icy {struct BoundaryEdge; struct Element;}

struct icy::BoundaryEdge : public std::pair<Node*,Node*>
{
    BoundaryEdge(Node *nd1, Node *nd2, Element *elem_) : elem(elem_)
    {
        if(nd1->locId > nd2->locId) std::swap(nd1,nd2);
        this->first = nd1;
        this->second = nd2;
    }
    BoundaryEdge() {}

    Element *elem = nullptr;
    bool isDeformable() const {return elem!=nullptr;}

    bool operator==(const std::pair<Node*,Node*> &other)
    {
        return ((first == other.first && second == other.second) || (first == other.second && second == other.first));
    }

    std::pair<Eigen::Vector2d,Eigen::Vector2d> offsetEdge(double offsetMag)
    {
        Eigen::Vector2d originalDirection = second->xt - first->xt;
        Eigen::Vector2d offsetV{originalDirection[1],-originalDirection[0]};
        offsetV.normalize();
        offsetV *= offsetMag;
        Eigen::Vector2d resultA = first->xt + offsetV;
        Eigen::Vector2d resultB = second->xt + offsetV;
        return {resultA,resultB};
    }
};

#endif // BOUNDARYEDGE_H
