#ifndef BOUNDARYEDGE_H
#define BOUNDARYEDGE_H

#include <utility>
#include "node.h"
#include "element.h"

namespace icy {struct BoundaryEdge; struct Element;}

struct icy::BoundaryEdge
{
    Element *elem = nullptr;
    BoundaryEdge(Node *nd1, Node *nd2, Element *elem_ = nullptr) : elem(elem_)
    {
        if(elem_!=nullptr)
        {
            // select the order of nd1, nd2, depending on which side is elem
            if(!elem_->isEdgeCW(nd1,nd2)) std::swap(nd1,nd2);

        }
        this->vertices.first = nd1;
        this->vertices.second = nd2;
    }

    BoundaryEdge() = delete;
    BoundaryEdge(const BoundaryEdge &other) = default;
    BoundaryEdge& operator=(const BoundaryEdge& other) = default;
    ~BoundaryEdge() = default;

    std::pair<Node*,Node*> vertices;
    bool isDeformable() const {return elem!=nullptr;}

    bool operator==(const icy::BoundaryEdge &other)
    {
        return ((vertices.first == other.vertices.first && vertices.second == other.vertices.second)||
                (vertices.first == other.vertices.second && vertices.second == other.vertices.first));
    }

    std::pair<Eigen::Vector2d,Eigen::Vector2d> offsetEdge(double offsetMag)
    {
        Eigen::Vector2d originalDirection = vertices.second->xt - vertices.first->xt;
        Eigen::Vector2d offsetV{originalDirection[1],-originalDirection[0]};
        offsetV.normalize();
        offsetV *= offsetMag;
        Eigen::Vector2d resultA = vertices.first->xt + offsetV;
        Eigen::Vector2d resultB = vertices.second->xt + offsetV;
        return {resultA,resultB};
    }
};

#endif // BOUNDARYEDGE_H
