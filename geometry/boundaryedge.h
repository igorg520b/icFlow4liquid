#ifndef BOUNDARYEDGE_H
#define BOUNDARYEDGE_H

#include <utility>
#include "node.h"
#include "element.h"

namespace icy {struct BoundaryEdge; struct Element;}

struct icy::BoundaryEdge
{
    Node* nds[2];
    Element *elem = nullptr;            // element to which this boundary belongs
    uint8_t edge_idx;                   // side in the element to which this boundary is attached
    uint8_t status;                     // color/type for visualization


    BoundaryEdge(Element *elem_, uint8_t edge_idx_, uint8_t status_ = 0) : elem(elem_), edge_idx(edge_idx_), status(status_)
    {
        UpdateNodes();
    }

    // infer nds[2] from elem and edge_idx
    void UpdateNodes()
    {
        uint8_t nd1Idx = (edge_idx+2)%3;
        uint8_t nd2Idx = (edge_idx+1)%3;
        MeshFragment *fr = elem->nds[0]->fragment;
        nds[0] = elem->nds[nd1Idx];
        nds[1] = elem->nds[nd2Idx];
    }

    BoundaryEdge(Node *nd1, Node *nd2, Element *elem_ = nullptr, uint8_t status_ = 0) : elem(elem_), status(status_)
    {
        if(elem_!=nullptr)
        {
            // select the order of nd1, nd2, depending on which side is elem
            if(!elem_->isEdgeCW(nd1,nd2)) std::swap(nd1,nd2);

        }
        nds[0] = nd1;
        nds[1] = nd2;
    }

    BoundaryEdge() = delete;
    BoundaryEdge(const BoundaryEdge &other) = default;
    BoundaryEdge& operator=(const BoundaryEdge& other) = default;
    ~BoundaryEdge() = default;

    bool isDeformable() const {return elem!=nullptr;}
};

#endif // BOUNDARYEDGE_H
