#include "edge.h"
#include "node.h"
#include "element.h"
#include "cohesivezone.h"

icy::Edge::Edge(icy::Node* nd0, icy::Node* nd1)
{
    if(nd0->locId > nd1->locId) std::swap(nd0,nd1);
    nds[0] = nd0;
    nds[1] = nd1;
    elems[0] = elems[1] = nullptr;
}

bool icy::Edge::containsCZ() const
{
    return ((elems[0] != nullptr && elems[0]->type == BaseElement::ElementType::CZ) ||
            (elems[1] != nullptr && elems[1]->type == BaseElement::ElementType::CZ));
}


void icy::Edge::AddElement(icy::BaseElement* be, short idx)
{
    if(be->type == BaseElement::ElementType::TElem)
    {
        Element *elem = dynamic_cast<Element*>(be);
        CohesiveZone *cz = static_cast<icy::CohesiveZone*>(elem->incident_elems[idx]);

        if(elem->isEdgeCW(nds[0],nds[1]) && elems[1]==nullptr)
        {
            elems[1] = elem;
            edge_in_elem_idx[1] = idx;
            if(cz != nullptr) elems[0] = cz;
        }
        else if(!elem->isEdgeCW(nds[0],nds[1]) && elems[0]==nullptr)
        {
            elems[0] = elem;
            edge_in_elem_idx[0] = idx;
            if(cz != nullptr) elems[1] = cz;
        }
        else throw std::runtime_error("AddElement: mesh topology error");
    }
    else throw std::runtime_error("AddElement: mesh topology error cz");
}

/*
icy::Node* icy::Edge::getOtherNode(const icy::Node* const nd) const
{
    if(nd==nds[0]) return nds[1];
    else if(nd==nds[1]) return nds[0];
    else throw std::runtime_error("icy::Edge::getOtherNode: node does not belong to the edge");
}

Eigen::Vector2d icy::Edge::getVec(const icy::Node* const from_node) const
{
    const Node* other = getOtherNode(from_node);
    return other->xt - from_node->xt;
}
*/

