#include "edge.h"
#include "node.h"
#include "element.h"

icy::Edge::Edge(icy::Node* nd0, icy::Node* nd1)
{
    if(nd0->locId > nd1->locId) std::swap(nd0,nd1);
    nds[0] = nd0;
    nds[1] = nd1;
    elems[0] = elems[1] = nullptr;
}

void icy::Edge::AddElement(icy::Element* elem, short idx)
{
    Eigen::Vector2d u = nds[1]->x_initial - nds[0]->x_initial;
    Eigen::Vector3d u3;
    u3 << u[0],u[1],0;

    icy::Node *opposite_node = elem->nds[idx];

    Eigen::Vector2d v = opposite_node->x_initial - nds[0]->x_initial;
    Eigen::Vector3d v3;
    v3 << v[0],v[1],0;

    bool elem1_isCCW = u3.cross(v3).z() > 0;

    if(elem1_isCCW && elems[0] == nullptr) { elems[0] = elem; edge_in_elem_idx[0] = idx; }
    else if(!elem1_isCCW && elems[1] == nullptr) { elems[1] = elem; edge_in_elem_idx[1] = idx; }
    else throw std::runtime_error("mesh topology error");
}


