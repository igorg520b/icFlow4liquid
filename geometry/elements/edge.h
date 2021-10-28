#if !defined(Q_MOC_RUN)
#ifndef EDGE_H
#define EDGE_H

#include <functional>
#include <Eigen/Core>

namespace icy {struct Edge; struct Node; struct Element;}

struct icy::Edge
{
    Edge() {}
    Edge(icy::Node* nd0, icy::Node* nd1);

    icy::Node* nds[2];
    icy::Element* elems[2];

    short edge_in_elem_idx[2];  // index in the range [0,3] indicating where the edge is located in the corresponding element
    bool isBoundary;    // belongs to only one element

    void AddElement(icy::Element* elem, short idx);
    Eigen::Vector2d getVec(const icy::Node* const from_node) const;  // edge as vector
    icy::Node* getOtherNode(const icy::Node* const nd) const;
};

#endif // EDGE_H
#endif
