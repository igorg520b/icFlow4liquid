#if !defined(Q_MOC_RUN)
#ifndef EDGE_H
#define EDGE_H

#include <functional>
#include <Eigen/Core>

namespace icy {struct Edge; struct Node; struct Element; class BaseElement;}

struct icy::Edge
{
    Edge() {}
    Edge(icy::Node* nd0, icy::Node* nd1);

    icy::Node* nds[2];
    icy::BaseElement* elems[2];

    short edge_in_elem_idx[2];  // index in the range [0,3] indicating where the edge is located in the corresponding element
    bool isBoundary() const {return (elems[0] == nullptr || elems[1] == nullptr);};    // belongs to only one element
    bool containsCZ() const;

    void AddElement(icy::BaseElement* elem, short idx);
};

#endif // EDGE_H
#endif
