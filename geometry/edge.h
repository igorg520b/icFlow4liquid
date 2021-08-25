#if !defined(Q_MOC_RUN)
#ifndef EDGE_H
#define EDGE_H

#include <functional>
#include "node.h"

namespace icy {struct Edge; struct Node; struct Element;}

struct icy::Edge
{
    Edge() {}
    Edge(icy::Node* nd0, icy::Node* nd1);

    icy::Node* nds[2];
    icy::Element* elems[2];

    short edge_in_elem_idx[2];  // index in the range [0,3] indicating where the edge is located in the corresponding element
    bool isBoundary;    // belongs to only one element
    bool toSplit = false; // for mesh splitting

    void AddElement(icy::Element* elem, short idx);
    Eigen::Vector2d getVec(const icy::Node* const from_node) const;  // edge as vector
    icy::Node* getOtherNode(const icy::Node* const nd) const;

/*
    // to be used if creating unordered_map
    bool operator==(const Edge &other) const
    {
        return (nds[0]==other.nds[0] && nds[1]==other.nds[1]);
    }
    struct Hasher
    {
        size_t operator()(const Edge & obj) const
        {
            return std::hash<int>()(icy::Node::make_key(obj.nds[0],obj.nds[1]));
        }
    };
*/
};

#endif // EDGE_H
#endif
