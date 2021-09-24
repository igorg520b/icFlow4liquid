#ifndef BOUNDARYEDGE_H
#define BOUNDARYEDGE_H

#include <utility>
#include "node.h"

namespace icy {class BoundaryEdge; class Element;}

struct icy::BoundaryEdge : public std::pair<Node*,Node*>
{
    Element *elem = nullptr;

    bool operator==(const BoundaryEdge &other)
    {
        return ((first == other.first && second == other.second) || (first == other.second && second == other.first));
    }
};

#endif // BOUNDARYEDGE_H
