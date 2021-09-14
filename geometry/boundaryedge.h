#ifndef BOUNDARYEDGE_H
#define BOUNDARYEDGE_H

#include "node.h"
#include <utility>

namespace icy {class BoundaryEdge;}

class icy::BoundaryEdge : public std::pair<Node*, Node*>
{
public:
    BoundaryEdge(Node* &nda, Node* &ndb) : std::pair<Node*,Node*>(nda,ndb), isActive(true) {};
    BoundaryEdge(Node* &nda, Node* &ndb, bool isInserted_) : std::pair<Node*,Node*>(nda,ndb), isInserted(isInserted_), isActive(false) {};
    uint64_t global_key() const {return Node::make_global_key(this->first,this->second);}
    uint64_t local_key() const {return Node::make_local_key(this->first,this->second);}
    bool isInserted = false;    // boundary was created via fracture
    bool isActive;       // boundary participates in collisions/interaction
    int VisualizationValue() const { return isInserted ? (isActive ? 2 : 1) : 0; }
};

#endif // BOUNDARYEDGE_H
