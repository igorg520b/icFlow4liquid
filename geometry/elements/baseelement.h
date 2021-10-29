#ifndef BASEELEMENT_H
#define BASEELEMENT_H

#include <cstdint>
//#include "node.h"
//#include "equationofmotionsolver.h"

namespace icy {class BaseElement; class Node;}

class icy::BaseElement
{
public:
    enum ElementType {BEdge, TElem, CZ};
    ElementType type;

    // for identifying disjoint regions and n-regions around crack tips
    bool traversed;
//    BaseElement **incident_begin = nullptr;
//    BaseElement **incident_end = nullptr;

    virtual ~BaseElement() = default;


//    virtual ElementType type() = 0;
//    Node *nds_begin, *nds_end;
//    virtual void AddToSparsityStructure(EquationOfMotionSolver &eq) const = 0;
};

#endif // BASEELEMENT_H
