#ifndef BASEELEMENT_H
#define BASEELEMENT_H

#include <cstdint>
#include "node.h"
#include "equationofmotionsolver.h"

namespace icy {class BaseElement; class Node;}

class icy::BaseElement
{
public:
    virtual uint8_t type() = 0;
    Node *nds_begin, *nds_end;

    virtual void AddToSparsityStructure(EquationOfMotionSolver &eq) const = 0;
};

#endif // BASEELEMENT_H
