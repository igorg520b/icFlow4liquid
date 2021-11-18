#ifndef BASEELEMENT_H
#define BASEELEMENT_H

#include <cstdint>

namespace icy {class BaseElement; class Node; class Element; }

class icy::BaseElement
{
public:
    enum ElementType {BEdge, TElem, CZ};
    ElementType type;

    // for identifying disjoint regions and n-regions around crack tips
    bool traversed;

    virtual ~BaseElement() = default;

    // all derived classes keep some sort of list of incident elems and need a way to substitute one Element for another
    virtual void ReplaceAdjacentElem(const Element* originalElem, Element* insertedElem, uint8_t idx) = 0;

    // if BaseElement is CZ or BoudnaryEdge, then UpdateNodes re-creates the list of nodes who which this element is conected
    virtual void UpdateNodes() {};
};

#endif // BASEELEMENT_H
