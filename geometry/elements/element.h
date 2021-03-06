#ifndef ELEMENT123_H
#define ELEMENT123_H

#include <tuple>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "baseelement.h"
#include "parameters_sim.h"
#include "node.h"
#include "edge.h"
#include "equationofmotionsolver.h"

namespace icy { struct Element; struct Node; }

struct icy::Element : public icy::BaseElement
{
    icy::Node* nds[3];          // initialized when the geometry is loaded or remeshed
    std::bitset<8> group;
    icy::BaseElement* incident_elems[3];    // nullptr or the element lying opposite of corresponding node
    Eigen::Matrix2d PiMultiplier;   // multiplicative plasticity

    Element() { type = ElementType::TElem; }
    ~Element() = default;
    Element& operator=(Element&) = delete;

    void Reset(void);
    void Initialize(Node *nd0, Node *nd1, Node *nd2);
    void PrecomputeInitialArea();

    void AddToSparsityStructure(EquationOfMotionSolver &eq) const;
    bool ComputeEquationEntries(EquationOfMotionSolver &eq, const SimParams &prms, double timeStep);
    void ComputeVisualizedVariables();  // Cauchy stress, Green strain, etc.
    bool PlasticDeformation(const SimParams &prms, double timeStep);  // true if plastic deformation occurred
    void RecalculatePiMultiplierFromDeformationGradient(Eigen::Matrix2d F_tilda);

    Eigen::Matrix2d CauchyStress, GreenStrain;
    double area_initial, area_current;
    double strain_energy_density;   // (not multiplied by the volume)
    double principal_stress1, principal_stress2, max_shear_stress, hydrostatic_stress;
    double volume_change, velocity_divergence;
    double plasticity_gamma, plasticity_tau_ratio;

    Eigen::Matrix2d Dm, DmInv;  // reference shape matrix
    Eigen::Matrix2d F;  // deformation gradient
    Eigen::Matrix2d P;  // First Piola-Kirchhoff stress tensor

    Eigen::Matrix2d leftCauchyGreenDeformationTensor, qualityMetricTensor;
    double quality_measure_Wicke;

    Eigen::Matrix2d getF_at_n() { return getDs_at_n()*DmInv*PiMultiplier; }
    Eigen::Matrix2d getDs_at_n() { return (Eigen::Matrix2d() << nds[0]->xn-nds[2]->xn, nds[1]->xn-nds[2]->xn).finished(); }

private:
    static Eigen::Matrix2d DDs[6]; // derivatives of Ds with respect to x1,y1,x2,y2,x3,y3
    static Eigen::Matrix<double,6,6> consistentMassMatrix;

    void computeDm() { Dm << nds[0]->x_initial-nds[2]->x_initial, nds[1]->x_initial-nds[2]->x_initial; DmInv = Dm.inverse(); }

// FRACTURE ALGORITHM
public:
    std::pair<Node*,Node*> SplitElem(Node *nd, Node *nd0, Node *nd1, double where); // split the element by inserting a node between nd0 and nd1

    bool isBoundary() {return std::any_of(std::begin(nds),std::end(nds),[](Node *nd){return nd->isBoundary;});}

    uint8_t getNodeIdx(const Node* nd) const;
    uint8_t getEdgeIdx(const Node *nd1, const Node *nd2) const;
    std::pair<Node*,Node*> CW_CCW_Node(const Node* nd) const;

    bool isBoundaryEdge(const uint8_t idx) const {return incident_elems[idx]->type != ElementType::TElem;}
    bool isOnBoundary(const Node* nd) const;
    bool isCWBoundary(const Node* nd) const;
    bool isCCWBoundary(const Node* nd) const;
    bool isEdgeCW(const Node *nd1, const Node *nd2) const; // true if nd1-nd2 is oriented clockwise
    bool containsEdge(const Node *nd1, const Node *nd2) const;

    bool containsNode(const Node* nd) const {return (nds[0]==nd || nds[1]==nd || nds[2]==nd);}
    Eigen::Vector2d getCenter() const {return (nds[0]->x_initial + nds[1]->x_initial + nds[2]->x_initial)/3.0;};
    icy::Node* getOppositeNode(Node* nd0, Node* nd1);
    void ReplaceNode(Node* replaceWhat, Node* replaceWith);
    void ReplaceIncidentElem(const BaseElement* which, BaseElement* withWhat);
    void ReplaceAdjacentElem(const Element* originalElem, Element* insertedElem, uint8_t idx) override;

private:
    constexpr static double threshold_area = 1e-7;
    Node* SplitBoundaryElem(Node *nd, Node *nd0, Node *nd1, double where);
    Node* SplitNonBoundaryElem(Node *nd, Node *nd0, Node *nd1, double where);
    std::pair<icy::Node*,icy::Node*> SplitElemWithCZ(Node *nd, Node *nd0, Node *nd1, double where);

    BaseElement* getIncidentElementOppositeToNode(Node* nd);

};

#endif // ELEMENT123_H
