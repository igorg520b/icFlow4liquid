#ifndef ELEMENT123_H
#define ELEMENT123_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "parameters_sim.h"
#include "node.h"
#include "edge.h"
#include "equationofmotionsolver.h"

namespace icy { struct Element; struct Node; }

struct icy::Element
{
    icy::Node* nds[3];          // initialized when the geometry is loaded or remeshed
    int group;
    std::size_t gmshTag;        // tag of the element in the original gmsh representation

    Eigen::Matrix2d PiMultiplier;   // multiplicative plasticity

    void Reset(void);
    void Initialize(Node *nd0, Node *nd1, Node *nd2);
    void PrecomputeInitialArea();

    void AddToSparsityStructure(EquationOfMotionSolver &eq) const;
    bool ComputeEquationEntries(EquationOfMotionSolver &eq, const SimParams &prms, double timeStep);
    void ComputeVisualizedVariables();  // Cauchy stress, Green strain, etc.
    bool PlasticDeformation(const SimParams &prms, double timeStep);  // true if plastic deformation occurred

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

private:
    static Eigen::Matrix2d DDs[6]; // derivatives of Ds with respect to x1,y1,x2,y2,x3,y3
    static Eigen::Matrix<double,6,6> consistentMassMatrix;


// FRACTURE ALGORITHM
public:
    icy::Element* incident_elems[3];    // nullptr or the element lying opposite of corresponding node
//    icy::Edge edges[3];
    unsigned traversal;     // used for identifying disjoint regions and n-regions around crack tips
    bool isBoundary(const short idx) const {return incident_elems[idx]==nullptr;}
    bool isOnBoundary(const Node* nd) const;
    bool isCWBoundary(const Node* nd) const;
    bool isCCWBoundary(const Node* nd) const;
    std::pair<Node*,Node*> CW_CCW_Node(const Node* nd) const;
//    const Edge& CWEdge(const Node* nd) const;   // clockwise edge
//    const Edge& CCWEdge(const Node* nd) const;  // counter-clockwise edge
//    const Edge& OppositeEdge(const Node* nd) const;

    bool containsNode(const Node* nd) const {return (nds[0]==nd || nds[1]==nd || nds[2]==nd);}
    Eigen::Vector2d getCenter() const {return (nds[0]->x_initial + nds[1]->x_initial + nds[2]->x_initial)/3.0;};
    void getIdxs(const icy::Node* nd, short &thisIdx, short &CWIdx, short &CCWIdx) const;
    Element* getAdjacentElementOppositeToNode(Node* nd);
    short getNodeIdx(const Node* nd) const;
    icy::Node* getOppositeNode(Node* nd0, Node* nd1);
    void ReplaceNode(Node* replaceWhat, Node* replaceWith);
    void DisconnectFromElem(Element* other);


private:
    constexpr static double threshold_area = 1e-7;
};

#endif // ELEMENT123_H
