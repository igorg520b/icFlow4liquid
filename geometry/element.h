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
    Element() { Reset();}


    icy::Node* nds[3];          // initialized when the geometry is loaded or remeshed
    int group;
    std::size_t gmshTag;        // tag of the element in the original gmsh representation

    Eigen::Matrix2d PiMultiplier;   // multiplicative plasticity
//    boost::container::small_vector<unsigned, 22> adj_elems;

    void Reset(void);
    void Reset(Node *nd0, Node *nd1, Node *nd2, std::size_t gmshTag);
    void PrecomputeInitialArea();

    void AddToSparsityStructure(EquationOfMotionSolver &eq);
    bool ComputeEquationEntries(EquationOfMotionSolver &eq, SimParams &prms, double timeStep);
    void ComputeVisualizedVariables();  // Cauchy stress, Green strain, etc.
    bool PlasticDeformation(SimParams &prms, double timeStep);  // true if plastic deformation occurred

    Eigen::Matrix2d CauchyStress, GreenStrain;
    double area_initial, area_current;
    double strain_energy_density;   // (not multiplied by the volume)
    double principal_stress1, principal_stress2, max_shear_stress, hydrostatic_stress;
    double volume_change, velocity_divergence;

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
    icy::Edge edges[3];

};

#endif // ELEMENT123_H
