#ifndef ELEMENT123_H
#define ELEMENT123_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "parameters_sim.h"
#include "node.h"
#include "equationofmotionsolver.h"

namespace icy { class Element; class Node; }

class icy::Element
{
public:

    icy::Node* nds[3];          // initialized when the geometry is loaded or remeshed
    Eigen::Matrix2d PiMultiplier;   // multiplicative plasticity
    int group;

    Element() { Reset();}
    void Reset(void);
    void Reset(Node *nd0, Node *nd1, Node *nd2);
    void PrecomputeInitialArea();

    void AddToSparsityStructure(EquationOfMotionSolver &eq);
    bool ComputeEquationEntries(EquationOfMotionSolver &eq, SimParams &prms, double timeStep);
    void ComputeVisualizedVariables();  // Cauchy stress, Green strain, etc.
    bool PlasticDeformation(SimParams &prms, double timeStep);  // true if plastic deformation occurred

    double strain_energy_density;   // (not multiplied by the volume)

    Eigen::Matrix2d CauchyStress, GreenStrain;
    double principal_stress1, principal_stress2, max_shear_stress, hydrostatic_stress;
    double volume_change, velocity_divergence;
    double area_initial, area_current;

private:
    static Eigen::Matrix2d DDs[6]; // derivatives of Ds with respect to x1,y1,x2,y2,x3,y3
    static Eigen::Matrix<double,6,6> consistentMassMatrix;

    Eigen::Matrix2d Dm, DmInv;  // reference shape matrix
    Eigen::Matrix2d F;  // deformation gradient
    Eigen::Matrix2d P;  // First Piola-Kirchhoff stress tensor
};

#endif // ELEMENT123_H
