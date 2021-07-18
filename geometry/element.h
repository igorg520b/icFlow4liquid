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
    icy::Element* adj_elems[3]; // nullptr if no adjacent element
    int group;
    bool fluid;
    double area_initial;

    Element() { Reset();}
    void Reset(void);
    
    void PrecomputeInitialArea();

    void AddToSparsityStructure(EquationOfMotionSolver &eq);
    bool ComputeEquationEntries(EquationOfMotionSolver &eq, SimParams &prms, double timeStep);
    void EvaluateVelocityDivergence();

    double strain_energy_density;   // (not multiplied by volume!)
    Eigen::Matrix<double, 6, 1> DE;    // energy gradient
    Eigen::Matrix<double, 6, 6> HE; // energy hessian

    Eigen::Matrix2d CauchyStress;
    double principal_stress1, principal_stress2, max_shear_stress;
    double volume_change, velocity_divergence;

private:
    void SpringModel(EquationOfMotionSolver &eq, SimParams &prms, double timeStep, Node *nd1, Node *nd2);
    bool NeoHookeanElasticity(EquationOfMotionSolver &eq, SimParams &prms, double h);

};

#endif // ELEMENT123_H
