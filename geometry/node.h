#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers
#ifndef NODE_H
#define NODE_H

#include <Eigen/Core>

#include "equationofmotionsolver.h"
#include "parameters_sim.h"

namespace icy { class Node; }

class icy::Node
{
public:
    Node();
    void Reset();

    int locId, globId, eqId;       // sequential number of a node; identificator in the equation of motion (if not pinned)
    bool pinned;
    double area;        // area that the node "represents", for applying various forces

    Eigen::Vector2d x_initial;  // initial configuration
    Eigen::Vector2d xn, vn;     // position and velocity at step n
    Eigen::Vector2d xt;         // tentative coordinates
    Eigen::Vector2d x_hat;

    Eigen::Vector2d intended_position; // when manipulating via GUI during a running simulation

    void ComputeEquationEntries(EquationOfMotionSolver &eq, SimParams &prms, double timeStep);

};

#endif // NODE_H
#endif // Q_MOC_RUN
