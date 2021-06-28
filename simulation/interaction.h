#ifndef INTERACTION_H
#define INTERACTION_H

#include<Eigen/Core>
#include "equationofmotionsolver.h"
#include "parameters_sim.h"

namespace icy { class Interaction; class Node;}

class icy::Interaction
{
public:
    Node *ndA, *ndB, *ndP;
    Eigen::Vector2d D;

    void AddToSparsityStructure(EquationOfMotionSolver &eq);
    void Evaluate(EquationOfMotionSolver &eq, SimParams &prms, double h);
    double static SegmentPointDistance(Eigen::Vector2d A, Eigen::Vector2d B, Eigen::Vector2d P, Eigen::Vector2d &D);

private:
    void static distance(Eigen::Vector2d (&p)[4], double &d, double &t, Eigen::Matrix<double,6,1> &Dd, Eigen::Matrix<double,6,6> &DDd);

    void static potential(double dHat, double d, Eigen::Matrix<double,6,1> &Dd, Eigen::Matrix<double,6,6> &DDd,
                          double &p, Eigen::Matrix<double,6,1> &Dp, Eigen::Matrix<double,6,6> &DDp);
};

#endif // INTERACTION_H
