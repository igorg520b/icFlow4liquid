#include <cmath>
#include "node.h"

void icy::Node::Reset()
{
    x_material = x_initial = xn = vn = xt = Eigen::Vector2d::Zero();
    eqId = locId = globId = -1;
    area = 0;
    pinned = false;
    spring_attached = 0;
    group.reset();
}

void icy::Node::Reset(int locId_, double x, double y)
{
    vn = Eigen::Vector2d::Zero();
    locId = locId_;
    x_initial << x,y;
    x_material = intended_position = xt = xn = x_initial;
    eqId = globId = -1;
    area = 0;
    pinned = false;
    spring_attached = 0;
    group.reset();
}



void icy::Node::ComputeEquationEntries(EquationOfMotionSolver &eq, SimParams &prms, double timeStep)
{
    if(this->eqId<0) return;

    double mass = area * prms.Density * prms.Thickness;
    //if(mass<1e-15) throw std::runtime_error("node's mass is zero");

    // for nodes that are not pinned, add the lumped mass matrix to the quadratic term of the equation
    Eigen::Matrix2d M_nd = Eigen::Matrix2d::Identity()*mass;

    Eigen::Vector2d lambda_n = xt-x_hat;
    Eigen::Vector2d linear_term = M_nd*lambda_n;

    double const_term = lambda_n.dot(M_nd*lambda_n)/2;
    eq.AddToEquation(const_term, linear_term, M_nd, eqId);
}

void icy::Node::AddSpringEntries(EquationOfMotionSolver &eq, SimParams &prms, double h, Eigen::Vector2d &spring)
{
    if(eqId<0 || spring_attached<1) return;
    double k = prms.YoungsModulus/1000;
    Eigen::Vector2d spr = spring_attachment_position+spring;
    double hsqk = h*h*k;

    Eigen::Vector2d fd = (spr-xt);

    double const_term = hsqk*fd.dot(fd)/2;
    Eigen::Vector2d linear_term = hsqk*(xt-spr);
    Eigen::Matrix2d quadratic_term = hsqk*Eigen::Matrix2d::Identity();
    eq.AddToEquation(const_term, linear_term, quadratic_term, eqId);
}
