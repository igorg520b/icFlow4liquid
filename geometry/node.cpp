#include <cmath>
#include "node.h"

icy::Node::Node()
{
    Reset();
}

void icy::Node::Reset()
{
    x_initial = xn = vn = xt = Eigen::Vector2d::Zero();
    eqId = locId = globId = -1;
    area = 0;
    pinned = false;
}

void icy::Node::ComputeEquationEntries(EquationOfMotionSolver &eq, SimParams &prms, double timeStep)
{
    if(this->eqId<0) return;

    double mass = area * prms.Density * prms.Thickness;

    // for nodes that are not pinned, add the lumped mass matrix to the quadratic term of the equation
    Eigen::Matrix2d M_nd = Eigen::Matrix2d::Identity()*mass;
    eq.AddToQ(eqId, eqId, M_nd);

    Eigen::Vector2d lambda_n = xt-x_hat;
    Eigen::Vector2d linear_term = M_nd*lambda_n;
    if(std::isnan(linear_term.x()) || std::isnan(linear_term.y()))
    {
        std::cout << "isnan node " << this->eqId << std::endl;
        throw std::runtime_error("isnan node");
    }

    eq.AddToC(eqId, linear_term);

    double const_term = lambda_n.dot(M_nd*lambda_n)/2;
    eq.AddToConstTerm(const_term);
}

