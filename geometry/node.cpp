#include <cmath>
#include "node.h"

void icy::Node::Reset()
{
    x_initial = xn = vn = xt = Eigen::Vector2d::Zero();
    eqId = locId = globId = -1;
    area = 0;
    pinned = false;
    spring_attached = 0;
    group.reset();
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

void icy::Node::AddSpringEntries(EquationOfMotionSolver &eq, SimParams &prms, double h, Eigen::Vector2d &spring)
{
    if(eqId<0 || spring_attached<1) return;
    double k = prms.YoungsModulus/1000;
    Eigen::Vector2d spr = spring_attachment_position+spring;
    double hsqk = h*h*k;

    Eigen::Vector2d fd;
    fd = (spr-xt);
    eq.AddToConstTerm(hsqk*fd.dot(fd)/2);

    eq.AddToC(eqId,hsqk*(xt-spr));

    Eigen::Matrix2d Hessian = hsqk*Eigen::Matrix2d::Identity();
    eq.AddToQ(eqId,eqId,Hessian);


}
