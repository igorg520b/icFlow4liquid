#include "interaction.h"
#include "node.h"

double icy::Interaction::SegmentPointDistance(Eigen::Vector2d A, Eigen::Vector2d B, Eigen::Vector2d P, Eigen::Vector2d &D)
{
    Eigen::Vector2d seg = B-A;
    Eigen::Vector2d v = P-A;
    double t = v.dot(seg)/seg.squaredNorm();
    t = std::clamp(t, 0.0, 1.0);
    D = A+seg*t;
    double dist = (D-P).norm();
    return dist;
}

void icy::Interaction::AddToSparsityStructure(EquationOfMotionSolver &eq) const
{
    int ids[] {ndA->eqId, ndB->eqId, ndP->eqId};
    eq.AddEntriesToStructure(std::begin(ids), std::end(ids));
}

void icy::Interaction::Evaluate(EquationOfMotionSolver &eq, SimParams &prms, double h)
{
    double dHat = prms.InteractionDistance;
    double k = prms.Kappa*h*h;

    Eigen::Vector2d pts[4];
    pts[0]=ndA->xt;
    pts[1]=ndB->xt;
    pts[2]=ndP->xt;

    double dist, proj_t;
    Eigen::Matrix<double,6,1> Dd;
    Eigen::Matrix<double,6,6> DDd;

    distance(pts, dist, proj_t, Dd, DDd);

    double p;
    Eigen::Matrix<double,6,1> Dp;
    Eigen::Matrix<double,6,6> DDp;
    if(pushing) potential_pushing(dHat, dist, Dd, DDd, p, Dp, DDp);
    else potential_pulling(dHat, dist, Dd, DDd, p, Dp, DDp);

    p*=k;
    Dp*=k;
    DDp*=k;

    // assemble the equation of motion
    eq.AddToEquation(Dp.data(), DDp.data(), {ndA->eqId,ndB->eqId,ndP->eqId});
}

void icy::Interaction::potential_pushing(double dHat, double d, Eigen::Matrix<double,6,1> &Dd, Eigen::Matrix<double,6,6> &DDd,
                      double &p, Eigen::Matrix<double,6,1> &Dp, Eigen::Matrix<double,6,6> &DDp)
{
    if(d<dHat)
    {
        double d_dHat = d-dHat;
        double d_dHat_sq = d_dHat*d_dHat;
        double dLog = std::log(d/dHat);
        p = -d_dHat_sq*dLog;
        double p_prime = -d_dHat_sq/d - 2*d_dHat*dLog;
        double p_double_prime = -4*d_dHat/d + d_dHat_sq/(d*d) - 2*dLog;
        Dp = p_prime*Dd;
        DDp = p_double_prime*Dd*Dd.transpose() + p_prime*DDd;
    }
    else
    {
        p=0;
        Dp=Eigen::Matrix<double,6,1>::Zero();
        DDp=Eigen::Matrix<double,6,6>::Zero();
    }
}

void icy::Interaction::potential_pulling(double dHat, double s, Eigen::Matrix<double,6,1> &Dd, Eigen::Matrix<double,6,6> &DDd,
                      double &p, Eigen::Matrix<double,6,1> &Dp, Eigen::Matrix<double,6,6> &DDp)
{
    constexpr double threshold = 0.9;
    if(s>dHat*threshold) s=dHat*threshold;
    double dLog = std::log((dHat-s)/dHat);
    double s_sq = s*s;
    p = -s_sq*dLog;
    double p_prime = s_sq/(dHat-s) -2*s*dLog;
    double p_double_prime = s_sq/((dHat-s)*(dHat-s)) +4*s/(dHat-s) - 2*dLog;
    Dp = p_prime*Dd;
    DDp = p_double_prime*Dd*Dd.transpose() + p_prime*DDd;

/*
    if(s<dHat)
    {


    }
    else
    {
        p=0;
        Dp=Eigen::Matrix<double,6,1>::Zero();
        DDp=Eigen::Matrix<double,6,6>::Zero();
    }
    */
}


void icy::Interaction::distance(Eigen::Vector2d (&p)[4], double &d, double &t,
    Eigen::Matrix<double,6,1> &Dd, Eigen::Matrix<double,6,6> &DDd)
{
    double x0 = p[0].x();
    double y0 = p[0].y();
    double x1 = p[1].x();
    double y1 = p[1].y();
    double x2 = p[2].x();
    double y2 = p[2].y();

    double t_num = (-x0 + x1)*(-x0 + x2) + (-y0 + y1)*(-y0 + y2); // (p2 - p0).(p1 - p0)
    double t_den = (-x0 + x1)*(-x0 + x1) + (-y0 + y1)*(-y0 + y1); // (p1 - p0).(p1 - p0)
    t = t_num/t_den;

    Eigen::Matrix<double,6,1> Dtt, Ddsq;
    Eigen::Matrix<double,6,6> DDtt, DDdsq;

    if(t>1)
    {
        t=1;
        Dtt = Eigen::Matrix<double,6,1>::Zero();
        DDtt = Eigen::Matrix<double,6,6>::Zero();
    }
    else if(t<0)
    {
        t=0;
        Dtt = Eigen::Matrix<double,6,1>::Zero();
        DDtt = Eigen::Matrix<double,6,6>::Zero();
    }
    else
    {
        Dtt << ((2*x0 - x1 - x2)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1)) +
               2*(-x0 + x1)*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)))/
             (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
               ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
            (2*(-y0 + y1)*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)) +
               ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(2*y0 - y1 - y2))/
             (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
               ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
            ((-x0 + x2)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1)) -
               2*(-x0 + x1)*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)))/
             (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
               ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
            (-2*(-y0 + y1)*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)) +
               ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y2))/
             (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
               ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
            (-x0 + x1)/((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1)),
            (-y0 + y1)/((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));

        DDtt << (2*(-(x1*x1*x1*x1) + x0*x0*x0*(x1 - x2) + x1*x1*x1*x2 -
                   3*x1*x2*((y0 - y1)*(y0 - y1)) +
                   3*x0*(x1*x1*x1 - x1*x1*x2 + x2*((y0 - y1)*(y0 - y1)) -
                      x1*(y0 - y1)*(y0 + y1 - 2*y2)) -
                   3*(x0*x0)*(x1*x1 - x1*x2 - (y0 - y1)*(y1 - y2)) + 3*(x1*x1)*(y0 - y1)*(y0 - y2) -
                   (y0 - y1)*(y0 - y1)*(y0 - y1)*(y1 - y2)))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (2*((2*x0 - x1 - x2)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y1) +
                   4*(-x0 + x1)*(-y0 + y1)*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)) +
                   (-x0 + x1)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(2*y0 - y1 - y2)))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (-(x0*x0*x0*x0) + x1*x1*x1*x1 - 2*(x1*x1*x1)*x2 + 2*(x0*x0*x0)*(x1 + x2) +
                 6*x1*x2*((y0 - y1)*(y0 - y1)) -
                 2*x0*(x1*x1*x1 - 3*(x1*x1)*x2 + 3*x2*((y0 - y1)*(y0 - y1)) -
                    3*x1*(y0 - y1)*(y0 + y1 - 2*y2)) - 6*(x0*x0)*(x1*x2 + (y0 - y1)*(y1 - y2)) +
                 (y0 - y1)*(y0 - y1)*(y0 - y1)*(y0 + y1 - 2*y2) - 6*(x1*x1)*(y0 - y1)*(y0 - y2))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (2*(-((2*x0 - x1 - x2)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y1)) -
                   4*(-x0 + x1)*(-y0 + y1)*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)) +
                   (-x0 + x1)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y2)))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (x0*x0 - 2*x0*x1 + x1*x1 - (y0 - y1)*(y0 - y1))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (2*(-x0 + x1)*(-y0 + y1))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (2*((2*x0 - x1 - x2)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y1) +
                   4*(-x0 + x1)*(-y0 + y1)*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)) +
                   (-x0 + x1)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(2*y0 - y1 - y2)))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (2*(x1*x1*x1*x1 - x1*x1*x1*x2 + x0*x0*x0*(-x1 + x2) + 3*x1*x2*((y0 - y1)*(y0 - y1)) -
                   3*x0*(x1*x1*x1 - x1*x1*x2 + x2*((y0 - y1)*(y0 - y1)) -
                      x1*(y0 - y1)*(y0 + y1 - 2*y2)) +
                   3*(x0*x0)*(x1*x1 - x1*x2 - (y0 - y1)*(y1 - y2)) - 3*(x1*x1)*(y0 - y1)*(y0 - y2) +
                   (y0 - y1)*(y0 - y1)*(y0 - y1)*(y1 - y2)))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (2*((-x0 + x2)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y1) -
                   4*(-x0 + x1)*(-y0 + y1)*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)) -
                   (-x0 + x1)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(2*y0 - y1 - y2)))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (x0*x0*x0*x0 - x1*x1*x1*x1 + 2*(x1*x1*x1)*x2 - 2*(x0*x0*x0)*(x1 + x2) -
                 6*x1*x2*((y0 - y1)*(y0 - y1)) +
                 2*x0*(x1*x1*x1 - 3*(x1*x1)*x2 + 3*x2*((y0 - y1)*(y0 - y1)) -
                    3*x1*(y0 - y1)*(y0 + y1 - 2*y2)) + 6*(x0*x0)*(x1*x2 + (y0 - y1)*(y1 - y2)) -
                 (y0 - y1)*(y0 - y1)*(y0 - y1)*(y0 + y1 - 2*y2) + 6*(x1*x1)*(y0 - y1)*(y0 - y2))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (2*(-x0 + x1)*(-y0 + y1))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (-(x0*x0) + 2*x0*x1 - x1*x1 + (y0 - y1)*(y0 - y1))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (-(x0*x0*x0*x0) + x1*x1*x1*x1 - 2*(x1*x1*x1)*x2 + 2*(x0*x0*x0)*(x1 + x2) +
                 6*x1*x2*((y0 - y1)*(y0 - y1)) -
                 2*x0*(x1*x1*x1 - 3*(x1*x1)*x2 + 3*x2*((y0 - y1)*(y0 - y1)) -
                    3*x1*(y0 - y1)*(y0 + y1 - 2*y2)) - 6*(x0*x0)*(x1*x2 + (y0 - y1)*(y1 - y2)) +
                 (y0 - y1)*(y0 - y1)*(y0 - y1)*(y0 + y1 - 2*y2) - 6*(x1*x1)*(y0 - y1)*(y0 - y2))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (2*((-x0 + x2)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y1) -
                   4*(-x0 + x1)*(-y0 + y1)*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)) -
                   (-x0 + x1)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(2*y0 - y1 - y2)))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (2*(-2*(-x0 + x1)*(-x0 + x2)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1)) +
                   4*((x0 - x1)*(x0 - x1))*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)) -
                   ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                    (x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2))))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (2*(-((-x0 + x2)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y1)) +
                   4*(-x0 + x1)*(-y0 + y1)*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)) -
                   (-x0 + x1)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y2)))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (-(x0*x0) + 2*x0*x1 - x1*x1 + (y0 - y1)*(y0 - y1))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (-2*(-x0 + x1)*(-y0 + y1))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (2*(-((2*x0 - x1 - x2)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y1)) -
                   4*(-x0 + x1)*(-y0 + y1)*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)) +
                   (-x0 + x1)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y2)))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (x0*x0*x0*x0 - x1*x1*x1*x1 + 2*(x1*x1*x1)*x2 - 2*(x0*x0*x0)*(x1 + x2) -
                 6*x1*x2*((y0 - y1)*(y0 - y1)) +
                 2*x0*(x1*x1*x1 - 3*(x1*x1)*x2 + 3*x2*((y0 - y1)*(y0 - y1)) -
                    3*x1*(y0 - y1)*(y0 + y1 - 2*y2)) + 6*(x0*x0)*(x1*x2 + (y0 - y1)*(y1 - y2)) -
                 (y0 - y1)*(y0 - y1)*(y0 - y1)*(y0 + y1 - 2*y2) + 6*(x1*x1)*(y0 - y1)*(y0 - y2))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (2*(-((-x0 + x2)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y1)) +
                   4*(-x0 + x1)*(-y0 + y1)*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)) -
                   (-x0 + x1)*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y2)))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (2*(-(((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                      (x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2))) +
                   4*((y0 - y1)*(y0 - y1))*(x0*x0 + x1*x2 - x0*(x1 + x2) + (y0 - y1)*(y0 - y2)) -
                   2*((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*(-y0 + y1)*(-y0 + y2)))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (-2*(-x0 + x1)*(-y0 + y1))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (x0*x0 - 2*x0*x1 + x1*x1 - (y0 - y1)*(y0 - y1))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (x0*x0 - 2*x0*x1 + x1*x1 - (y0 - y1)*(y0 - y1))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (2*(-x0 + x1)*(-y0 + y1))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (-(x0*x0) + 2*x0*x1 - x1*x1 + (y0 - y1)*(y0 - y1))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (-2*(-x0 + x1)*(-y0 + y1))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),0,0,
              (2*(-x0 + x1)*(-y0 + y1))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (-(x0*x0) + 2*x0*x1 - x1*x1 + (y0 - y1)*(y0 - y1))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),
              (-2*(-x0 + x1)*(-y0 + y1))/
               (((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))*
                 ((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1))),
              (x0*x0 - 2*x0*x1 + x1*x1 - (y0 - y1)*(y0 - y1))/
               ((x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))*
                 (x0*x0 - 2*x0*x1 + x1*x1 + (y0 - y1)*(y0 - y1))),0,0;
    }

    double dsq = ((1 - t)*x0 + t*x1 - x2)*((1 - t)*x0 + t*x1 - x2) +  ((1 - t)*y0 + t*y1 - y2)*((1 - t)*y0 + t*y1 - y2);

    Ddsq << 2*(y0 - y1)*(-y0 + t*(y0 - y1) + y2)*Dtt(0) +
            2*(-x0 + t*(x0 - x1) + x2)*(-1 + t + (x0 - x1)*Dtt(0)),
           2*(x0 - x1)*(-x0 + t*(x0 - x1) + x2)*Dtt(1) +
            2*(-y0 + t*(y0 - y1) + y2)*(-1 + t + (y0 - y1)*Dtt(1)),
           2*(y0 - y1)*(-y0 + t*(y0 - y1) + y2)*Dtt(2) +
            2*(x0 + t*(-x0 + x1) - x2)*(t + (-x0 + x1)*Dtt(2)),
           2*(x0 - x1)*(-x0 + t*(x0 - x1) + x2)*Dtt(3) +
            2*(y0 + t*(-y0 + y1) - y2)*(t + (-y0 + y1)*Dtt(3)),
           2*(y0 - y1)*(-y0 + t*(y0 - y1) + y2)*Dtt(4) +
            2*(-x0 + t*(x0 - x1) + x2)*(1 + (x0 - x1)*Dtt(4)),
           2*(x0 - x1)*(-x0 + t*(x0 - x1) + x2)*Dtt(5) +
            2*(-y0 + t*(y0 - y1) + y2)*(1 + (y0 - y1)*Dtt(5));

    DDdsq << 2*((y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(0,0) +
                (y0 - y1)*(y0 - y1)*(Dtt(0)*Dtt(0)) +
                (-x0 + t*(x0 - x1) + x2)*((x0 - x1)*DDtt(0,0) + 2*Dtt(0)) +
                (-1 + t + (x0 - x1)*Dtt(0))*(-1 + t + (x0 - x1)*Dtt(0))),
             2*((-y0 + t*(y0 - y1) + y2)*((y0 - y1)*DDtt(0,1) + Dtt(0)) +
                (x0 - x1)*(-1 + t + (x0 - x1)*Dtt(0))*Dtt(1) +
                (-x0 + t*(x0 - x1) + x2)*((x0 - x1)*DDtt(0,1) + Dtt(1)) +
                (y0 - y1)*Dtt(0)*(-1 + t + (y0 - y1)*Dtt(1))),
             2*((y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(0,2) + (y0 - y1)*(y0 - y1)*Dtt(0)*Dtt(2) +
                (-x0 + t*(x0 - x1) + x2)*((x0 - x1)*DDtt(0,2) - Dtt(0) + Dtt(2)) -
                (-1 + t + (x0 - x1)*Dtt(0))*(t + (-x0 + x1)*Dtt(2))),
             2*((y0 + t*(-y0 + y1) - y2)*((-y0 + y1)*DDtt(0,3) + Dtt(0)) +
                (x0 - x1)*(-1 + t + (x0 - x1)*Dtt(0))*Dtt(3) +
                (-x0 + t*(x0 - x1) + x2)*((x0 - x1)*DDtt(0,3) + Dtt(3)) +
                (-y0 + y1)*Dtt(0)*(t + (-y0 + y1)*Dtt(3))),
             2*((y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(0,4) + (y0 - y1)*(y0 - y1)*Dtt(0)*Dtt(4) +
                (-x0 + t*(x0 - x1) + x2)*((x0 - x1)*DDtt(0,4) + Dtt(4)) +
                (-1 + t + (x0 - x1)*Dtt(0))*(1 + (x0 - x1)*Dtt(4))),
             2*((y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(0,5) +
                (x0 - x1)*(-1 + t + (x0 - x1)*Dtt(0))*Dtt(5) +
                (-x0 + t*(x0 - x1) + x2)*((x0 - x1)*DDtt(0,5) + Dtt(5)) +
                (y0 - y1)*Dtt(0)*(1 + (y0 - y1)*Dtt(5))),
             2*((-y0 + t*(y0 - y1) + y2)*((y0 - y1)*DDtt(0,1) + Dtt(0)) +
                (x0 - x1)*(-1 + t + (x0 - x1)*Dtt(0))*Dtt(1) +
                (-x0 + t*(x0 - x1) + x2)*((x0 - x1)*DDtt(0,1) + Dtt(1)) +
                (y0 - y1)*Dtt(0)*(-1 + t + (y0 - y1)*Dtt(1))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(1,1) + (x0 - x1)*(x0 - x1)*(Dtt(1)*Dtt(1)) +
                (-y0 + t*(y0 - y1) + y2)*((y0 - y1)*DDtt(1,1) + 2*Dtt(1)) +
                (-1 + t + (y0 - y1)*Dtt(1))*(-1 + t + (y0 - y1)*Dtt(1))),
             2*((x0 + t*(-x0 + x1) - x2)*((-x0 + x1)*DDtt(1,2) + Dtt(1)) +
                (y0 - y1)*(-1 + t + (y0 - y1)*Dtt(1))*Dtt(2) +
                (-y0 + t*(y0 - y1) + y2)*((y0 - y1)*DDtt(1,2) + Dtt(2)) +
                (-x0 + x1)*Dtt(1)*(t + (-x0 + x1)*Dtt(2))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(1,3) + (x0 - x1)*(x0 - x1)*Dtt(1)*Dtt(3) +
                (-y0 + t*(y0 - y1) + y2)*((y0 - y1)*DDtt(1,3) - Dtt(1) + Dtt(3)) -
                (-1 + t + (y0 - y1)*Dtt(1))*(t + (-y0 + y1)*Dtt(3))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(1,4) +
                (y0 - y1)*(-1 + t + (y0 - y1)*Dtt(1))*Dtt(4) +
                (-y0 + t*(y0 - y1) + y2)*((y0 - y1)*DDtt(1,4) + Dtt(4)) +
                (x0 - x1)*Dtt(1)*(1 + (x0 - x1)*Dtt(4))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(1,5) + (x0 - x1)*(x0 - x1)*Dtt(1)*Dtt(5) +
                (-y0 + t*(y0 - y1) + y2)*((y0 - y1)*DDtt(1,5) + Dtt(5)) +
                (-1 + t + (y0 - y1)*Dtt(1))*(1 + (y0 - y1)*Dtt(5))),
             2*((y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(0,2) + (y0 - y1)*(y0 - y1)*Dtt(0)*Dtt(2) +
                (-x0 + t*(x0 - x1) + x2)*((x0 - x1)*DDtt(0,2) - Dtt(0) + Dtt(2)) -
                (-1 + t + (x0 - x1)*Dtt(0))*(t + (-x0 + x1)*Dtt(2))),
             2*((x0 + t*(-x0 + x1) - x2)*((-x0 + x1)*DDtt(1,2) + Dtt(1)) +
                (y0 - y1)*(-1 + t + (y0 - y1)*Dtt(1))*Dtt(2) +
                (-y0 + t*(y0 - y1) + y2)*((y0 - y1)*DDtt(1,2) + Dtt(2)) +
                (-x0 + x1)*Dtt(1)*(t + (-x0 + x1)*Dtt(2))),
             2*((y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(2,2) +
                (-x0 + t*(x0 - x1) + x2)*((x0 - x1)*DDtt(2,2) - 2*Dtt(2)) +
                (y0 - y1)*(y0 - y1)*(Dtt(2)*Dtt(2)) + (t + (-x0 + x1)*Dtt(2))*(t + (-x0 + x1)*Dtt(2)))
              ,2*((y0 + t*(-y0 + y1) - y2)*((-y0 + y1)*DDtt(2,3) + Dtt(2)) +
                (-x0 + x1)*(t + (-x0 + x1)*Dtt(2))*Dtt(3) +
                (x0 + t*(-x0 + x1) - x2)*((-x0 + x1)*DDtt(2,3) + Dtt(3)) +
                (-y0 + y1)*Dtt(2)*(t + (-y0 + y1)*Dtt(3))),
             2*((y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(2,4) + (y0 - y1)*(y0 - y1)*Dtt(2)*Dtt(4) +
                (x0 + t*(-x0 + x1) - x2)*((-x0 + x1)*DDtt(2,4) + Dtt(4)) +
                (t + (-x0 + x1)*Dtt(2))*(-1 + (-x0 + x1)*Dtt(4))),
             2*((y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(2,5) +
                (-x0 + x1)*(t + (-x0 + x1)*Dtt(2))*Dtt(5) +
                (x0 + t*(-x0 + x1) - x2)*((-x0 + x1)*DDtt(2,5) + Dtt(5)) +
                (y0 - y1)*Dtt(2)*(1 + (y0 - y1)*Dtt(5))),
             2*((y0 + t*(-y0 + y1) - y2)*((-y0 + y1)*DDtt(0,3) + Dtt(0)) +
                (x0 - x1)*(-1 + t + (x0 - x1)*Dtt(0))*Dtt(3) +
                (-x0 + t*(x0 - x1) + x2)*((x0 - x1)*DDtt(0,3) + Dtt(3)) +
                (-y0 + y1)*Dtt(0)*(t + (-y0 + y1)*Dtt(3))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(1,3) + (x0 - x1)*(x0 - x1)*Dtt(1)*Dtt(3) +
                (-y0 + t*(y0 - y1) + y2)*((y0 - y1)*DDtt(1,3) - Dtt(1) + Dtt(3)) -
                (-1 + t + (y0 - y1)*Dtt(1))*(t + (-y0 + y1)*Dtt(3))),
             2*((y0 + t*(-y0 + y1) - y2)*((-y0 + y1)*DDtt(2,3) + Dtt(2)) +
                (-x0 + x1)*(t + (-x0 + x1)*Dtt(2))*Dtt(3) +
                (x0 + t*(-x0 + x1) - x2)*((-x0 + x1)*DDtt(2,3) + Dtt(3)) +
                (-y0 + y1)*Dtt(2)*(t + (-y0 + y1)*Dtt(3))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(3,3) +
                (-y0 + t*(y0 - y1) + y2)*((y0 - y1)*DDtt(3,3) - 2*Dtt(3)) +
                (x0 - x1)*(x0 - x1)*(Dtt(3)*Dtt(3)) + (t + (-y0 + y1)*Dtt(3))*(t + (-y0 + y1)*Dtt(3)))
              ,2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(3,4) +
                (-y0 + y1)*(t + (-y0 + y1)*Dtt(3))*Dtt(4) +
                (y0 + t*(-y0 + y1) - y2)*((-y0 + y1)*DDtt(3,4) + Dtt(4)) +
                (x0 - x1)*Dtt(3)*(1 + (x0 - x1)*Dtt(4))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(3,5) + (x0 - x1)*(x0 - x1)*Dtt(3)*Dtt(5) +
                (y0 + t*(-y0 + y1) - y2)*((-y0 + y1)*DDtt(3,5) + Dtt(5)) +
                (t + (-y0 + y1)*Dtt(3))*(-1 + (-y0 + y1)*Dtt(5))),
             2*((y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(0,4) + (y0 - y1)*(y0 - y1)*Dtt(0)*Dtt(4) +
                (-x0 + t*(x0 - x1) + x2)*((x0 - x1)*DDtt(0,4) + Dtt(4)) +
                (-1 + t + (x0 - x1)*Dtt(0))*(1 + (x0 - x1)*Dtt(4))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(1,4) +
                (y0 - y1)*(-1 + t + (y0 - y1)*Dtt(1))*Dtt(4) +
                (-y0 + t*(y0 - y1) + y2)*((y0 - y1)*DDtt(1,4) + Dtt(4)) +
                (x0 - x1)*Dtt(1)*(1 + (x0 - x1)*Dtt(4))),
             2*((y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(2,4) + (y0 - y1)*(y0 - y1)*Dtt(2)*Dtt(4) +
                (x0 + t*(-x0 + x1) - x2)*((-x0 + x1)*DDtt(2,4) + Dtt(4)) +
                (t + (-x0 + x1)*Dtt(2))*(-1 + (-x0 + x1)*Dtt(4))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(3,4) +
                (-y0 + y1)*(t + (-y0 + y1)*Dtt(3))*Dtt(4) +
                (y0 + t*(-y0 + y1) - y2)*((-y0 + y1)*DDtt(3,4) + Dtt(4)) +
                (x0 - x1)*Dtt(3)*(1 + (x0 - x1)*Dtt(4))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(4,4) +
                (y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(4,4) + (y0 - y1)*(y0 - y1)*(Dtt(4)*Dtt(4)) +
                (1 + (x0 - x1)*Dtt(4))*(1 + (x0 - x1)*Dtt(4))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(4,5) +
                (y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(4,5) +
                (x0 - x1)*(1 + (x0 - x1)*Dtt(4))*Dtt(5) + (y0 - y1)*Dtt(4)*(1 + (y0 - y1)*Dtt(5))),
             2*((y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(0,5) +
                (x0 - x1)*(-1 + t + (x0 - x1)*Dtt(0))*Dtt(5) +
                (-x0 + t*(x0 - x1) + x2)*((x0 - x1)*DDtt(0,5) + Dtt(5)) +
                (y0 - y1)*Dtt(0)*(1 + (y0 - y1)*Dtt(5))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(1,5) + (x0 - x1)*(x0 - x1)*Dtt(1)*Dtt(5) +
                (-y0 + t*(y0 - y1) + y2)*((y0 - y1)*DDtt(1,5) + Dtt(5)) +
                (-1 + t + (y0 - y1)*Dtt(1))*(1 + (y0 - y1)*Dtt(5))),
             2*((y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(2,5) +
                (-x0 + x1)*(t + (-x0 + x1)*Dtt(2))*Dtt(5) +
                (x0 + t*(-x0 + x1) - x2)*((-x0 + x1)*DDtt(2,5) + Dtt(5)) +
                (y0 - y1)*Dtt(2)*(1 + (y0 - y1)*Dtt(5))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(3,5) + (x0 - x1)*(x0 - x1)*Dtt(3)*Dtt(5) +
                (y0 + t*(-y0 + y1) - y2)*((-y0 + y1)*DDtt(3,5) + Dtt(5)) +
                (t + (-y0 + y1)*Dtt(3))*(-1 + (-y0 + y1)*Dtt(5))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(4,5) +
                (y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(4,5) +
                (x0 - x1)*(1 + (x0 - x1)*Dtt(4))*Dtt(5) + (y0 - y1)*Dtt(4)*(1 + (y0 - y1)*Dtt(5))),
             2*((x0 - x1)*(-x0 + t*(x0 - x1) + x2)*DDtt(5,5) +
                (y0 - y1)*(-y0 + t*(y0 - y1) + y2)*DDtt(5,5) + (x0 - x1)*(x0 - x1)*(Dtt(5)*Dtt(5)) +
                (1 + (y0 - y1)*Dtt(5))*(1 + (y0 - y1)*Dtt(5)));

    d = sqrt(dsq);
    if(d!=0)
    {
        Dd = Ddsq/(2*d);
        DDd = DDdsq/(2*d)-Ddsq*Ddsq.transpose()/(4*d*dsq);
    }
    else
    {
        Dd = Eigen::Matrix<double,6,1>::Zero();
        DDd = Eigen::Matrix<double,6,6>::Zero();
    }

    p[3]=p[0]*(1-t)+p[1]*t;
}
