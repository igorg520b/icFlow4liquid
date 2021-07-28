#include "element.h"
#include "node.h"

#include <cstdlib>
#include <algorithm>

#include <iostream>
#include <vector>
#include <utility>
#include <cmath>


Eigen::Matrix2d icy::Element::DDs[6] = {
    (Eigen::Matrix2d() << 1, 0, 0, 0).finished(),
    (Eigen::Matrix2d() << 0, 0, 1, 0).finished(),
    (Eigen::Matrix2d() << 0, 1, 0, 0).finished(),
    (Eigen::Matrix2d() << 0, 0, 0, 1).finished(),
    (Eigen::Matrix2d() << -1, -1, 0, 0).finished(),
    (Eigen::Matrix2d() << 0, 0, -1, -1).finished(),
};


void icy::Element::Reset(void)
{
    nds[0] = nds[1] = nds[2] = nullptr;
    area_initial = area_current = 0;
    group = -1;
    fluid = false;
    PiMultiplier = Eigen::Matrix2d::Identity();
}

void icy::Element::PrecomputeInitialArea()
{
    Eigen::Matrix2d J;
    J << nds[0]->x_initial-nds[2]->x_initial, nds[1]->x_initial-nds[2]->x_initial;
    area_initial = area_current = J.determinant()/2;
    if(area_initial==0) throw std::runtime_error("element's initial area is zero");
}


void icy::Element::AddToSparsityStructure(EquationOfMotionSolver &eq)
{
    // register the positions of non-zero entries with the EquationOfMotionSolver
    eq.AddElementToStructure(nds[0]->eqId, nds[1]->eqId);
    eq.AddElementToStructure(nds[0]->eqId, nds[2]->eqId);
    eq.AddElementToStructure(nds[1]->eqId, nds[2]->eqId);
}


bool icy::Element::ComputeEquationEntries(EquationOfMotionSolver &eq, SimParams &prms, double h)
{
    // NeoHookeanElasticity
    double lambda = prms.lambda;
    double mu = prms.mu;

    Eigen::Matrix2d Dm, Dm_inv, Ds, F, Finv, FT, FinvT;

    // reference shape matrix
    Dm << nds[0]->x_initial-nds[2]->x_initial, nds[1]->x_initial-nds[2]->x_initial;
    double W = prms.Thickness*Dm.determinant()/2;   // element's initial "volume"
    Dm_inv = Dm.inverse();

    // deformed shape matrix
    Ds << nds[0]->xt-nds[2]->xt, nds[1]->xt-nds[2]->xt;
    if(Ds.determinant()<=0) return false; // mesh is inverted

    F = Ds*Dm_inv*PiMultiplier;    // deformation gradient (multiplied by a coefficient)

    double J = F.determinant();     // represents the change of volume in comparison with the reference
    volume_change = J;
    FT = F.transpose();
    Finv = F.inverse();
    FinvT = Finv.transpose();

    double log_J = log(J);
    strain_energy_density = (mu/2.0)*((F*FT).trace()-2.0)-mu*log_J+(lambda/2.0)*log_J*log_J;

    // First Piola-Kirchhoff stress tensor
    Eigen::Matrix2d P = F*mu + FinvT*(lambda*log_J-mu);

    // forces on nodes 1 and 2 (inverted sign)
    Eigen::Matrix2d H = W*P*Dm_inv.transpose();

    // energy gradient with respect to x1,y1,x2,y2,x3,y3
    Eigen::Matrix<double, 6, 1> DE;    // energy gradient
    DE << H(0,0), H(1,0), H(0,1), H(1,1), -(H(0,0)+H(0,1)),-(H(1,0)+H(1,1));

    Eigen::Matrix<double, 6, 6> HE; // energy Hessian, 6x6 symmetric
    for(int i=0;i<6;i++)
    {
        Eigen::Matrix2d DF_i = DDs[i]*Dm_inv; // derivative of F with respect to x_i
        Eigen::Matrix2d dP = mu*DF_i + (mu-lambda*log_J)*FinvT*DF_i.transpose()*FinvT + lambda*(Finv*DF_i).trace()*FinvT;
        Eigen::Matrix2d dH = W*dP*Dm_inv.transpose();
        HE(0,i) = dH(0,0);
        HE(1,i) = dH(1,0);
        HE(2,i) = dH(0,1);
        HE(3,i) = dH(1,1);
        HE(4,i) = -dH(0,0)-dH(0,1);
        HE(5,i) = -dH(1,0)-dH(1,1);
    }


    // assemble the equation of motion
    double hsq = h*h;
    for(int i=0;i<3;i++)
    {
        int row = nds[i]->eqId;
        Eigen::Vector2d locDE = DE.block(i*2,0,2,1)*hsq;
        eq.AddToC(row, locDE);
        for(int j=0;j<3;j++)
        {
            int col = nds[j]->eqId;
            Eigen::Matrix2d locHE = HE.block(i*2,j*2,2,2)*hsq;
            eq.AddToQ(row, col, locHE);
        }
    }
    eq.AddToConstTerm(strain_energy_density*W*hsq);


    // compute various variables, mostly for visualization
    CauchyStress = F*P.transpose()/J;   // symmetric up to roundoff error

    double sx = CauchyStress(0,0);
    double sy = CauchyStress(1,1);
    double tauxy = (CauchyStress(0,1)+CauchyStress(1,0))/2;
    CauchyStress(0,1) = CauchyStress(1,0) = tauxy;

    max_shear_stress = sqrt((sx-sy)*(sx-sy)/4 + tauxy*tauxy);
    principal_stress1 = (sx+sy)/2 + max_shear_stress;
    principal_stress2 = (sx+sy)/2 - max_shear_stress;
    hydrostatic_stress = CauchyStress.trace()/2;
    GreenStrain = (F.transpose()*F - Eigen::Matrix2d::Identity())*0.5;

    return true;


}




void icy::Element::EvaluateVelocityDivergence()
{
    Eigen::Matrix2d D, DinvT, DDot;

    // deformed shape matrix
    D << nds[1]->xn-nds[0]->xn, nds[2]->xn-nds[0]->xn;
    DinvT = D.inverse().transpose();

    // velocity matrix
    DDot << nds[1]->vn-nds[0]->vn, nds[2]->vn-nds[0]->vn;
    velocity_divergence = DinvT.cwiseProduct(DDot).sum();
}

void icy::Element::PlasticDeformation(SimParams &prms, double timeStep)
{
    double stressNorm = CauchyStress.norm();
    double tau = prms.PlasticYieldThreshold;
    double gamma = timeStep*prms.PlasticFlowRate*((stressNorm-tau)/stressNorm);
    gamma = std::clamp(gamma, 0.0, 1.0);
    if(gamma < 1e-10) return;

    Eigen::Matrix2d Dm, Dm_inv, Ds, F, V, S_hat;

    // reference shape matrix
    Dm << nds[0]->x_initial-nds[2]->x_initial, nds[1]->x_initial-nds[2]->x_initial;
    Dm_inv = Dm.inverse();

    // deformed shape matrix
    Ds << nds[0]->xn-nds[2]->xn, nds[1]->xn-nds[2]->xn;
    area_current = Ds.determinant()/2;
    if(area_current<=0) throw std::runtime_error("PlasticDeformation: inverted mesh");

    F = Ds*Dm_inv*PiMultiplier;    // deformation gradient (multiplied by a coefficient)
    Eigen::JacobiSVD<Eigen::Matrix2d> svd(F, Eigen::ComputeFullV);

    Eigen::Vector2d SVDvals = svd.singularValues();
    double detS_sqRoot = 1;//sqrt(SVDvals[0]*SVDvals[1]);
    double val1 = std::pow(SVDvals[0]/detS_sqRoot,-gamma);
    double val2 = std::pow(SVDvals[1]/detS_sqRoot,-gamma);
    double sqRoot = sqrt(val1*val2);
    S_hat << val1/sqRoot, 0, 0, val2/sqRoot;
    V = svd.matrixV();
    PiMultiplier *= V*S_hat*V.transpose();
}


