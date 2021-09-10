#include "element.h"
#include "node.h"

#include <cstdlib>
#include <algorithm>

//#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <tuple>

#include "spdlog/spdlog.h"


Eigen::Matrix2d icy::Element::DDs[6] = {
    (Eigen::Matrix2d() << 1, 0, 0, 0).finished(),
    (Eigen::Matrix2d() << 0, 0, 1, 0).finished(),
    (Eigen::Matrix2d() << 0, 1, 0, 0).finished(),
    (Eigen::Matrix2d() << 0, 0, 0, 1).finished(),
    (Eigen::Matrix2d() << -1, -1, 0, 0).finished(),
    (Eigen::Matrix2d() << 0, 0, -1, -1).finished(),
};

Eigen::Matrix<double,6,6> icy::Element::consistentMassMatrix =
        (Eigen::Matrix<double,6,6>() <<
         2,0,1,0,1,0,
         0,2,0,1,0,1,
         1,0,2,0,1,0,
         0,1,0,2,0,1,
         1,0,1,0,2,0,
         0,1,0,1,0,2
         ).finished()*(1.0/12.0);


void icy::Element::Reset(void)
{
    nds[0] = nds[1] = nds[2] = nullptr;
    incident_elems[0] = incident_elems[1] = incident_elems[2] = nullptr;
    area_initial = area_current = 0;
    group = -1;
    PiMultiplier = Eigen::Matrix2d::Identity();
    gmshTag = 0; // cannot stay zero
    quality_measure_Wicke = 1;
    leftCauchyGreenDeformationTensor = Eigen::Matrix2d::Identity();
}

void icy::Element::Initialize(Node *nd0, Node *nd1, Node *nd2)
{
    nds[0] = nd0;
    nds[1] = nd1;
    nds[2] = nd2;
}

void icy::Element::PrecomputeInitialArea()
{
    auto get_angle = [](Eigen::Vector2d u, Eigen::Vector2d v)
    { return (180.0/M_PI)*acos(std::clamp((double)u.normalized().dot(v.normalized()),-1.0,1.0)); };

    // reference shape matrix
    Dm << nds[0]->x_initial-nds[2]->x_initial, nds[1]->x_initial-nds[2]->x_initial;
    area_initial = area_current = Dm.determinant()/2;

    if(area_initial==0) throw std::runtime_error("element's initial area is zero");
    else if(area_initial < 0)
    {
        std::swap(nds[0],nds[1]);
        Dm << nds[0]->x_initial-nds[2]->x_initial, nds[1]->x_initial-nds[2]->x_initial;
        area_initial = area_current = Dm.determinant()/2;
    }
    DmInv = Dm.inverse();


    double angles[3];
    for(int i=0;i<3;i++) angles[i] = get_angle(nds[i]->x_initial-nds[(i+1)%3]->x_initial,
            nds[i]->x_initial-nds[(i+2)%3]->x_initial);
    constexpr double angle_thredhold = 3;
    if(angles[0]<angle_thredhold || angles[1]<angle_thredhold || angles[2]<angle_thredhold || area_initial < threshold_area)
    {
        spdlog::critical("Area {0} of elem {1}-{2}-{3}\nangles: {4:.1f}; {5:.1f}; {6:.1f}",
                         area_initial, nds[0]->locId,nds[1]->locId,nds[2]->locId,angles[0],angles[1],angles[2]);
        //throw std::runtime_error("PrecomputeInitialArea(): degenerate element created");
    }

}


void icy::Element::AddToSparsityStructure(EquationOfMotionSolver &eq) const
{
    eq.AddEntriesToStructure(nds[0]->eqId, nds[1]->eqId, nds[2]->eqId);
}


bool icy::Element::ComputeEquationEntries(EquationOfMotionSolver &eq, const SimParams &prms, double h)
{
    // NeoHookeanElasticity
    double lambda = prms.lambda;
    double mu = prms.mu;

    Eigen::Matrix2d Ds, Finv, FT, FinvT;

    double W = prms.Thickness*area_initial;   // element's initial "volume"

    // deformed shape matrix
    Ds << nds[0]->xt-nds[2]->xt, nds[1]->xt-nds[2]->xt;
    if(Ds.determinant()<=0) return false; // mesh is inverted

    F = Ds*DmInv*PiMultiplier;    // deformation gradient (multiplied by a coefficient)
    FT = F.transpose();
    Finv = F.inverse();
    FinvT = Finv.transpose();
    double J = F.determinant();     // represents the change of volume in comparison with the reference
    volume_change = J;

    double log_J = log(J);
    strain_energy_density = (mu/2.0)*((F*FT).trace()-2.0)-mu*log_J+(lambda/2.0)*log_J*log_J;

    // First Piola-Kirchhoff stress tensor
    P = F*mu + FinvT*(lambda*log_J-mu);

    // forces on nodes 1 and 2 (inverted sign)
    Eigen::Matrix2d H = W*P*DmInv.transpose();

    // energy gradient with respect to x1,y1,x2,y2,x3,y3
    Eigen::Matrix<double, 6, 1> DE;    // energy gradient
    DE << H(0,0), H(1,0), H(0,1), H(1,1), -(H(0,0)+H(0,1)),-(H(1,0)+H(1,1));

    Eigen::Matrix<double, 6, 6> HE; // energy Hessian, 6x6 symmetric
    for(int i=0;i<6;i++)
    {
        Eigen::Matrix2d DF_i = DDs[i]*DmInv; // derivative of F with respect to x_i
        Eigen::Matrix2d dP = mu*DF_i + (mu-lambda*log_J)*FinvT*DF_i.transpose()*FinvT + lambda*(Finv*DF_i).trace()*FinvT;
        Eigen::Matrix2d dH = W*dP*DmInv.transpose();
        HE(0,i) = dH(0,0);
        HE(1,i) = dH(1,0);
        HE(2,i) = dH(0,1);
        HE(3,i) = dH(1,1);
        HE(4,i) = -dH(0,0)-dH(0,1);
        HE(5,i) = -dH(1,0)-dH(1,1);
    }

    double hsq = h*h;   // squared time step
    DE*=hsq;
    HE*=hsq;
    double constTerm = strain_energy_density*W*hsq;

    // Apply body forces via consistent mass matrix
    double massMatrixMultiplier = area_initial * prms.Density * prms.Thickness;
    Eigen::Matrix<double,6,6> massMatrix = consistentMassMatrix*massMatrixMultiplier;
    HE += massMatrix;

    Eigen::Matrix<double,6,1> xt, x_hat, lambda_n, linear_term_mass;
    xt << nds[0]->xt[0],nds[0]->xt[1],nds[1]->xt[0],nds[1]->xt[1],nds[2]->xt[0],nds[2]->xt[1];
    x_hat << nds[0]->x_hat[0],nds[0]->x_hat[1],nds[1]->x_hat[0],nds[1]->x_hat[1],nds[2]->x_hat[0],nds[2]->x_hat[1];
    lambda_n = xt-x_hat;
    linear_term_mass = massMatrix*lambda_n;
    DE+=linear_term_mass;

    double const_term_mass = lambda_n.dot(massMatrix*lambda_n)/2;
    constTerm += const_term_mass;

    // distribute to the equation of motion
    eq.AddToEquation(constTerm, DE, HE, {nds[0]->eqId,nds[1]->eqId,nds[2]->eqId});

    return true;
}

void icy::Element::ComputeVisualizedVariables()
{
    // compute various variables, mostly for visualization
    CauchyStress = F*P.transpose()/volume_change;   // symmetric up to the roundoff error

    double sx = CauchyStress(0,0);
    double sy = CauchyStress(1,1);
    double tauxy = (CauchyStress(0,1)+CauchyStress(1,0))/2;
    CauchyStress(0,1) = CauchyStress(1,0) = tauxy;

    max_shear_stress = sqrt((sx-sy)*(sx-sy)/4 + tauxy*tauxy);
    principal_stress1 = (sx+sy)/2 + max_shear_stress;
    principal_stress2 = (sx+sy)/2 - max_shear_stress;
    hydrostatic_stress = CauchyStress.trace()/2;
    GreenStrain = (F.transpose()*F - Eigen::Matrix2d::Identity())/2;

    // velocity divergence
    Eigen::Matrix2d D, DinvT, DDot;

    // deformed shape matrix
    D << nds[1]->xn-nds[0]->xn, nds[2]->xn-nds[0]->xn;
    DinvT = D.inverse().transpose();

    // velocity matrix
    DDot << nds[1]->vn-nds[0]->vn, nds[2]->vn-nds[0]->vn;
    velocity_divergence = DinvT.cwiseProduct(DDot).sum();

    //  quality measures
    double V = D.determinant()/2;
    double e0sq = (nds[1]->xn-nds[2]->xn).squaredNorm();
    double e1sq = (nds[0]->xn-nds[2]->xn).squaredNorm();
    double e2sq = (nds[1]->xn-nds[0]->xn).squaredNorm();
    double e0 = sqrt(e0sq);
    double e1 = sqrt(e1sq);
    double e2 = sqrt(e2sq);
    double l_harm = 3.0/(1.0/e0 + 1.0/e1 + 1.0/e2);
    double l_rms_sq = (1.0/3.0)*(e0sq+e1sq+e2sq);

    constexpr double coeff = 2.30940107*1.2;
    quality_measure_Wicke = coeff*V*l_harm/(l_rms_sq*sqrt(l_rms_sq));
    quality_measure_Wicke = std::min(1.0,quality_measure_Wicke);
    leftCauchyGreenDeformationTensor = F*F.transpose();
}

bool icy::Element::PlasticDeformation(const SimParams &prms, double timeStep)
{
    constexpr double epsilon = 1e-5;
    double stressNorm = CauchyStress.norm();
    double tau = prms.PlasticYieldThreshold;

    plasticity_tau_ratio = (stressNorm-tau)/stressNorm;
    plasticity_gamma = timeStep*prms.PlasticFlowRate*((stressNorm-tau)/stressNorm);
    plasticity_gamma = std::clamp(plasticity_gamma, 0.0, 1.0);
    if(plasticity_gamma < epsilon) return false;

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
    double val1 = std::pow(SVDvals[0]/detS_sqRoot,-plasticity_gamma);
    double val2 = std::pow(SVDvals[1]/detS_sqRoot,-plasticity_gamma);
    double sqRoot = sqrt(val1*val2);
    S_hat << val1/sqRoot, 0, 0, val2/sqRoot;
    V = svd.matrixV();
    PiMultiplier *= V*S_hat*V.transpose();
    return true;
}


// FRACTURE ALGORITHM
void icy::Element::getIdxs(const icy::Node*nd, short &thisIdx, short &CWIdx, short &CCWIdx) const
{
    if(nd==nds[0]) thisIdx=0;
    else if(nd==nds[1]) thisIdx=1;
    else if(nd==nds[2]) thisIdx=2;
    else throw std::runtime_error("getIdxs: node does not belong to the element");
    CWIdx = (thisIdx+1)%3;
    CCWIdx = (thisIdx+2)%3;
}

std::pair<icy::Node*,icy::Node*> icy::Element::CW_CCW_Node(const Node* nd) const
{
    short idx = getNodeIdx(nd);
    return {nds[(idx+1)%3],nds[(idx+2)%3]};
}

bool icy::Element::isOnBoundary(const Node* nd) const
{
    short idx = getNodeIdx(nd);
    short cw_idx = (idx+1)%3;
    short ccw_idx = (idx+2)%3;
    return isBoundary(cw_idx) || isBoundary(ccw_idx);
}

bool icy::Element::isCWBoundary(const Node* nd) const
{
    short idx = getNodeIdx(nd);
    short cw_idx = (idx+2)%3;
    return isBoundary(cw_idx);
}

bool icy::Element::isCCWBoundary(const Node* nd) const
{
    short idx = getNodeIdx(nd);
    short ccw_idx = (idx+1)%3;
    return isBoundary(ccw_idx);
}



icy::Element* icy::Element::getAdjacentElementOppositeToNode(Node *nd)
{
    if(nds[0]==nd) return incident_elems[0];
    else if(nds[1]==nd)  return incident_elems[1];
    else if(nds[2]==nd)  return incident_elems[2];
    else throw std::runtime_error("getAdjacentElementOppositeToNode");
}

short icy::Element::getNodeIdx(const Node *nd) const
{
    if(nds[0]==nd) return 0;
    else if(nds[1]==nd) return 1;
    else if(nds[2]==nd) return 2;
    else
    {
        spdlog::critical("getNodeIdx: trying to obtain index of the node {} in element {}-{}-{}",
                         nd->locId,nds[0]->locId,nds[1]->locId,nds[2]->locId);
        throw std::runtime_error("getNodeIdx");
    }
}

icy::Node* icy::Element::getOppositeNode(Node *nd0, Node* nd1)
{
    for(int i=0;i<3;i++)
    {
        int idx_next = (i+1)%3;
        if((nds[i] == nd0 && nds[idx_next] == nd1)||
                (nds[i] == nd1 && nds[idx_next] == nd0))
            return nds[(i+2)%3];
    }

    spdlog::critical("getOppositeNode: trying to find edge {}-{} in element {}-{}-{}",
                     nd0->locId,nd1->locId,nds[0]->locId,nds[1]->locId,nds[2]->locId);
    throw std::runtime_error("getOppositeNode: opposite node not found");
}

void icy::Element::ReplaceNode(Node *replaceWhat, Node *replaceWith)
{
    //spdlog::info("Replacing Node {} with {} in elem {}-{}-{}",
    //             replaceWhat->locId,replaceWith->locId, nds[0]->locId,nds[1]->locId,nds[2]->locId);
    if(nds[0] == replaceWhat) nds[0] = replaceWith;
    else if(nds[1] == replaceWhat) nds[1] = replaceWith;
    else if(nds[2] == replaceWhat) nds[2] = replaceWith;
    else throw std::runtime_error("icy::Element::ReplaceNode: replaced node not found");
    PrecomputeInitialArea();
    //spdlog::info("Replaced Node {} with {}; elem {}-{}-{}",
    //             replaceWhat->locId,replaceWith->locId, nds[0]->locId,nds[1]->locId,nds[2]->locId);
}

void icy::Element::DisconnectFromElem(Element* other)
{
    if(incident_elems[0]==other) incident_elems[0]=nullptr;
    else if(incident_elems[1]==other) incident_elems[1]=nullptr;
    else if(incident_elems[2]==other) incident_elems[2]=nullptr;
    else throw std::runtime_error("DisconnectFromElem: incident elem not found");
}

void icy::Element::ReplaceIncidentElem(Element* which, Element* withWhat)
{
    if(incident_elems[0]==which) incident_elems[0]=withWhat;
    else if(incident_elems[1]==which) incident_elems[1]=withWhat;
    else if(incident_elems[2]==which) incident_elems[2]=withWhat;
    else throw std::runtime_error("ReplaceIncidentElem: incident elem not found");
}

