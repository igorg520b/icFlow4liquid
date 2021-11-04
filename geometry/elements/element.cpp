#include "element.h"
#include "boundaryedge.h"
#include "node.h"
#include "mesh.h"

#include <cstdlib>
#include <algorithm>

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
    PiMultiplier = Eigen::Matrix2d::Identity();
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
    computeDm();
    area_initial = area_current = Dm.determinant()/2;

    if(area_initial==0) throw std::runtime_error("element's initial area is zero");
    else if(area_initial < 0)
    {
        std::swap(nds[0],nds[1]);
        computeDm();
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
    int idxs[] {nds[0]->eqId, nds[1]->eqId, nds[2]->eqId};
    eq.AddEntriesToStructure(std::begin(idxs),std::end(idxs));
}


bool icy::Element::ComputeEquationEntries(EquationOfMotionSolver &eq, const SimParams &prms, double h)
{
    // NeoHookeanElasticity
    double lambda = prms.lambda;
    double mu = prms.mu;
    double W = prms.Thickness*area_initial;   // element's initial "volume"

    // deformed shape matrix
    Eigen::Matrix2d Ds;
    Ds << (nds[0]->xt - nds[2]->xt), (nds[1]->xt - nds[2]->xt);
    if(Ds.determinant()<=0) return false; // mesh is inverted

    F = Ds*DmInv*PiMultiplier;    // deformation gradient (multiplied by a coefficient)
    Eigen::Matrix2d FT = F.transpose();
    Eigen::Matrix2d Finv = F.inverse();
    Eigen::Matrix2d FinvT = Finv.transpose();
    double J = volume_change = F.determinant();     // represents the change of volume in comparison with the reference

    double log_J = log(J);
    strain_energy_density = (mu/2.0)*((F*FT).trace()-2.0) - mu*log_J + (lambda/2.0)*log_J*log_J;

    // First Piola-Kirchhoff stress tensor
    P = F*mu + FinvT*(lambda*log_J - mu);

    // forces on nodes 1 and 2 (inverted sign)
    Eigen::Matrix2d H = W * P * DmInv.transpose();

    // energy gradient with respect to x1,y1,x2,y2,x3,y3
    Eigen::Matrix<double,6,1> DE; // energy gradient
    DE << H(0,0), H(1,0), H(0,1), H(1,1), -(H(0,0) + H(0,1)), -(H(1,0) + H(1,1));

    Eigen::Matrix < double, 6, 6 > HE; // energy Hessian, 6x6 symmetric
    for (int i=0;i<6;i++)
    {
        Eigen::Matrix2d DF_i = DDs[i] * DmInv; // derivative of F with respect to x_i
        Eigen::Matrix2d dP = mu*DF_i + (mu - lambda*log_J)*FinvT*DF_i.transpose()*FinvT + lambda*(Finv*DF_i).trace()*FinvT;
        Eigen::Matrix2d dH = W*dP*DmInv.transpose();
        HE(0,i) = dH(0,0);
        HE(1,i) = dH(1,0);
        HE(2,i) = dH(0,1);
        HE(3,i) = dH(1,1);
        HE(4,i) = -dH(0,0) - dH(0,1);
        HE(5,i) = -dH(1,0) - dH(1,1);
    }

    double hsq = h*h; // squared time step
    DE *= hsq;
    HE *= hsq;

    // Apply body forces via consistent mass matrix
    Eigen::Matrix<double,6,6> massMatrix = consistentMassMatrix*(area_initial*prms.Thickness*prms.Density);
    HE += massMatrix;

    Eigen::Matrix<double,6,1> xt_xhat;
    xt_xhat << (nds[0]->xt - nds[0]->x_hat), (nds[1]->xt - nds[1]->x_hat), (nds[2]->xt - nds[2]->x_hat);
    DE += massMatrix*(xt_xhat);

    // distribute to the equation of motion
    eq.AddToEquation(DE.data(), HE.data(), {nds[0]->eqId,nds[1]->eqId,nds[2]->eqId});

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

    // calculate velocity vector divergence
    // velocity matrix (same as Ds, but with velocities)
    Eigen::Matrix2d Ds, DDot;
    Ds = getDs_at_n();
    DDot << (nds[0]->vn - nds[2]->vn), (nds[1]->vn - nds[2]->vn);
    velocity_divergence = Ds.inverse().transpose().cwiseProduct(DDot).sum();

    //  quality measures
    double V = Ds.determinant()/2;
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




// GEOMETRICAL HELPER FUNCTIONS

uint8_t icy::Element::getNodeIdx(const Node *nd) const
{
    if(nds[0]==nd) return 0;
    else if(nds[1]==nd) return 1;
    else if(nds[2]==nd) return 2;
    else
    {
        spdlog::info("getNodeIdx: node not found; *nd = {}", (void*)nd);
        spdlog::info("getNodeIdx: nd->locId {}", nd->locId);
        spdlog::info("getNodeIdx: elem nds {} - {} - {}", (void*)nds[0],(void*)nds[1],(void*)nds[2]);
        spdlog::info("getNodeIdx: elem nds locid {} - {} - {}", nds[0]->locId,nds[1]->locId,nds[2]->locId);
        throw std::runtime_error("getNodeIdx");
    }
}

void icy::Element::getIdxs(const icy::Node*nd, uint8_t &thisIdx, uint8_t &CWIdx, uint8_t &CCWIdx) const
{
    thisIdx = getNodeIdx(nd);
    CWIdx = (thisIdx+1)%3;
    CCWIdx = (thisIdx+2)%3;
}

std::pair<icy::Node*,icy::Node*> icy::Element::CW_CCW_Node(const Node* nd) const
{
    uint8_t idx = getNodeIdx(nd);
    return {nds[(idx+1)%3],nds[(idx+2)%3]};
}

icy::Node* icy::Element::CW_Node(const Node* nd) const
{
    uint8_t idx = getNodeIdx(nd);
    return nds[(idx+1)%3];
}

icy::Node* icy::Element::CCW_Node(const Node* nd) const
{
    uint8_t idx = getNodeIdx(nd);
    return nds[(idx+2)%3];
}

bool icy::Element::isEdgeCW(const Node *nd1, const Node *nd2) const
{
    for(int i=0;i<3;i++)
    {
        if(nds[i%3]==nd1 && nds[(i+2)%3]==nd2) return true;
        if(nds[i%3]==nd1 && nds[(i+1)%3]==nd2) return false;
    }
    throw std::runtime_error("icy::Element::isEdgeCW: edge not found");
}

bool icy::Element::containsEdge(const Node *nd1, const Node *nd2) const
{
    for(uint8_t i=0;i<3;i++)
        if((nds[(i+1)%3]==nd1 && nds[(i+2)%3]==nd2) || (nds[(i+2)%3]==nd1 && nds[(i+1)%3]==nd2))
            return true;
    return false;
}


uint8_t icy::Element::getEdgeIdx(const Node *nd1, const Node *nd2) const
{
    for(uint8_t i=0;i<3;i++)
        if((nds[(i+1)%3]==nd1 && nds[(i+2)%3]==nd2) || (nds[(i+2)%3]==nd1 && nds[(i+1)%3]==nd2))
            return i;
    throw std::runtime_error("icy::Element::getEdgeIdx: edge not found");
}


bool icy::Element::isOnBoundary(const Node* nd) const
{
    uint8_t idx = getNodeIdx(nd);
    uint8_t cw_idx = (idx+1)%3;
    uint8_t ccw_idx = (idx+2)%3;
    return isBoundaryEdge(cw_idx) || isBoundaryEdge(ccw_idx);
}

bool icy::Element::isCWBoundary(const Node* nd) const
{
    uint8_t idx = getNodeIdx(nd);
    uint8_t cw_idx = (idx+2)%3;
    bool result = isBoundaryEdge(cw_idx);
    return result;
}

bool icy::Element::isCCWBoundary(const Node* nd) const
{
    uint8_t idx = getNodeIdx(nd);
    uint8_t ccw_idx = (idx+1)%3;
    return isBoundaryEdge(ccw_idx);
}

icy::BaseElement* icy::Element::getIncidentElementOppositeToNode(Node *nd)
{
    return incident_elems[getNodeIdx(nd)];
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




// FRACTURE ALGORITHM

void icy::Element::ReplaceNode(Node *replaceWhat, Node *replaceWith)
{
    uint8_t nd_idx = getNodeIdx(replaceWhat);
    nds[nd_idx] = replaceWith;

    PrecomputeInitialArea();

    uint8_t cw_idx = (nd_idx+1)%3;
    uint8_t ccw_idx = (nd_idx+2)%3;

    // update incident elements after replacing the node
    for(uint8_t idx : {cw_idx,ccw_idx})
    {
        if(incident_elems[idx]->type == ElementType::BEdge)
            dynamic_cast<BoundaryEdge*>(incident_elems[idx])->UpdateNodes();
        else if(incident_elems[idx]->type == ElementType::CZ)
            dynamic_cast<CohesiveZone*>(incident_elems[idx])->UpdateNodes();
    }

}


void icy::Element::ReplaceIncidentElem(BaseElement* which, BaseElement* withWhat)
{
    if(incident_elems[0] == which) incident_elems[0] = withWhat;
    else if(incident_elems[1] == which) incident_elems[1] = withWhat;
    else if(incident_elems[2] == which) incident_elems[2] = withWhat;
    else throw std::runtime_error("ReplaceIncidentElem: incident elem not found");
}

void icy::Element::RecalculatePiMultiplierFromDeformationGradient(Eigen::Matrix2d F_tilda)
{
    PiMultiplier = Dm * getDs_at_n().inverse() * F_tilda;
}



std::pair<icy::Node*,icy::Node*> icy::Element::SplitElem(Node *nd, Node *nd0, Node *nd1, double where)
{
    BaseElement* incident_elem = getIncidentElementOppositeToNode(nd);

    if(incident_elem->type == ElementType::BEdge)
        return {SplitBoundaryElem(nd, nd0, nd1, where),nullptr};
    else if(incident_elem->type == ElementType::TElem)
        return {SplitNonBoundaryElem(nd, nd0, nd1, where),nullptr};
    else if(incident_elem->type == ElementType::CZ)
        return SplitElemWithCZ(nd, nd0, nd1, where);
    else throw std::runtime_error("SplitElem: unknown incident element type");
}


icy::Node* icy::Element::SplitBoundaryElem(Node *nd, Node *nd0, Node *nd1, double where)
{
    MeshFragment *fragment = nd->fragment;

    uint8_t ndIdx = getNodeIdx(nd);
    uint8_t nd0Idx = getNodeIdx(nd0);
    uint8_t nd1Idx = getNodeIdx(nd1);

    Eigen::Matrix2d F_orig = getF_at_n();     // save the deformation gradient


    // insert element
    Element *insertedElem = fragment->AddElement();
    nd->adj_elems.push_back(insertedElem);

    // insert the node between nd0 and nd1; initialize its coordinates; connect to adjacent elements
    Node *split = fragment->AddNode();
    split->InitializeLERP(nd0, nd1, where);
    split->isBoundary = true;
    split->adj_elems.push_back(this);
    split->adj_elems.push_back(insertedElem);

    // modify the original element
    nds[nd1Idx] = split;

    // initialize the inserted element's nodes
    insertedElem->nds[ndIdx] = nd;
    insertedElem->nds[nd1Idx] = nd1;
    insertedElem->nds[nd0Idx] = split;

    // if the original element had a boundary at nd0Idx, the inserted element now takes that boundary
    if(incident_elems[nd0Idx]->type == ElementType::BEdge)
        static_cast<BoundaryEdge*>(incident_elems[nd0Idx])->Initialize(insertedElem,nd0Idx,3);
    else if(incident_elems[nd0Idx]->type == ElementType::TElem)
        static_cast<Element*>(incident_elems[nd0Idx])->ReplaceIncidentElem(this,insertedElem);
    else if(incident_elems[nd0Idx]->type == ElementType::CZ)
        static_cast<CohesiveZone*>(incident_elems[nd0Idx])->ReplaceAdjacentElem(this, insertedElem, nd0Idx);

    // add the boundary that has just appeared
    // automatically initialize the inserted element's adjacency data (insertedElem->incident_elems[ndIdx])
    fragment->AddBoundary(insertedElem,ndIdx,3);
    dynamic_cast<BoundaryEdge*>(incident_elems[ndIdx])->UpdateNodes();

    insertedElem->incident_elems[nd0Idx] = incident_elems[nd0Idx];
    incident_elems[nd0Idx] = insertedElem;
    insertedElem->incident_elems[nd1Idx] = this;

    // from node "nd1", disconnect the original element and replace it with the inserted element
    nd1->ReplaceAdjacentElement(this, insertedElem);

    // compute the new area and reference shape matrix
    this->PrecomputeInitialArea();
    insertedElem->PrecomputeInitialArea();

    // re-evaluate PiMultiplier on both elements to maintain consistent plasticity
    this->RecalculatePiMultiplierFromDeformationGradient(F_orig);
    insertedElem->RecalculatePiMultiplierFromDeformationGradient(F_orig);

    // "fix" the fan for the node, whose element was just replaced
    nd1->PrepareFan();
    return split;
}

icy::Node* icy::Element::SplitNonBoundaryElem(Node *nd, Node *nd0, Node *nd1, double where)
{
    MeshFragment *fragment = nd->fragment;

    icy::Element *adjElem = dynamic_cast<icy::Element*>(getIncidentElementOppositeToNode(nd));
    if(adjElem == nullptr) throw std::runtime_error("icy::Element::SplitNonBoundaryElem dynamic cast issue");

    uint8_t ndIdx_orig = getNodeIdx(nd);
    uint8_t nd0Idx_orig = getNodeIdx(nd0);
    uint8_t nd1Idx_orig = getNodeIdx(nd1);

    // preserve deformation gradient
    Eigen::Matrix2d F_orig = getF_at_n();
    Eigen::Matrix2d F_adj = adjElem->getF_at_n();

    Node *oppositeNode = adjElem->getOppositeNode(nd0, nd1);
    uint8_t nd0Idx_adj = adjElem->getNodeIdx(nd0);
    uint8_t nd1Idx_adj = adjElem->getNodeIdx(nd1);
    uint8_t oppIdx_adj = adjElem->getNodeIdx(oppositeNode);

    // insert "main" element
    Element *insertedElem = fragment->AddElement();
    nd->adj_elems.push_back(insertedElem);

    // insert "adjacent" element
    Element *insertedElem_adj = fragment->AddElement();

    // insert the "split" node between nd0 and nd1
    Node *split = fragment->AddNode();
    split->InitializeLERP(nd0, nd1, where);
    split->adj_elems.insert(split->adj_elems.end(),{this,insertedElem,adjElem,insertedElem_adj});
    split->isBoundary = false;

    // modify the original element
    nds[nd1Idx_orig] = split;

    nd1->ReplaceAdjacentElement(this,insertedElem);

    // initialize the inserted "main" element
    insertedElem->nds[ndIdx_orig] = nd;
    insertedElem->nds[nd1Idx_orig] = nd1;
    insertedElem->nds[nd0Idx_orig] = split;

    // if the original element had a boundary at nd0Idx, the inserted element now takes that boundary
    if(incident_elems[nd0Idx_orig]->type == ElementType::BEdge)
        static_cast<BoundaryEdge*>(incident_elems[nd0Idx_orig])->Initialize(insertedElem,nd0Idx_orig);
    else if(incident_elems[nd0Idx_orig]->type == ElementType::TElem)
        static_cast<Element*>(incident_elems[nd0Idx_orig])->ReplaceIncidentElem(this,insertedElem);
    else if(incident_elems[nd0Idx_orig]->type == ElementType::CZ)
        static_cast<CohesiveZone*>(incident_elems[nd0Idx_orig])->ReplaceAdjacentElem(this, insertedElem, nd0Idx_orig);

    insertedElem->incident_elems[ndIdx_orig] = insertedElem_adj;
    insertedElem->incident_elems[nd0Idx_orig] = incident_elems[nd0Idx_orig];
    insertedElem->incident_elems[nd1Idx_orig] = this;
    incident_elems[nd0Idx_orig] = insertedElem;

    // similarly, modify the existing adjacent element
    adjElem->nds[nd1Idx_adj] = split;
    insertedElem_adj->nds[oppIdx_adj] = oppositeNode;
    insertedElem_adj->nds[nd1Idx_adj] = nd1;
    insertedElem_adj->nds[nd0Idx_adj] = split;

    // if the original element had a boundary at nd0Idx, the inserted element now takes that boundary
    if(adjElem->incident_elems[nd0Idx_adj]->type == ElementType::BEdge)
        static_cast<BoundaryEdge*>(adjElem->incident_elems[nd0Idx_adj])->Initialize(insertedElem_adj,nd0Idx_adj);
    else if(adjElem->incident_elems[nd0Idx_adj]->type == ElementType::TElem)
        static_cast<Element*>(adjElem->incident_elems[nd0Idx_adj])->ReplaceIncidentElem(adjElem,insertedElem_adj);
    else if(adjElem->incident_elems[nd0Idx_adj]->type == ElementType::CZ)
        static_cast<CohesiveZone*>(adjElem->incident_elems[nd0Idx_adj])->ReplaceAdjacentElem(adjElem, insertedElem_adj, nd0Idx_adj);

    insertedElem_adj->incident_elems[oppIdx_adj] = insertedElem;
    insertedElem_adj->incident_elems[nd1Idx_adj] = adjElem;
    insertedElem_adj->incident_elems[nd0Idx_adj] = adjElem->incident_elems[nd0Idx_adj];
    adjElem->incident_elems[nd0Idx_adj] = insertedElem_adj;

    oppositeNode->adj_elems.push_back(insertedElem_adj);
    nd1->ReplaceAdjacentElement(adjElem,insertedElem_adj);

    this->PrecomputeInitialArea();
    insertedElem->PrecomputeInitialArea();
    adjElem->PrecomputeInitialArea();
    insertedElem_adj->PrecomputeInitialArea();

    // "fix" palsticity on all four elements
    this->RecalculatePiMultiplierFromDeformationGradient(F_orig);
    insertedElem->RecalculatePiMultiplierFromDeformationGradient(F_orig);
    adjElem->RecalculatePiMultiplierFromDeformationGradient(F_adj);
    insertedElem_adj->RecalculatePiMultiplierFromDeformationGradient(F_adj);

    oppositeNode->PrepareFan();
    nd1->PrepareFan();
    return split;
}

std::pair<icy::Node*,icy::Node*> icy::Element::SplitElemWithCZ(Node *nd, Node *nd0, Node *nd1, double where)
{
    icy::CohesiveZone *adjCZ = dynamic_cast<icy::CohesiveZone*>(getIncidentElementOppositeToNode(nd));
    if(adjCZ == nullptr) throw std::runtime_error("icy::Element::SplitElemWithCZ dynamic cast issue");

    Mesh *m = adjCZ->parentMesh;
    icy::Element *adjElem = adjCZ->getOtherElem(this);
    Node *nd2 = adjCZ->getOtherNode(nd0);
    Node *nd3 = adjCZ->getOtherNode(nd1);

    // for testing: remove CZ, proceed as boundary elem

    adjCZ->Disconnect();
//    auto iter = std::find(m->allCZs.begin(),m->allCZs.end(),adjCZ);
//    if(iter == m->allCZs.end()) throw std::runtime_error("icy::Element::SplitElemWithCZ: CZ not found");
//    m->allCZs.erase(iter);

    m->allCZs.erase(std::remove(m->allCZs.begin(), m->allCZs.end(), adjCZ), m->allCZs.end());
    icy::MeshFragment::CZFactory.release(adjCZ);






    // ===================
    uint8_t ndIdx_orig = getNodeIdx(nd);
    uint8_t nd0Idx_orig = getNodeIdx(nd0);
    uint8_t nd1Idx_orig = getNodeIdx(nd1);

    // preserve deformation gradient
    Eigen::Matrix2d F_orig = getF_at_n();


    icy::Node* split_orig;
//    split_orig = SplitBoundaryElem(nd,nd0,nd1,where);

        MeshFragment *fragment_orig = nd->fragment;

        // insert element
        Element *insertedElem = fragment_orig->AddElement();
        nd->adj_elems.push_back(insertedElem);

        // insert the node between nd0 and nd1; initialize its coordinates; connect to adjacent elements
        split_orig = fragment_orig->AddNode();
        split_orig->InitializeLERP(nd0, nd1, where);
        split_orig->isBoundary = true;
        split_orig->adj_elems.push_back(this);
        split_orig->adj_elems.push_back(insertedElem);

        // modify the original element
        nds[nd1Idx_orig] = split_orig;

        // initialize the inserted element's nodes
        insertedElem->nds[ndIdx_orig] = nd;
        insertedElem->nds[nd1Idx_orig] = nd1;
        insertedElem->nds[nd0Idx_orig] = split_orig;

        // if the original element had a boundary at nd0Idx, the inserted element now takes that boundary
        if(incident_elems[nd0Idx_orig]->type == ElementType::BEdge)
            static_cast<BoundaryEdge*>(incident_elems[nd0Idx_orig])->Initialize(insertedElem,nd0Idx_orig,3);
        else if(incident_elems[nd0Idx_orig]->type == ElementType::TElem)
            static_cast<Element*>(incident_elems[nd0Idx_orig])->ReplaceIncidentElem(this,insertedElem);
        else if(incident_elems[nd0Idx_orig]->type == ElementType::CZ)
            static_cast<CohesiveZone*>(incident_elems[nd0Idx_orig])->ReplaceAdjacentElem(this, insertedElem, nd0Idx_orig);

        // add the boundary that has just appeared
        // automatically initialize the inserted element's adjacency data (insertedElem->incident_elems[ndIdx])
        fragment_orig->AddBoundary(insertedElem,ndIdx_orig,3);
        dynamic_cast<BoundaryEdge*>(incident_elems[ndIdx_orig])->UpdateNodes();

        insertedElem->incident_elems[nd0Idx_orig] = incident_elems[nd0Idx_orig];
        incident_elems[nd0Idx_orig] = insertedElem;
        insertedElem->incident_elems[nd1Idx_orig] = this;

        // from node "nd1", disconnect the original element and replace it with the inserted element
        nd1->ReplaceAdjacentElement(this, insertedElem);

        // compute the new area and reference shape matrix
        this->PrecomputeInitialArea();
        insertedElem->PrecomputeInitialArea();

        // re-evaluate PiMultiplier on both elements to maintain consistent plasticity
        this->RecalculatePiMultiplierFromDeformationGradient(F_orig);
        insertedElem->RecalculatePiMultiplierFromDeformationGradient(F_orig);

        // "fix" the fan for the node, whose element was just replaced
        nd1->PrepareFan();



    Eigen::Matrix2d F_adj = adjElem->getF_at_n();


    Node *oppositeNode = adjElem->getOppositeNode(nd2, nd3);
    uint8_t nd2Idx_adj = adjElem->getNodeIdx(nd2);
    uint8_t nd3Idx_adj = adjElem->getNodeIdx(nd3);
    uint8_t oppIdx_adj = adjElem->getNodeIdx(oppositeNode);

    MeshFragment *fragment = nd2->fragment;

    // insert "adjacent" element
    Element *insertedElem_adj = fragment->AddElement();

    // insert the "split" node between nd0 and nd1
    Node *split_adj = fragment->AddNode();


    split_adj->InitializeLERP(nd2, nd3, where);
    split_adj->adj_elems.insert(split_adj->adj_elems.end(),{adjElem,insertedElem_adj});
    split_adj->isBoundary = true;
    split_adj->weakening_direction = (split_orig->xn - nd->xn).normalized();
    split_adj->isCrackTip = true;

    adjElem->nds[nd3Idx_adj] = split_adj;
    static_cast<BoundaryEdge*>(adjElem->incident_elems[oppIdx_adj])->UpdateNodes();

    insertedElem_adj->nds[oppIdx_adj] = oppositeNode;
    insertedElem_adj->nds[nd3Idx_adj] = nd3;
    insertedElem_adj->nds[nd2Idx_adj] = split_adj;


    // if the original element had a boundary at nd0Idx, the inserted element now takes that boundary
    if(adjElem->incident_elems[nd2Idx_adj]->type == ElementType::BEdge)
        static_cast<BoundaryEdge*>(adjElem->incident_elems[nd2Idx_adj])->Initialize(insertedElem_adj,nd2Idx_adj);
    else if(adjElem->incident_elems[nd2Idx_adj]->type == ElementType::TElem)
        static_cast<Element*>(adjElem->incident_elems[nd2Idx_adj])->ReplaceIncidentElem(adjElem,insertedElem_adj);
    else if(adjElem->incident_elems[nd2Idx_adj]->type == ElementType::CZ)
        static_cast<CohesiveZone*>(adjElem->incident_elems[nd2Idx_adj])->ReplaceAdjacentElem(adjElem, insertedElem_adj, nd2Idx_adj);


    fragment->AddBoundary(insertedElem_adj,oppIdx_adj,3);
    dynamic_cast<BoundaryEdge*>(insertedElem_adj->incident_elems[oppIdx_adj])->UpdateNodes();

    insertedElem_adj->incident_elems[nd3Idx_adj] = adjElem;
    insertedElem_adj->incident_elems[nd2Idx_adj] = adjElem->incident_elems[nd2Idx_adj];
    adjElem->incident_elems[nd2Idx_adj] = insertedElem_adj;

    oppositeNode->adj_elems.push_back(insertedElem_adj);
    nd3->ReplaceAdjacentElement(adjElem,insertedElem_adj);




    adjElem->PrecomputeInitialArea();
    insertedElem_adj->PrecomputeInitialArea();

    adjElem->RecalculatePiMultiplierFromDeformationGradient(F_adj);
    insertedElem_adj->RecalculatePiMultiplierFromDeformationGradient(F_adj);

    oppositeNode->PrepareFan();
    nd3->PrepareFan();
    split_adj->PrepareFan();


    // remove the old boundaries
    nd->fragment->RemoveBoundary(this,ndIdx_orig);
    nd2->fragment->RemoveBoundary(adjElem,oppIdx_adj);
    CohesiveZone *cz1 = m->AddCZ();
    cz1->Initialize(this, ndIdx_orig, adjElem, oppIdx_adj);


    nd->fragment->RemoveBoundary(insertedElem,ndIdx_orig);
    nd2->fragment->RemoveBoundary(insertedElem_adj,oppIdx_adj);
    CohesiveZone *cz2 = m->AddCZ();
    cz2->Initialize(insertedElem,ndIdx_orig, insertedElem_adj,oppIdx_adj);


    return {split_orig,split_adj};

}
