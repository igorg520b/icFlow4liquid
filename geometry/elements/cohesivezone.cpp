#include <iomanip>
#include "cohesivezone.h"
#include "element.h"
#include "node.h"
#include "meshfragment.h"
#include <spdlog/spdlog.h>

const Eigen::Matrix<double,8,2> icy::CohesiveZone::B[nQPts] = {
    (Eigen::Matrix<double,8,2>() <<
         0.9305681557970262,0,
     0,0.9305681557970262,
     0.06943184420297371,0,
     0,0.06943184420297371,
     -0.9305681557970262,-0,
     -0,-0.9305681557970262,
     -0.06943184420297371,-0,
     -0,-0.06943184420297371).finished(),
    (Eigen::Matrix<double,8,2>() <<
         0.6699905217924281,0,
     0,0.6699905217924281,
     0.3300094782075719,0,
     0,0.3300094782075719,
     -0.6699905217924281,-0,
     -0,-0.6699905217924281,
     -0.3300094782075719,-0,
     -0,-0.3300094782075719).finished(),
    (Eigen::Matrix<double,8,2>() <<
         0.3300094782075719,0,
     0,0.3300094782075719,
     0.6699905217924281,0,
     0,0.6699905217924281,
     -0.3300094782075719,-0,
     -0,-0.3300094782075719,
     -0.6699905217924281,-0,
     -0,-0.6699905217924281).finished(),
    (Eigen::Matrix<double,8,2>() <<
         0.06943184420297371,0,
     0,0.06943184420297371,
     0.9305681557970262,0,
     0,0.9305681557970262,
     -0.06943184420297371,-0,
     -0,-0.06943184420297371,
     -0.9305681557970262,-0,
     -0,-0.9305681557970262).finished(),
    };


icy::CohesiveZone::CohesiveZone()
{
    type = ElementType::CZ;
}


void icy::CohesiveZone::Initialize(Element *elem0, uint8_t edgeIdx0, Element *elem1, uint8_t edgeIdx1)
{
    for(int i=0;i<nQPts;i++)
        pmax[i] = tmax[i] = 0;
    isActive = true;
    isDamaged = false;

    elems2[0] = elem0;
    elems2[1] = elem1;
    edgeIds[0] = edgeIdx0;
    edgeIds[1] = edgeIdx1;

    elem0->incident_elems[edgeIdx0] = this;
    elem1->incident_elems[edgeIdx1] = this;
    GetNodes();
}

void icy::CohesiveZone::GetNodes()
{
    nds[0] = elems2[0]->nds[(edgeIds[0]+2)%3];
    nds[1] = elems2[0]->nds[(edgeIds[0]+1)%3];
    nds[2] = elems2[1]->nds[(edgeIds[1]+1)%3];
    nds[3] = elems2[1]->nds[(edgeIds[1]+2)%3];
}

void icy::CohesiveZone::Disconnect()
{
    MeshFragment *fr0 = elems2[0]->nds[0]->fragment;
    MeshFragment *fr1 = elems2[1]->nds[0]->fragment;

    fr0->AddBoundary(elems2[0],edgeIds[0],3);
    fr1->AddBoundary(elems2[1],edgeIds[1],2);

    spdlog::info("disconnecting CZ: elem {} side {}; elem {} side {}", (void*)elems2[0], edgeIds[0], (void*)elems2[1], edgeIds[1]);
}




void icy::CohesiveZone::AddToSparsityStructure(EquationOfMotionSolver &eq)
{
    if(!isActive) return;
    GetNodes();
    int idxs[] {nds[0]->eqId, nds[1]->eqId, nds[2]->eqId, nds[3]->eqId};
    eq.AddEntriesToStructure(std::begin(idxs),std::end(idxs));
}

bool icy::CohesiveZone::ComputeEquationEntries(EquationOfMotionSolver &eq, const SimParams &prms, double h)
{
    if(!isActive) return true;
    GetNodes();

    Eigen::Vector2d xc[4];  // coordinates of the cz nodes
    for(int i=0;i<4;i++)
        xc[i] = nds[i]->xt;

    // center
    Eigen::Vector2d center = (xc[0]+xc[1]+xc[2]+xc[3])*0.25;
    for(int i=0;i<4;i++) xc[i] -= center;

    Eigen::Vector2d dir;
    dir = ((xc[1]-xc[0])+(xc[3]-xc[2]))*0.5;

    double cz_area = dir.norm()*prms.Thickness; // "length" of the zone times an assumed thickness of the 2D material
    dir.normalize();
    Eigen::Matrix2d r;  // aligns midplane with x-axis
    r << dir.y(), dir.x(),
        -dir.x(), dir.y();

    Eigen::Vector2d xr[4];  // rotated cz nodes
    for(int i=0;i<4;i++) xr[i] = r*xc[i];

    Eigen::Vector2d deltaA = xr[2]-xr[0];   // deltaA.x() - tangential opening; deltaA.y() - normal opening
    Eigen::Vector2d deltaB = xr[3]-xr[1];

    bool contact_gp[nQPts] = {};
    bool failed_gp[nQPts] = {};

    Eigen::Matrix<double,8,1> DE;       // energy gradient
    Eigen::Matrix<double,8,8> HE;       // energy Hessian
    DE.setZero();
    HE.setZero();

    // iterate over QPs
    for(int qp=0;qp<nQPts;qp++)
    {
        const double qp_coord = quadraturePoints[qp];
        const double qp_weight = quadratureWeights[qp]/2;

        const double N1 = (1-qp_coord)/2;   // shape functions (their values at current QP)
        const double N2 = (1+qp_coord)/2;

        Eigen::Vector2d openingDisplacement = deltaA*N1 + deltaB*N2;

        double dt = openingDisplacement.x();
        double opt = std::abs(dt);
        double opn = openingDisplacement.y();
        double opt_sign;
        if(opt < epsilon_abs) opt_sign = 0;
        else if(dt < 0) opt_sign = -1;
        else opt_sign = 1;

        tentative_pmax[qp] = pmax[qp];
        tentative_tmax[qp] = tmax[qp];

        double Tn, Tt, Dnn, Dtt, Dnt, Dtn;

        PPR_cohesive_zone_formulation(prms, opn, opt,
                                      contact_gp[qp], failed_gp[qp],
                                      tentative_pmax[qp], tentative_tmax[qp], Tn, Tt, Dnn, Dtt, Dnt, Dtn);

        Eigen::Vector2d T(Tt*opt_sign,Tn);
        Eigen::Matrix2d DT;
        DT << Dtt*opt_sign*opt_sign, Dtn*opt_sign, Dtn*opt_sign, Dnn;

        DE -= (B[qp]*T)*(cz_area*qp_weight);
        HE -= (B[qp]*DT*B[qp].transpose())*(cz_area*qp_weight);

        // ensure that CZ is not openend too quickly
        double delta_n = tentative_pmax[qp] - pmax[qp];
        double delta_t = tentative_tmax[qp] - tmax[qp];
        constexpr double coeff = 0.35;   // max progression coeff

        if(pmax[qp]<prms.cz_nThreshold && delta_n > prms.cz_nThreshold*coeff) return false;
        if(pmax[qp]>=prms.cz_nThreshold && delta_n > prms.cz_del_n*coeff) return false;
        if(tmax[qp] < prms.cz_tThreshold && delta_t > prms.cz_tThreshold*coeff) return false;
        if(tmax[qp] >= prms.cz_tThreshold && delta_t > prms.cz_del_t*coeff) return false;

    }

    // rotate back to the initial reference frace
    Eigen::Matrix2d rT = r.transpose();
    Eigen::Matrix<double,8,8> RT;
    Eigen::Matrix2d z;
    z.setZero();
    RT << rT,z, z, z,
          z, rT,z, z,
          z, z, rT,z,
          z, z, z, rT;

    DE=RT*DE;
    HE=RT*HE;

    double hsq = h*h;
    DE *= hsq;
    HE *= hsq;

    tentative_pmax_final = *std::max_element(std::begin(tentative_pmax),std::end(tentative_pmax));
    tentative_tmax_final = *std::max_element(std::begin(tentative_tmax),std::end(tentative_tmax));

    tentative_failed = std::any_of(std::begin(failed_gp),std::end(failed_gp),[](bool k){return k;});
    tentative_contact = std::any_of(std::begin(contact_gp),std::end(contact_gp),[](bool k){return k;});

    tentative_damaged = false;
    if(!tentative_failed)
    {
        for(int i=0;i<nQPts;i++)
            if(tentative_pmax[i] >= prms.cz_nThreshold || tmax[i] >= prms.cz_tThreshold)
            { tentative_damaged = true; break; }
    }

    eq.AddToEquation(DE.data(), HE.data(), {nds[0]->eqId,nds[1]->eqId,nds[2]->eqId,nds[3]->eqId});
    return true;
}

void icy::CohesiveZone::AcceptValues()
{
    if(!isActive) return;
    isActive = !tentative_failed;
    for(int i=0;i<nQPts;i++)
    {
        pmax[i] = tentative_pmax[i];
        tmax[i] = tentative_tmax[i];
    }

    if(tentative_damaged) isDamaged = true;
}

double icy::CohesiveZone::Tn_(const SimParams &prms, const double Dn, const double Dt)
{
    const double &deln = prms.cz_del_n;
    const double &delt = prms.cz_del_t;
    const double &p_m = prms.cz_p_m;
    const double &p_n = prms.cz_p_n;
    const double &alpha = prms.cz_alpha;
    const double &beta = prms.cz_beta;
    const double &gam_n = prms.cz_gam_n;
    const double &gam_t = prms.cz_gam_t;
    const double &pMtn = prms.cz_pMtn;

    double Dndn = Dn / deln;
    double Dtdt = Dt / delt;
    double expr2 = p_m / alpha + Dndn;
    double pr1 = gam_n / deln;
    double pr2 = (p_m * pow(1 - Dndn, alpha) * pow(expr2, p_m - 1)) -
                 (alpha * pow(1 - Dndn, alpha - 1) * pow(expr2, p_m));
    double pr3 = gam_t * pow(1 - Dtdt, beta) * pow(p_n / beta + Dtdt, p_n) + pMtn;
    return pr1 * pr2 * pr3;
}

double icy::CohesiveZone::Tt_(const SimParams &prms, const double Dn, const double Dt)
{
    const double &deln = prms.cz_del_n;
    const double &delt = prms.cz_del_t;
    const double &p_m = prms.cz_p_m;
    const double &p_n = prms.cz_p_n;
    const double &alpha = prms.cz_alpha;
    const double &beta = prms.cz_beta;
    const double &gam_n = prms.cz_gam_n;
    const double &gam_t = prms.cz_gam_t;
    const double &pMnt = prms.cz_pMnt;

    double Dndn = Dn / deln;
    double Dtdt = Dt / delt;
    double expr1 = 1 - Dtdt;
    double expr2 = p_n / beta + Dtdt;
    double pr1 = gam_t / delt;
    double pr2 = p_n * pow(expr1, beta) * pow(expr2, p_n - 1) - beta * pow(expr1, beta - 1) * pow(expr2, p_n);
    double pr3 = gam_n * pow(1 - Dndn, alpha) * pow(p_m / alpha + Dndn, p_m) + pMnt;
    return pr1 * pr2 * pr3;
}

double icy::CohesiveZone::Dnn_(const SimParams &prms, const double opn, const double opt)
{
    const double &deln = prms.cz_del_n;
    const double &delt = prms.cz_del_t;
    const double &p_m = prms.cz_p_m;
    const double &p_n = prms.cz_p_n;
    const double &alpha = prms.cz_alpha;
    const double &beta = prms.cz_beta;
    const double &gam_n = prms.cz_gam_n;
    const double &gam_t = prms.cz_gam_t;
    const double &pMtn = prms.cz_pMtn;

    double coeff = gam_n / (deln * deln);
    double expr1 = (p_m * p_m - p_m) * pow(1.0 - (opn / deln), alpha) * pow((p_m / alpha) + (opn / deln), p_m - 2.0);
    double expr2 = (alpha * alpha - alpha) * pow(1.0 - (opn / deln), alpha - 2.0) * pow((p_m / alpha) + (opn / deln), p_m);
    double expr3 = 2.0 * alpha * p_m * pow(1.0 - (opn / deln), alpha - 1.0) * pow((p_m / alpha) + (opn / deln), p_m - 1.0);
    double expr4 = gam_t * pow((1.0 - (opt / delt)), beta) * pow(((p_n / beta) + (opt / delt)), p_n) + pMtn;
    double result = coeff * (expr1 + expr2 - expr3) * expr4;
    return result;
}

double icy::CohesiveZone::Dtt_(const SimParams &prms, const double opn, const double opt)
{
    const double &deln = prms.cz_del_n;
    const double &delt = prms.cz_del_t;
    const double &p_m = prms.cz_p_m;
    const double &p_n = prms.cz_p_n;
    const double &alpha = prms.cz_alpha;
    const double &beta = prms.cz_beta;
    const double &gam_n = prms.cz_gam_n;
    const double &gam_t = prms.cz_gam_t;
    const double &pMnt = prms.cz_pMnt;

    double coeff = gam_t / (delt * delt);
    double expr1 = (p_n * p_n - p_n) * pow(1.0 - (opt / delt), beta) * pow((p_n / beta) + (opt / delt), p_n - 2.0);
    double expr2 = (beta * beta - beta) * pow(1.0 - (opt / delt), beta - 2.0) * pow((p_n / beta) + (opt / delt), p_n);
    double expr3 = 2.0 * beta * p_n * pow(1.0 - (opt / delt), beta - 1.0) * pow((p_n / beta) + (opt / delt), p_n - 1.0);
    double expr4 = gam_n * pow(1.0 - (opn / deln), alpha) * pow((p_m / alpha) + (opn / deln), p_m) + pMnt;
    double result = coeff * (expr1 + expr2 - expr3) * expr4;
    return result;
}

double icy::CohesiveZone::Dnt_(const SimParams &prms, const double opn, const double opt)
{
    const double &deln = prms.cz_del_n;
    const double &delt = prms.cz_del_t;
    const double &p_m = prms.cz_p_m;
    const double &p_n = prms.cz_p_n;
    const double &alpha = prms.cz_alpha;
    const double &beta = prms.cz_beta;
    const double &gam_n = prms.cz_gam_n;
    const double &gam_t = prms.cz_gam_t;

    double coeff = gam_n * gam_t / (deln * delt);
    double expr1 = p_m * pow(1.0 - (opn / deln), alpha) * pow((p_m / alpha) + (opn / deln), p_m - 1.0);
    double expr2 = alpha * pow(1.0 - (opn / deln), alpha - 1.0) * pow((p_m / alpha) + (opn / deln), p_m);
    double expr3 = p_n * pow(1.0 - (opt / delt), beta) * pow((p_n / beta) + (opt / delt), p_n - 1.0);
    double expr4 = beta * pow(1.0 - (opt / delt), beta - 1.0) * pow((p_n / beta) + (opt / delt), p_n);
    double result = coeff * (expr1 - expr2) * (expr3 - expr4);
    return result;
}

void icy::CohesiveZone::PPR_cohesive_zone_formulation(
    const SimParams &prms,
    const double opn, const double opt,
    bool &cz_contact, bool &cz_failed,
    double &pmax, double &tmax,
    double &Tn, double &Tt, double &Dnn,
    double &Dtt, double &Dnt, double &Dtn)
{
    const double &deln = prms.cz_del_n;
    const double &delt = prms.cz_del_t;
    const double &sigma_max = prms.cz_sigma_max;
    const double &tau_max = prms.cz_tau_max;
    const double &lambda_n = prms.cz_lambda_n;
    const double &lambda_t = prms.cz_lambda_t;

    Tn = Tt = Dnn = Dtt = Dnt = Dtn = 0;
    if (opn > deln || opt > delt)
    {
        cz_contact = false;
        cz_failed = true;
        return;
    }
    cz_contact = (opn < 0);
    cz_failed = false;
    const double threshold_tangential = tau_max * epsilon_fail_traction;
    const double threshold_normal = sigma_max * epsilon_fail_traction;

    if (cz_contact)
    {
        Dnt = 0;
        if (pmax != 0)
        {
            double peakTn = Tn_(prms, pmax, tmax);
            Tn = peakTn * opn / pmax;
            Dnn = peakTn / pmax;
        }
        else
        {
            Dnn = Dnn_(prms, 0, tmax);
            Tn = Dnn * opn;
        }

        Tt = Tt_(prms, 0, opt);
        if (Tt >= epsilon && !(opt > delt * lambda_t / 5 && Tt < threshold_tangential))
        {
            if (opt >= tmax)
            {
                // tangential softening
                tmax = opt;
                Dtt = Dtt_(prms, 0, opt);
            }
            else
            {
                // unload/reload
                double peakTt = Tt_(prms, 0, tmax);
                Tt = peakTt * opt / tmax;
                Dtt = peakTt / tmax;
            }

        }
        else
        {
            // cz failed in tangential direction while in contact
            Tt = Dtt = Dnt = 0;
            Tn = Dnn = 0;
            cz_failed = true;
        }
    }
    else
    {
        // not in contact
        Tt = Tt_(prms, opn, opt);
        Tn = Tn_(prms, opn, opt);
        if (Tt >= epsilon && Tn >= epsilon &&
            !(opt > delt * lambda_t / 5 && Tt < threshold_tangential) &&
            !(opn > deln * lambda_n / 5 && Tn < threshold_normal))
        {
            // tangential component
            bool tsoft = (opt >= tmax);
            bool nsoft = (opn >= pmax);
            if (tsoft && nsoft)
            {
                // tangential and normal softening
                tmax = opt;
                pmax = opn;
                Dnn = Dnn_(prms, opn, opt);
                Dnt = Dnt_(prms, opn, opt);
                Dtt = Dtt_(prms, opn, opt);
            }
            else if (tsoft && !nsoft)
            {
                Dnt = 0;
                if (pmax != 0)
                {
                    double peakTn = Tn_(prms, pmax, tmax);
                    Tn = peakTn * opn / pmax;
                    Dnn = peakTn / pmax;
                }
                else
                {
                    Tn = 0; Dnn = Dnn_(prms, 0, tmax);
                }

                // normal unload/reload
                tmax = opt;
                Tt = Tt_(prms, pmax, opt);
                Dtt = Dtt_(prms, pmax, opt);
            }
            else if (!tsoft && nsoft)
            {
                Dnt = 0;
                if (tmax != 0)
                {
                    double peakTt = Tt_(prms, pmax, tmax);
                    Tt = peakTt * opt / tmax;
                    Dtt = peakTt / tmax;
                }
                else
                {
                    Tt = 0; Dtt = Dtt_(prms, pmax, 0);
                }

                pmax = opn;
                Tn = Tn_(prms, pmax, tmax);
                Dnn = Dnn_(prms, pmax, tmax);

            }
            else
            {
                Dnt = 0;
                // reloading in both tangential and normal
                double peakTn = Tn_(prms, pmax, tmax);
                if (pmax != 0)
                {
                    Tn = peakTn * opn / pmax;
                    Dnn = peakTn / pmax;
                }
                else
                {
                    Tn = 0; Dnn = Dnn_(prms, 0, tmax);
                }

                if (tmax != 0)
                {
                    double peakTt = Tt_(prms, pmax, tmax);
                    Tt = peakTt * opt / tmax;
                    Dtt = peakTt / tmax;
                }
                else
                {
                    Tt = 0; Dtt = Dtt_(prms, pmax, 0);
                }
            }
        }
        else
        {
            cz_failed = true;
            Tn = Tt = Dnn = Dtt = Dnt = 0;
        }
    }
    Dtn = Dnt;
}

void icy::CohesiveZone::CalculateAndPrintBMatrix()
{
    std::cout << "icy::CohesiveZone::CalculateAndPrintBMatrix()\n";
    for(int i=0;i<nQPts;i++)
    {
        Eigen::Matrix<double,8,2> m;
        m.setZero();
        double xi = quadraturePoints[i];
        const double N1 = (1-xi)/2;   // shape functions (their values at current QP)
        const double N2 = (1+xi)/2;
        m << Eigen::Matrix2d::Identity()*N1, Eigen::Matrix2d::Identity()*N2,
            -Eigen::Matrix2d::Identity()*N1, -Eigen::Matrix2d::Identity()*N2;
        std::cout << std::setprecision(16);
        std::cout << "B[" << i << "]\n";
        for(int k=0;k<8;k++)
        {
            for(int l=0;l<2;l++)
                std::cout << m(k,l) << ',';
            std::cout << '\n';
        }
    }
}

