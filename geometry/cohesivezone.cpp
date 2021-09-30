#include "cohesivezone.h"
#include "node.h"

void icy::CohesiveZone::Reset()
{
    pmax[0] = pmax[1] = 0;
    tmax[0] = tmax[1] = 0;
    isActive = true;
    avgDn = avgDt = avgTn = avgTt = 0; // average traction-separations for subsequent analysis
    maxAvgDn = maxAvgDt = 0;
}

void icy::CohesiveZone::Initialize(Node *nd1a, Node *nd2a, Node *nd1b, Node *nd2b)
{

}


void icy::CohesiveZone::AddToSparsityStructure(EquationOfMotionSolver &eq) const
{

}

bool icy::CohesiveZone::ComputeEquationEntries(EquationOfMotionSolver &eq, const SimParams &prms, double h)
{
    Eigen::Vector2d xc[4];  // coordinates of the cz nodes
    for(int i=0;i<4;i++) xc[i] = nds[i]->xt;

    // "midplane" between the two boundaries of the cz
    Eigen::Vector2d mp[2] = {(xc[0]+xc[2])/2, (xc[1]+xc[3])/2};

    // center
    Eigen::Vector2d center = (mp[0]+mp[1])/2;
    for(int i=0;i<4;i++) xc[i] -= center;
    for(int i=0;i<2;i++) mp[i] -= center;

    Eigen::Vector2d dir = (mp[1]-mp[0]).normalized();
    Eigen::Matrix2d R;  // aligns midplane with x-axis
    R << dir.x(), dir.y(), -dir.x(), dir.y();
    Eigen::Matrix2d RT = R.transpose();

    Eigen::Vector2d xr[4];  // rotated cz nodes
    for(int i=0;i<4;i++) xr[i] = R*xc[i];
    Eigen::Vector2d deltaA = xr[2]-xr[0];   // deltaA.x() - tangential opening; deltaA.y() - normal opening
    Eigen::Vector2d deltaB = xr[3]-xr[1];

    // from https://en.wikipedia.org/wiki/Gaussian_quadrature
    constexpr double quadraturePoints[nQPts] {-0.8611363115940526,-0.3399810435848563,0.3399810435848563,0.8611363115940526};
    constexpr double quadratureWeights[nQPts] {0.3478548451374539,0.6521451548625461,0.6521451548625461,0.3478548451374539};

    bool contact_gp[nQPts] = {};
    bool failed_gp[nQPts] = {};

    Eigen::Matrix<double,8,1> DE = Eigen::Matrix<double,8,1>::Zero();    // energy gradient
    Eigen::Matrix<double,8,8> HE = Eigen::Matrix<double,8,8>::Zero();       // energy Hessian

    // iterate over QPs
    for(int qp=0;qp<nQPts;qp++)
    {
        const double qp_coord = quadraturePoints[qp];
        const double qp_weight = quadratureWeights[qp];

        const double N1 = (1-qp_coord)/2;   // shape functions (their values at current QP)
        const double N2 = (1+qp_coord)/2;

        Eigen::Vector2d openingDisplacement = deltaA*N1 + deltaB*N2;
        double dt = openingDisplacement.x();
        double opt = std::abs(dt);
        double opn = openingDisplacement.y();
        double opt_sign;
        if(dt < epsilon) opt_sign = 0;
        else if(dt < 0) opt_sign = -1;
        else opt_sign = 1;

        tentative_pmax[qp] = pmax[qp];
        tentative_tmax[qp] = tmax[qp];

        double Tn, Tt, Dnn, Dtt, Dnt, Dtn;

        PPR_cohesive_zone_formulation(prms,
            opn, opt,
            contact_gp[qp], failed_gp[qp],
            tentative_pmax[qp], tentative_tmax[qp], Tn, Tt, Dnn, Dtt, Dnt, Dtn);

        Eigen::Vector2d T(Tt*opt_sign,Tn);
        Eigen::Matrix2d DT;
        DT << Dtt, Dtn*opt_sign, Dtn*opt_sign, Dnn;

    }
}

void icy::CohesiveZone::AcceptValues()
{

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
