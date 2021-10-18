// Park, Kyoungsoo, and Glaucio H. Paulino. "Computational implementation of
// the PPR potential-based cohesive model in ABAQUS: Educational perspective."
// Engineering fracture mechanics 93 (2012): 239-262.


#ifndef COHESIVEZONE_H
#define COHESIVEZONE_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "equationofmotionsolver.h"
#include "parameters_sim.h"
// #include "node.h"

namespace icy {struct CohesiveZone; struct Node; struct Element; }

struct icy::CohesiveZone
{

    constexpr static int nQPts = 4;     // number of quadrature points
    Node *nds[4];                       // two nodes per side: Side_A_Node_1, Side_A_Node_2, Side_B_Node_1, Side_B_Node_2
    Element *elems[2];                  // each CZ connects two elements
    double pmax[nQPts];                 // max normal separation reached
    double tmax[nQPts];                 // max tangential separation reached
    bool isActive;                      // if not active, CZ has failed
    bool isDamaged;
    //double avgDn, avgDt, avgTn, avgTt;  // average traction-separations for subsequent analysis
    //double maxAvgDn, maxAvgDt;

    void Reset();
    void Initialize(Node *nd1a, Node *nd2a, Node *nd1b, Node *nd2b);
    void Disconnect();
    void AddToSparsityStructure(EquationOfMotionSolver &eq) const;
    bool ComputeEquationEntries(EquationOfMotionSolver &eq, const SimParams &prms, double timeStep);
    void AcceptValues();
    static void CalculateAndPrintBMatrix(); // used to calculate B[]

private:
    bool tentative_contact, tentative_failed, tentative_damaged;
    //double tentative_avgDn, tentative_avgDt, tentative_avgTn, tentative_avgTt;
    double tentative_pmax_final, tentative_tmax_final;
    double tentative_pmax[nQPts], tentative_tmax[nQPts];

    constexpr static double epsilon = -1e-9;
    constexpr static double epsilon_abs = 1e-9;
    constexpr static double epsilon_fail_traction = 0.05; // if traction is <5% of max, CZ fails
    // from https://en.wikipedia.org/wiki/Gaussian_quadrature
    constexpr static double quadraturePoints[nQPts] {-0.8611363115940526,-0.3399810435848563,0.3399810435848563,0.8611363115940526};
    constexpr static double quadratureWeights[nQPts] {0.3478548451374539,0.6521451548625461,0.6521451548625461,0.3478548451374539};
    const static Eigen::Matrix<double,8,2> B[nQPts];

    static double Tn_(const SimParams &prms, const double Dn, const double Dt);
    static double Tt_(const SimParams &prms, const double Dn, const double Dt);
    static double Dnn_(const SimParams &prms, const double opn, const double opt);
    static double Dtt_(const SimParams &prms, const double opn, const double opt);
    static double Dnt_(const SimParams &prms, const double opn, const double opt);

    static void PPR_cohesive_zone_formulation(
        const SimParams &prms,
        const double opn, const double opt,
        bool &cz_contact, bool &cz_failed,
        double &pmax, double &tmax,
        double &Tn, double &Tt, double &Dnn,
        double &Dtt, double &Dnt, double &Dtn);

};

#endif // COHESIVEZONE_H
