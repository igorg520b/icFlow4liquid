// Park, Kyoungsoo, and Glaucio H. Paulino. "Computational implementation of
// the PPR potential-based cohesive model in ABAQUS: Educational perspective."
// Engineering fracture mechanics 93 (2012): 239-262.


#ifndef COHESIVEZONE_H
#define COHESIVEZONE_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include "equationofmotionsolver.h"
#include "parameters_sim.h"
#include "node.h"

namespace icy {struct CohesiveZone; struct Node; struct Element; }

struct icy::CohesiveZone
{
    Node *nds[4];           // two nodes per side: Side_A_Node_1, Side_A_Node_2, Side_B_Node_1, Side_B_Node_2
    Element *elems[2];      // each CZ connects two elements
    double pmax[2];         // max normal separation reached
    double tmax[2];         // max tangential separation reached
    bool isActive;          // if not active, CZ has failed
    double avgDn, avgDt, avgTn, avgTt; // average traction-separations for subsequent analysis
    double maxAvgDn, maxAvgDt;

    void Reset();
    void Initialize(Node *nd1a, Node *nd2a, Node *nd1b, Node *nd2b);
    void AddToSparsityStructure(EquationOfMotionSolver &eq) const;
    bool ComputeEquationEntries(EquationOfMotionSolver &eq, const SimParams &prms);
    void AcceptValues();

private:
    bool tentative_contact, tentative_failed, tentative_damaged;
    double tentative_avgDn, tentative_avgDt, tentative_avgTn, tentative_avgTt;
    double tentative_pmax_final, tentative_tmax_final;
    double tentative_pmax[2], tentative_tmax[2];

    constexpr static double epsilon = -1e-9;
    constexpr static double epsilon_fail_traction = 0.05; // if traction is <5% of max, CZ fails


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