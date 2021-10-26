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
    Element *elems2[2];                  // each CZ connects two elements
    uint8_t edgeIds[2];                 // side in the corresponding element to which CZ connects

    // two nodes per side: Side_A_Node_1, Side_A_Node_2, Side_B_Node_1, Side_B_Node_2
    // NOTE: this is populated at a later stage from elems[] and edgeIds[]
    Node *nds[4];


    double pmax[nQPts];                 // max normal separation reached
    double tmax[nQPts];                 // max tangential separation reached
    bool isActive;                      // if not active, CZ has failed
    bool isDamaged;
    bool traversed;                 // used for identifying disjoint regions and n-regions around crack tips


    void Reset();
    void Initialize(Node *nd1a, Node *nd2a, Node *nd1b, Node *nd2b);

    // first element/edge must be CW, sedond will be CCW
    void Initialize(Element *elem0, uint8_t edgeIdx0, Element *elem1, uint8_t edgeIdx1);
    void Disconnect();
    void AddToSparsityStructure(EquationOfMotionSolver &eq);
    bool ComputeEquationEntries(EquationOfMotionSolver &eq, const SimParams &prms, double timeStep);
    void AcceptValues();
    static void CalculateAndPrintBMatrix(); // used to calculate B[]
    void GetNodes(); // convert elems[2], edgeIds[2] into nodes[4]

private:
    bool tentative_contact, tentative_failed, tentative_damaged;
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