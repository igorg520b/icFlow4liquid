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
    Node *nds[4];           // two nodes per side
    Element *elems[2];      // each CZ connects two elements
    double pmax[2];         // max normal separation reached
    double tmax[2];         // max tangential separation reached
    bool isActive;          // if not active, CZ has failed
    double avgDn, avgDt, avgTn, avgTt; // average traction-separations for subsequent analysis
    double maxAvgDn, maxAvgDt;

    void Reset();
    void Initialize(Node *nd1a, Node *nd2a, Node *nd1b, Node *nd2b);
    void AddToSparsityStructure(EquationOfMotionSolver &eq) const;
    bool ComputeEquationEntries(EquationOfMotionSolver &eq, const SimParams &prms, double timeStep);
    void AcceptValues();

private:
    bool tentative_contact, tentative_failed, tentative_damaged;
    double tentative_avgDn, tentative_avgDt, tentative_avgTn, tentative_avgTt;
    double tentative_pmax, tentative_tmax;
    double tentative_pmax[2], tentative_tmax[2];

};

#endif // COHESIVEZONE_H
