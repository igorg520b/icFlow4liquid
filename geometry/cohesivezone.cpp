#include "cohesivezone.h"

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

bool icy::CohesiveZone::ComputeEquationEntries(EquationOfMotionSolver &eq, const SimParams &prms, double timeStep)
{

}

void icy::CohesiveZone::AcceptValues()
{

}
