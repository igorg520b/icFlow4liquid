#ifndef INTERACTION_H
#define INTERACTION_H

#include<Eigen/Core>
#include <functional>
#include <node.h>
#include "equationofmotionsolver.h"
#include "parameters_sim.h"

namespace icy { class Interaction; struct Node;}

class icy::Interaction
{
public:
    Node *ndA, *ndB, *ndP;
    Eigen::Vector2d D;

    Interaction() {};
    Interaction(Node* A, Node* B, Node* P, Eigen::Vector2d D_) : ndA(A), ndB(B), ndP(P), D(D_) {};
    void AddToSparsityStructure(EquationOfMotionSolver &eq) const;
    void Evaluate(EquationOfMotionSolver &eq, SimParams &prms, double h);
    double static SegmentPointDistance(Eigen::Vector2d A, Eigen::Vector2d B, Eigen::Vector2d P, Eigen::Vector2d &D);

private:
    void static distance(Eigen::Vector2d (&p)[4], double &d, double &t, Eigen::Matrix<double,6,1> &Dd, Eigen::Matrix<double,6,6> &DDd);

    void static potential(double dHat, double d, Eigen::Matrix<double,6,1> &Dd, Eigen::Matrix<double,6,6> &DDd,
                          double &p, Eigen::Matrix<double,6,1> &Dp, Eigen::Matrix<double,6,6> &DDp);

    // for unordered_map
public:
    bool operator==(const Interaction& other) const
    {
        return (ndP==other.ndP && ((ndA==other.ndA && ndB==other.ndB)||(ndA==other.ndB && ndB==other.ndA)));
    }
};

namespace std
{
    template <>
    struct hash<icy::Interaction>
    {
        size_t operator()(const icy::Interaction& k) const
        {
            // Compute individual hash values for two data members and combine them using XOR and bit shifting
            uint64_t c = ((uint64_t)k.ndA->globId << 32 || k.ndB->globId) ^ (k.ndP->globId << 24);
            return std::hash<uint64_t>{}(c);
        }
    };
}

#endif // INTERACTION_H
