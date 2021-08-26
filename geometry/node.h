#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers
#ifndef NODE_H
#define NODE_H

#include <bitset>
#include <Eigen/Core>
#include <boost/container/small_vector.hpp>
#include "equationofmotionsolver.h"
#include "parameters_sim.h"
#include "edge.h"

namespace icy { struct Node; struct Element;}

struct icy::Node
{
    Node(){ Reset();}
    void Reset();
    void Reset(int locId, double x, double y);
    void Reset(Node *other) { Reset(other->locId, other->x_initial.x(), other->x_initial.y()); }

    int locId, globId, eqId, indId;    // id in fragment; id in mesh; id in freenode list; id in movable boundary
    std::size_t gmshTag;
    bool pinned;                // the position of the node is set externally
    double area;                // area that the node represents, for applying body forces
    std::bitset<8> group;       // for meshing/initialization

    Eigen::Vector2d x_initial;  // initial configuration
    Eigen::Vector2d xn, vn;     // position and velocity at step n
    Eigen::Vector2d xt;         // tentative coordinates
    Eigen::Vector2d x_hat;      // refer to Minchen Li's paper
    Eigen::Vector2d intended_position; // when manipulating via GUI during a running simulation

    // add spring forces when manipulating the deformable object via GUI
    void AddSpringEntries(EquationOfMotionSolver &eq, SimParams &prms, double timeStep, Eigen::Vector2d &spring);
    double spring_attached;
    Eigen::Vector2d spring_attachment_position;


// FRACTURE MODEL
    struct Sector
    {
        icy::Element *face;
        icy::Node* nd[2];
        double centerAngle; // angle from the node to the center of the adjacent element
        double angle0, angle1;
        Eigen::Vector2d u_normalized, v_normalized, u_p, v_p;
        Eigen::Vector2d t0, t1;
    };

    // separation stress
    struct SepStressResult
    {
        double angle_fwd, angle_bwd;
        icy::Element* faces[2];
        Eigen::Vector2d traction[2];
        Eigen::Vector2d tn, tn_p;
        double phi[2], theta[2];
        double trac_normal, trac_tangential;
        icy::Edge e[4];
        icy::Edge e_opposite[2]; // edges that lie opposite of the center node
    };

    boost::container::small_vector<icy::Element*, 8> adj_elems;
    boost::container::small_vector<icy::Node::Sector,8> fan;
    double fan_angle_span;  // assigned in UpdateFan();
    bool isCrackTip;
    bool isBoundary;
    SepStressResult result_with_max_traction;
    Eigen::Vector2d dir;
    Eigen::Vector2d weakening_direction;    // used when isCrackTip==true
    double max_normal_traction;

    void PrepareFan();  // performed when topology changes
    void PrintoutFan(); // for testing
    void ComputeFanVariables(SimParams &prms);

    static uint64_t make_key(Node *nd0, Node *nd1); // return unique id for a segment defined by two nodes

private:
    void UpdateFan();   // performed when tentative displacements and stress distribution change; invoked from ComputeFanVariables()
    double NormalTraction(double angle_fwd, double weakening_coeff) const;
    void EvaluateTractions(double angle_fwd, SepStressResult &ssr, const double weakening_coeff) const;
};

#endif // NODE_H
#endif // Q_MOC_RUN
