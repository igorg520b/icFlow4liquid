#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers
#ifndef NODE_H
#define NODE_H

#include <bitset>
#include <Eigen/Core>
#include <boost/container/small_vector.hpp>
#include "equationofmotionsolver.h"
#include "parameters_sim.h"
#include "edge.h"

namespace icy { struct Node; struct Element; struct CohesiveZone; class MeshFragment;}

struct icy::Node
{
    Node() { Reset(); }
    ~Node() = default;
    Node& operator=(Node&) = delete;

    void Reset();
    void Initialize(double x, double y);
    void Initialize(const Node *other);
    void InitializeLERP(const Node *nd0, const Node *nd1, double f);    // linear interpolaiton between two other nodes


    int locId, globId, eqId, indId;    // id in fragment; id in mesh; id in freenode list; id in movable boundary
    bool pinned;                // the position of the node is set externally
    double area;                // area that the node represents, for applying body forces
    std::bitset<8> group;       // for meshing/initialization
    bool isBoundary;

    Eigen::Vector2d x_initial;  // initial configuration
    Eigen::Vector2d xn, vn;     // position and velocity at step n
    Eigen::Vector2d xt;         // tentative coordinates
    Eigen::Vector2d x_hat;      // refer to Minchen Li's paper
    Eigen::Vector2d intended_position; // when manipulating via GUI during a running simulation

    // add spring forces when manipulating the deformable object via GUI
    void AddSpringEntries(EquationOfMotionSolver &eq, const SimParams &prms, double timeStep, Eigen::Vector2d &spring);
    double spring_attached;
    Eigen::Vector2d spring_attachment_position;


// FRACTURE MODEL
    MeshFragment *fragment;         // mesh fragment to which the node belongs
    unsigned short traversal;

    struct Sector
    {
        icy::Element *face;
        icy::Node* nd[2];
        double centerAngle; // angle from the node to the center of the adjacent element
        double angle0, angle1;
        Eigen::Vector2d u_normalized, v_normalized, u_p, v_p;   // u is CW, v is CCW
        Eigen::Vector2d t0, t1;
        bool operator<(const Sector& other) const {return centerAngle < other.centerAngle;}
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
        double angle0[2], angle1[2];
        double sectorSpan(const int idx) const {return theta[idx]+phi[idx];}
    };

    boost::container::small_vector<icy::Element*, 8> adj_elems;
    boost::container::small_vector<icy::Node::Sector,8> fan;
    double fan_angle_span;  // assigned in UpdateFan();
    bool isCrackTip;
    SepStressResult result_with_max_traction;
    Eigen::Vector2d dir;
    Eigen::Vector2d weakening_direction;    // used when isCrackTip==true
    double max_normal_traction;
    double time_loaded_above_threshold;

    void CreateUnrotatedFan();
    void PrepareFan();  // performed when topology changes
    void PrintoutFan(); // for testing
    void ComputeFanVariables(const SimParams &prms);
    void ReplaceAdjacentElement(Element *originalElem, Element *replacement);

    static uint64_t make_local_key(Node *nd0, Node *nd1); // return unique id for a segment defined by two nodes
    static uint64_t make_global_key(Node *nd0, Node *nd1); // return unique id for a segment defined by two nodes

private:
    void UpdateFan();   // performed when tentative displacements and stress distribution change; invoked from ComputeFanVariables()
    double NormalTraction(double angle_fwd, double weakening_coeff) const;
    void EvaluateTractions(double angle_fwd, SepStressResult &ssr, const double weakening_coeff) const;
};

#endif // NODE_H
#endif // Q_MOC_RUN
