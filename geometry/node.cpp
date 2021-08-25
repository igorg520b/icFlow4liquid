#include "node.h"
#include "element.h"
#include <algorithm>
#include <cmath>
#include "spdlog/spdlog.h"

void icy::Node::Reset()
{
    x_initial = xn = vn = xt = Eigen::Vector2d::Zero();
    eqId = locId = globId = -1;
    area = 0;
    pinned = false;
    spring_attached = 0;
    group.reset();
    crack_tip = false;
}

void icy::Node::Reset(int locId_, double x, double y)
{
    Reset();
    locId = locId_;
    x_initial << x,y;
    intended_position = xt = xn = x_initial;
}

void icy::Node::AddSpringEntries(EquationOfMotionSolver &eq, SimParams &prms, double h, Eigen::Vector2d &spring)
{
    if(eqId<0 || spring_attached<1) return;
    double k = prms.YoungsModulus/300;
    Eigen::Vector2d spr = spring_attachment_position+spring;
    double hsqk = h*h*k;

    Eigen::Vector2d fd = (spr-xt);

    double const_term = hsqk*fd.dot(fd)/2;
    Eigen::Vector2d linear_term = hsqk*(xt-spr);
    Eigen::Matrix2d quadratic_term = hsqk*Eigen::Matrix2d::Identity();
    eq.AddToEquation(const_term, linear_term, quadratic_term, eqId);
}


void icy::Node::PrepareFan()
{
    if(adj_elems.size()==0) throw std::runtime_error("PrepareFan: disconnected node");

    fan.clear();
    fan.resize(adj_elems.size());
    area = 0;
    isBoundary = false;
    for(unsigned k=0;k<adj_elems.size();k++)
    {
        icy::Element *elem = adj_elems[k];
        if(!elem->containsNode(this))  throw std::runtime_error("PrepareFan: mesh topology error");
        area += elem->area_initial/3;

        Sector &s = fan[k];
        s.face = elem;
        Eigen::Vector2d tcv = elem->getCenter() - x_initial;
        s.centerAngle = atan2(tcv.y(), tcv.x());

        short thisIdx, CWIdx, CCWIdx;
        elem->getIdxs(this, thisIdx, CWIdx, CCWIdx);

        s.nd[0] = elem->nds[CWIdx];
        s.nd[1] = elem->nds[CCWIdx];

        const Edge &e0 = elem->CWEdge(this);
        const Edge &e1 = elem->CCWEdge(this);
        if(e0.isBoundary || e1.isBoundary) isBoundary = true;
    }

    std::sort(fan.begin(), fan.end(), [](const Sector &f0, const Sector &f1) {return f0.centerAngle < f1.centerAngle; });

    // if boundary, then ensure that sectors start with a boundary element and end with a boundary element
    if(isBoundary)
    {
        // find the fan element with the border on the CW direction
        auto cw_boundary = std::find_if(fan.begin(), fan.end(), [this](const Sector &f){return f.face->CWEdge(this).isBoundary;});
        if(cw_boundary == fan.end())
        {
            PrintoutFan();
            throw std::runtime_error("cw boundary not found");
        }
        else
        {
            std::rotate(fan.begin(), cw_boundary, fan.end());
        }
    }

    // assert that the nodes of the fan connect
    for(std::size_t i = 0;i<fan.size()-1;i++)
    {
        if(fan[i].nd[1] != fan[i+1].nd[0])
        {
            spdlog::critical("fan nodes are not contiguous for node {}", locId);
            PrintoutFan();
            throw std::runtime_error("fan nodes are not contiguous");
        }
    }
}

void icy::Node::UpdateFan()
{
    // assumes that u and v are normalized
    auto get_angle = [](Eigen::Vector2d u, Eigen::Vector2d v)
    { return acos(std::clamp((double)u.dot(v),-1.0,1.0)); };

    fan_angle_span = 0;

    for(Sector &f : fan)
    {
        // TODO: this can be simplified by using icy::Element::getIdxs()
        f.u_normalized = f.face->CWEdge(this).getVec(this).normalized();
        f.v_normalized = f.face->CCWEdge(this).getVec(this).normalized();

        f.angle0 = fan_angle_span;
        fan_angle_span += get_angle(f.u_normalized,f.v_normalized);
        f.angle1 = fan_angle_span;

        f.u_p << -f.u_normalized.y(), f.u_normalized.x();
        f.v_p << -f.v_normalized.y(), f.v_normalized.x();

        f.t0 = f.face->CauchyStress * f.u_p;
        f.t1 = f.face->CauchyStress * f.v_p;
    }
}


void icy::Node::PrintoutFan()
{
    spdlog::info("Printing fan for node {}; crack_tip: {}; isBoundary{}", locId, crack_tip, isBoundary);
    spdlog::info("fan.size {}; adj_elems.size {}", fan.size(), adj_elems.size());
    for(Sector &s : fan)
        spdlog::info("│ {0:>4}-{1:<4} │ {2: >4}-{3:0^4}-{4: <4} │ ",
                     s.nd[0]->locId, s.nd[1]->locId,
                     s.face->nds[0]->locId,s.face->nds[1]->locId,s.face->nds[2]->locId);
}


uint64_t icy::Node::make_key(Node *nd0, Node *nd1)
{
    int nd0idx = nd0->locId;
    int nd1idx = nd1->locId;

    if(nd0idx > nd1idx) std::swap(nd0idx, nd1idx);
    return ((uint64_t)nd0idx << 32) | nd1idx;
}
