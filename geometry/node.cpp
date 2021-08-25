//#include <cmath>
#include "node.h"
#include "spdlog/spdlog.h"
#include "element.h"

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
    throw std::runtime_error("not implemented");
}

void icy::Node::UpdateFan()
{
    throw std::runtime_error("not implemented");
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


/*

void icy::Node::PrepareFan2()
{
    Eigen::Vector3d nd_vec = x_initial;

    unsigned nElems = adjacent_elems.size();
    if(nElems == 0)
    {
        qDebug() << "node " << this->locId;
        throw std::runtime_error("disconnected node");
    }
    fan.clear();
    fan.reserve(nElems);
    area = 0;
    for(unsigned k=0;k<nElems;k++)
    {
        icy::Element *elem = adjacent_elems[k];
        if(!elem->ContainsNode(this))  throw std::runtime_error("PrepareFan2: mesh topology error 0");
        area += elem->area_initial/3;

        Sector s;
        s.face = elem;
        Eigen::Vector3d tcv = elem->getCenter() - nd_vec;
        s.centerAngle = atan2(tcv.y(), tcv.x());

        short thisIdx, CWIdx, CCWIdx;
        elem->getIdxs(this, thisIdx, CWIdx, CCWIdx);

        s.nd[0] = elem->nds[CWIdx];
        s.nd[1] = elem->nds[CCWIdx];

        // note that the indices are swapped
        Edge e0 = elem->CWEdge(this);
        Edge e1 = elem->CCWEdge(this);
        fan.push_back(s);

        if(e0.isBoundary || e1.isBoundary) isBoundary = true;
    }

    std::sort(fan.begin(), fan.end(),
              [](const Sector &f0, const Sector &f1)
    {return f0.centerAngle < f1.centerAngle; });

    if(isBoundary) // assert means that PrepareFan2 is not called from Fix_X
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
            std::cout << "\n\n\nfan nodes are not contiguous " << locId << std::endl;
            PrintoutFan();
            throw std::runtime_error("fan nodes are not contiguous");
        }
//        if(fan[i].e[1].nds[0] != fan[i+1].e[0].nds[0] || fan[i].e[1].nds[1] != fan[i+1].e[0].nds[1])
//            throw std::runtime_error("edges not shared");
    }
}




void icy::Node::InitializeFan()
{
    auto get_angle = [](Eigen::Vector2f u, Eigen::Vector2f v)
    {
        double dot = u.dot(v);
//        double dot = u.dot(v)/(u.norm()*v.norm());
        if(dot > 1) dot = 1.0;
        else if(dot < -1.0) dot = -1.0;
        return acos(dot);
    };

    fan_angle_span = 0;

    for(Sector &f : fan)
    {
        // TODO: Simplify these two statements
        f.u_normalized = f.face->CWEdge(this).getVec(this).normalized();
        f.v_normalized = f.face->CCWEdge(this).getVec(this).normalized();

        f.angle0 = fan_angle_span;
        fan_angle_span += get_angle(f.u_normalized,f.v_normalized);
        f.angle1 = fan_angle_span;

        f.u_p << -f.u_normalized.y(), f.u_normalized.x();
        f.v_p << -f.v_normalized.y(), f.v_normalized.x();

        f.t0_top << f.face->str_top * f.u_p;
        f.t1_top << f.face->str_top * f.v_p;
        f.t0_bottom << f.face->str_bottom * f.u_p;
        f.t1_bottom << f.face->str_bottom * f.v_p;
    }
}

*/


uint64_t icy::Node::make_key(Node *nd0, Node *nd1)
{
    int nd0idx = nd0->locId;
    int nd1idx = nd1->locId;

    if(nd0idx > nd1idx) std::swap(nd0idx, nd1idx);
    return ((uint64_t)nd0idx << 32) | nd1idx;
}
