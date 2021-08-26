#include "node.h"
#include "element.h"
#include <algorithm>
#include <cmath>
#include "spdlog/spdlog.h"
#include "boost/math/tools/minima.hpp"

void icy::Node::Reset()
{
    x_initial = xn = vn = xt = Eigen::Vector2d::Zero();
    eqId = locId = globId = -1;
    area = 0;
    pinned = false;
    spring_attached = 0;
    group.reset();
    isCrackTip = false;
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
    spdlog::info("Printing fan for node {}; isCrackTip: {}; isBoundary{}", locId, isCrackTip, isBoundary);
    spdlog::info("fan.size {}; adj_elems.size {}", fan.size(), adj_elems.size());
    for(Sector &s : fan)
        spdlog::info("│ {0:>4}-{1:<4} │ {2: >4}-{3:0^4}-{4: <4} │ ",
                     s.nd[0]->locId, s.nd[1]->locId,
                     s.face->nds[0]->locId,s.face->nds[1]->locId,s.face->nds[2]->locId);
}

void icy::Node::ComputeFanVariables(SimParams &prms)
{
    UpdateFan();

    dir = Eigen::Vector2d::Zero();
    max_normal_traction = 0;
    unsigned nFan = fan.size();

    double weakening_coeff = prms.FractureWeakeningCoeff;

    unsigned gridPts = isBoundary ? nFan+1 : nFan;

    float grid_results[gridPts];
    for(unsigned i=0; i<nFan; i++)
    {
        grid_results[i] = normal_traction(fan[i].angle0, weakening_coeff);
        if(std::isnan(grid_results[i])) throw std::runtime_error("traction is nan");
    }
    if(isBoundary) {
        grid_results[nFan] = normal_traction(fan[nFan-1].angle1, weakening_coeff);
        if(std::isnan(grid_results[nFan])) throw std::runtime_error("traction is nan");
    }

    float *highest_grid_pt = std::max_element(grid_results, &grid_results[gridPts]);
    unsigned idx = std::distance(grid_results, highest_grid_pt);

    // reject if the grid max is low
    if(*highest_grid_pt < prms.normal_traction_threshold*prms.cutoff_coefficient) return;

    // sectors
    int sector1, sector2;

    if(isBoundary && (idx == 0 || idx==gridPts-1))
    {
        sector1 = idx == 0 ? 0 : gridPts-2;
        sector2 = -1;
    }
    else
    {
        sector1 = idx;
        sector2 = (idx-1+nFan)%nFan;
    }

    int bits = std::numeric_limits<float>::digits/2;

    boost::uintmax_t max_iter = 15;
    auto [fracture_angle, max1] = boost::math::tools::brent_find_minima(
                    [=](double x){return -normal_traction(x, weakening_coeff);},
        fan[sector1].angle0, fan[sector1].angle1, bits, max_iter);
    max_normal_traction = -max1;

    if(sector2 > -1)
    {
        max_iter = 15;
        auto [fracture_angle2, max2] = boost::math::tools::brent_find_minima(
                        [=](double x){return -normal_traction(x, weakening_coeff);},
            fan[sector2].angle0, fan[sector2].angle1, bits, max_iter);
        max2 = -max2;
        if(max2 > max_normal_traction) fracture_angle = fracture_angle2;
    }

    evaluate_tractions(fracture_angle, result_with_max_traction, weakening_coeff);
    if(result_with_max_traction.faces[0]==result_with_max_traction.faces[1])
        throw std::runtime_error("evaluate_tractions: face0==face1");
    if(!result_with_max_traction.faces[0]->ContainsNode(this))
        throw std::runtime_error("ComputeFanVariablesAlt: mesh topology error 0");
    if(result_with_max_traction.faces[1]!= nullptr && !result_with_max_traction.faces[1]->ContainsNode(this))
        throw std::runtime_error("ComputeFanVariablesAlt: mesh topology error 1");
    max_normal_traction = result_with_max_traction.trac_normal_max;
    dir = result_with_max_traction.tn;

    const float threshold_angle = fan_angle_span*0.1;
    if(isBoundary && (fracture_angle < threshold_angle ||
                      fracture_angle > fan_angle_span-threshold_angle || fan_angle_span < M_PI/2))
    {max_normal_traction=0; return;}
}

double icy::Node::NormalTraction(double angle_fwd, double weakening_coeff) const
{
    SepStressResult tmpSsr;
    EvaluateTractions(angle_fwd, tmpSsr, weakening_coeff);
    return tmpSsr.trac_normal;
}

void icy::Node::EvaluateTractions(double angle_fwd, SepStressResult &ssr, const double weakening_coeff) const
{
    ssr.traction_top[0] = ssr.traction_top[1] = Eigen::Vector2f::Zero();
    ssr.traction_bottom[0] = ssr.traction_bottom[1] = Eigen::Vector2f::Zero();
    ssr.faces[0] = ssr.faces[1] = nullptr;

    if(angle_fwd == fan_angle_span) angle_fwd -= 1e-4;
    ssr.angle_fwd = angle_fwd;

    double angle_bwd = angle_fwd+fan_angle_span/2;
    if (angle_bwd >= fan_angle_span) angle_bwd -= fan_angle_span;
    ssr.angle_bwd = angle_bwd;

    // integrate traction
    int sector = (isBoundary || angle_fwd < angle_bwd) ? 0 : 1;

    std::size_t nFans = fan.size();

    for (std::size_t f=0; f < nFans; f++)
    {
        const Sector &fp = fan[f];

        if (angle_fwd >= fp.angle0 && angle_fwd < fp.angle1)
        {
            ssr.faces[0] = fp.face;
            ssr.e[0] = fp.face->CWEdge(this);
            ssr.e[1] = fp.face->CCWEdge(this);
            ssr.e_opposite[0] = fp.face->OppositeEdge(this);

            float phi = ssr.phi[0] = angle_fwd - fp.angle0;
            ssr.theta[0] = fp.angle1 - angle_fwd;

            float ratio = phi/(fp.angle1-fp.angle0);
            ssr.tn = (fp.u_normalized*(1-ratio) + fp.v_normalized*ratio).normalized();
            ssr.tn_p = (fp.u_p*(1-ratio) + fp.v_p*ratio).normalized(); // perpendicular to tn
            //ssr.tn_p = normal_n.cross(ssr.tn).normalized();
            Eigen::Vector2f tmult_top = fp.face->str_top * ssr.tn_p;
            Eigen::Vector2f tmult_bottom = fp.face->str_bottom * ssr.tn_p;

            ssr.traction_top[sector] += tmult_top - fp.t0_top;
            ssr.traction_bottom[sector] += tmult_bottom - fp.t0_bottom;
            sector = 1-sector;
            ssr.traction_top[sector] += fp.t1_top - tmult_top;
            ssr.traction_bottom[sector] += fp.t1_bottom - tmult_bottom;
        }
        else if (!isBoundary && angle_bwd >= fp.angle0 && angle_bwd < fp.angle1)
        {
            ssr.faces[1] = fp.face;
            ssr.e[2] = fp.face->CWEdge(this);
            ssr.e[3] = fp.face->CCWEdge(this);
            ssr.e_opposite[1] = fp.face->OppositeEdge(this);

            float phi = ssr.phi[1] = angle_bwd - fp.angle0;
            ssr.theta[1] = fp.angle1 - angle_bwd;

            float ratio = phi/(fp.angle1-fp.angle0);
            Eigen::Vector2f tn_p = (fp.u_p*(1-ratio) + fp.v_p*ratio).normalized(); // perpendicular to tn

            Eigen::Vector2f tmult_top = fp.face->str_top * tn_p;
            Eigen::Vector2f tmult_bottom = fp.face->str_bottom * tn_p;

            ssr.traction_top[sector] += tmult_top - fp.t0_top;
            ssr.traction_bottom[sector] += tmult_bottom - fp.t0_bottom;
            sector = 1-sector;
            ssr.traction_top[sector] += fp.t1_top - tmult_top;
            ssr.traction_bottom[sector] += fp.t1_bottom - tmult_bottom;
        }
        else
        {
            ssr.traction_top[sector] += fp.t1_top - fp.t0_top;
            ssr.traction_bottom[sector] += fp.t1_bottom - fp.t0_bottom;
        }
    }   // nFans

    float t0_tangential_top = ssr.traction_top[0].dot(ssr.tn);
    float t1_tangential_top = ssr.traction_top[1].dot(ssr.tn);
    float t0_normal_top = ssr.tn_p.dot(ssr.traction_top[0]);
    float t1_normal_top = -ssr.tn_p.dot(ssr.traction_top[1]);
    ssr.trac_normal_top = t0_normal_top + t1_normal_top;
    ssr.trac_tangential_top = t0_tangential_top - t1_tangential_top;

    float t0_tangential_bottom = ssr.traction_bottom[0].dot(ssr.tn);
    float t1_tangential_bottom = ssr.traction_bottom[1].dot(ssr.tn);
    float t0_normal_bottom = ssr.tn_p.dot(ssr.traction_bottom[0]);
    float t1_normal_bottom = -ssr.tn_p.dot(ssr.traction_bottom[1]);
    ssr.trac_normal_bottom = t0_normal_bottom + t1_normal_bottom;
    ssr.trac_tangential_bottom = t0_tangential_bottom - t1_tangential_bottom;

    if(!isBoundary)
    {
        ssr.trac_normal_bottom /= 2;
        ssr.trac_tangential_bottom /= 2;
        ssr.trac_normal_top /= 2;
        ssr.trac_tangential_top /= 2;
    }

    if(isCrackTip)
    {
        float coeff = ((1-weakening_coeff)+(weakening_coeff)*pow((weakening_direction.dot(ssr.tn)+1)/2, 5));
        ssr.trac_normal_bottom*=coeff;
        ssr.trac_normal_top*=coeff;
    }

    ssr.trac_normal_max = std::max(ssr.trac_normal_top, ssr.trac_normal_bottom);
}


uint64_t icy::Node::make_key(Node *nd0, Node *nd1)
{
    int nd0idx = nd0->locId;
    int nd1idx = nd1->locId;

    if(nd0idx > nd1idx) std::swap(nd0idx, nd1idx);
    return ((uint64_t)nd0idx << 32) | nd1idx;
}
