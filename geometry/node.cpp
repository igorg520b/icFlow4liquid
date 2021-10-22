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
    isBoundary = isCrackTip = false;
    time_loaded_above_threshold = 0;
    adj_elems.clear();
    adj_czs.clear();
}

void icy::Node::Initialize(double x, double y)
{
    x_initial << x,y;
    intended_position = xt = xn = x_initial;
}

void icy::Node::Initialize(const Node *nd)
{
    x_initial = nd->x_initial;
    xt = nd->xt;
    xn = nd->xn;
    vn = nd->vn;
}


void icy::Node::AddSpringEntries(EquationOfMotionSolver &eq, const SimParams &prms, double h, Eigen::Vector2d &spring)
{
    if(eqId<0 || spring_attached<1) return;
    double k = prms.YoungsModulus/300;
    Eigen::Vector2d spr = spring_attachment_position+spring;
    double hsqk = h*h*k;

    Eigen::Vector2d fd = (spr-xt);

//    double const_term = hsqk*fd.dot(fd)/2;
    Eigen::Vector2d linear_term = -hsqk*fd;
    Eigen::Matrix2d quadratic_term = hsqk*Eigen::Matrix2d::Identity();
    eq.AddToEquation(linear_term.data(), quadratic_term.data(), {eqId});
}


void icy::Node::CreateUnrotatedFan()
{
    if(adj_elems.size() == 0) throw std::runtime_error("Node::CreateUnorderedFan: disconnected node");
    fan.clear();
    area = 0;
    for(Element *elem : adj_elems)
    {
        area += elem->area_initial/3;
        Node::Sector s;
        s.face = elem;
        Eigen::Vector2d tcv = elem->getCenter() - x_initial;
        s.centerAngle = atan2(tcv.y(), tcv.x());
        uint8_t thisIdx, CWIdx, CCWIdx;
        elem->getIdxs(this, thisIdx, CWIdx, CCWIdx);
        s.nd[0] = elem->nds[CWIdx];
        s.nd[1] = elem->nds[CCWIdx];
        fan.push_back(s);
    }
    std::sort(fan.begin(), fan.end());  // order by angle of the element, counter-clockwise
}


void icy::Node::PrepareFan()
{
    CreateUnrotatedFan();
    isBoundary = std::any_of(fan.begin(),fan.end(),[this](Sector &s){return s.face->isOnBoundary(this);});
/*
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

        if(elem->isOnBoundary(this)) isBoundary = true;
    }

    std::sort(fan.begin(), fan.end());
*/
    // if boundary, then ensure that sectors start with a boundary element and end with a boundary element
    if(isBoundary)
    {
        // find the fan element with the border on the CW direction
        auto cw_boundary = std::find_if(fan.begin(), fan.end(), [this](const Sector &f){return f.face->isCWBoundary(this);});
        if(cw_boundary == fan.end()) { PrintoutFan(); throw std::runtime_error("cw boundary not found"); }
        std::rotate(fan.begin(), cw_boundary, fan.end());
        if(!fan.back().face->isCCWBoundary(this)) { PrintoutFan(); throw std::runtime_error("PrepareFan(): no CCW boundary"); }

        // assign CWBoundaryNode and CCWBoundaryNode
        this->CWBoundaryNode = fan.front().face->CW_Node(this);
        this->CCWBoundaryNode = fan.back().face->CCW_Node(this);
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
        Node *cw_node, *ccw_node;
        std::tie(cw_node,ccw_node) = f.face->CW_CCW_Node(this);

        f.u_normalized = (cw_node->xt-this->xt).normalized();
        f.v_normalized = (ccw_node->xt-this->xt).normalized();

        f.angle0 = fan_angle_span;
        fan_angle_span += get_angle(f.u_normalized,f.v_normalized);
        f.angle1 = fan_angle_span;

        f.u_p << -f.u_normalized[1], f.u_normalized[0];
        f.v_p << -f.v_normalized[1], f.v_normalized[0];

        f.t0 = f.face->CauchyStress * f.u_p;
        f.t1 = f.face->CauchyStress * f.v_p;
    }
}


void icy::Node::PrintoutFan()
{
    spdlog::info("Printing fan for node {}; isCrackTip: {}; isBoundary: {}", locId, isCrackTip, isBoundary);
    spdlog::info("fan.size {}; adj_elems.size {}", fan.size(), adj_elems.size());
    spdlog::info("fan_angle_span: {}", fan_angle_span);

    spdlog::info("┌ {0: ^9} ┬ {1: ^14} ┬ {2: ^14} ┬ {3: ^6} ┬ {4: ^10} ┬ {5: ^8} ┬ {6: ^8}",
                 "nd1-nd2", "nds 0-1-2", "angle0-angle1", "cAngle", "area", "CWB", "CCWB");

    for(Sector &s : fan)
        spdlog::info("│ {0:>4}-{1:<4} │ {2: >4}-{3: >4}-{4: <4} │ {5:>6.4f}-{6:<6.4f} │ {7: ^6.4f} | {8:^6.4e} | {9:^8} | {10:^8}",
                     s.nd[0]->locId, s.nd[1]->locId,
                     s.face->nds[0]->locId,s.face->nds[1]->locId,s.face->nds[2]->locId,
                s.angle0, s.angle1, s.centerAngle, s.face->area_initial, s.face->isCWBoundary(this), s.face->isCCWBoundary(this));
}

void icy::Node::ComputeFanVariables(const SimParams &prms)
{
    double dont_split_nearly_degenerate_elems = prms.FractureAngleThreshold*M_PI/180;

    if(fan.size()==0) throw std::runtime_error("invoking ComputeFanVariables on a Node without elements");
    dir = Eigen::Vector2d::Zero();
    max_normal_traction = 0;
    UpdateFan();
    unsigned nFan = fan.size();
    if(nFan==1 || fan_angle_span < dont_split_nearly_degenerate_elems) return;

    double weakening_coeff = prms.FractureWeakeningCoeff;
    unsigned gridPts = isBoundary ? nFan+1 : nFan;

    double grid_results[gridPts];
    for(unsigned i=0; i<nFan; i++)
    {
        grid_results[i] = NormalTraction(fan[i].angle0, weakening_coeff);
        if(std::isnan(grid_results[i])) throw std::runtime_error("ComputeFanVariables: traction is NaN");
    }
    if(isBoundary)
    {
        grid_results[nFan] = NormalTraction(fan[nFan-1].angle1, weakening_coeff);
        if(std::isnan(grid_results[nFan])) throw std::runtime_error("ComputeFanVariables: traction is NaN");
    }

    double *highest_grid_pt = std::max_element(grid_results, &grid_results[gridPts]);
    unsigned idx = std::distance(grid_results, highest_grid_pt);

    // reject if the grid max is low
    if(*highest_grid_pt < prms.FractureTractionThreshold*0.4) return;

    // sectors
    int sector1, sector2;

    if(isBoundary && (idx==0 || idx==gridPts-1))
    {
        sector1 = idx==0 ? 0 : gridPts-2;
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
                    [&](double x){return -NormalTraction(x, weakening_coeff);},
        fan[sector1].angle0, fan[sector1].angle1, bits, max_iter);
    max_normal_traction = -max1;

    if(sector2 > -1)
    {
        max_iter = 15;
        auto [fracture_angle2, max2] = boost::math::tools::brent_find_minima(
                        [&](double x){return -NormalTraction(x, weakening_coeff);},
            fan[sector2].angle0, fan[sector2].angle1, bits, max_iter);
        max2 = -max2;
        if(max2 > max_normal_traction) fracture_angle = fracture_angle2;
    }

    EvaluateTractions(fracture_angle, result_with_max_traction, weakening_coeff);

    if(result_with_max_traction.faces[0]==result_with_max_traction.faces[1] || result_with_max_traction.faces[0]==nullptr)
    {
        spdlog::critical("evaluate_tractions: face0=={}; face1=={}",
                         (void*)result_with_max_traction.faces[0],
                (void*)result_with_max_traction.faces[1]);
        spdlog::critical("fracture_angle: {}; fan_angle_span {}",fracture_angle, fan_angle_span);
        PrintoutFan();
        throw std::runtime_error("evaluate_tractions: face0==face1");
    }

    if(!result_with_max_traction.faces[0]->containsNode(this))
    {
        spdlog::critical("ComputeFanVariables: mesh topology error 0");
        throw std::runtime_error("ComputeFanVariables: mesh topology error 0");
    }

    if(result_with_max_traction.faces[1] != nullptr && !result_with_max_traction.faces[1]->containsNode(this))
    {
        spdlog::critical("ComputeFanVariables: mesh topology error 1");
        throw std::runtime_error("ComputeFanVariables: mesh topology error 1");
    }

    max_normal_traction = result_with_max_traction.trac_normal;
    dir = result_with_max_traction.tn;

    // don't break nearly-degenerate sectors
    double span0 = result_with_max_traction.sectorSpan(0);
    if(result_with_max_traction.phi[0] < dont_split_nearly_degenerate_elems ||
            (result_with_max_traction.faces[0]->area_initial < prms.FractureAreaThreshold &&
             result_with_max_traction.phi[0] < result_with_max_traction.theta[0]))
    {
            result_with_max_traction.phi[0] = 0;
            result_with_max_traction.theta[0] = span0;
            fracture_angle = result_with_max_traction.angle0[0];
    }
    else if(result_with_max_traction.theta[0] < dont_split_nearly_degenerate_elems ||
            (result_with_max_traction.faces[0]->area_initial < prms.FractureAreaThreshold &&
             result_with_max_traction.phi[0] > result_with_max_traction.theta[0]))
    {
        result_with_max_traction.theta[0] = 0;
        result_with_max_traction.phi[0] = span0;
        fracture_angle = result_with_max_traction.angle1[0];
    }

    const double threshold_angle = dont_split_nearly_degenerate_elems;
    if(isBoundary && (fracture_angle < threshold_angle ||
                      fracture_angle > fan_angle_span-threshold_angle))
    {max_normal_traction=0; return;}

    if(!isBoundary)
    {
        double fracture_angle_bwd;
        // similar snap on the other side
        // TODO: there is a narrow case that may result in invalid topology
        double span1 = result_with_max_traction.sectorSpan(1);
        if(result_with_max_traction.phi[1] < dont_split_nearly_degenerate_elems ||
                (result_with_max_traction.faces[1]->area_initial < prms.FractureAreaThreshold &&
                 result_with_max_traction.phi[1] < result_with_max_traction.theta[1]))
        {
                result_with_max_traction.phi[1] = 0;
                result_with_max_traction.theta[1] = span1;
                fracture_angle_bwd = result_with_max_traction.angle0[1];
                if(fracture_angle_bwd == fracture_angle) {max_normal_traction=0; return;}
        }
        else if(result_with_max_traction.theta[1] < dont_split_nearly_degenerate_elems ||
                (result_with_max_traction.faces[1]->area_initial < prms.FractureAreaThreshold &&
                 result_with_max_traction.phi[1] > result_with_max_traction.theta[1]))
        {
            result_with_max_traction.theta[1] = 0;
            result_with_max_traction.phi[1] = span1;
            fracture_angle_bwd = result_with_max_traction.angle1[1];
            if(fracture_angle_bwd == fracture_angle) {max_normal_traction=0; return;}
        }
    }
}

double icy::Node::NormalTraction(double angle_fwd, double weakening_coeff) const
{
    SepStressResult tmpSsr;
    EvaluateTractions(angle_fwd, tmpSsr, weakening_coeff);
    return tmpSsr.trac_normal;
}

void icy::Node::EvaluateTractions(double angle_fwd, SepStressResult &ssr, const double weakening_coeff) const
{
    ssr.traction[0] = ssr.traction[1] = Eigen::Vector2d::Zero();
    ssr.faces[0] = ssr.faces[1] = nullptr;

    if(angle_fwd == fan_angle_span) angle_fwd = std::max(0.0,angle_fwd-1e-11);
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
//            ssr.e[0] = fp.face->CWEdge(this);
//            ssr.e[1] = fp.face->CCWEdge(this);
//            ssr.e_opposite[0] = fp.face->OppositeEdge(this);

            double phi = ssr.phi[0] = angle_fwd - fp.angle0;
            ssr.theta[0] = fp.angle1 - angle_fwd;
            ssr.angle0[0] = fp.angle0;
            ssr.angle1[0] = fp.angle1;

            double ratio = phi/(fp.angle1-fp.angle0);
            ssr.tn = (fp.u_normalized*(1-ratio) + fp.v_normalized*ratio).normalized();
            ssr.tn_p = (fp.u_p*(1-ratio) + fp.v_p*ratio).normalized(); // perpendicular to tn
            Eigen::Vector2d tmult = fp.face->CauchyStress * ssr.tn_p;

            ssr.traction[sector] += tmult - fp.t0;
            sector = 1-sector;
            ssr.traction[sector] += fp.t1 - tmult;
        }
        else if (!isBoundary && angle_bwd >= fp.angle0 && angle_bwd < fp.angle1)
        {
            ssr.faces[1] = fp.face;
//            ssr.e[2] = fp.face->CWEdge(this);
//            ssr.e[3] = fp.face->CCWEdge(this);
//            ssr.e_opposite[1] = fp.face->OppositeEdge(this);

            float phi = ssr.phi[1] = angle_bwd - fp.angle0;
            ssr.theta[1] = fp.angle1 - angle_bwd;
            ssr.angle0[1] = fp.angle0;
            ssr.angle1[1] = fp.angle1;

            double ratio = phi/(fp.angle1-fp.angle0);
            Eigen::Vector2d tn_p = (fp.u_p*(1-ratio) + fp.v_p*ratio).normalized(); // perpendicular to tn

            Eigen::Vector2d tmult = fp.face->CauchyStress * tn_p;

            ssr.traction[sector] += tmult - fp.t0;
            sector = 1-sector;
            ssr.traction[sector] += fp.t1 - tmult;
        }
        else
        {
            ssr.traction[sector] += fp.t1 - fp.t0;
        }
    }   // nFans

    double t0_tangential = ssr.traction[0].dot(ssr.tn);
    double t1_tangential = ssr.traction[1].dot(ssr.tn);
    double t0_normal = ssr.tn_p.dot(ssr.traction[0]);
    double t1_normal = -ssr.tn_p.dot(ssr.traction[1]);
    ssr.trac_normal = t0_normal + t1_normal;
    ssr.trac_tangential = t0_tangential - t1_tangential;

    if(!isBoundary)
    {
        ssr.trac_normal /= 2;
        ssr.trac_tangential /= 2;
    }

    if(isCrackTip)
    {
        double coeff = ((1-weakening_coeff)+(weakening_coeff)*pow((weakening_direction.dot(ssr.tn)+1)/2, 5));
        ssr.trac_normal*=coeff;
    }
}

void icy::Node::InitializeLERP(const Node *nd0, const Node *nd1, double f)
{
    x_initial = nd0->x_initial*f + nd1->x_initial*(1-f);
    xt = nd0->xt*f + nd1->xt*(1-f);
    xn = nd0->xn*f + nd1->xn*(1-f);
    vn = nd0->vn*f + nd1->vn*(1-f);
}


uint64_t icy::Node::make_local_key(Node *nd0, Node *nd1)
{
    int nd0idx = nd0->locId;
    int nd1idx = nd1->locId;
    if(nd0idx > nd1idx) std::swap(nd0idx, nd1idx);
    return ((uint64_t)nd0idx << 32) | nd1idx;
}


uint64_t icy::Node::make_global_key(Node *nd0, Node *nd1)
{
    int nd0idx = nd0->globId;
    int nd1idx = nd1->globId;
    if(nd0idx > nd1idx) std::swap(nd0idx, nd1idx);
    return ((uint64_t)nd0idx << 32) | nd1idx;
}

void icy::Node::ReplaceAdjacentElement(Element *originalElem, Element *replacement)
{
    auto iter = std::find(adj_elems.begin(),adj_elems.end(),originalElem);
    if(iter == adj_elems.end())
        throw std::runtime_error("icy::Node::ReplaceAdjacentElement: can't find the original element to replace");
    else *iter = replacement;
}

