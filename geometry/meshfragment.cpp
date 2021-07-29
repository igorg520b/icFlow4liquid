#include <unordered_map>
#include <unordered_set>
//#include <hdf5.h>
#include <H5Cpp.h>
#include <gmsh.h>
#include "meshfragment.h"

icy::ConcurrentPool<icy::BVHN> icy::MeshFragment::BVHNLeafFactory(50000);
icy::ConcurrentPool<icy::Node> icy::MeshFragment::NodeFactory(50000);
icy::ConcurrentPool<icy::Element> icy::MeshFragment::ElementFactory(50000);


icy::MeshFragment::~MeshFragment()
{
    // release objects back to their pools
    BVHNLeafFactory.release(leaves_for_ccd);
    BVHNLeafFactory.release(leaves_for_contact);
    NodeFactory.release(nodes);
    NodeFactory.release(nodes_tmp);
    ElementFactory.release(elems);
    ElementFactory.release(elems_tmp);
}

void icy::MeshFragment::GenerateSpecialBrick(double ElementSize)
{
    deformable = true;

    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::model::add("block1");

    double width = 2;
    double height = 1;

    int point1 = gmsh::model::occ::addPoint(-width/2, 1e-10, 0, 1.0);
    int point2 = gmsh::model::occ::addPoint(-width/2, height, 0, 1.0);
    int point3 = gmsh::model::occ::addPoint(width/2, height*1.0, 0, 1.0);
    int point4 = gmsh::model::occ::addPoint(width/2, 1e-10, 0, 1.0);

    int line1 = gmsh::model::occ::addLine(point1, point2);
    int line2 = gmsh::model::occ::addLine(point2, point3);
    int line3 = gmsh::model::occ::addLine(point3, point4);
    int line4 = gmsh::model::occ::addLine(point4, point1);

    int loopTag1 = gmsh::model::occ::addCurveLoop({line1, line2, line3, line4});

    double offset = 1.2*ElementSize;
    int point11 = gmsh::model::occ::addPoint(-width/2+offset, offset, 0, 1.0);
    int point12 = gmsh::model::occ::addPoint(-width/2+offset, height-offset, 0, 1.0);
    int point13 = gmsh::model::occ::addPoint(width/2-offset, height-offset, 0, 1.0);
    int point14 = gmsh::model::occ::addPoint(width/2-offset, offset, 0, 1.0);

    int line11 = gmsh::model::occ::addLine(point11, point12);
    int line12 = gmsh::model::occ::addLine(point12, point13);
    int line13 = gmsh::model::occ::addLine(point13, point14);
    int line14 = gmsh::model::occ::addLine(point14, point11);

    int loopTag2 = gmsh::model::occ::addCurveLoop({line11, line12, line13, line14});

    double radius = 0.15;
    int ellipseTag = gmsh::model::occ::addEllipse(0, height/2, 0, radius, radius*0.999);
    int loopTag3 = gmsh::model::occ::addCurveLoop({ellipseTag});

    int surfaceTag1 = gmsh::model::occ::addPlaneSurface({loopTag1, loopTag2});
    int surfaceTag2 = gmsh::model::occ::addPlaneSurface({loopTag2, loopTag3});
    int surfaceTag3 = gmsh::model::occ::addPlaneSurface({loopTag3});


    gmsh::model::occ::synchronize();

    // physical groups
    unsigned nGroups = 5;
    int groups[nGroups];
    groups[0] = gmsh::model::addPhysicalGroup(1, {line1, line2, line3, line4});
    groups[1] = gmsh::model::addPhysicalGroup(1, {line11, line12, line13, line14});
    groups[2] = gmsh::model::addPhysicalGroup(2, {surfaceTag1});
    groups[3] = gmsh::model::addPhysicalGroup(2, {surfaceTag2,surfaceTag3});
    groups[4] = gmsh::model::addPhysicalGroup(2, {surfaceTag1,surfaceTag2,surfaceTag3});   // everything

    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);
    gmsh::model::mesh::generate(2);



    // GET NODES
//    freeNodeCount = 0;
    std::unordered_map<std::size_t, int> mtags; // nodeTag -> sequential position in nodes[]
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords;

    // nodes of the "inner boundary" are placed first
    gmsh::model::mesh::getNodesForPhysicalGroup(1, groups[1], nodeTags, nodeCoords);
    nInnerBoundaryNodes = nodeTags.size();
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)) continue;  // node was already added
        int sequential_id = nodes.size();
        mtags[tag] = sequential_id;
        icy::Node* nd = NodeFactory.take();
        nd->Reset(sequential_id, nodeCoords[i*3+0], nodeCoords[i*3+1]);
        nodes.push_back(nd);
    }

    // place the rest of the nodes of the "first group", i.e. elastic layer
    nodeTags.clear();
    nodeCoords.clear();
    gmsh::model::mesh::getNodesForPhysicalGroup(2, groups[2], nodeTags, nodeCoords);
    nFirstGroupNodes = nodeTags.size();
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)) continue;  // node was already added
        int sequential_id = nodes.size();
        mtags[tag] = sequential_id;
        icy::Node* nd = NodeFactory.take();
        nd->Reset(sequential_id, nodeCoords[i*3+0], nodeCoords[i*3+1]);
        nodes.push_back(nd);
    }
    std::cout << "\nStep1 nodes added:" << nodes.size() << std::endl;

    // place the rest of the nodes
    nodeTags.clear();
    nodeCoords.clear();
    gmsh::model::mesh::getNodesForPhysicalGroup(2, groups[3], nodeTags, nodeCoords);
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)) continue;  // node was already added

        int sequential_id = nodes.size();
        mtags[tag] = sequential_id;
        icy::Node* nd = NodeFactory.take();
        nd->Reset(sequential_id, nodeCoords[i*3+0], nodeCoords[i*3+1]);
        nodes.push_back(nd);
    }


    // mark the "outer" and "inner" boundaries
    for(int i=0;i<2;i++)
    {
        gmsh::model::mesh::getNodesForPhysicalGroup(1, groups[i], nodeTags, nodeCoords);
        for(const std::size_t tag : nodeTags) nodes[mtags.at(tag)]->group.set(i);
    }

    std::cout << "special brick: pricessing elements" << std::endl;
    // GET ELEMENTS
    if(elems.size()!=0) throw std::runtime_error("elems.size != 0");
    std::vector<int> surfaces = {surfaceTag1, surfaceTag2, surfaceTag3};

    for(int k=0;k<3;k++)
    {
        std::vector<std::size_t> trisTags, nodeTagsInTris;
        gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris, surfaces[k]);

        for(std::size_t i=0;i<trisTags.size();i++)
        {
            icy::Element *elem = ElementFactory.take();
            elem->Reset();
            for(int j=0;j<3;j++) elem->nds[j] = nodes[mtags.at(nodeTagsInTris[i*3+j])];
            elem->PrecomputeInitialArea();
            if(elem->area_initial < 0)
            {
                std::swap(elem->nds[0], elem->nds[1]);
                elem->PrecomputeInitialArea();
                if(elem->area_initial < 0) throw std::runtime_error("icy::Mesh::Reset - error");
            }
            for(int j=0;j<3;j++) elem->nds[j]->area += elem->area_initial/3;
            elem->group=k;
//            if(k==2) elem->PiMultiplier << 1, 0, 0, 0.5;
            elems.push_back(elem);

        }
        if(k==0) nFirstGroupElems = elems.size();
    }



    // BOUNDARIRES - EDGES

    std::vector<std::size_t> edgeTags, nodeTagsInEdges;
    gmsh::model::mesh::getElementsByType(1, edgeTags, nodeTagsInEdges);

    boundary_edges.clear();
    inner_boundary_edges.clear();
    for(std::size_t i=0;i<nodeTagsInEdges.size()/2;i++)
    {
        int idx1 = mtags.at(nodeTagsInEdges[i*2+0]);
        int idx2 = mtags.at(nodeTagsInEdges[i*2+1]);
        Node *nd1, *nd2;
        if(idx1<idx2)
        {
            nd1 = nodes[idx1];
            nd2 = nodes[idx2];
        }
        else
        {
            nd1 = nodes[idx2];
            nd2 = nodes[idx1];
        }
        if(nd1->group.test(0) && nd2->group.test(0)) boundary_edges.emplace_back(nd1,nd2);
        else if(nd1->group.test(1) && nd2->group.test(1)) inner_boundary_edges.emplace_back(idx1,idx2);
    }

    gmsh::clear();

    std::cout << "\nSpecial Brick:\nTotal elems: " << elems.size();
    std::cout << "\nFirst-group elems: " << nFirstGroupElems;
    std::cout << "\nTotal nodes: " << nodes.size();
    std::cout << "\nFirst-group nodes: " << nFirstGroupNodes;
    std::cout << "\nInner boundary: " << inner_boundary_edges.size() << std::endl;
    std::cout << "brick generated" << std::endl;
}

void icy::MeshFragment::GenerateBrick(double ElementSize)
{
    deformable = true;

    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    double width = 2;
    double height = 1;
    int point1 = gmsh::model::occ::addPoint(-width/2, 1e-10, 0, 1.0);
    int point2 = gmsh::model::occ::addPoint(-width/2, height, 0, 1.0);
    int point3 = gmsh::model::occ::addPoint(width/2, height*1.1, 0, 1.0);
    int point4 = gmsh::model::occ::addPoint(width/2, 1e-10, 0, 1.0);

    int line1 = gmsh::model::occ::addLine(point1, point2);
    int line2 = gmsh::model::occ::addLine(point2, point3);
    int line3 = gmsh::model::occ::addLine(point3, point4);
    int line4 = gmsh::model::occ::addLine(point4, point1);

    std::vector<int> curveTags;
    curveTags.push_back(line1);
    curveTags.push_back(line2);
    curveTags.push_back(line3);
    curveTags.push_back(line4);
    int loopTag = gmsh::model::occ::addCurveLoop(curveTags);

    std::vector<int> loops;
    loops.push_back(loopTag);
    gmsh::model::occ::addPlaneSurface(loops);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);

    GetFromGmsh();
}

void icy::MeshFragment::GenerateContainer(double ElementSize, double offset)
{
    std::cout << "\nGenerateContainer" << std::endl;;

    deformable = false;

    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::clear();
    gmsh::model::add("block1");

    double width = 2;
    double height = 1;
    int point1 = gmsh::model::occ::addPoint(-width/2-offset, height*1.5, 0, 1.0);
    int point2 = gmsh::model::occ::addPoint(-width/2-offset, -offset, 0, 1.0);
    int point3 = gmsh::model::occ::addPoint(width/2+offset, -offset, 0, 1.0);
    int point4 = gmsh::model::occ::addPoint(width/2+offset, height*1.5, 0, 1.0);

    int line1 = gmsh::model::occ::addLine(point1, point2);
    int line2 = gmsh::model::occ::addLine(point2, point3);
    int line3 = gmsh::model::occ::addLine(point3, point4);

    std::vector<int> curveTags;
    curveTags.push_back(line1);
    curveTags.push_back(line2);
    curveTags.push_back(line3);
    gmsh::model::occ::addWire(curveTags);

    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);

    GetFromGmsh();
}

void icy::MeshFragment::GenerateIndenter(double ElementSize)
{
    std::cout << "\nGenerateIndenter\n";

    deformable = false;

    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::clear();
    gmsh::model::add("block1");

    double height = 1;
    double radius = 0.15;
//    int point1 = gmsh::model::occ::addPoint(0, height+radius*1.1, 0, 1.0);

    gmsh::model::occ::addEllipse(0, height+radius*1.1, 0, radius, radius/2);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);
    GetFromGmsh();

}

void icy::MeshFragment::GenerateBall(double x, double y, double r1, double r2, double ElementSize)
{
    deformable = true;

    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    int diskTag = gmsh::model::occ::addDisk(x, y, 0, r1, r2);
    gmsh::vectorpair vp;
    vp.push_back(std::make_pair(2,diskTag));
    gmsh::model::occ::rotate(vp, x,y,0,0,0,1,M_PI/2);
    gmsh::model::occ::synchronize();
    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);
    GetFromGmsh();
}


void icy::MeshFragment::GetFromGmsh()
{
    if(elems.size() || nodes.size() || boundary_edges.size()) throw std::runtime_error("GetFromGmsh(): elem/node/edge arrays not empty");
//    freeNodeCount = 0;

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, int> mtags; // nodeTag -> sequential position in nodes[]

    if(deformable)
    {
        gmsh::model::mesh::generate(2);

        // GET NODES
        gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);

        // set the size of the resulting nodes array
        for(unsigned i=0;i<nodeTags.size();i++)
        {
            std::size_t tag = nodeTags[i];
            if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");
            unsigned sequential_idx = nodes.size();
            mtags[tag] = sequential_idx;
            double x = nodeCoords[i*3+0];
            double y = nodeCoords[i*3+1];

            icy::Node* nd = NodeFactory.take();
            nodes.push_back(nd);
            nd->Reset(sequential_idx, x, y);
            // if(y==0 || !deformable) nd->pinned=true;
        }

        // GET ELEMENTS
        std::vector<std::size_t> trisTags, nodeTagsInTris;
        gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

        for(std::size_t i=0;i<trisTags.size();i++)
        {
            icy::Element *elem = ElementFactory.take();
            elem->Reset();
            for(int j=0;j<3;j++) elem->nds[j] = nodes[mtags.at(nodeTagsInTris[i*3+j])];
            elem->PrecomputeInitialArea();
            if(elem->area_initial < 0)
            {
                std::swap(elem->nds[0], elem->nds[1]);
                elem->PrecomputeInitialArea();
                if(elem->area_initial < 0) throw std::runtime_error("icy::Mesh::Reset - error");
            }
            for(int j=0;j<3;j++) elem->nds[j]->area += elem->area_initial/3;
            elems.push_back(elem);
        }

    }
    else
    {
        gmsh::model::mesh::generate(1);

        // GET NODES
        gmsh::model::mesh::getNodesByElementType(1, nodeTags, nodeCoords, parametricCoords);

        // set the size of the resulting nodes array
        for(unsigned i=0;i<nodeTags.size();i++)
        {
            std::size_t tag = nodeTags[i];
            if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in boundary");
            unsigned sequential_id = nodes.size();
            mtags[tag] = sequential_id;
            double x = nodeCoords[i*3+0];
            double y = nodeCoords[i*3+1];

            icy::Node* nd = NodeFactory.take();
            nd->Reset(sequential_id, nodeCoords[i*3+0], nodeCoords[i*3+1]);
            nodes.push_back(nd);
            nd->pinned=true;
        }


        // BOUNDARIRES - EDGES

        std::vector<std::size_t> edgeTags, nodeTagsInEdges;
        gmsh::model::mesh::getElementsByType(1, edgeTags, nodeTagsInEdges);

        for(std::size_t i=0;i<edgeTags.size();i++)
        {
            int idx1 = mtags.at(nodeTagsInEdges[i*2+0]);
            int idx2 = mtags.at(nodeTagsInEdges[i*2+1]);
            if(idx1<idx2)
                boundary_edges.emplace_back(nodes[idx1],nodes[idx2]);
            else
                boundary_edges.emplace_back(nodes[idx2],nodes[idx1]);
        }
    }

    gmsh::clear();
}

void icy::MeshFragment::GenerateLeafs(unsigned edge_idx)
{
    root_ccd.isLeaf = root_contact.isLeaf = false;
    root_contact.test_self_collision = root_ccd.test_self_collision = false; // this->deformable;

    BVHNLeafFactory.release(leaves_for_ccd);
    BVHNLeafFactory.release(leaves_for_contact);

    leaves_for_ccd.reserve(boundary_edges.size());
    leaves_for_contact.reserve(boundary_edges.size());

    for(auto edge : boundary_edges)
    {
        Node *nd1 = edge.first;
        Node *nd2 = edge.second;
        BVHN *leaf_ccd = BVHNLeafFactory.take();
        BVHN *leaf_contact = BVHNLeafFactory.take();

        leaf_contact->feature_idx = leaf_ccd->feature_idx = edge_idx++;
        leaf_contact->test_self_collision = leaf_ccd->test_self_collision = false;
        leaf_contact->isLeaf = leaf_ccd->isLeaf = true;

        leaves_for_ccd.push_back(leaf_ccd);
        leaves_for_contact.push_back(leaf_contact);
    }
    qDebug() << "icy::MeshFragment::GenerateLeafs() " << leaves_for_ccd.size() << " " << leaves_for_contact.size();
}



void icy::MeshFragment::RemeshSpecialBrick(double ElementSize)
{
    std::cout << "\nicy::MeshFragment::RemeshSpecialBrick(double ElementSize)" << std::endl;
    NodeFactory.release(nodes_tmp);
    ElementFactory.release(elems_tmp);

    for(unsigned i=0;i<nFirstGroupNodes;i++)
    {
        Node *nd = NodeFactory.take();
        nd->Reset(nodes[i]);
        nodes_tmp.push_back(nd);
    }
    std::cout << "first group nodes count " << nodes_tmp.size() << std::endl;

    for(unsigned i=0;i<nFirstGroupElems;i++)
    {
        Element *elem = ElementFactory.take();
        elem->Reset();
        elem->group=1;
        for(int j=0;j<3;j++) elem->nds[j]=nodes_tmp[elems[i]->nds[j]->locId];
        elem->PrecomputeInitialArea();
        for(int j=0;j<3;j++) elem->nds[j]->area += elem->area_initial/3;
        elems_tmp.push_back(elem);
    }

    // re-create the boundary in gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("submesh1");

    std::unordered_map<int, unsigned> btags; // gmsh tag -> index in nodes_tmp
    for(unsigned i=0;i<nInnerBoundaryNodes;i++)
    {
        int tag = gmsh::model::occ::addPoint(nodes_tmp[i]->x_initial.x(), nodes_tmp[i]->x_initial.y(), 0, 0, i+1);
        if(tag != i+1) throw std::runtime_error("tag != i+1");
        btags[tag] = i;
    }

    std::vector<int> ltags; // tags of boundary lines
    for(auto edge : inner_boundary_edges)
    {
        if(edge.first >= nInnerBoundaryNodes || edge.second >= nInnerBoundaryNodes) throw std::out_of_range("inner boudnary");
        int lineTag = gmsh::model::occ::addLine(edge.first+1,edge.second+1);
        ltags.push_back(lineTag);
    }

    int innerLoopTag = gmsh::model::occ::addCurveLoop(ltags);
    int surfaceTag = gmsh::model::occ::addPlaneSurface({innerLoopTag});
    gmsh::model::occ::synchronize();

    for(int lineTag : ltags) gmsh::model::mesh::setTransfiniteCurve(lineTag,0);


    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);
    gmsh::model::mesh::generate(2);

    // get from GMSH -> nodes_tmp, elems_tmp

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, unsigned> mtags; // nodeTag -> sequential position in nodes[]

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);
    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        double x = nodeCoords[i*3+0];
        double y = nodeCoords[i*3+1];

        if(btags.count(tag))
        {
            // node already exists in nodes_tmp
            unsigned idx = btags.at(tag);
            icy::Node* nd = nodes_tmp[btags.at(tag)];
            mtags[tag] = idx;
            if(nd->x_initial.x()!=x || nd->x_initial.y()!=y) throw std::runtime_error("some tags don't match after remeshing");
        }
        else
        {
            if(mtags.count(tag)>0) continue;
            icy::Node* nd = NodeFactory.take();
            unsigned idx = nodes_tmp.size();
            nd->Reset(idx, x, y);
            mtags[tag] = idx;
            nodes_tmp.push_back(nd);
        }

    }
//    for(Node *nd : nodes_tmp) nd->area = 0;
    std::cout << "\n node count after remeshing " << nodes_tmp.size() << std::endl;
    std::cout << "\n nodeTags " << nodeTags.size() << std::endl;

    // GET ELEMENTS
    std::vector<std::size_t> trisTags, nodeTagsInTris;
    gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

    for(std::size_t i=0;i<trisTags.size();i++)
    {

        icy::Element *elem = ElementFactory.take();
        elem->Reset();
        elem->group = 2;
        for(int j=0;j<3;j++)
        {
            unsigned idx0 = nodeTagsInTris[i*3+j];
            unsigned idx = mtags.at(idx0);
            if(idx>=nodes_tmp.size()) throw std::runtime_error("index error");
            elem->nds[j] = nodes_tmp[idx];
        }
        elem->PrecomputeInitialArea();
        if(elem->area_initial < 0)
        {
            std::swap(elem->nds[0], elem->nds[1]);
            elem->PrecomputeInitialArea();
            if(elem->area_initial <= 0) throw std::runtime_error("icy::Mesh::Reset - error");
        }
        for(int j=0;j<3;j++) elem->nds[j]->area += elem->area_initial/3;
        elems_tmp.push_back(elem);
    }

    std::cout << "new elem count " << elems_tmp.size() << std::endl;

    boundary_edges_tmp.clear();
    for(unsigned i=0;i<boundary_edges.size();i++)
        boundary_edges_tmp.emplace_back(nodes_tmp[boundary_edges[i].first->locId],nodes_tmp[boundary_edges[i].second->locId]);

    //gmsh::fltk::run();
    gmsh::clear();
    Swap();

    // double-check

    std::unordered_set<Node*> connected_nds;
    for(Element *e : elems)
        for(unsigned j=0;j<3;j++) connected_nds.insert(e->nds[j]);

    if(connected_nds.size() != nodes.size()) throw std::runtime_error("remeshing error");
}

void icy::MeshFragment::Swap()
{
    std::cout << "SWAP" << std::endl;
    std::vector<icy::Node*> n_tmp(nodes);
    std::vector<icy::Element*> e_tmp(elems);

    nodes.resize(nodes_tmp.size());
    elems.resize(elems_tmp.size());
    std::copy(nodes_tmp.begin(),nodes_tmp.end(),nodes.begin());
    std::copy(elems_tmp.begin(),elems_tmp.end(),elems.begin());
    nodes_tmp.resize(n_tmp.size());
    elems_tmp.resize(e_tmp.size());
    std::copy(n_tmp.begin(),n_tmp.end(),nodes_tmp.begin());
    std::copy(e_tmp.begin(),e_tmp.end(),elems_tmp.begin());

    std::vector<std::pair<icy::Node*,icy::Node*>> be_tmp(boundary_edges);
    boundary_edges.resize(boundary_edges_tmp.size());
    std::copy(boundary_edges_tmp.begin(),boundary_edges_tmp.end(),boundary_edges.begin());
    boundary_edges_tmp.resize(be_tmp.size());
    std::copy(be_tmp.begin(),be_tmp.end(),boundary_edges_tmp.begin());

}

void icy::MeshFragment::SaveFragment(std::string fileName)
{
    constexpr unsigned ElemDataFields = 4;
    H5::IntType datatype_int(H5::PredType::NATIVE_INT);
    H5::FloatType datatype_double(H5::PredType::NATIVE_DOUBLE);

    H5::H5File file(fileName, H5F_ACC_TRUNC);

    // NEW MESH (no displacements, velocity, elem data)
    double *nodes_buffer = new double[nodes.size()*2];
    for(unsigned i=0;i<nodes.size();i++)
    {
        Node* nd = nodes[i];
        nodes_buffer[i*2+0] = nd->x_initial.x();
        nodes_buffer[i*2+1] = nd->x_initial.y();
    }
    hsize_t dimsf_nodes_new[2] = {nodes.size(), 2};
    H5::DataSpace dataspace_nodes_new(2, dimsf_nodes_new);
    H5::DataSet dataset_nodes_new = file.createDataSet("Nodes_New", datatype_double, dataspace_nodes_new);
    dataset_nodes_new.write(nodes_buffer, H5::PredType::NATIVE_DOUBLE);

    int *elems_buffer_nodes = new int[elems.size()*3];
    for(unsigned i=0;i<elems.size();i++)
    {
        Element* e = elems[i];
        int idx0 = e->nds[0]->locId;
        int idx1 = e->nds[1]->locId;
        int idx2 = e->nds[2]->locId;
        if(idx0>=nodes.size() || idx1>=nodes.size() || idx2>=nodes.size()) throw std::out_of_range("nodal index error");
        elems_buffer_nodes[i*3+0] = e->nds[0]->locId;
        elems_buffer_nodes[i*3+1] = e->nds[1]->locId;
        elems_buffer_nodes[i*3+2] = e->nds[2]->locId;
    }

    hsize_t dimsf_elems_nodes_new[2] = {elems.size(), 3};
    H5::DataSpace dataspace_elems_nodes_new(2, dimsf_elems_nodes_new);
    H5::DataSet dataset_elems_nodes = file.createDataSet("Elements_New", datatype_int, dataspace_elems_nodes_new);

    dataset_elems_nodes.write(elems_buffer_nodes, H5::PredType::NATIVE_INT);

    delete[] nodes_buffer;
    delete[] elems_buffer_nodes;



    // OLD MESH (includes nodal displacements, velocities, and per-element data)
    nodes_buffer = new double[nodes_tmp.size()*6];
    for(unsigned i=0;i<nodes_tmp.size();i++)
    {
        Node* nd = nodes_tmp[i];
        if(nd->x_initial.x()==0 && nd->x_initial.y()==0) std::cout << "node " << i << " is zero\n";
        nodes_buffer[i*6+0] = nd->x_initial.x();
        nodes_buffer[i*6+1] = nd->x_initial.y();
        nodes_buffer[i*6+2] = nd->xn.x();
        nodes_buffer[i*6+3] = nd->xn.y();
        nodes_buffer[i*6+4] = nd->vn.x();
        nodes_buffer[i*6+5] = nd->vn.y();
    }
    hsize_t dimsf_nodes_old[2] = {nodes_tmp.size(), 6};
    H5::DataSpace dataspace_nodes_old(2, dimsf_nodes_old);
    H5::DataSet dataset_old = file.createDataSet("Nodes_Old", datatype_double, dataspace_nodes_old);
    dataset_old.write(nodes_buffer, H5::PredType::NATIVE_DOUBLE);

    double *elems_buffer_data = new double[elems_tmp.size()*ElemDataFields];
    elems_buffer_nodes = new int[elems_tmp.size()*3];
    for(unsigned i=0;i<elems_tmp.size();i++)
    {
        Element* e = elems_tmp[i];
        elems_buffer_nodes[i*3+0] = e->nds[0]->locId;
        elems_buffer_nodes[i*3+1] = e->nds[1]->locId;
        elems_buffer_nodes[i*3+2] = e->nds[2]->locId;
        elems_buffer_data[i*ElemDataFields+0] = e->hydrostatic_stress;
        elems_buffer_data[i*ElemDataFields+1] = e->GreenStrain(0,0);
        elems_buffer_data[i*ElemDataFields+2] = e->GreenStrain(1,1);
        elems_buffer_data[i*ElemDataFields+3] = e->GreenStrain(1,0);
    }

    hsize_t dimsf_elems_nodes_old[2] = {elems_tmp.size(), 3};
    hsize_t dimsf_elems_data_old[2] = {elems_tmp.size(), ElemDataFields};
    H5::DataSpace dataspace_elems_nodes_old(2, dimsf_elems_nodes_old);
    H5::DataSpace dataspace_elems_data_old(2, dimsf_elems_data_old);

    H5::DataSet dataset_elems_nodes_old = file.createDataSet("Elements_Old", datatype_int, dataspace_elems_nodes_old);
    H5::DataSet dataset_elems_data_old = file.createDataSet("Elements_Old_Data", datatype_double, dataspace_elems_data_old);

    dataset_elems_nodes_old.write(elems_buffer_nodes, H5::PredType::NATIVE_INT);
    dataset_elems_data_old.write(elems_buffer_data, H5::PredType::NATIVE_DOUBLE);

    delete[] nodes_buffer;
    delete[] elems_buffer_nodes;
    delete[] elems_buffer_data;

}


