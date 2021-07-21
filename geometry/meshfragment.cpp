#include <unordered_map>
#include <hdf5.h>
#include <gmsh.h>
#include "meshfragment.h"

icy::ConcurrentPool<icy::BVHN> icy::MeshFragment::BVHNLeafFactory(50000);
icy::ConcurrentPool<icy::Node> icy::MeshFragment::NodeFactory(50000);
icy::ConcurrentPool<icy::Element> icy::MeshFragment::ElementFactory(50000);


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
    int point3 = gmsh::model::occ::addPoint(width/2, height*1.1, 0, 1.0);
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
    freeNodeCount = 0;
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
        freeNodeCount++;
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
        freeNodeCount++;
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
        freeNodeCount++;
    }


    // mark the "outer" and "inner" boundaries
    for(int i=0;i<2;i++)
    {
        gmsh::model::mesh::getNodesForPhysicalGroup(1, groups[i], nodeTags, nodeCoords);
        for(const std::size_t tag : nodeTags) nodes[mtags[tag]]->group.set(i);
    }


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
            for(int j=0;j<3;j++) elem->nds[j] = nodes[mtags[nodeTagsInTris[i*3+j]]];
            elem->PrecomputeInitialArea();
            if(elem->area_initial < 0)
            {
                std::swap(elem->nds[0], elem->nds[1]);
                elem->PrecomputeInitialArea();
                if(elem->area_initial < 0) throw std::runtime_error("icy::Mesh::Reset - error");
            }
            for(int j=0;j<3;j++) elem->nds[j]->area += elem->area_initial/3;
            elem->group=k;
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
        int idx1 = mtags[nodeTagsInEdges[i*2+0]];
        int idx2 = mtags[nodeTagsInEdges[i*2+1]];
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
    std::cout << "\nicy::MeshFragment::GenerateContainer\n";

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
    std::cout << "\nicy::MeshFragment::GenerateIndenter\n";

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
    freeNodeCount = 0;

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, int> mtags; // nodeTag -> sequential position in nodes[]

    if(deformable)
    {
        gmsh::model::mesh::generate(2);

        // GET NODES
        gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);

        // set the size of the resulting nodes array
        nodes.resize(nodeTags.size());
        for(unsigned i=0;i<nodeTags.size();i++)
        {
            std::size_t tag = nodeTags[i];
            mtags[tag] = i;
            double x = nodeCoords[i*3+0];
            double y = nodeCoords[i*3+1];

            icy::Node* nd = NodeFactory.take();
            nodes[i] = nd;
            nd->Reset(i, nodeCoords[i*3+0], nodeCoords[i*3+1]);
            if(y==0 || !deformable) nd->pinned=true;
            else freeNodeCount++;
        }

        // GET ELEMENTS
        std::vector<std::size_t> trisTags, nodeTagsInTris;
        gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

        for(std::size_t i=0;i<trisTags.size();i++)
        {
            icy::Element *elem = ElementFactory.take();
            elem->Reset();
            for(int j=0;j<3;j++) elem->nds[j] = nodes[mtags[nodeTagsInTris[i*3+j]]];
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
        nodes.resize(nodeTags.size());
        for(unsigned i=0;i<nodeTags.size();i++)
        {
            std::size_t tag = nodeTags[i];
            mtags[tag] = i;
            double x = nodeCoords[i*3+0];
            double y = nodeCoords[i*3+1];

            icy::Node* nd = NodeFactory.take();
            nd->Reset(i, nodeCoords[i*3+0], nodeCoords[i*3+1]);
            nodes[i] = nd;
            nd->pinned=true;
        }


        // BOUNDARIRES - EDGES

        std::vector<std::size_t> edgeTags, nodeTagsInEdges;
        gmsh::model::mesh::getElementsByType(1, edgeTags, nodeTagsInEdges);

        for(std::size_t i=0;i<edgeTags.size();i++)
        {
            int idx1 = mtags[nodeTagsInEdges[i*2+0]];
            int idx2 = mtags[nodeTagsInEdges[i*2+1]];
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

    leafs_for_ccd.clear();
    leafs_for_contact.clear();
    leafs_for_ccd.reserve(boundary_edges.size());
    leafs_for_contact.reserve(boundary_edges.size());

    for(auto edge : boundary_edges)
    {
        Node *nd1 = edge.first;
        Node *nd2 = edge.second;
        BVHN *leaf_ccd = BVHNLeafFactory.take();
        BVHN *leaf_contact = BVHNLeafFactory.take();

        leaf_contact->feature_idx = leaf_ccd->feature_idx = edge_idx++;
        leaf_contact->test_self_collision = leaf_ccd->test_self_collision = false;
        leaf_contact->isLeaf = leaf_ccd->isLeaf = true;

        leafs_for_ccd.push_back(leaf_ccd);
        leafs_for_contact.push_back(leaf_contact);
    }
    qDebug() << "icy::MeshFragment::GenerateLeafs() " << leafs_for_ccd.size() << " " << leafs_for_contact.size();
}

void icy::MeshFragment::SaveFragment(std::string fileName)
{

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

    for(unsigned i=0;i<nFirstGroupElems;i++)
    {
        Element *elem = ElementFactory.take();
        elem->Reset();
        elem->group=1;
        for(int j=0;j<3;j++) elem->nds[j]=nodes_tmp[elems[i]->nds[j]->locId];
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
    std::unordered_map<std::size_t, int> mtags; // nodeTag -> sequential position in nodes[]


    // GET NODES
    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);

    std::cout << "GET NODES" << std::endl;
    // set the size of the resulting nodes array
    nodes.resize(nodeTags.size());
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        mtags[tag] = i;

        double x = nodeCoords[i*3+0];
        double y = nodeCoords[i*3+1];
        if(!btags.count(tag))
        {
            icy::Node* nd = NodeFactory.take();
            nd->Reset(i, x, y);
            nodes_tmp.push_back(nd);
        }
        else
        {
            // node already exists in nodes_tmp
            icy::Node* nd = nodes_tmp[btags.at(tag)];
            if(nd->x_initial.x()!=x || nd->x_initial.y()!=y) throw std::runtime_error("some tags don't match after remeshing");
        }

    }

    freeNodeCount = nodes_tmp.size();

    std::cout << "GET ELEMENTS" << std::endl;
    // GET ELEMENTS
    std::vector<std::size_t> trisTags, nodeTagsInTris;
    gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

    for(std::size_t i=0;i<trisTags.size();i++)
    {
        std::cout << "elem" << i << std::endl;

        icy::Element *elem = ElementFactory.take();
        elem->Reset();
        for(int j=0;j<3;j++) elem->nds[j] = nodes_tmp[mtags[nodeTagsInTris[i*3+j]]];
        elem->PrecomputeInitialArea();
        if(elem->area_initial < 0)
        {
            std::swap(elem->nds[0], elem->nds[1]);
            elem->PrecomputeInitialArea();
            if(elem->area_initial < 0) throw std::runtime_error("icy::Mesh::Reset - error");
        }
        for(int j=0;j<3;j++) elem->nds[j]->area += elem->area_initial/3;
        elems_tmp.push_back(elem);
    }


    //gmsh::fltk::run();
    gmsh::clear();
    Swap();

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
}


//
