#include "meshfragment.h"
#include <gmsh.h>
#include <unordered_set>

icy::ConcurrentPool<icy::BVHN> icy::MeshFragment::BVHNLeafFactory(50000);


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

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", ElementSize);

    GetFromGmsh();
}

void icy::MeshFragment::GenerateContainer(double ElementSize, double offset)
{
    deformable = false;

    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
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

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", ElementSize);

    GetFromGmsh();
}


void icy::MeshFragment::GenerateIndenter(double ElementSize)
{
    deformable = false;

    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    double height = 1;
    double radius = 0.15;
//    int point1 = gmsh::model::occ::addPoint(0, height+radius*1.1, 0, 1.0);

    gmsh::model::occ::addEllipse(0, height+radius*1.1, 0, radius, radius/2);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", ElementSize);

    GetFromGmsh();
}

void icy::MeshFragment::GenerateCup(double ElementSize)
{
    deformable = false;

    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    double width = 2;
    double height = 1;
    int point0 = gmsh::model::occ::addPoint(-width/2, height*2, 0, 1.0);
    int point1 = gmsh::model::occ::addPoint(-width/2, height, 0, 1.0);
    int point2 = gmsh::model::occ::addPoint(-width/20, 0, 0, 1.0);
    int point3 = gmsh::model::occ::addPoint(width/20, 0, 0, 1.0);
    int point4 = gmsh::model::occ::addPoint(width/2, height, 0, 1.0);
    int point5 = gmsh::model::occ::addPoint(width/2, height*2, 0, 1.0);

    int line0 = gmsh::model::occ::addLine(point0, point1);
    int line1 = gmsh::model::occ::addLine(point1, point2);
    int line2 = gmsh::model::occ::addLine(point2, point3);
    int line3 = gmsh::model::occ::addLine(point3, point4);
    int line4 = gmsh::model::occ::addLine(point4, point5);
    int line5 = gmsh::model::occ::addLine(point5, point0);

    std::vector<int> curveTags;
    curveTags.push_back(line0);
    curveTags.push_back(line1);
    curveTags.push_back(line2);
    curveTags.push_back(line3);
    curveTags.push_back(line4);
    curveTags.push_back(line5);
    gmsh::model::occ::addWire(curveTags);

    gmsh::model::occ::synchronize();
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", ElementSize*10);

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
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", ElementSize);
    GetFromGmsh();
}


void icy::MeshFragment::GetFromGmsh()
{
    gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary",0);
    gmsh::model::mesh::generate(deformable ? 2 : 1);


    // retrieve the result
    elems.clear();
    nodes.clear();
    boundary_nodes.clear();
    boundary_edges.clear();
    freeNodeCount = 0;


    // nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords);

    // elems
    std::vector<std::size_t> trisTags, nodeTagsInTris;
    gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

    // boundary
    std::vector<std::size_t> edgeTags, nodeTagsInEdges;
    gmsh::model::mesh::getElementsByType(1, edgeTags, nodeTagsInEdges);

    // nodeTags, nodeCoords, nodeTagsInTris => nodes, elems
    std::map<std::size_t, int> nodeTagsMap1; // nodeTag -> its sequential position in nodeTag
    for(std::size_t i=0;i<nodeTags.size();i++) nodeTagsMap1[nodeTags[i]] = i;

    std::unordered_set<std::size_t> tagSet; // only keep nodes from the tris and edges
    for(std::size_t &tag : nodeTagsInTris) tagSet.insert(tag);
    for(std::size_t &tag : nodeTagsInEdges) tagSet.insert(tag);

    // set the size of the resulting nodes array
    nodes.resize(tagSet.size());

    std::map<std::size_t, int> mtags; // nodeTag -> sequential position
    int count = 0;

    for(const std::size_t &tag : tagSet)
    {
        int idx1 = nodeTagsMap1[tag];
        double x = nodeCoords[idx1*3+0];
        double y = nodeCoords[idx1*3+1];

        icy::Node &nd = nodes[count];
        nd.Reset();
        nd.locId = count;
        nd.x_initial << x, y;
        nd.intended_position = nd.xt = nd.xn = nd.x_initial;
        if(y==0 || !deformable) nd.pinned=true;
        else freeNodeCount++;
        mtags[tag] = count;
        count++;
    }

    // resulting elements array
    elems.resize(nodeTagsInTris.size()/3);

    for(std::size_t i=0;i<nodeTagsInTris.size()/3;i++)
    {
        icy::Element &elem = elems[i];
        for(int j=0;j<3;j++) elem.nds[j] = &(nodes[mtags[nodeTagsInTris[i*3+j]]]);
        elem.PrecomputeInitialArea();
        if(elem.area_initial == 0) throw std::runtime_error("icy::Mesh::Reset - element's area is zero");
        if(elem.area_initial < 0)
        {
            for(int j=0;j<3;j++) elem.nds[2-j] = &(nodes[mtags[nodeTagsInTris[i*3+j]]]);
            elem.PrecomputeInitialArea();
            if(elem.area_initial < 0) throw std::runtime_error("icy::Mesh::Reset - error");
        }
        for(int j=0;j<3;j++) elem.nds[j]->area += elem.area_initial/3;
    }

    std::unordered_set<int> set_nds; // set of boundary nodes
    boundary_edges.resize(nodeTagsInEdges.size()/2);
    for(std::size_t i=0;i<nodeTagsInEdges.size()/2;i++)
    {
        int idx1 = mtags[nodeTagsInEdges[i*2+0]];
        int idx2 = mtags[nodeTagsInEdges[i*2+1]];
        set_nds.insert(idx1);
        set_nds.insert(idx2);
        Node *nd1, *nd2;
        if(idx1<idx2)
        {
            nd1 = &nodes[idx1];
            nd2 = &nodes[idx2];
        }
        else
        {
            nd1 = &nodes[idx2];
            nd2 = &nodes[idx1];
        }
        boundary_edges[i]=std::make_pair(nd1,nd2);
    }

    boundary_nodes.resize(set_nds.size());
    std::copy(set_nds.begin(), set_nds.end(), boundary_nodes.begin());
    gmsh::clear();
}

void icy::MeshFragment::GenerateLeafs(unsigned edge_idx)
{
    root_ccd.isLeaf = false;
    root_ccd.test_self_collision = false;
    root_contact.isLeaf = false;
    root_contact.test_self_collision = false;

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




/*

void MeshFragment::GenerateSelfCollisionTest(double ElementSize)
{
    deformable = true;

    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    GetFromGmsh();
}

void MeshFragment::GenerateCircle(double x, double y, double r, double ElementSize)
{
    deformable = true;


    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    GetFromGmsh();

}


void MeshFragment::GenerateCup(double ElementSize)
{
    deformable = false;

    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");


    GetFromGmsh();
}
*/
