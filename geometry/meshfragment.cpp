#include <unordered_map>
#include <unordered_set>
#include <H5Cpp.h>
#include <gmsh.h>
#include "meshfragment.h"
#include "tbb/concurrent_unordered_map.h"
#include "edge.h"
#include "mesh.h"

icy::ConcurrentPool<icy::BVHN> icy::MeshFragment::BVHNLeafFactory(nPreallocate);
icy::ConcurrentPool<icy::Node> icy::MeshFragment::NodeFactory(nPreallocate);
icy::ConcurrentPool<icy::Element> icy::MeshFragment::ElementFactory(nPreallocate);
icy::ConcurrentPool<icy::BoundaryEdge> icy::MeshFragment::BoundaryEdgeFactory(nPreallocate);
icy::ConcurrentPool<icy::CohesiveZone> icy::MeshFragment::CZFactory(nPreallocate);

icy::MeshFragment::~MeshFragment()
{
    // release objects back to their pools
    BVHNLeafFactory.release(leaves_for_ccd);
    NodeFactory.release(nodes);
    ElementFactory.release(elems);
}

icy::Node* icy::MeshFragment::AddNode()
{
    Node* result = NodeFactory.take();
    result->Reset();
    result->fragment = this;
    result->locId = (int)nodes.size();
    nodes.push_back(result);
    result->globId = (int)parentMesh->allNodes.size();
    parentMesh->allNodes.push_back(result);
    return result;
}

icy::Element* icy::MeshFragment::AddElement()
{
    Element* elem = ElementFactory.take();
    elem->Reset();
    elems.push_back(elem);
    parentMesh->allElems.push_back(elem);
    return elem;
}

icy::CohesiveZone* icy::MeshFragment::AddCZ()
{
    CohesiveZone *cz = CZFactory.take();
    czs.push_back(cz);
    parentMesh->allCZs.push_back(cz);
    return cz;
}



void icy::MeshFragment::GenerateBrick(double ElementSize, double width, double height)
{
    isDeformable = true;

    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    int point1 = gmsh::model::occ::addPoint(-width/2, 0, 0, 1.0);
    int point2 = gmsh::model::occ::addPoint(-width/2, height, 0, 1.0);
    int point3 = gmsh::model::occ::addPoint(width/2, height, 0, 1.0);
    int point4 = gmsh::model::occ::addPoint(width/2, 0, 0, 1.0);

    int line1 = gmsh::model::occ::addLine(point1, point2);
    int line2 = gmsh::model::occ::addLine(point2, point3);
    int line3 = gmsh::model::occ::addLine(point3, point4);
    int line4 = gmsh::model::occ::addLine(point4, point1);

    int loopTag = gmsh::model::occ::addCurveLoop({line1,line2,line3,line4});

    std::vector<int> loops;
    loops.push_back(loopTag);
    gmsh::model::occ::addPlaneSurface(loops);
    gmsh::model::occ::synchronize();

    int groupTag0 = gmsh::model::addPhysicalGroup(1, {line1,line2,line3,line4});  // whole boundary
    int groupTag1 = gmsh::model::addPhysicalGroup(1, {line3});  // right boundary
    int groupTag2 = gmsh::model::addPhysicalGroup(1, {line1});  // left boundary
    int groupTag3 = gmsh::model::addPhysicalGroup(0, {point2, point4});
    std::array<int,3> dim1groups {groupTag0,groupTag1,groupTag2};

    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    mtags.clear();
    gmsh::model::mesh::generate(2);

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);

    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");

        Node *nd = AddNode();
        mtags[tag] = nd->locId;
        nd->Initialize(nodeCoords[i*3+0], nodeCoords[i*3+1]);
    }

    // mark the nodes of the special boundary
    for(unsigned i=0;i<dim1groups.size();i++)
    {
        nodeTags.clear();
        nodeCoords.clear();
        gmsh::model::mesh::getNodesForPhysicalGroup(1, dim1groups[i], nodeTags, nodeCoords);
        for(unsigned j=0;j<nodeTags.size();j++) nodes[mtags[nodeTags[j]]]->group.set(i);
    }
    for(Node *nd : nodes) nd->isBoundary = nd->group.test(0);

    // dim 0 boundary
    nodeTags.clear();
    nodeCoords.clear();
    gmsh::model::mesh::getNodesForPhysicalGroup(0, groupTag3, nodeTags, nodeCoords);
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        icy::Node* nd = nodes[mtags[nodeTags[i]]];
        nd->group.set(3);
    }

    // GET ELEMENTS
    std::vector<std::size_t> trisTags, nodeTagsInTris;
    gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

    for(std::size_t i=0;i<trisTags.size();i++)
    {
        icy::Element *elem = AddElement();
        elem->Initialize(nodes[mtags.at(nodeTagsInTris[i*3+0])],
                nodes[mtags.at(nodeTagsInTris[i*3+1])],
                nodes[mtags.at(nodeTagsInTris[i*3+2])]);
    }

    // BOUNDARIRES - EDGES
    std::vector<std::size_t> edgeTags, nodeTagsInEdges;
    gmsh::model::mesh::getElementsByType(1, edgeTags, nodeTagsInEdges);

    for(std::size_t i=0;i<edgeTags.size();i++)
    {
        icy::Node *nd0 = nodes[mtags.at(nodeTagsInEdges[i*2+0])];
        icy::Node *nd1 = nodes[mtags.at(nodeTagsInEdges[i*2+1])];
        if(nd0->group.test(1) && nd1->group.test(1)) special_boundary.emplace_back(nd0,nd1);
    }

    PostMeshingEvaluations();
    gmsh::clear();
    ConnectIncidentElements();
}

void icy::MeshFragment::GenerateSelfCollisionBrick(double ElementSize, double width, double height)
{
    isDeformable = true;

    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    int point1 = gmsh::model::occ::addPoint(-width/2, 0, 0, 1.0);
    int point2 = gmsh::model::occ::addPoint(-width/2, height, 0, 1.0);
    int point3 = gmsh::model::occ::addPoint(-width/50, height, 0, 1.0);
    int point4 = gmsh::model::occ::addPoint(0, height/2, 0, 1.0);
    int point5 = gmsh::model::occ::addPoint(width/50, height, 0, 1.0);
    int point6 = gmsh::model::occ::addPoint(width/2, height*1.1, 0, 1.0);
    int point7 = gmsh::model::occ::addPoint(width/2, 0, 0, 1.0);

    int line1 = gmsh::model::occ::addLine(point1, point2);
    int line2 = gmsh::model::occ::addLine(point2, point3);
    int line3 = gmsh::model::occ::addLine(point3, point4);
    int line4 = gmsh::model::occ::addLine(point4, point5);
    int line5 = gmsh::model::occ::addLine(point5, point6);
    int line6 = gmsh::model::occ::addLine(point6, point7);
    int line7 = gmsh::model::occ::addLine(point7, point1);

    int loopTag = gmsh::model::occ::addCurveLoop({line1,line2,line3,line4,line5,line6,line7});

    std::vector<int> loops;
    loops.push_back(loopTag);
    gmsh::model::occ::addPlaneSurface(loops);
    gmsh::model::occ::synchronize();

    int groupTag0 = gmsh::model::addPhysicalGroup(1, {line1,line2,line3,line4,line5,line6,line7});
    int groupTag1 = gmsh::model::addPhysicalGroup(1, {line6});
    int groupTag2 = gmsh::model::addPhysicalGroup(1, {line1});
    std::array<int,3> dim1groups {groupTag0,groupTag1,groupTag2};

    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    mtags.clear();
    gmsh::model::mesh::generate(2);

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);

    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // don't duplicate nodes
        Node *nd = AddNode();
        mtags[tag] = nd->locId;
        nd->Initialize(nodeCoords[i*3+0], nodeCoords[i*3+1]);
    }

    // mark the nodes of the special boundary
    // nodes of the "inner boundary" are placed first
    for(unsigned i=0;i<dim1groups.size();i++)
    {
        nodeTags.clear();
        nodeCoords.clear();
        gmsh::model::mesh::getNodesForPhysicalGroup(1, dim1groups[i], nodeTags, nodeCoords);
        for(unsigned j=0;j<nodeTags.size();j++) nodes[mtags[nodeTags[j]]]->group.set(i);
    }
    for(Node *nd : nodes)
    {
        nd->isBoundary = nd->group.test(0);
        if(nd->group.test(1) || nd->group.test((2))) nd->pinned = true;
    }

    // GET ELEMENTS
    std::vector<std::size_t> trisTags, nodeTagsInTris;
    gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

    for(std::size_t i=0;i<trisTags.size();i++)
    {
        icy::Element *elem = AddElement();
        elem->Initialize(nodes[mtags.at(nodeTagsInTris[i*3+0])],
                nodes[mtags.at(nodeTagsInTris[i*3+1])],
                nodes[mtags.at(nodeTagsInTris[i*3+2])]);
    }

    // BOUNDARIRES - EDGES
    std::vector<std::size_t> edgeTags, nodeTagsInEdges;
    gmsh::model::mesh::getElementsByType(1, edgeTags, nodeTagsInEdges);

    for(std::size_t i=0;i<edgeTags.size();i++)
    {
        icy::Node *nd0 = nodes[mtags.at(nodeTagsInEdges[i*2+0])];
        icy::Node *nd1 = nodes[mtags.at(nodeTagsInEdges[i*2+1])];
        if(nd0->group.test(1) && nd1->group.test(1)) special_boundary.emplace_back(nd0,nd1);
    }

    PostMeshingEvaluations();
    gmsh::clear();
    ConnectIncidentElements();
}


void icy::MeshFragment::GenerateCZBrick(double ElementSize, double width, double height)
{
    isDeformable = true;

    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    std::vector<std::pair<double,double>> ptCoords {
        {0,0},                      // 0
        {-width/2,0},               // 1
        {-width/2,height},          // 2
        {-width/50,height},         // 3
        {0,height/2},               // 4
        {width/50,height},          // 5
        {width/2,height*1.1},       // 6
        {width/2,0}                 // 7
        };

    std::vector<int> ptTags, lts;
    for(auto p : ptCoords) ptTags.push_back(gmsh::model::occ::addPoint(p.first,p.second,0,1.0));

    for(unsigned i=0;i<ptTags.size();i++) lts.push_back(gmsh::model::occ::addLine(ptTags[i],ptTags[(i+1)%ptTags.size()]));
    int czLineTag = gmsh::model::occ::addLine(ptTags[4],ptTags[0]);

//    int loop1 = gmsh::model::occ::addCurveLoop({lts[0],lts[1],lts[2],lts[3],lts[4],lts[5],lts[6],lts[7]});
    int loop2 = gmsh::model::occ::addCurveLoop({lts[0],lts[1],lts[2],lts[3],czLineTag});
    int loop3 = gmsh::model::occ::addCurveLoop({lts[4],lts[5],lts[6],lts[7],-czLineTag});


    gmsh::model::occ::addPlaneSurface({loop2});
    gmsh::model::occ::addPlaneSurface({loop3});
    gmsh::model::occ::synchronize();

//    gmsh::model::mesh::embed(1, {czLineTag},2, surfaceTag);

    int groupTag0 = gmsh::model::addPhysicalGroup(1, {lts[0],lts[1],lts[2],lts[3],lts[4],lts[5],lts[6],lts[7]});
    int groupTag1 = gmsh::model::addPhysicalGroup(1, {lts[6]});
    int groupTag2 = gmsh::model::addPhysicalGroup(1, {lts[1]});
    int groupTag3 = gmsh::model::addPhysicalGroup(1, {czLineTag});
    std::vector<int> dim1groups {groupTag0,groupTag1,groupTag2,groupTag3};

    int groupTag6 = gmsh::model::addPhysicalGroup(2, {loop2});
    int groupTag7 = gmsh::model::addPhysicalGroup(2, {loop3});

    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    mtags.clear();
    gmsh::model::mesh::generate(2);

    // GET NODES
    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);

    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // don't duplicate nodes
        Node *nd = AddNode();
        mtags[tag] = nd->locId;
        nd->Initialize(nodeCoords[i*3+0], nodeCoords[i*3+1]);
    }

    // mark the nodes of the special boundary
    // nodes of the "inner boundary" are placed first
    for(unsigned i=0;i<dim1groups.size();i++)
    {
        nodeTags.clear();
        nodeCoords.clear();
        gmsh::model::mesh::getNodesForPhysicalGroup(1, dim1groups[i], nodeTags, nodeCoords);
        for(unsigned j=0;j<nodeTags.size();j++) nodes[mtags[nodeTags[j]]]->group.set(i);
    }

    std::unordered_map<Node*,Node*> splitNodes; // for CZ insertion
    for(Node *nd : nodes)
    {
        nd->isBoundary = nd->group.test(0);
        if(nd->group.test(1) || nd->group.test(2)) nd->pinned = true;
        if(nd->group.test(3))
        {
            Node *nd2 = AddNode();
            nd2->Initialize(nd);
            splitNodes[nd] = nd2;
        }
    }

    // left group
    nodeTags.clear(); nodeCoords.clear();
    gmsh::model::mesh::getNodesForPhysicalGroup(2, groupTag6, nodeTags, nodeCoords);
    for(unsigned j=0;j<nodeTags.size();j++) nodes[mtags[nodeTags[j]]]->group.set(6);

    // right group
    nodeTags.clear(); nodeCoords.clear();
    gmsh::model::mesh::getNodesForPhysicalGroup(2, groupTag7, nodeTags, nodeCoords);
    for(unsigned j=0;j<nodeTags.size();j++) nodes[mtags[nodeTags[j]]]->group.set(7);


    // GET ELEMENTS
    std::vector<std::size_t> trisTags, nodeTagsInTris;
    gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

    for(std::size_t i=0;i<trisTags.size();i++)
    {
        icy::Element *elem = AddElement();
        elem->Initialize(nodes[mtags.at(nodeTagsInTris[i*3+0])],
                         nodes[mtags.at(nodeTagsInTris[i*3+1])],
                         nodes[mtags.at(nodeTagsInTris[i*3+2])]);
    }

    // define element groups and replace nodes with splits
    for(Element *e : elems)
    {
        if(e->nds[0]->group.test(6) && e->nds[1]->group.test(6) &&e->nds[2]->group.test(6)) e->group.set(0);
        if(e->nds[0]->group.test(7) && e->nds[1]->group.test(7) &&e->nds[2]->group.test(7)) e->group.set(1);

        for(int i=0;i<3;i++)
            if(e->group.test(0) && e->nds[i]->group.test(3)) e->nds[i] = splitNodes.at(e->nds[i]);
    }

    PostMeshingEvaluations();

    // Edges and cohesive zones
    std::vector<std::size_t> edgeTags, nodeTagsInEdges;
    gmsh::model::mesh::getElementsByType(1, edgeTags, nodeTagsInEdges);

    for(std::size_t i=0;i<edgeTags.size();i++)
    {
        icy::Node *nd0 = nodes[mtags.at(nodeTagsInEdges[i*2+0])];
        icy::Node *nd1 = nodes[mtags.at(nodeTagsInEdges[i*2+1])];
        if(nd0->group.test(1) && nd1->group.test(1)) special_boundary.emplace_back(nd0,nd1);

        // insert cohesive zones
        if(nd0->group.test(3) && nd1->group.test(3))
        {
            CohesiveZone *cz = AddCZ();
            Node *nd2 = splitNodes.at(nd0);
            Node *nd3 = splitNodes.at(nd1);
            auto iter1 = std::find_if(elems.begin(),elems.end(),
                                   [nd0,nd1](Element *elem){return elem->containsEdge(nd0,nd1);});

            auto iter2 = std::find_if(elems.begin(),elems.end(),
                                   [nd2,nd3](Element *elem){return elem->containsEdge(nd2,nd3);});
            if(iter1==elems.end() || iter2==elems.end()) throw std::runtime_error("GenerateCZBrick");

            Element *e1 = *iter1;
            Element *e2 = *iter2;

            uint8_t edgeIdx1 = e1->getEdgeIdx(nd0,nd1);
            uint8_t edgeIdx2 = e2->getEdgeIdx(nd2,nd3);

            if(e1->isEdgeCW(nd0,nd1)) cz->Initialize(e1,edgeIdx1,e2,edgeIdx2);
            else cz->Initialize(e2,edgeIdx2,e1,edgeIdx1);
        }
    }

    gmsh::clear();
    ConnectIncidentElements();
    spdlog::info("GenerateCZBrick done; czs {}",czs.size());
    for(auto *cz : czs) cz->GetNodes();

}

void icy::MeshFragment::GenerateContainer(double ElementSize, double offset)
{
    std::cout << "\nGenerateContainer" << std::endl;;

    isDeformable = false;

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

void icy::MeshFragment::GenerateIndenter(double ElementSize, double cx, double cy, double radius, double aspect)
{
    isDeformable = false;
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::clear();
    gmsh::model::add("indenter1");
    gmsh::model::occ::addEllipse(cx, cy, 0, radius, radius/aspect);
    gmsh::model::occ::synchronize();
    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);
    GetFromGmsh();
}


void icy::MeshFragment::GetFromGmsh()
{
    if(elems.size() || nodes.size()) throw std::runtime_error("GetFromGmsh(): elem/node/edge arrays not empty");

    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    mtags.clear();

    if(isDeformable)
    {
        gmsh::model::mesh::generate(2);

        // GET NODES
        gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);

        // set the size of the resulting nodes array
        for(unsigned i=0;i<nodeTags.size();i++)
        {
            std::size_t tag = nodeTags[i];
            if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");

            Node *nd = AddNode();
            mtags[tag] = nd->locId;
            nd->Initialize(nodeCoords[i*3+0], nodeCoords[i*3+1]);
        }

        // GET ELEMENTS
        std::vector<std::size_t> trisTags, nodeTagsInTris;
        gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

        for(std::size_t i=0;i<trisTags.size();i++)
        {
            icy::Element *elem = AddElement();
            elem->Initialize(nodes[mtags.at(nodeTagsInTris[i*3+0])],
                    nodes[mtags.at(nodeTagsInTris[i*3+1])],
                    nodes[mtags.at(nodeTagsInTris[i*3+2])]);
        }

        PostMeshingEvaluations();
        ConnectIncidentElements();
    }
    else
    {
        // not deformable - only generate the boundary
        gmsh::model::mesh::generate(1);

        // GET NODES
        gmsh::model::mesh::getNodesByElementType(1, nodeTags, nodeCoords, parametricCoords);

        // set the size of the resulting nodes array
        for(unsigned i=0;i<nodeTags.size();i++)
        {
            std::size_t tag = nodeTags[i];
            if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in boundary");

            Node *nd = AddNode();
            mtags[tag] = nd->locId;
            nd->Initialize(nodeCoords[i*3+0], nodeCoords[i*3+1]);
            nd->pinned = true;
            nd->isBoundary = true;
        }

        // BOUNDARIRES - EDGES
        std::vector<std::size_t> edgeTags, nodeTagsInEdges;
        gmsh::model::mesh::getElementsByType(1, edgeTags, nodeTagsInEdges);
        boundaryEdges.reserve(edgeTags.size());
        for(std::size_t i=0;i<edgeTags.size();i++)
        {
            Node *nd0 = nodes[mtags.at(nodeTagsInEdges[i*2+0])];
            Node *nd1 = nodes[mtags.at(nodeTagsInEdges[i*2+1])];
            BoundaryEdge *be = BoundaryEdgeFactory.take();
            be->Initialize(nd0,nd1);
            boundaryEdges.push_back(be);
        }
    }


    gmsh::clear();
}


void icy::MeshFragment::CreateLeaves()
{
    root_ccd.boundaryEdge = nullptr;
    root_ccd.test_self_collision = this->isDeformable;
    root_ccd.isLeaf = false;

    BVHNLeafFactory.release(leaves_for_ccd);

    leaves_for_ccd.reserve(boundaryEdges.size());

    for(const BoundaryEdge *boundary : boundaryEdges)
    {
        BVHN *leaf_ccd = BVHNLeafFactory.take();
        leaves_for_ccd.push_back(leaf_ccd);
        leaf_ccd->test_self_collision = false;
        leaf_ccd->boundaryEdge = boundary;
        leaf_ccd->isLeaf = true;
    }
}







void icy::MeshFragment::SaveFragment(std::string fileName)
{
    /*
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
        elems_buffer_data[i*ElemDataFields+0] = e->PiMultiplier(0,0);
        elems_buffer_data[i*ElemDataFields+1] = e->PiMultiplier(0,1);
        elems_buffer_data[i*ElemDataFields+2] = e->PiMultiplier(1,0);
        elems_buffer_data[i*ElemDataFields+3] = e->PiMultiplier(1,1);
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
    */
}


void icy::MeshFragment::PostMeshingEvaluations()
{
    unsigned nElems = elems.size();
    unsigned nNodes = nodes.size();

#pragma omp parallel for
    for(unsigned i=0;i<nElems;i++) elems[i]->PrecomputeInitialArea();

    for(unsigned i=0;i<nNodes;i++) nodes[i]->area = 0;

    for(unsigned i=0;i<nElems;i++)
    {
        icy::Element *elem = elems[i];
        for(int j=0;j<3;j++)
        {
            icy::Node *nd = elem->nds[j];
            nd->area += elem->area_initial/3;
        }
    }
}



void icy::MeshFragment::ConnectIncidentElements()
{
    // clear adj_elems vector and isBoundary flag
#pragma omp parallel for
    for(std::size_t i=0;i<nodes.size();i++)
    {
        icy::Node* nd = nodes[i];
        nd->adj_elems.clear();
    }

    // edges_map will hold all edges and their connected elements
    std::unordered_map<uint64_t, Edge> edges_map;

    // associate edges with one or two adjacent elements
    // for simplicity, this part of the algorithm is single-threaded

    for(icy::Element *elem : elems)
    {
        for(int i=0;i<3;i++)
        {
            //elem->incident_elems[i] = nullptr;
            elem->nds[i]->adj_elems.push_back(elem);

            // insert element's edges into edges_map
            Node *nd0 = elem->nds[i];
            Node *nd1 = elem->nds[(i+1)%3];
            uint64_t key = Node::make_local_key(nd0,nd1);

            auto result = edges_map.try_emplace(key, nd0, nd1);
            Edge &e = result.first->second;
            e.AddElement(elem, (i+2)%3);
        }
    }


    for(const auto &kvpair : edges_map)
    {
        const icy::Edge &e = kvpair.second;
        short idx0 = e.edge_in_elem_idx[0];
        short idx1 = e.edge_in_elem_idx[1];

        if(e.containsCZ()) continue;

        icy::Element *elem_of_edge0 = static_cast<Element*>(e.elems[0]);
        icy::Element *elem_of_edge1 = static_cast<Element*>(e.elems[1]);
        if(elem_of_edge0 == nullptr && elem_of_edge1 == nullptr) throw std::runtime_error("CreateEdges(): disconnected edge");

        // ensure that each element "knows" its incident elements (remains nullptr if element is on the boundary)
        if(!e.isBoundary())
        {
            elem_of_edge0->incident_elems[idx0] = elem_of_edge1;
            elem_of_edge1->incident_elems[idx1] = elem_of_edge0;
        }
        else
        {
            if(e.elems[0] != nullptr)
                AddBoundary(elem_of_edge0,e.edge_in_elem_idx[0]);
            else
                AddBoundary(elem_of_edge1,e.edge_in_elem_idx[1]);
        }
    }

    //ConvertBoundaryEdges();
#pragma omp parallel for
    for(std::size_t i=0;i<nodes.size();i++) nodes[i]->PrepareFan();

    for(auto elem : elems)
        for(int k=0;k<3;k++)
            if(elem->incident_elems[k]==nullptr) throw std::runtime_error("ConnectIncidentElements assertion fail");
}

void icy::MeshFragment::AddBoundary(Element *elem, uint8_t edge_idx, uint8_t status)
{
    BoundaryEdge *be = BoundaryEdgeFactory.take()->Initialize(elem,edge_idx,status);
    be->elem->incident_elems[be->edge_idx] = be;
    spdlog::info("boundaryEdges before: {}",boundaryEdges.size());
    boundaryEdges.push_back(be);
    spdlog::info("boundaryEdges after: {}",boundaryEdges.size());

}
