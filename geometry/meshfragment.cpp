#include "meshfragment.h"
#include <gmsh.h>
#include <unordered_set>

icy::ConcurrentPool<icy::BVHN> icy::MeshFragment::BVHNLeafFactory(50000);




void icy::MeshFragment::GenerateSpecialBrick(double ElementSize)
{
    std::cout << "\nicy::MeshFragment::GenerateSpecialBrick\n";
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
    unsigned nGroups = 1;
    int groups[nGroups];
    groups[0] = gmsh::model::addPhysicalGroup(1, {line1, line2, line3, line4});
//    groups[1] = gmsh::model::addPhysicalGroup(2, {loopTag2});
//    groups[2] = gmsh::model::addPhysicalGroup(2, {loopTag3});

    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);
    gmsh::model::mesh::generate(2);



    // GET NODES
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);

    // set the size of the resulting nodes array
    nodes.resize(nodeTags.size());
    freeNodeCount = 0;
    std::map<std::size_t, int> mtags; // nodeTag -> sequential position in nodes[]
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        mtags[tag] = i;
        double x = nodeCoords[i*3+0];
        double y = nodeCoords[i*3+1];

        icy::Node &nd = nodes[i];
        nd.Reset();
        nd.locId = i;
        nd.x_initial << x, y;
        nd.intended_position = nd.xt = nd.xn = nd.x_initial;
        if(y==0 || !deformable) nd.pinned=true;
        else freeNodeCount++;
    }


    // physical groups of nodes
    for(int i=0;i<nGroups;i++)
    {
        std::vector<std::size_t> nodeTagsGroup;
        gmsh::model::mesh::getNodesForPhysicalGroup(1, groups[i], nodeTagsGroup, nodeCoords);
        for(const std::size_t tag : nodeTagsGroup) nodes[mtags[tag]].group.set(i);
    }


    // GET ELEMENTS
    elems.clear();
    std::vector<int> surfaces = {surfaceTag1, surfaceTag2, surfaceTag3};

    for(int j=0;j<3;j++)
    {
        std::vector<std::size_t> trisTags, nodeTagsInTris;
//        trisTags.clear();nodeTagsInTris.clear();
        gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris, surfaces[j]);

        for(std::size_t i=0;i<trisTags.size();i++)
        {
            icy::Element elem;// = elems[i];
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
            //        for(int j=0;j<3;j++)
            //            if(elem.nds[0]->group.test(j) && elem.nds[1]->group.test(j) && elem.nds[2]->group.test(j))
            //                elem.group=j;
            elem.group=j;
            elems.push_back(elem);

        }
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
            nd1 = &nodes[idx1];
            nd2 = &nodes[idx2];
        }
        else
        {
            nd1 = &nodes[idx2];
            nd2 = &nodes[idx1];
        }
        if(nd1->group.test(0) && nd2->group.test(0)) boundary_edges.emplace_back(nd1,nd2);
    }

    gmsh::clear();
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

    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);

    GetFromGmsh();
}

void icy::MeshFragment::GenerateSpecialIndenter(double ElementSize)
{
    std::cout << "\nGenerateSpecialIndenter\n";
    deformable = false;

    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    double height = 1;
    double radius = 0.15;
//    int point1 = gmsh::model::occ::addPoint(0, height+radius*1.1, 0, 1.0);

    gmsh::model::addDiscreteEntity(0, 1);
    gmsh::model::addDiscreteEntity(0, 2);
    gmsh::model::addDiscreteEntity(0, 3);

    gmsh::model::setCoordinates(1, 0, height+0.01, 0);
    gmsh::model::setCoordinates(2, -radius, height+radius, 0);
    gmsh::model::setCoordinates(3, radius, height+radius, 0);

    gmsh::model::addDiscreteEntity(1, 1, {1,2});
    gmsh::model::addDiscreteEntity(1, 2, {2,3});
    gmsh::model::addDiscreteEntity(1, 3, {3,1});

    gmsh::model::addDiscreteEntity(2, 1, {1,2,3});

    gmsh::model::mesh::createGeometry();

    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);


    gmsh::model::mesh::generate(1);


    // retrieve the result
    elems.clear();
    nodes.clear();
    boundary_edges.clear();
    freeNodeCount = 0;


    // nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    gmsh::model::mesh::getNodesByElementType(1, nodeTags, nodeCoords, parametricCoords);

    // boundary
    std::vector<std::size_t> edgeTags, nodeTagsInEdges;
    gmsh::model::mesh::getElementsByType(1, edgeTags, nodeTagsInEdges);

    // nodeTags, nodeCoords, nodeTagsInTris => nodes, elems
    std::map<std::size_t, int> nodeTagsMap1; // nodeTag -> its sequential position in nodeTag
    for(std::size_t i=0;i<nodeTags.size();i++) nodeTagsMap1[nodeTags[i]] = i;

    // set the size of the resulting nodes array
    nodes.resize(nodeTags.size());

    std::map<std::size_t, int> mtags; // nodeTag -> sequential position
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        mtags[tag] = i;

        int idx1 = nodeTagsMap1[tag];
        double x = nodeCoords[idx1*3+0];
        double y = nodeCoords[idx1*3+1];

        icy::Node &nd = nodes[i];
        nd.Reset();
        nd.locId = i;
        nd.x_initial << x, y;
        nd.intended_position = nd.xt = nd.xn = nd.x_initial;
        if(y==0 || !deformable) nd.pinned=true;
        else freeNodeCount++;
    }

    boundary_edges.resize(nodeTagsInEdges.size()/2);

    for(std::size_t i=0;i<nodeTagsInEdges.size()/2;i++)
    {
        int idx1 = mtags[nodeTagsInEdges[i*2+0]];
        int idx2 = mtags[nodeTagsInEdges[i*2+1]];
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

    std::vector<std::pair<int,int>> dimTags;
    gmsh::model::getEntities(dimTags,-1);
    std::cout << "PRINTING entities FOR THE INDENTER\n";
    for(std::pair<int,int> &val : dimTags)
        std::cout << "DIM " << val.first << "; TAG " << val.second << std::endl;

    std::vector<std::pair<int, int> > boundary_dimtags;
    gmsh::model::getBoundary({{1, 1}}, boundary_dimtags, false, false);
    std::vector<int> boundary_tags, complement_tags;
    std::cout << "PRINTING BOUNDARY_DIMTAGS FOR {1,1}\n";
      for(auto e : boundary_dimtags) {
          std::cout << "DIM " << e.first << "; TAG " << e.second << std::endl;
      }


    gmsh::clear();
}


void icy::MeshFragment::GenerateIndenter(double ElementSize)
{
    std::cout << "\nicy::MeshFragment::GenerateIndenter\n";

    deformable = false;

    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    double height = 1;
    double radius = 0.15;
//    int point1 = gmsh::model::occ::addPoint(0, height+radius*1.1, 0, 1.0);

    gmsh::model::occ::addEllipse(0, height+radius*1.1, 0, radius, radius/2);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize);

    gmsh::model::mesh::generate(1);

    // retrieve the result
    elems.clear();
    nodes.clear();
    boundary_edges.clear();
    freeNodeCount = 0;


    // nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    gmsh::model::mesh::getNodesByElementType(1, nodeTags, nodeCoords, parametricCoords);

    // boundary
    std::vector<std::size_t> edgeTags, nodeTagsInEdges;
    gmsh::model::mesh::getElementsByType(1, edgeTags, nodeTagsInEdges);

    // nodeTags, nodeCoords, nodeTagsInTris => nodes, elems
    std::map<std::size_t, int> nodeTagsMap1; // nodeTag -> its sequential position in nodeTag
    for(std::size_t i=0;i<nodeTags.size();i++) nodeTagsMap1[nodeTags[i]] = i;

    // set the size of the resulting nodes array
    nodes.resize(nodeTags.size());

    std::map<std::size_t, int> mtags; // nodeTag -> sequential position
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        mtags[tag] = i;

        int idx1 = nodeTagsMap1[tag];
        double x = nodeCoords[idx1*3+0];
        double y = nodeCoords[idx1*3+1];

        icy::Node &nd = nodes[i];
        nd.Reset();
        nd.locId = i;
        nd.x_initial << x, y;
        nd.intended_position = nd.xt = nd.xn = nd.x_initial;
        if(y==0 || !deformable) nd.pinned=true;
        else freeNodeCount++;
    }

    boundary_edges.resize(nodeTagsInEdges.size()/2);

    for(std::size_t i=0;i<nodeTagsInEdges.size()/2;i++)
    {
        int idx1 = mtags[nodeTagsInEdges[i*2+0]];
        int idx2 = mtags[nodeTagsInEdges[i*2+1]];
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

    std::vector<std::pair<int,int>> dimTags;
    gmsh::model::getEntities(dimTags,-1);
    std::cout << "PRINTING entities FOR THE INDENTER\n";
    for(std::pair<int,int> &val : dimTags)
        std::cout << "DIM " << val.first << "; TAG " << val.second << std::endl;

    std::vector<std::pair<int, int> > boundary_dimtags;
    gmsh::model::getBoundary({{1, 1}}, boundary_dimtags, false, false);
    std::vector<int> boundary_tags, complement_tags;
    std::cout << "PRINTING BOUNDARY_DIMTAGS FOR {1,1}\n";
      for(auto e : boundary_dimtags) {
          std::cout << "DIM " << e.first << "; TAG " << e.second << std::endl;
      }

    // PRINT OUT

      // Get all the elementary entities in the model, as a vector of (dimension,
        // tag) pairs:
        std::vector<std::pair<int, int> > entities;
        gmsh::model::getEntities(entities);

        for(auto e : entities) {
          // Dimension and tag of the entity:
          int dim = e.first, tag = e.second;

          // Mesh data is made of `elements' (points, lines, triangles, ...), defined
          // by an ordered list of their `nodes'. Elements and nodes are identified by
          // `tags' as well (strictly positive identification numbers), and are stored
          // ("classified") in the model entity they discretize. Tags for elements and
          // nodes are globally unique (and not only per dimension, like entities).

          // A model entity of dimension 0 (a geometrical point) will contain a mesh
          // element of type point, as well as a mesh node. A model curve will contain
          // line elements as well as its interior nodes, while its boundary nodes
          // will be stored in the bounding model points. A model surface will contain
          // triangular and/or quadrangular elements and all the nodes not classified
          // on its boundary or on its embedded entities. A model volume will contain
          // tetrahedra, hexahedra, etc. and all the nodes not classified on its
          // boundary or on its embedded entities.

          // Get the mesh nodes for the entity (dim, tag):
          std::vector<std::size_t> nodeTags;
          std::vector<double> nodeCoords, nodeParams;
          gmsh::model::mesh::getNodes(nodeTags, nodeCoords, nodeParams, dim, tag);

          // Get the mesh elements for the entity (dim, tag):
          std::vector<int> elemTypes;
          std::vector<std::vector<std::size_t> > elemTags, elemNodeTags;
          gmsh::model::mesh::getElements(elemTypes, elemTags, elemNodeTags, dim, tag);

          // Elements can also be obtained by type, by using `getElementTypes()'
          // followed by `getElementsByType()'.

          // Let's print a summary of the information available on the entity and its
          // mesh.

          // * Type of the entity:
          std::string type;
          gmsh::model::getType(dim, tag, type);
          std::string name;
          gmsh::model::getEntityName(dim, tag, name);
          if(name.size()) name += " ";
          std::cout << "Entity " << name << "(" << dim << "," << tag << ") of type "
                    << type << "\n";

          // * Number of mesh nodes and elements:
          int numElem = 0;
          for(auto &tags : elemTags) numElem += tags.size();
          std::cout << " - Mesh has " << nodeTags.size() << " nodes and " << numElem
                    << " elements\n";
      /*
          // * Upward and downward adjacencies:
          std::vector<int> up, down;
          gmsh::model::getAdjacencies(dim, tag, up, down);
          if(up.size()) {
            std::cout << " - Upward adjacencies: ";
            for(auto e : up) std::cout << e << " ";
            std::cout << "\n";
          }
          if(down.size()) {
            std::cout << " - Downward adjacencies: ";
            for(auto e : down) std::cout << e << " ";
            std::cout << "\n";
          }
      */
          // * Does the entity belong to physical groups?
          std::vector<int> physicalTags;
          gmsh::model::getPhysicalGroupsForEntity(dim, tag, physicalTags);
          if(physicalTags.size()) {
            std::cout << " - Physical group: ";
            for(auto physTag : physicalTags) {
              std::string n;
              gmsh::model::getPhysicalName(dim, physTag, n);
              if(n.size()) n += " ";
              std::cout << n << "(" << dim << ", " << physTag << ") ";
            }
            std::cout << "\n";
          }

          // * Is the entity a partition entity? If so, what is its parent entity?
          std::vector<int> partitions;
          gmsh::model::getPartitions(dim, tag, partitions);
          if(partitions.size()) {
            std::cout << " - Partition tags:";
            for(auto part : partitions) std::cout << " " << part;
            int parentDim, parentTag;
            gmsh::model::getParent(dim, tag, parentDim, parentTag);
            std::cout << " - parent entity (" << parentDim << "," << parentTag
                      << ")\n";
          }

          // * List all types of elements making up the mesh of the entity:
          for(auto elemType : elemTypes) {
            std::string name;
            int d, order, numv, numpv;
            std::vector<double> param;
            gmsh::model::mesh::getElementProperties(elemType, name, d, order, numv,
                                                    param, numpv);
            std::cout << " - Element type: " << name << ", order " << order << "\n";
            std::cout << "   with " << numv << " nodes in param coord: (";
            for(auto p : param) std::cout << p << " ";
            std::cout << ")\n";
          }
        }


    gmsh::clear();
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
    gmsh::option::setNumber("Mesh.MeshSizeMax", ElementSize*10);

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
//    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary",0);
    gmsh::model::mesh::generate(deformable ? 2 : 1);


    // retrieve the result
    elems.clear();
    nodes.clear();
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
