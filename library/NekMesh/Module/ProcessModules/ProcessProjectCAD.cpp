////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessProjectCAD.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
////////////////////////////////////////////////////////////////////////////////

#include "ProcessProjectCAD.h"
#include <NekMesh/MeshElements/Element.h>

#include <NekMesh/CADSystem/CADCurve.h>

#include <LibUtilities/Foundations/ManagerAccess.h>

#include <boost/algorithm/string.hpp>

using namespace std;
using namespace Nektar::NekMesh;

namespace Nektar::NekMesh
{

const NekDouble prismU1[6] = {-1.0, 1.0, 1.0, -1.0, -1.0, -1.0};
const NekDouble prismV1[6] = {-1.0, -1.0, 1.0, 1.0, -1.0, 1.0};
const NekDouble prismW1[6] = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0};

const int pyramidV0[4] = {0, 1, 2, 3};
const int pyramidV1[4] = {1, 2, 3, 0};
const int pyramidV2[4] = {3, 0, 1, 2};
const int pyramidV3[4] = {4, 4, 4, 4};

ModuleKey ProcessProjectCAD::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "projectcad"), ProcessProjectCAD::create,
        "Projects mesh to CAD");

ProcessProjectCAD::ProcessProjectCAD(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["file"]  = ConfigOption(false, "", "CAD file");
    m_config["order"] = ConfigOption(false, "4", "Enforce a polynomial order");
    m_config["surfopti"] = ConfigOption(false, "1", "Run HO-Surface Module");
    m_config["varopti"] =
        ConfigOption(false, "0", "Run the Variational Optmiser");
    m_config["cLength"] =
        ConfigOption(false, "1", "Characteristic Length CAD-Reconstruction");
    m_config["tolv1"] =
        ConfigOption(false, "1e-7",
                     "(Optional) min distance of initial Vertex to CADSurface");
    m_config["tolv2"] = ConfigOption(
        false, "1e-5",
        " (Optional) max distance of initial Vertex to CADSurface");
}

ProcessProjectCAD::~ProcessProjectCAD()
{
}

bool ProcessProjectCAD::FindAndProject(
    bgi::rtree<boxI, bgi::quadratic<16>> &rtree, std::array<NekDouble, 3> &in,
    [[maybe_unused]] int &surf)
{
    point q(in[0], in[1], in[2]);
    vector<boxI> result;
    rtree.query(bgi::intersects(q), back_inserter(result));

    if (result.size() == 0)
    {
        // along a projecting edge the node is too far from any surface boxes
        // this is hardly surprising but rare, return false and linearise
        m_log(VERBOSE) << "FindAndProject() - Edge node not in any of the "
                          "Bounding boxes - No projection! "
                       << endl;
        return false;
    }
    else
    {
        cout << "FindAndProject() - Edge node found in " << result.size()
             << " bounding boxes" << endl;
    }

    int minsurf       = 0;
    NekDouble minDist = numeric_limits<double>::max();
    NekDouble dist;

    for (int j = 0; j < result.size(); j++)
    {
        m_mesh->m_cad->GetSurf(result[j].second)->locuv(in, dist);

        if (dist < minDist)
        {
            minDist = dist;
            minsurf = result[j].second;
        }
    }

    auto uv = m_mesh->m_cad->GetSurf(minsurf)->locuv(in, dist);

    in = m_mesh->m_cad->GetSurf(minsurf)->P(uv);

    return true;
}

bool ProcessProjectCAD::IsNotValid(vector<ElementSharedPtr> &els)
{
    // short algebraic method to figure out the vailidy of elements
    // test the volume of tetrahedrons constituted by three sibling edges
    for (int i = 0; i < els.size(); i++)
    {
        if (els[i]->GetShapeType() == LibUtilities::ePrism)
        {
            vector<NodeSharedPtr> ns = els[i]->GetVertexList();
            for (int j = 0; j < 6; j++)
            {
                NekDouble a2 = 0.5 * (1 + prismU1[j]);
                NekDouble b1 = 0.5 * (1 - prismV1[j]);
                NekDouble b2 = 0.5 * (1 + prismV1[j]);
                NekDouble c2 = 0.5 * (1 + prismW1[j]);
                NekDouble d  = 0.5 * (prismU1[j] + prismW1[j]);

                std::array<NekDouble, 9> jac;

                jac[0] = -0.5 * b1 * ns[0]->m_x + 0.5 * b1 * ns[1]->m_x +
                         0.5 * b2 * ns[2]->m_x - 0.5 * b2 * ns[3]->m_x;
                jac[1] = -0.5 * b1 * ns[0]->m_y + 0.5 * b1 * ns[1]->m_y +
                         0.5 * b2 * ns[2]->m_y - 0.5 * b2 * ns[3]->m_y;
                jac[2] = -0.5 * b1 * ns[0]->m_z + 0.5 * b1 * ns[1]->m_z +
                         0.5 * b2 * ns[2]->m_z - 0.5 * b2 * ns[3]->m_z;

                jac[3] = 0.5 * d * ns[0]->m_x - 0.5 * a2 * ns[1]->m_x +
                         0.5 * a2 * ns[2]->m_x - 0.5 * d * ns[3]->m_x -
                         0.5 * c2 * ns[4]->m_x + 0.5 * c2 * ns[5]->m_x;
                jac[4] = 0.5 * d * ns[0]->m_y - 0.5 * a2 * ns[1]->m_y +
                         0.5 * a2 * ns[2]->m_y - 0.5 * d * ns[3]->m_y -
                         0.5 * c2 * ns[4]->m_y + 0.5 * c2 * ns[5]->m_y;
                jac[5] = 0.5 * d * ns[0]->m_z - 0.5 * a2 * ns[1]->m_z +
                         0.5 * a2 * ns[2]->m_z - 0.5 * d * ns[3]->m_z -
                         0.5 * c2 * ns[4]->m_z + 0.5 * c2 * ns[5]->m_z;

                jac[6] = -0.5 * b1 * ns[0]->m_x - 0.5 * b2 * ns[3]->m_x +
                         0.5 * b1 * ns[4]->m_x + 0.5 * b2 * ns[5]->m_x;
                jac[7] = -0.5 * b1 * ns[0]->m_y - 0.5 * b2 * ns[3]->m_y +
                         0.5 * b1 * ns[4]->m_y + 0.5 * b2 * ns[5]->m_y;
                jac[8] = -0.5 * b1 * ns[0]->m_z - 0.5 * b2 * ns[3]->m_z +
                         0.5 * b1 * ns[4]->m_z + 0.5 * b2 * ns[5]->m_z;

                NekDouble jc = jac[0] * (jac[4] * jac[8] - jac[5] * jac[7]) -
                               jac[3] * (jac[1] * jac[8] - jac[2] * jac[7]) +
                               jac[6] * (jac[1] * jac[5] - jac[2] * jac[4]);

                if (jc < NekConstants::kNekZeroTol)
                {
                    return true;
                }
            }
        }
        else if (els[i]->GetShapeType() == LibUtilities::ePyramid)
        {
            vector<NodeSharedPtr> ns = els[i]->GetVertexList();
            for (int j = 0; j < 4; j++)
            {
                std::array<NekDouble, 9> jac;

                jac[0] = 0.5 * (ns[pyramidV1[j]]->m_x - ns[pyramidV0[j]]->m_x);
                jac[1] = 0.5 * (ns[pyramidV1[j]]->m_y - ns[pyramidV0[j]]->m_y);
                jac[2] = 0.5 * (ns[pyramidV1[j]]->m_z - ns[pyramidV0[j]]->m_z);
                jac[3] = 0.5 * (ns[pyramidV2[j]]->m_x - ns[pyramidV0[j]]->m_x);
                jac[4] = 0.5 * (ns[pyramidV2[j]]->m_y - ns[pyramidV0[j]]->m_y);
                jac[5] = 0.5 * (ns[pyramidV2[j]]->m_z - ns[pyramidV0[j]]->m_z);
                jac[6] = 0.5 * (ns[pyramidV3[j]]->m_x - ns[pyramidV0[j]]->m_x);
                jac[7] = 0.5 * (ns[pyramidV3[j]]->m_y - ns[pyramidV0[j]]->m_y);
                jac[8] = 0.5 * (ns[pyramidV3[j]]->m_z - ns[pyramidV0[j]]->m_z);

                NekDouble jc = jac[0] * (jac[4] * jac[8] - jac[5] * jac[7]) -
                               jac[3] * (jac[1] * jac[8] - jac[2] * jac[7]) +
                               jac[6] * (jac[1] * jac[5] - jac[2] * jac[4]);

                if (jc < NekConstants::kNekZeroTol)
                {
                    return true;
                }
            }
        }
        else if (els[i]->GetShapeType() == LibUtilities::eTetrahedron)
        {
            vector<NodeSharedPtr> ns = els[i]->GetVertexList();
            std::array<NekDouble, 9> jac;

            jac[0] = 0.5 * (ns[1]->m_x - ns[0]->m_x);
            jac[1] = 0.5 * (ns[1]->m_y - ns[0]->m_y);
            jac[2] = 0.5 * (ns[1]->m_z - ns[0]->m_z);
            jac[3] = 0.5 * (ns[2]->m_x - ns[0]->m_x);
            jac[4] = 0.5 * (ns[2]->m_y - ns[0]->m_y);
            jac[5] = 0.5 * (ns[2]->m_z - ns[0]->m_z);
            jac[6] = 0.5 * (ns[3]->m_x - ns[0]->m_x);
            jac[7] = 0.5 * (ns[3]->m_y - ns[0]->m_y);
            jac[8] = 0.5 * (ns[3]->m_z - ns[0]->m_z);

            NekDouble jc = jac[0] * (jac[4] * jac[8] - jac[5] * jac[7]) -
                           jac[3] * (jac[1] * jac[8] - jac[2] * jac[7]) +
                           jac[6] * (jac[1] * jac[5] - jac[2] * jac[4]);

            if (jc < NekConstants::kNekZeroTol)
            {
                return true;
            }
        }
        else if (els[i]->GetShapeType() == LibUtilities::eHexahedron)
        {
            // Hexes are not checked !
            NekDouble jc = 1.0;
            if (jc < NekConstants::kNekZeroTol)
            {
                return true;
            }
        }
        else
        {
            m_log(FATAL) << "Only prisms, pyramids and tetrahedra supported."
                         << endl;
        }
    }

    return false;
}

void ProcessProjectCAD::Process()
{
    m_log(VERBOSE) << "Projecting CAD back onto linear mesh." << endl;

    m_log(WARNING) << "ProcessAssignCAD: Warning: This module is designed for "
                   << "use with Star-CCM+ meshes only; it also requires that "
                   << "the star mesh was created in a certain way." << endl;

    if (!m_config["order"].beenSet)
    {
        m_log(VERBOSE) << "Mesh order not set: will assume order 4" << endl;
    }

    // Projection Order
    int order         = m_config["order"].as<int>();
    m_mesh->m_nummode = order + 1;

    // Tolerances for vertex association
    NekDouble tolv1, tolv2;
    tolv1 = m_config["tolv1"].as<NekDouble>();
    tolv2 = m_config["tolv2"].as<NekDouble>();

    // Characteristic Length for CAD Reconstruction (auto calculation of tolv1
    // and tolv2)
    NekDouble cLength = 1.0;
    if (m_config["cLength"].beenSet)
    {
        cLength = m_config["cLength"].as<NekDouble>();
        tolv1 *= cLength;
        tolv2 *= cLength;
    }

    // 1. Load CAD instance of the CAD model
    std::string filename = m_config["file"].as<string>();
    LoadCAD(filename);

    // 2. Create Bounding boxes of the CAD surfaces into a k-d tree
    NekDouble scale = cLength;
    bgi::rtree<boxI, bgi::quadratic<16>> rtree;
    bgi::rtree<boxI, bgi::quadratic<16>> rtreeCurve;
    bgi::rtree<boxI, bgi::quadratic<16>> rtreeNode;
    CreateBoundingBoxes(rtree, rtreeCurve, rtreeNode, scale);

    m_log(VERBOSE) << "Bounding Boxes Surf/Curv/Vertex= " << rtree.size() << " "
                   << rtreeCurve.size() << " " << rtreeNode.size() << endl;
    // 3. Auxilaries ( can be moved to Module.cpp)
    // SurfNodes , surfNodeToEl, minConEdge
    Auxilaries();

    // 4.  Link Surface Vertices to CAD and Project them to the closest CAD
    LinkVertexToCAD(m_mesh, false, lockedNodes, tolv1, tolv2, rtree, rtreeCurve,
                    rtreeNode);

    // 5. clear the associations with CAD surfaces
    // necessary since the projection of the edges will change the surface uv
    // and the association will be wrong to some surfaces that were closed
    // beforehand
    for (auto vertex = surfNodes.begin(); vertex != surfNodes.end(); vertex++)
    {
        (*vertex)->ClearCADSurfs();
    }

    // 6. Update the secondary tolerances on already projected nodes and do the
    //  final Linking Vertex - CAD Surface / Curve
    tolv1 = 1e-9, tolv2 = 1e-8;
    LinkVertexToCAD(m_mesh, true, lockedNodes, tolv1, tolv2, rtree, rtreeCurve,
                    rtreeNode);

    // // 7. Associate Edges to CAD
    LinkEdgeToCAD(surfEdges, tolv1);

    // // 8. Associate Faces to CAD
    LinkFaceToCAD();

    // Project the Edges to CAD that
    // ProjectEdges(surfEdges, order, rtree);

    Diagnostics();

    ////**** HOSurface ****////
    int m_surfopti         = m_config["surfopti"].as<bool>();
    ModuleSharedPtr module = GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "hosurface"), m_mesh);
    module->SetLogger(m_log);
    if (m_surfopti == 1)
    {
        //        module->RegisterConfig("opti","");
        module->RegisterConfig("third_party", "");

        try
        {
            module->SetDefaults();
            module->Process();
        }
        catch (runtime_error &e)
        {
            m_log(WARNING)
                << "High-order surface meshing has failed with message:"
                << endl;
            m_log(WARNING) << e.what() << endl;
            m_log(WARNING) << "The mesh will be written as normal but the "
                           << "incomplete surface will remain faceted" << endl;
            return;
        }
    }

    Diagnostics();
    // ExportCAD();

    m_log(VERBOSE) << "HO-Surface CAD complete." << endl;
}

void ProcessProjectCAD::LoadCAD(std::string filename)
{
    m_log(VERBOSE) << "Start Loading the CAD file " << endl;
    ModuleSharedPtr module = GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "loadcad"), m_mesh);
    module->RegisterConfig("filename", filename);
    module->SetDefaults();
    module->Process();
    m_log(VERBOSE) << "CAD loaded succesfully!" << endl;
}

void ProcessProjectCAD::Auxilaries()
{
    // find nodes on the surface
    // find unique nodes on the surface
    for (int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        ElementSharedPtr el      = m_mesh->m_element[2][i];
        vector<NodeSharedPtr> ns = el->GetVertexList();
        for (int j = 0; j < ns.size(); j++)
        {
            surfNodes.insert(ns[j]);
        }
    }

    // link surface nodes to their 3D element
    for (int i = 0; i < m_mesh->m_element[3].size(); i++)
    {
        if (m_mesh->m_element[3][i]->HasBoundaryLinks())
        {
            vector<NodeSharedPtr> ns = m_mesh->m_element[3][i]->GetVertexList();
            for (int j = 0; j < ns.size(); j++)
            {
                if (surfNodes.count(ns[j]) > 0)
                {
                    surfNodeToEl[ns[j]].push_back(m_mesh->m_element[3][i]);
                }
            }
        }
    }

    // Calculate the min edge length of the elements
    // Calculate MinConEdge
    CalculateMinEdgeLength();

    // make edges of surface mesh unique
    ClearElementLinks();
    // EdgeSet surfEdges;
    vector<ElementSharedPtr> &elmt = m_mesh->m_element[2];
    map<int, int> surfIdToLoc;
    for (int i = 0; i < elmt.size(); i++)
    {
        surfIdToLoc.insert(make_pair(elmt[i]->GetId(), i));
        for (int j = 0; j < elmt[i]->GetEdgeCount(); ++j)
        {
            pair<EdgeSet::iterator, bool> testIns;
            EdgeSharedPtr ed = elmt[i]->GetEdge(j);
            testIns          = surfEdges.insert(ed);

            if (testIns.second)
            {
                EdgeSharedPtr ed2 = *testIns.first;
                ed2->m_elLink.push_back(
                    pair<ElementSharedPtr, int>(elmt[i], j));
            }
            else
            {
                EdgeSharedPtr e2 = *(testIns.first);
                elmt[i]->SetEdge(j, e2);

                // Update edge to element map.
                e2->m_elLink.push_back(pair<ElementSharedPtr, int>(elmt[i], j));
            }
        }
    }
}

void ProcessProjectCAD::CreateBoundingBoxes(
    bgi::rtree<boxI, bgi::quadratic<16>> &rtreeSurf,
    bgi::rtree<boxI, bgi::quadratic<16>> &rtreeCurve,
    bgi::rtree<boxI, bgi::quadratic<16>> &rtreeNode, NekDouble scale)
{
    // CAD Surfs
    vector<boxI> boxes;
    for (int i = 1; i <= m_mesh->m_cad->GetNumSurf(); i++)
    {
        // m_log(VERBOSE).Progress(i, m_mesh->m_cad->GetNumSurf(),
        //                         "building surface bboxes", i - 1);
        auto bx = m_mesh->m_cad->GetSurf(i)->BoundingBox(scale);
        boxes.push_back(make_pair(
            box(point(bx[0], bx[1], bx[2]), point(bx[3], bx[4], bx[5])), i));

        m_log(VERBOSE) << " Boundaing box = " << i << endl;
        m_log(VERBOSE) << bx[0] << " " << bx[1] << " " << bx[2] << endl;
        m_log(VERBOSE) << bx[3] << " " << bx[4] << " " << bx[5] << endl;
    }

    m_log(VERBOSE).Newline();
    m_log(VERBOSE) << "Building Surf admin data structures." << endl;
    rtreeSurf.insert(boxes.begin(), boxes.end());
    m_log(VERBOSE) << "Bounding Box ." << endl;

    // CAD Curves
    boxes.clear();
    for (int i = 1; i <= m_mesh->m_cad->GetNumCurve(); i++)
    {
        auto bx = m_mesh->m_cad->GetCurve(i)->BoundingBox(scale);

        boxes.push_back(make_pair(
            box(point(bx[0], bx[1], bx[2]), point(bx[3], bx[4], bx[5])), i));
    }
    rtreeCurve.insert(boxes.begin(), boxes.end());

    // CAD Vertices
    boxes.clear();
    NekDouble tol = 1e-8 * scale; // 1e-8 * scale
    for (int i = 1; i <= m_mesh->m_cad->GetNumVerts(); i++)
    {
        auto vert                     = m_mesh->m_cad->GetVert(i);
        std::array<NekDouble, 3> locT = vert->GetLoc();
        boxes.push_back(
            make_pair(box(point(locT[0] - tol, locT[1] - tol, locT[2] - tol),
                          point(locT[0] + tol, locT[1] + tol, locT[2] + tol)),
                      i));
    }
    rtreeNode.insert(boxes.begin(), boxes.end());
}

void ProcessProjectCAD::CalculateMinEdgeLength()
{
    // link the surface node to a value for the shortest connecting edge to it
    for (int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        ElementSharedPtr el      = m_mesh->m_element[2][i];
        vector<NodeSharedPtr> ns = el->GetVertexList();
        NekDouble l1             = ns[0]->Distance(ns[1]);
        NekDouble l2             = ns[1]->Distance(ns[2]);
        NekDouble l3             = ns[2]->Distance(ns[0]);

        if (minConEdge.count(ns[0]))
        {
            NekDouble l       = minConEdge[ns[0]];
            minConEdge[ns[0]] = min(l, min(l1, l3));
        }
        else
        {
            minConEdge.insert(make_pair(ns[0], min(l1, l3)));
        }
        if (minConEdge.count(ns[1]))
        {
            NekDouble l       = minConEdge[ns[1]];
            minConEdge[ns[1]] = min(l, min(l1, l1));
        }
        else
        {
            minConEdge.insert(make_pair(ns[1], min(l1, l1)));
        }
        if (minConEdge.count(ns[2]))
        {
            NekDouble l       = minConEdge[ns[2]];
            minConEdge[ns[2]] = min(l, min(l2, l3));
        }
        else
        {
            minConEdge.insert(make_pair(ns[2], min(l2, l3)));
        }
    }
}

void ProcessProjectCAD::LinkVertexToCAD(
    NekMesh::MeshSharedPtr &m_mesh, bool CADCurve, NodeSet &lockedNodes,
    NekDouble tolv1, NekDouble tolv2,
    bgi::rtree<boxI, bgi::quadratic<16>> &rtree,
    bgi::rtree<boxI, bgi::quadratic<16>> &rtreeCurve,
    bgi::rtree<boxI, bgi::quadratic<16>> &rtreeNode)
{
    map<int, vector<int>> finds;

    m_log(VERBOSE) << "Searching tree." << endl;

    NekDouble maxNodeCor = 0;

    // find nodes surface and parametric location
    int ct = 0;
    for (auto i = surfNodes.begin(); i != surfNodes.end(); i++, ct++)
    {
        // m_log(VERBOSE).Progress(ct, surfNodes.size(), "projecting verts",
        //                         ct - 1);

        point q((*i)->m_x, (*i)->m_y, (*i)->m_z);
        vector<boxI> result;
        rtree.query(bgi::intersects(q), back_inserter(result));

        if (result.size() == 0)
        {
            // Vertex is too far from any surface bounding boxes
            m_log(WARNING)
                << "Vertex  " << (*i)
                << " is not in any boundin boxes. 1. "
                   "Make sure you use the correct STEP file. 2. The problem is "
                   "likely in the linear mesh -> try to refine this region.  "
                << endl;
            continue;
        }

        // Vertex Tolerances to the CAD Surface - 0.5 * min edge length + [min
        // tol max tol]
        NekDouble tol = minConEdge[*i] * 0.5;
        tol           = min(tol, tolv1);
        tol           = max(tol, tolv2);

        vector<int> distId;
        vector<NekDouble> distList;
        // sort the surfaces by distance to the node
        for (int j = 0; j < result.size(); j++)
        {
            NekDouble dist;
            m_mesh->m_cad->GetSurf(result[j].second)
                ->locuv((*i)->GetLoc(), dist);
            distList.push_back(dist);
            distId.push_back(result[j].second);
        }

        bool repeat = true;
        while (repeat)
        {
            repeat = false;
            for (int j = 0; j < distId.size() - 1; j++)
            {
                if (distList[j + 1] < distList[j])
                {
                    repeat = true;
                    swap(distList[j + 1], distList[j]);
                    swap(distId[j + 1], distId[j]);
                }
            }
        }

        int pos = 0;
        for (int j = 0; j < distId.size(); j++)
        {
            if (distList[j] < tol)
            {
                pos++;
            }
        }

        distId.resize(pos);

        finds[pos].push_back(0);
        // if the node is not close to any surface lock it (no CAD given + no
        // projection for its edges)
        if (pos == 0)
        {
            lockedNodes.insert(*i);
            m_log(WARNING) << "surface " << distList[0] << " unknown "
                           << "(tolerance: " << tol << ")" << endl;
        }
        else
        {
            NekDouble shift;
            bool st = false;
            for (int j = 0; j < distId.size(); j++)
            {
                if (distList[j] > tol)
                {
                    continue;
                }
                if (m_mesh->m_cad->GetSurf(distId[j])->IsPlanar())
                {
                    // continue;
                }

                shift                     = distList[j];
                NekDouble dist            = 0;
                CADSurfSharedPtr s        = m_mesh->m_cad->GetSurf(distId[j]);
                auto l                    = (*i)->GetLoc();
                [[maybe_unused]] auto uvt = s->locuv(l, dist);

                NekDouble tmpX = (*i)->m_x;
                NekDouble tmpY = (*i)->m_y;
                NekDouble tmpZ = (*i)->m_z;

                (*i)->m_x = s->P(uvt)[0];
                (*i)->m_y = s->P(uvt)[1];
                (*i)->m_z = s->P(uvt)[2];

                if (ProcessProjectCAD::IsNotValid(surfNodeToEl[*i]))
                {
                    (*i)->m_x = tmpX;
                    (*i)->m_y = tmpY;
                    (*i)->m_z = tmpZ;

                    m_log(VERBOSE) << "Element not valid after vertex ";
                    m_log(VERBOSE)
                        << "projection reset it and lock the vertex" << endl;
                    break;
                }

                st = true;
                break;
            }

            if (!st)
            {
                lockedNodes.insert(*i);
                continue;
            }

            for (int j = 0; j < distId.size(); j++)
            {
                if (distList[j] > tol)
                {
                    continue;
                }
                if (m_mesh->m_cad->GetSurf(distId[j])->IsPlanar())
                {
                    // continue;
                }

                CADSurfSharedPtr s = m_mesh->m_cad->GetSurf(distId[j]);
                NekDouble dist     = 0;
                auto loc           = (*i)->GetLoc();
                auto uv            = s->locuv(loc, dist);
                (*i)->SetCADSurf(s, uv);
            }
            maxNodeCor = max(maxNodeCor, shift);
        }
    }

    // for the vertices with multiple CAD Surfaces, check CADCurves with the
    // bounding boxes Intersect with CAD Curves rtree
    if (CADCurve)
    {
        for (auto vertex : surfNodes)
        {
            if (vertex->GetCADSurfs().size() > 2)
            {
                point point((vertex)->m_x, (vertex)->m_y, (vertex)->m_z);
                vector<boxI> result;
                rtreeNode.query(bgi::intersects(point), back_inserter(result));

                if (result.size() == 1)
                {
                    // Single CAD Vertex
                    vertex->SetCADVertex(m_mesh->m_cad->GetVert(
                        result[0].second)); // No Adj Curves Assigned to the V !
                }
                else
                {
                    // Multiple CAD Vertices
                    m_log(WARNING) << "Multiple CAD Vertices found for the "
                                   << "vertex " << vertex << endl;
                }
            }
            else if(vertex->GetCADSurfs().size() == 2)
            {
                point point((vertex)->m_x, (vertex)->m_y, (vertex)->m_z);
                vector<boxI> result;
                rtreeCurve.query(bgi::intersects(point), back_inserter(result));

                if(result.size()==1)
                {
                    // Single CAD Curve
                    CADCurveSharedPtr CADCurve_t = m_mesh->m_cad->GetCurve(result[0].second);
                    NekDouble t0, t1, dist0, dist1;
                    NekDouble tmin, tmax;
                    CADCurve_t->GetBounds(tmin, tmax);
    
                    dist0 = CADCurve_t->loct(vertex->GetLoc(), t0, tmin, tmax);
                    if(dist0 < tolv1)
                    {
                        vertex->SetCADCurve(CADCurve_t, t0 ); 
                    }
                }
                else
                {
                    m_log(WARNING) << "Multiple CAD Curves found for the "
                                   << "vertex " << vertex << endl;
                }
                
            }
        }
    }

    m_log(VERBOSE) << "  - max surface Vertex correction " << maxNodeCor
                   << endl;
    m_log(VERBOSE) << "  - lockedNodes N= " << lockedNodes.size() << endl;
}

void ProcessProjectCAD::LinkEdgeToCAD(EdgeSet &surfEdges, NekDouble tolv1)
{
    for (auto edge : surfEdges)
    {
        NodeSharedPtr v1 = edge->m_n1;
        NodeSharedPtr v2 = edge->m_n2;

        if (lockedNodes.count(v1) || lockedNodes.count(v2))
        {
            continue;
        }

        vector<int> cmn =
            IntersectCADSurf(v1->GetCADSurfs(), v2->GetCADSurfs());

        if (cmn.size() == 0)
        {
            // no CAD surface found for the edge (CASE3)
            // cout << "Case 3 edge association v1 = " << v1 << "  = v2 " << v2
            //<< endl;
            continue;
        }

        if (cmn.size() == 1)
        {
            // Clearly Edge is on a single CAD surface (internal)
            edge->m_parentCAD = m_mesh->m_cad->GetSurf(cmn[0]);
        }
        else if (cmn.size() == 2)
        {
            // Could be CAD-curve or CAD-surface (CASE2)
            vector<CADSurfSharedPtr> v1_CAD = v1->GetCADSurfs();
            vector<CADSurfSharedPtr> v2_CAD = v2->GetCADSurfs();

            // Associate to the closest CAD Surface

            // 1.Create vi1 , vi2
            vector<int> CADCurves_uv1;
            vector<int> CADCurves_uv2;

            CADSurfSharedPtr EdgeSurf1 = m_mesh->m_cad->GetSurf(cmn[0]);
            CADSurfSharedPtr EdgeSurf2 = m_mesh->m_cad->GetSurf(cmn[1]);

            for (auto EdgeLoop : EdgeSurf1->GetEdges())
            {
                for (auto CADCurve : EdgeLoop->edges)
                {
                    CADCurves_uv1.push_back(CADCurve->GetId());
                }
            }
            for (auto EdgeLoop : EdgeSurf2->GetEdges())
            {
                for (auto CADCurve : EdgeLoop->edges)
                {
                    CADCurves_uv2.push_back(CADCurve->GetId());
                }
            }

            sort(CADCurves_uv1.begin(), CADCurves_uv1.end());
            sort(CADCurves_uv2.begin(), CADCurves_uv2.end());

            vector<int> commonCADCurves;
            set_intersection(CADCurves_uv1.begin(), CADCurves_uv1.end(),
                             CADCurves_uv2.begin(), CADCurves_uv2.end(),
                             back_inserter(commonCADCurves));

            NekDouble tolDist = tolv1; // POSSIBLE PROBLEM !!!!
            if (commonCADCurves.size() == 1)
            {
                CADCurveSharedPtr CADCurve_t =
                    m_mesh->m_cad->GetCurve(commonCADCurves[0]);

                NekDouble t0, t1, dist0, dist1;
                NekDouble tmin, tmax;
                CADCurve_t->GetBounds(tmin, tmax);

                dist0 = CADCurve_t->loct(edge->m_n1->GetLoc(), t0, tmin, tmax);
                dist1 = CADCurve_t->loct(edge->m_n2->GetLoc(), t1, tmin, tmax);

                if ((dist0 < tolDist) && (dist1 < tolDist))
                {
                    edge->m_n1->SetCADCurve(CADCurve_t, t0);
                    edge->m_n2->SetCADCurve(CADCurve_t, t1);
                    edge->m_parentCAD = CADCurve_t;
                }
                else
                {
                    m_log(VERBOSE)
                        << " dist > distol (cmnCADCurve.size=1) = " << dist0
                        << " " << dist1 << endl;
                }
                /*
                                        if((dist0 < tolDist) && (dist1 <
                    tolDist) )
                                        {
                                            //edge->m_n1->SetCADCurve(CADCurve_t,t0);
                                            edge->m_n2->SetCADCurve(CADCurve_t,t1);

                                            // IF USE GetMinDist -> Does
                    not give the correct distance to the vertex in case t
                    outside of curve [t_min t_max] bonds

                                            // vertex 1
                                            if(t0 >= tmin - tolDist &&
                    t0 <= tmax + tolDist )
                                            {
                                                edge->m_n1->SetCADCurve(CADCurve_t,t0);
                                            }
                                            else if( (t0 > tmin) && (t0
                    < tmax + tolDist))
                                            {
                                                // this is necessary to
                    accomodate for Vertices generated in StarCCM which
                    are slightly further away
                                                // Alternatively one can
                    project all vertices on the m_log(WARNING) <<
                    "Parametric m_n1 on CADCurve has been manually to
                    curve t_max " <<  endl ;
                                                edge->m_n1->SetCADCurve(CADCurve_t,tmax);
                                            }
                                            else if( (t0 > tmin-
                    tolDist) && (t0 < tmax ))
                                            {
                                                // this is necessary to
                    accomodate for Vertices generated in StarCCM which
                    are slightly further away
                                                // Alternatively one can
                    project all vertices on the m_log(WARNING) <<
                    "Parametric m_n1 on CADCurve has been manually to
                    curve tmin " <<  endl ;
                                                edge->m_n1->SetCADCurve(CADCurve_t,tmin);
                                            }
                                            else
                                            {
                                                    edge->m_parentCAD =
                    NULL ; m_log(WARNING) <<  "Parametric m_n1 on
                    CADCurve cmn=1 has been too far away from the Curve
                    Bounds " <<  endl ;
                                            }

                                            // // vertex 2
                                            if(t1 >= tmin - tolDist &&
                    t1 <= tmax + tolDist )
                                            {
                                                edge->m_n2->SetCADCurve(CADCurve_t,t1);
                                            }
                                            else if( (t1 > tmin) && (t1
                    < tmax + tolDist))
                                            {
                                                // this is necessary to
                    accomodate for Vertices generated in StarCCM which
                    are slightly further away
                                                // Alternatively one can
                    project all vertices on the m_log(WARNING) <<
                    "Parametric m_n2 on CADCurve has been manually to
                    curve t_max " <<  endl ;
                                                edge->m_n2->SetCADCurve(CADCurve_t,tmax);
                                            }
                                            else if( (t1 > tmin-
                    tolDist) && (t1 < tmax ))
                                            {
                                                // this is necessary to
                    accomodate for Vertices generated in StarCCM which
                    are slightly further away
                                                // Alternatively one can
                    project all vertices on the m_log(WARNING) <<
                    "Parametric m_n2 on CADCurve has been manually to
                    curve tmin " <<  endl ;
                                                edge->m_n2->SetCADCurve(CADCurve_t,tmin);
                                            }
                                            else
                                            {
                                                edge->m_parentCAD = NULL
                    ; m_log(WARNING) <<  "Parametric m_n2 on CADCurve
                    cmn=1 has been too far away from the Curve Bounds "
                    <<  endl ;
                                            }
                                        }
                                        else
                                        {
                                            cout << "dist of m_n1 to the
                    2 CAD Curves " << dist0 << " " << dist1 << endl ;
                                            m_log(WARNING) <<
                    "CommonCADCurves =1 - The edge vertices are too far
                    away from the CADCurve dist = " << dist0 << " " <<
                    dist1 << endl ;
                                        }
                */
            }
            else if (commonCADCurves.size() == 2)
            {
                m_log(WARNING)
                    << " CASE commonCADCurves.size()==2 for comn.size()  = 2 "
                    << endl;
                // // << endl ;
                // Array<OneD, NekDouble> xyz1 = edge->m_n1->GetLoc();
                // Array<OneD, NekDouble> xyz2 = edge->m_n2->GetLoc();

                // CADCurveSharedPtr Curve_t0 =
                //     m_mesh->m_cad->GetCurve(commonCADCurves[0]);
                // CADCurveSharedPtr Curve_t1 =
                //     m_mesh->m_cad->GetCurve(commonCADCurves[1]);

                // NekDouble dist01, dist02, dist11, dist12, m_n1_t0,
                //     m_n1_t1, m_n2_t0, m_n2_t1;
                // // cout << "Initialization m_n1/2 t0/1 " << m_n1_t0 << "
                // // " << m_n2_t0 << " " <<  m_n1_t1 <<" " <<  m_n2_t1 <<
                // // endl;

                // // NekDouble tmin0, tmax0, tmin1, tmax1;
                // Array<OneD, NekDouble> tmin_max0 =
                //     Curve_t0->GetBounds();
                // Array<OneD, NekDouble> tmin_max1 =
                //     Curve_t1->GetBounds();
                // // cout << "tmin_max0 =  " << tmin_max0[0] << " " <<
                // // tmin_max0[1] << endl ; cout << "tmin_max1 =  " <<
                // // tmin_max1[0] << " " << tmin_max1[1] << endl ;

                // dist01 = Curve_t0->loct(xyz1, m_n1_t0, tmin_max0[0],
                //                         tmin_max0[1]);
                // dist02 = Curve_t0->loct(xyz2, m_n2_t0, tmin_max0[0],
                //                         tmin_max0[1]);

                // dist11 = Curve_t1->loct(xyz1, m_n1_t1, tmin_max1[0],
                //                         tmin_max1[1]);
                // dist12 = Curve_t1->loct(xyz2, m_n2_t1, tmin_max1[0],
                //                         tmin_max1[1]);

                // // cout << "After Update : " << m_n1_t0 << " " <<
                // // m_n2_t0 << " " << m_n1_t1 << " " << m_n2_t1 << endl ;

                // // Just for V&V
                // // NekDouble dist01_GetMinDistance=
                // // Curve_t0->GetMinDistance(xyz1) ; NekDouble
                // // dist02_GetMinDistance= Curve_t0->GetMinDistance(xyz2)
                // // ; NekDouble dist11_GetMinDistance=
                // // Curve_t1->GetMinDistance(xyz1) ; NekDouble
                // // dist12_GetMinDistance= Curve_t1->GetMinDistance(xyz2)
                // // ;

                // // m_log(VERBOSE) <<" dist 01 "  << dist01 << " dist02 "
                // // << dist02 << endl ; m_log(VERBOSE) << " GetMin "<<  "
                // // dist 01 "  << dist01_GetMinDistance << " dist02 " <<
                // // dist02_GetMinDistance << endl ; m_log(VERBOSE) <<"
                // // dist 11 "  << dist11 << " dist12 " << dist12 << endl
                // // ; m_log(VERBOSE) << " GetMin "<<  " dist 01 "  <<
                // // dist11_GetMinDistance << " dist02 " <<
                // // dist12_GetMinDistance << endl ;

                // // if( (dist1 < tol)  ||  (dist0 < tol)  )
                // NekDouble dist0 = (dist01 + dist02) / 2.0;
                // NekDouble dist1 = (dist11 + dist12) / 2.0;
                // if ((dist0 < tolDist || dist1 < tolDist))
                // {
                //     if (dist0 < dist1 &&
                //         (m_n1_t0 > tmin_max0[0] - tolDist) &&
                //         (m_n1_t0 < tmin_max0[1] + tolDist) &&
                //         (m_n2_t0 > tmin_max0[0] - tolDist) &&
                //         (m_n2_t0 < tmin_max0[1] + tolDist))
                //     {
                //         // m_log(VERBOSE)  << "CAD1" << endl ;
                //         edge->m_parentCAD = Curve_t0;

                //         edge->m_n1->SetCADCurve(Curve_t0, m_n1_t0);
                //         edge->m_n2->SetCADCurve(Curve_t0, m_n2_t0);
                //         // cout << "CADCurve 0 m_n1_t0= " <<
                //         // edge->m_n1->GetCADCurveInfo(Curve_t0->GetId())
                //         // << endl ; cout << "CADCurve 0 m_n2_t0= " <<
                //         // edge->m_n2->GetCADCurveInfo(Curve_t0->GetId())
                //         // << endl ;
                //     }
                //     else if (dist1 < dist0 &&
                //                 (m_n1_t1 > tmin_max1[0] - tolDist) &&
                //                 (m_n1_t1 < tmin_max1[1] + tolDist) &&
                //                 (m_n2_t1 > tmin_max1[0] - tolDist) &&
                //                 (m_n2_t1 < tmin_max1[1] + tolDist))
                //     {
                //         // m_log(VERBOSE) << "CAD2" << endl ;
                //         edge->m_parentCAD = Curve_t1;
                //         edge->m_n1->SetCADCurve(Curve_t1, m_n1_t1);
                //         edge->m_n2->SetCADCurve(Curve_t1, m_n2_t1);
                //         // cout << "CADCurve 1 m_n1_t1= " <<
                //         // edge->m_n1->GetCADCurveInfo(Curve_t1->GetId())
                //         // << endl ; cout << "CADCurve 1 m_n2_t1 " <<
                //         // edge->m_n2->GetCADCurveInfo(Curve_t1->GetId())
                //         // << endl ;
                //     }
                //     else
                //     {
                //         m_log(WARNING)
                //             << "loct circle bug cmn.size()=2 " << endl;
                //     }
                // }
                // else
                // {
                //     m_log(WARNING) << " dist 01 " << dist01
                //                     << " dist02 " << dist02 << endl;
                //     m_log(WARNING) << " dist 11 " << dist11
                //                     << " dist12 " << dist12 << endl;

                //     cout << "dist of m_n1 to the 2 CAD Curves " << dist0
                //             << " " << dist1 << endl;
                //     m_log(WARNING)
                //         << "Cannot find close CADCurve in CornerCase "
                //             "of Edge with 2 common curves"
                //         << endl;
                // }
                // // cout << "end edge commonCADCurves.size()==2 " << endl
                // // ;
            }
            else
            {
                // Associate to the closest CAD Surface
                m_log(WARNING)
                    << "too many common CADcurves for Edge association  "
                       "(cmn>2) will use the element to associate. "
                    << endl;
            }
        }
        else
        {
            // Associate to the closest CAD Surface
            m_log(WARNING) << "too many common surfaces for Edge association  "
                              "(cmn>2) will use the element to associate. "
                           << endl;
        }
    }
}

vector<int> ProcessProjectCAD::IntersectCADSurf(
    vector<CADSurfSharedPtr> v1_CADs, vector<CADSurfSharedPtr> v2_CADs)
{
    vector<CADSurfSharedPtr> v1 = v1_CADs;
    vector<CADSurfSharedPtr> v2 = v2_CADs;

    vector<int> vi1, vi2;
    for (size_t j = 0; j < v1.size(); ++j)
    {
        vi1.push_back(v1[j]->GetId());
    }
    for (size_t j = 0; j < v2.size(); ++j)
    {
        vi2.push_back(v2[j]->GetId());
    }

    sort(vi1.begin(), vi1.end());
    sort(vi2.begin(), vi2.end());

    vector<int> cmn;
    set_intersection(vi1.begin(), vi1.end(), vi2.begin(), vi2.end(),
                     back_inserter(cmn));

    return cmn;
}

void ProcessProjectCAD::ProjectEdges(
    EdgeSet &surfEdges, int order, bgi::rtree<boxI, bgi::quadratic<16>> &rtree)
{
    m_log(VERBOSE) << " Projecting Edges to CAD (CASE3)" << endl;
    // Project the Edges to CAD

    LibUtilities::PointsKey ekey(order + 1,
                                 LibUtilities::eGaussLobattoLegendre);
    Array<OneD, NekDouble> gll;
    LibUtilities::PointsManager()[ekey]->GetPoints(gll);

    // make surface edges high-order
    int cnt = 0;
    for (auto i = surfEdges.begin(); i != surfEdges.end(); i++)
    {
        // IF the edge is already associated with a CAD surface, HOSurf will do
        // the curving job
        if ((*i)->m_parentCAD)
        {
            continue;
        }
        if (lockedNodes.count((*i)->m_n1) || lockedNodes.count((*i)->m_n2))
        {
            continue;
        }
        cnt++;

        vector<CADSurfSharedPtr> v1 = (*i)->m_n1->GetCADSurfs();
        vector<CADSurfSharedPtr> v2 = (*i)->m_n2->GetCADSurfs();

        vector<int> vi1, vi2, vi1vi2;
        for (size_t j = 0; j < v1.size(); ++j)
        {
            vi1.push_back(v1[j]->GetId());
        }
        for (size_t j = 0; j < v2.size(); ++j)
        {
            vi2.push_back(v2[j]->GetId());
        }

        sort(vi1.begin(), vi1.end());
        sort(vi2.begin(), vi2.end());

        vector<int> cmn;
        set_intersection(vi1.begin(), vi1.end(), vi2.begin(), vi2.end(),
                         back_inserter(cmn));

        (*i)->m_curveType = LibUtilities::eGaussLobattoLegendre;

        cout << "cmn.size() = " << cmn.size() << endl;

        if (cmn.size() == 1 || cmn.size() == 2)
        {
            for (int j = 0; j < cmn.size(); j++)
            {
                if (m_mesh->m_cad->GetSurf(cmn[j])->IsPlanar())
                {
                    // if its planar dont care
                    continue;
                }

                auto uvb = (*i)->m_n1->GetCADSurfInfo(cmn[j]);
                auto uve = (*i)->m_n2->GetCADSurfInfo(cmn[j]);

                // can compare the loction of the projection to the
                // corresponding position of the straight sided edge
                // if the two differ by more than the length of the edge
                // something has gone wrong
                NekDouble len = (*i)->m_n1->Distance((*i)->m_n2);

                for (int k = 1; k < order + 1 - 1; k++)
                {
                    std::array<NekDouble, 2> uv = {
                        uvb[0] * (1.0 - gll[k]) / 2.0 +
                            uve[0] * (1.0 + gll[k]) / 2.0,
                        uvb[1] * (1.0 - gll[k]) / 2.0 +
                            uve[1] * (1.0 + gll[k]) / 2.0};
                    auto loc = m_mesh->m_cad->GetSurf(cmn[j])->P(uv);

                    std::array<NekDouble, 3> locT;
                    locT[0] = (*i)->m_n1->m_x * (1.0 - gll[k]) / 2.0 +
                              (*i)->m_n2->m_x * (1.0 + gll[k]) / 2.0;
                    locT[1] = (*i)->m_n1->m_y * (1.0 - gll[k]) / 2.0 +
                              (*i)->m_n2->m_y * (1.0 + gll[k]) / 2.0;
                    locT[2] = (*i)->m_n1->m_z * (1.0 - gll[k]) / 2.0 +
                              (*i)->m_n2->m_z * (1.0 + gll[k]) / 2.0;

                    NekDouble d = sqrt((locT[0] - loc[0]) * (locT[0] - loc[0]) +
                                       (locT[1] - loc[1]) * (locT[1] - loc[1]) +
                                       (locT[2] - loc[2]) * (locT[2] - loc[2]));

                    if (d > len)
                    {
                        (*i)->m_edgeNodes.clear();
                        break;
                    }

                    NodeSharedPtr nn = std::shared_ptr<Node>(
                        new Node(0, loc[0], loc[1], loc[2]));

                    (*i)->m_edgeNodes.push_back(nn);
                }

                if ((*i)->m_edgeNodes.size() != 0)
                {
                    // it suceeded on this surface so skip the other possibility
                    break;
                }
            }
        }
        else if (cmn.size() == 0)
        {
            // projection, if the projection requires more than two surfaces
            // including the edge nodes, then,  in theory projection shouldnt be
            // used

            vi1vi2.insert(vi1vi2.end(), vi1.begin(), vi1.end());
            vi1vi2.insert(vi1vi2.end(), vi2.begin(), vi2.end());

            cout << vi1vi2.size() << " CAD surfaces found for mn1 and mn2"
                 << endl;
            cout << (*i)->m_n1 << "    v2 " << (*i)->m_n2 << endl;
            set<int> sused;
            for (int k = 1; k < order + 1 - 1; k++)
            {
                std::array<NekDouble, 3> locT = {
                    (*i)->m_n1->m_x * (1.0 - gll[k]) / 2.0 +
                        (*i)->m_n2->m_x * (1.0 + gll[k]) / 2.0,
                    (*i)->m_n1->m_y * (1.0 - gll[k]) / 2.0 +
                        (*i)->m_n2->m_y * (1.0 + gll[k]) / 2.0,
                    (*i)->m_n1->m_z * (1.0 - gll[k]) / 2.0 +
                        (*i)->m_n2->m_z * (1.0 + gll[k]) / 2.0};

                int s;
                if (!FindAndProject(rtree, locT, s))
                {
                    (*i)->m_edgeNodes.clear();
                    m_log(VERBOSE) << "failed to find CAD" << endl;
                    break;
                }

                sused.insert(s);

                if (sused.size() > 2)
                {
                    m_log(WARNING) << "found too many CAD " << endl;
                    (*i)->m_edgeNodes.clear();
                    break;
                }

                NodeSharedPtr nn = std::shared_ptr<Node>(
                    new Node(0, locT[0], locT[1], locT[2]));

                (*i)->m_edgeNodes.push_back(nn);
            }
        }
    }
    cout << "Cnt for projected edges = " << cnt << endl;
}

void ProcessProjectCAD::Diagnostics()
{
    m_log(VERBOSE) << endl
                   << " ---------Diagnostics---------     " << endl
                   << endl;

    if (true)
    {
        long int counterElementsNoCAD = 0;
        long int counterEdgesNoCAD    = 0;
        long int counterVerticesNoCAD = 0;

        m_log(VERBOSE) << "         CAD Stats       " << endl;
        m_log(VERBOSE) << "CAD Vertices =           "
                       << m_mesh->m_cad->GetNumVerts() << endl;
        m_log(VERBOSE) << "CAD Curves =             "
                       << m_mesh->m_cad->GetNumCurve() << endl;
        m_log(VERBOSE) << "CAD Surfaces =           "
                       << m_mesh->m_cad->GetNumSurf() << endl;

        m_log(VERBOSE) << "        Mesh Stats       " << endl;
        m_log(VERBOSE) << "Mesh Surface Vertices =  " << surfNodes.size()
                       << endl;
        m_log(VERBOSE) << "Mesh Surface Edges =     " << surfEdges.size()
                       << endl;
        m_log(VERBOSE) << "Mesh Surface Face =      "
                       << m_mesh->m_element[2].size() << endl;

        // Elements without CAD
        EdgeSet NoCADEdges;
        for (auto element : m_mesh->m_element[2])
        {
            if (element->m_parentCAD == NULL)
            {
                counterElementsNoCAD++;
                // m_log(VERBOSE) << "Element No CAD Vertex1= " <<
                // element->GetVertex(0)->m_x <<" " <<
                // element->GetVertex(0)->m_y << " " <<
                // element->GetVertex(0)->m_z  << endl ;
            }
        }

        // Edges without CAD
        int cntEdgeCurve = 0, cntEdgeSurf = 0;
        for (auto edge : surfEdges)
        {
            if (edge->m_parentCAD == NULL)
            {
                counterEdgesNoCAD++;
            }
            else if (edge->m_parentCAD->GetType() == CADType::eCurve)
            {
                cntEdgeCurve++;
            }
            else
            {
                cntEdgeSurf++;
            }
        }

        // Vertices CAD
        int cntCAD2 = 0, cntCAD1 = 0, cntCAD3orMore = 0;
        int cntCADVertices = 0 , cntCADCurves = 0;
        for (auto vertex : surfNodes)
        {
            if (vertex->GetCADSurfs().size() == 0 &&
                vertex->GetCADCurves().size() == 0)
            {
                counterVerticesNoCAD++;
            }

            if (vertex->GetCADSurfs().size() == 1)
            {
                cntCAD1++;
            }

            if (vertex->GetCADSurfs().size() == 2)
            {
                cntCAD2++;
            }

            if (vertex->GetCADSurfs().size() > 2)
            {
                cntCAD3orMore++;
            }

            if (vertex->GetCADVertex() != NULL)
            {
                cntCADVertices++;
            }
            if(vertex->GetCADCurves().size() > 0)
            {
                cntCADCurves++;
            }
        }

        // Surface Edges without CAD

        m_log(WARNING) << "Vertices No CAD (Includes PlanarSurf) N= "
                       << counterVerticesNoCAD << endl;
        // m_log(WARNING) << "Planar Vertices =                        " <<
        // planarCnt << endl;
        m_log(WARNING) << "Edges without CADObject               N= "
                       << counterEdgesNoCAD << endl;
        m_log(WARNING) << "Faces without CADObject               N= "
                       << counterElementsNoCAD << endl;

        m_log(VERBOSE) << "Vertices 1 CADSurf                    N= " << cntCAD1
                       << endl;
        m_log(VERBOSE) << "Vertices 2 CADSurf                    N= " << cntCAD2
                       << endl;
        m_log(VERBOSE) << "Vertices CAD3 or more Surf            N= "
                       << cntCAD3orMore << endl;
        m_log(VERBOSE) << "Vertices- CADVertex (also have CADSu) N= "
                       << cntCADVertices << endl;
        m_log(VERBOSE) << "Vertices- CADCurve                    N= "
                       << cntCADCurves << endl;
        m_log(VERBOSE) << "Edges CADCurve                        N= "
                       << cntEdgeCurve << endl;
        m_log(VERBOSE) << "Edges CADSurf                         N= "
                       << cntEdgeSurf << endl;
    }
}

void ProcessProjectCAD::LinkFaceToCAD()
{
    for (auto element : m_mesh->m_element[2])
    {
        vector<NodeSharedPtr> vertices = element->GetVertexList();
        vector<EdgeSharedPtr> edges    = element->GetEdgeList();

        // CASE1 - all Edges internal to same CADSurf
        bool internal = true;
        for (int i = 1; i < edges.size(); i++)
        {
            if (edges[i]->m_parentCAD != edges[i - 1]->m_parentCAD)
            {
                internal = false;
                break;
            }
        }
        if (internal)
        {
            element->m_parentCAD = edges[0]->m_parentCAD;
            continue;
        }

        // CASE2 and CASE3 - 2 or more CADSurfs or None (Use vertice CAD)
        vector<vector<int>> cmn;
        for (int i = 1; i < vertices.size(); i++)
        {
            vector<int> cmn_i = IntersectCADSurf(
                vertices[i]->GetCADSurfs(), vertices[i - 1]->GetCADSurfs());
            if (cmn_i.size() > 0)
            {
                cmn.push_back(cmn_i);
            }
        }

        if (cmn.size() == 0)
        {
            // CASE 3
            // cout << "face cmn.size() = 0 " << endl;
            continue;
        }

        std::vector<int> commonCAD = cmn[0];
        for (auto cmn_i : cmn)
        {
            std::sort(cmn_i.begin(), cmn_i.end());
            std::vector<int> temp;
            std::set_intersection(commonCAD.begin(), commonCAD.end(),
                                  cmn_i.begin(), cmn_i.end(),
                                  std::back_inserter(temp));
            commonCAD = temp;
        }

        // commonCAD CADSurf found
        if (commonCAD.size() == 1)
        {
            // Internal element based on the
            element->m_parentCAD = m_mesh->m_cad->GetSurf(commonCAD[0]);
            for (auto edge : edges)
            {
                if (edge->m_parentCAD == NULL)
                {
                    edge->m_parentCAD = m_mesh->m_cad->GetSurf(commonCAD[0]);
                }
            }
        }
        else
        {
            // CASE 2 - 2 or more CADSurfs
            // Project the face
            cout << " face cmn.size() = " << commonCAD.size() << endl;
        }
    }

    // // Identify CAD object for every boundary element based on the common
    // // CADobjects of the member vertices
    // for (auto el : m_mesh->m_element[2])
    // {
    //     // cout << "New Element ID = "  << endl ;
    //     //  .1 Identify the common CADs
    //     vector<int> commonSurfacesEL;
    //     vector<int> CADObjectIDs0, CADObjectIDs1;
    //     for (auto CADObject : el->GetVertex(0)->GetCADSurfs())
    //     {
    //         CADObjectIDs0.push_back(CADObject->GetId());
    //         // cout << "CAD0 =" << CADObject->GetId() << endl ;
    //     }

    //     for (auto CADObject : el->GetVertex(1)->GetCADSurfs())
    //     {
    //         CADObjectIDs1.push_back(CADObject->GetId());
    //         // cout << "CAD1 =" << CADObject->GetId() << endl ;
    //     }

    //     sort(CADObjectIDs0.begin(), CADObjectIDs0.end());
    //     sort(CADObjectIDs1.begin(), CADObjectIDs1.end());

    //     set_intersection(CADObjectIDs0.begin(), CADObjectIDs0.end(),
    //                      CADObjectIDs1.begin(), CADObjectIDs1.end(),
    //                      back_inserter(commonSurfacesEL));

    //     // 1.2 Intersect with the next objects
    //     for (int i = 2; i < el->GetVertexCount(); i++)
    //     {
    //         vector<int> VertexCADSurfIDs;
    //         for (auto CADSurf_i : el->GetVertex(i)->GetCADSurfs())
    //         {
    //             VertexCADSurfIDs.push_back(CADSurf_i->GetId());
    //             // cout << "CADi =" << CADSurf_i->GetId() << endl ;
    //         }

    //         sort(VertexCADSurfIDs.begin(), VertexCADSurfIDs.end());
    //         sort(commonSurfacesEL.begin(), commonSurfacesEL.end());

    //         vector<int> cmn;
    //         set_intersection(VertexCADSurfIDs.begin(),
    //         VertexCADSurfIDs.end(),
    //                          commonSurfacesEL.begin(),
    //                          commonSurfacesEL.end(), back_inserter(cmn));
    //         commonSurfacesEL = cmn;
    //         // m_log(VERBOSE) << " commonSurfacesEL.size =  " <<
    //         // commonSurfacesEL.size() << endl ;
    //     }

    //     // m_log(WARNING) << " commonSurfacesEL.size() = " <<
    //     // commonSurfacesEL.size() << endl ;

    //     // .2 Are there common CAD objects (if not ProjectCAD legacy)
    //     if (commonSurfacesEL.size() == 1)
    //     {
    //         // Internal or Corner to 1 NURB element!
    //         // Workflow - assign

    //         // .2.1 Ideally - we have internal to the CADSurf vertex so
    //         assign
    //         // it to the CADSurf associated to it
    //         for (auto vertex : el->GetVertexList())
    //         {
    //             // If one of the vertices has only 1 CAD Surf
    //             // -> hence internal vertex ->
    //             // 2D element is part of this surface
    //             // V&VTestCase 2 Semisphere

    //             // cout << " vertex->GetCADSurf().size()   =
    //             // "<<vertex->GetCADSurfs().size()  << endl ;
    //             if (vertex->GetCADSurfs().size() == 1)
    //             {
    //                 CADSurfSharedPtr CADSurf = vertex->GetCADSurfs()[0];
    //                 el->m_parentCAD          = CADSurf;
    //             }
    //         }

    //         // .2.2 If no internal vertices (EX proximity or TE with only 1
    //         Quad
    //         // element - this becomes CASE2 always)
    //         if (!el->m_parentCAD)
    //         {
    //             if (el->GetShapeType() != 3)
    //             {
    //                 // if a rectangular face - it should always have an
    //                 internal
    //                 // to the face vertex
    //                 m_log(VERBOSE)
    //                     << " Face Shape Type = " << el->GetShapeType()
    //                     << "  EdgeList.size()= " << el->GetEdgeList().size()
    //                     << endl;
    //                 ;
    //                 m_log(WARNING) << " NO Face internal vertices for a Quad
    //                 "
    //                                   "element - EX:Proximity May need some "
    //                                   "further testing/development? "
    //                                << endl;
    //             }
    //             el->m_parentCAD =
    //             m_mesh->m_cad->GetSurf(commonSurfacesEL[0]);
    //         }

    //         // .2.3 Assign !!CADCurve!! not CADSURF!!! to Edges that are on
    //         the
    //         // corner between two NNURBS -> its vertices have more than one
    //         // "common" CADSurf TestCase 1 - RING + TestCase2 -
    //         for (auto edge : el->GetEdgeList())
    //         {
    //             if (edge->m_parentCAD)
    //             {
    //                 continue;
    //             }
    //             // .2.3.1 Identified the common CADSurf
    //             if (lockedNodes.count(edge->m_n1) ||
    //                 lockedNodes.count(edge->m_n2))
    //             {
    //                 continue;
    //             }
    //             vector<CADSurfSharedPtr> v1 = edge->m_n1->GetCADSurfs();
    //             vector<CADSurfSharedPtr> v2 = edge->m_n2->GetCADSurfs();

    //             vector<int> vi1, vi2;
    //             for (size_t j = 0; j < v1.size(); ++j)
    //             {
    //                 vi1.push_back(v1[j]->GetId());
    //             }
    //             for (size_t j = 0; j < v2.size(); ++j)
    //             {
    //                 vi2.push_back(v2[j]->GetId());
    //             }

    //             sort(vi1.begin(), vi1.end());
    //             sort(vi2.begin(), vi2.end());

    //             vector<int> commonSurfacesEdge;
    //             set_intersection(vi1.begin(), vi1.end(), vi2.begin(),
    //             vi2.end(),
    //                              back_inserter(commonSurfacesEdge));

    //             // .2.3.2 NB! Main CAD Curve reconstruction algorithm
    //             if (commonSurfacesEdge.size() == 1)
    //             {
    //                 // if the edge is internal to the CAD Surf -> use the CAD
    //                 // Surf as parentCAD
    //                 edge->m_parentCAD = el->m_parentCAD;
    //             }
    //             else if (commonSurfacesEdge.size() == 2) // this might be 3 ?
    //             {
    //                 // Corner Edge -> !!!CAD CURVE procedure!!! INSTEAD of
    //                 // CADSurf to avoid #bug 1(TestCASE 1-2 for junction
    //                 betwen
    //                 // planar and NURBS surfaces)
    //                 if (edge->m_parentCAD)
    //                 {
    //                     // cout << "Corner Edge has a parentCAD -> continue "
    //                     <<
    //                     // endl ;
    //                     continue;
    //                 }
    //                 // m_log(WARNING) << "CORNER EDGE -> Curve Reconstruction
    //                 "
    //                 // << endl ; m_log(VERBOSE) <<
    //                 //
    //                 "m_mesh->m_cad->GetSurf(commonSurfacesEdge[0])->GetEdges().size()
    //                 // - EdgeLoops = " <<
    //                 //
    //                 m_mesh->m_cad->GetSurf(commonSurfacesEdge[0])->GetEdges().size()
    //                 // << endl ; m_log(VERBOSE) <<
    //                 //
    //                 "m_mesh->m_cad->GetSurf(commonSurfacesEdge[1])->GetEdges().size()
    //                 // - EdgeLoops " <<
    //                 //
    //                 m_mesh->m_cad->GetSurf(commonSurfacesEdge[1])->GetEdges().size()
    //                 // << endl ;

    //                 // m_log(VERBOSE) << "CADCurves per Surf1-> EdgeLoop[0] "
    //                 <<
    //                 //
    //                 m_mesh->m_cad->GetSurf(commonSurfacesEdge[0])->GetEdges()[0]->edges.size()
    //                 // << endl ; m_log(VERBOSE) << "CADCurves per Surf2->
    //                 // EdgeLoop[0] " <<
    //                 //
    //                 m_mesh->m_cad->GetSurf(commonSurfacesEdge[1])->GetEdges()[0]->edges.size()
    //                 // << endl ;

    //                 // 1.Create vi1 , vi2
    //                 vector<int> CADCurves_uv1;
    //                 vector<int> CADCurves_uv2;

    //                 CADSurfSharedPtr EdgeSurf1 =
    //                     m_mesh->m_cad->GetSurf(commonSurfacesEdge[0]);
    //                 CADSurfSharedPtr EdgeSurf2 =
    //                     m_mesh->m_cad->GetSurf(commonSurfacesEdge[1]);

    //                 for (auto EdgeLoop : EdgeSurf1->GetEdges())
    //                 {
    //                     for (auto CADCurve : EdgeLoop->edges)
    //                     {
    //                         CADCurves_uv1.push_back(CADCurve->GetId());
    //                         // cout << "Curves1 = " << CADCurve->GetId() <<
    //                         endl
    //                         // ;
    //                     }
    //                 }
    //                 for (auto EdgeLoop : EdgeSurf2->GetEdges())
    //                 {
    //                     for (auto CADCurve : EdgeLoop->edges)
    //                     {
    //                         CADCurves_uv2.push_back(CADCurve->GetId());
    //                         // cout << "Curves2 = " << CADCurve->GetId() <<
    //                         endl
    //                         // ;
    //                     }
    //                 }

    //                 sort(CADCurves_uv1.begin(), CADCurves_uv1.end());
    //                 sort(CADCurves_uv2.begin(), CADCurves_uv2.end());

    //                 vector<int> commonCADCurves;
    //                 set_intersection(CADCurves_uv1.begin(),
    //                 CADCurves_uv1.end(),
    //                                  CADCurves_uv2.begin(),
    //                                  CADCurves_uv2.end(),
    //                                  back_inserter(commonCADCurves));

    //                 NekDouble tolDist = 5e-5; // POSSIBLE PROBLEM !!!!
    //                 if (commonCADCurves.size() == 1)
    //                 {
    //                     // m_log(WARNING) << " CASE
    //                     commonCADCurves.size()==1"
    //                     // << endl ;
    //                     //  m_log(VERBOSE) << "edge id = " << edge->m_id <<
    //                     endl
    //                     //  ; m_log(VERBOSE) << "Vertex1  = " <<
    //                     //  edge->m_n1->GetLoc()[0] << ' ' <<
    //                     //  edge->m_n1->GetLoc()[1]  << " "<<
    //                     //  edge->m_n1->GetLoc()[2]  << endl ; m_log(VERBOSE)
    //                     <<
    //                     //  "VErtex2  = " << edge->m_n2->GetLoc()[0] << ' '
    //                     <<
    //                     //  edge->m_n2->GetLoc()[1]  << " "<<
    //                     //  edge->m_n2->GetLoc()[2]  << endl ;

    //                     // m_log(VERBOSE) << "CAD curve = " <<
    //                     // commonCADCurves[0] << endl ;

    //                     CADCurveSharedPtr CADCurve_t =
    //                         m_mesh->m_cad->GetCurve(commonCADCurves[0]);

    //                     // cout << "CAD VErtex 0 = "<<
    //                     // CADCurve_t->GetVertex()[0]->GetLoc()[0] << " "
    //                     // <<CADCurve_t->GetVertex()[0]->GetLoc()[1]  << " "
    //                     <<
    //                     // CADCurve_t->GetVertex()[0]->GetLoc()[2]  << endl ;
    //                     // cout << "CAD VErtex 1 = "<<
    //                     // CADCurve_t->GetVertex()[1]->GetLoc()[0] << " "
    //                     // <<CADCurve_t->GetVertex()[1]->GetLoc()[1]  << " "
    //                     <<
    //                     // CADCurve_t->GetVertex()[1]->GetLoc()[2]  << endl ;
    //                     edge->m_parentCAD = CADCurve_t;

    //                     NekDouble t0, t1, dist0, dist1;
    //                     NekDouble tmin, tmax;
    //                     CADCurve_t->GetBounds(tmin, tmax);

    //                     dist0 = CADCurve_t->loct(edge->m_n1->GetLoc(), t0,
    //                     tmin,
    //                                              tmax);
    //                     dist1 = CADCurve_t->loct(edge->m_n2->GetLoc(), t1,
    //                     tmin,
    //                                              tmax);
    //                     boost::ignore_unused(dist0, dist1);

    //                     // Array<OneD, NekDouble> xyz = edge->m_n1->GetLoc()
    //                     ;
    //                     // NekDouble dist0_GetMinDistance=
    //                     // CADCurve_t->GetMinDistance(xyz) ; xyz =
    //                     // edge->m_n2->GetLoc() ; NekDouble
    //                     // dist1_GetMinDistance=
    //                     CADCurve_t->GetMinDistance(xyz)
    //                     // ;
    //                     // cout << "dist0 = " << dist0 << " dist1= " << dist1
    //                     <<
    //                     // endl ;
    //                     // cout << "dist0 = " << dist0_GetMinDistance << "
    //                     // dist1= " << dist1_GetMinDistance << "GetMinDist"<<
    //                     // endl ; if(edge->m_id==796)
    //                     // {
    //                     //      m_log(WARNING) << "t0 "  << t0 << " t1 " <<t1
    //                     <<
    //                     //      endl ;
    //                     // }
    //                     // cout << "t0 = " << setprecision(10)  << t0 << "
    //                     t1= "
    //                     // << t1 << endl ; cout << "CAD tmin= " <<
    //                     // setprecision(10)  << tmin << " tmax= "<<  tmax <<
    //                     // endl ; NekDouble Tol_parametric =
    //                     // abs(CADCurve_t->GetTotLength())*1e-6;

    //                     /*
    //                                             if((dist0 < tolDist) &&
    //                                             (dist1 <
    //                        tolDist) )
    //                                             {
    //                                                 //edge->m_n1->SetCADCurve(CADCurve_t,t0);
    //                                                 edge->m_n2->SetCADCurve(CADCurve_t,t1);

    //                                                 // IF USE GetMinDist ->
    //                                                 Does
    //                        not give the correct distance to the vertex in
    //                        case t outside of curve [t_min t_max] bonds

    //                                                 // vertex 1
    //                                                 if(t0 >= tmin - tolDist
    //                                                 &&
    //                        t0 <= tmax + tolDist )
    //                                                 {
    //                                                     edge->m_n1->SetCADCurve(CADCurve_t,t0);
    //                                                 }
    //                                                 else if( (t0 > tmin) &&
    //                                                 (t0
    //                        < tmax + tolDist))
    //                                                 {
    //                                                     // this is necessary
    //                                                     to
    //                        accomodate for Vertices generated in StarCCM which
    //                        are slightly further away
    //                                                     // Alternatively one
    //                                                     can
    //                        project all vertices on the m_log(WARNING) <<
    //                        "Parametric m_n1 on CADCurve has been manually to
    //                        curve t_max " <<  endl ;
    //                                                     edge->m_n1->SetCADCurve(CADCurve_t,tmax);
    //                                                 }
    //                                                 else if( (t0 > tmin-
    //                        tolDist) && (t0 < tmax ))
    //                                                 {
    //                                                     // this is necessary
    //                                                     to
    //                        accomodate for Vertices generated in StarCCM which
    //                        are slightly further away
    //                                                     // Alternatively one
    //                                                     can
    //                        project all vertices on the m_log(WARNING) <<
    //                        "Parametric m_n1 on CADCurve has been manually to
    //                        curve tmin " <<  endl ;
    //                                                     edge->m_n1->SetCADCurve(CADCurve_t,tmin);
    //                                                 }
    //                                                 else
    //                                                 {
    //                                                      edge->m_parentCAD =
    //                        NULL ; m_log(WARNING) <<  "Parametric m_n1 on
    //                        CADCurve cmn=1 has been too far away from the
    //                        Curve Bounds " <<  endl ;
    //                                                 }

    //                                                 // // vertex 2
    //                                                 if(t1 >= tmin - tolDist
    //                                                 &&
    //                        t1 <= tmax + tolDist )
    //                                                 {
    //                                                     edge->m_n2->SetCADCurve(CADCurve_t,t1);
    //                                                 }
    //                                                 else if( (t1 > tmin) &&
    //                                                 (t1
    //                        < tmax + tolDist))
    //                                                 {
    //                                                     // this is necessary
    //                                                     to
    //                        accomodate for Vertices generated in StarCCM which
    //                        are slightly further away
    //                                                     // Alternatively one
    //                                                     can
    //                        project all vertices on the m_log(WARNING) <<
    //                        "Parametric m_n2 on CADCurve has been manually to
    //                        curve t_max " <<  endl ;
    //                                                     edge->m_n2->SetCADCurve(CADCurve_t,tmax);
    //                                                 }
    //                                                 else if( (t1 > tmin-
    //                        tolDist) && (t1 < tmax ))
    //                                                 {
    //                                                     // this is necessary
    //                                                     to
    //                        accomodate for Vertices generated in StarCCM which
    //                        are slightly further away
    //                                                     // Alternatively one
    //                                                     can
    //                        project all vertices on the m_log(WARNING) <<
    //                        "Parametric m_n2 on CADCurve has been manually to
    //                        curve tmin " <<  endl ;
    //                                                     edge->m_n2->SetCADCurve(CADCurve_t,tmin);
    //                                                 }
    //                                                 else
    //                                                 {
    //                                                     edge->m_parentCAD =
    //                                                     NULL
    //                        ; m_log(WARNING) <<  "Parametric m_n2 on CADCurve
    //                        cmn=1 has been too far away from the Curve Bounds
    //                        "
    //                        <<  endl ;
    //                                                 }
    //                                             }
    //                                             else
    //                                             {
    //                                                 cout << "dist of m_n1 to
    //                                                 the
    //                        2 CAD Curves " << dist0 << " " << dist1 << endl ;
    //                                                 m_log(WARNING) <<
    //                        "CommonCADCurves =1 - The edge vertices are too
    //                        far away from the CADCurve dist = " << dist0 << "
    //                        " << dist1 << endl ;
    //                                             }
    //                     */
    //                 }
    //                 else if (commonCADCurves.size() == 2)
    //                 {
    //                     // m_log(WARNING) << " CASE
    //                     commonCADCurves.size()==2"
    //                     // << endl ;
    //                     Array<OneD, NekDouble> xyz1 = edge->m_n1->GetLoc();
    //                     Array<OneD, NekDouble> xyz2 = edge->m_n2->GetLoc();

    //                     CADCurveSharedPtr Curve_t0 =
    //                         m_mesh->m_cad->GetCurve(commonCADCurves[0]);
    //                     CADCurveSharedPtr Curve_t1 =
    //                         m_mesh->m_cad->GetCurve(commonCADCurves[1]);

    //                     NekDouble dist01, dist02, dist11, dist12, m_n1_t0,
    //                         m_n1_t1, m_n2_t0, m_n2_t1;
    //                     // cout << "Initialization m_n1/2 t0/1 " << m_n1_t0
    //                     << "
    //                     // " << m_n2_t0 << " " <<  m_n1_t1 <<" " <<  m_n2_t1
    //                     <<
    //                     // endl;

    //                     // NekDouble tmin0, tmax0, tmin1, tmax1;
    //                     Array<OneD, NekDouble> tmin_max0 =
    //                         Curve_t0->GetBounds();
    //                     Array<OneD, NekDouble> tmin_max1 =
    //                         Curve_t1->GetBounds();
    //                     // cout << "tmin_max0 =  " << tmin_max0[0] << " " <<
    //                     // tmin_max0[1] << endl ; cout << "tmin_max1 =  " <<
    //                     // tmin_max1[0] << " " << tmin_max1[1] << endl ;

    //                     dist01 = Curve_t0->loct(xyz1, m_n1_t0, tmin_max0[0],
    //                                             tmin_max0[1]);
    //                     dist02 = Curve_t0->loct(xyz2, m_n2_t0, tmin_max0[0],
    //                                             tmin_max0[1]);

    //                     dist11 = Curve_t1->loct(xyz1, m_n1_t1, tmin_max1[0],
    //                                             tmin_max1[1]);
    //                     dist12 = Curve_t1->loct(xyz2, m_n2_t1, tmin_max1[0],
    //                                             tmin_max1[1]);

    //                     // cout << "After Update : " << m_n1_t0 << " " <<
    //                     // m_n2_t0 << " " << m_n1_t1 << " " << m_n2_t1 <<
    //                     endl ;

    //                     // Just for V&V
    //                     // NekDouble dist01_GetMinDistance=
    //                     // Curve_t0->GetMinDistance(xyz1) ; NekDouble
    //                     // dist02_GetMinDistance=
    //                     Curve_t0->GetMinDistance(xyz2)
    //                     // ; NekDouble dist11_GetMinDistance=
    //                     // Curve_t1->GetMinDistance(xyz1) ; NekDouble
    //                     // dist12_GetMinDistance=
    //                     Curve_t1->GetMinDistance(xyz2)
    //                     // ;

    //                     // m_log(VERBOSE) <<" dist 01 "  << dist01 << "
    //                     dist02 "
    //                     // << dist02 << endl ; m_log(VERBOSE) << " GetMin "<<
    //                     "
    //                     // dist 01 "  << dist01_GetMinDistance << " dist02 "
    //                     <<
    //                     // dist02_GetMinDistance << endl ; m_log(VERBOSE) <<"
    //                     // dist 11 "  << dist11 << " dist12 " << dist12 <<
    //                     endl
    //                     // ; m_log(VERBOSE) << " GetMin "<<  " dist 01 "  <<
    //                     // dist11_GetMinDistance << " dist02 " <<
    //                     // dist12_GetMinDistance << endl ;

    //                     // if( (dist1 < tol)  ||  (dist0 < tol)  )
    //                     NekDouble dist0 = (dist01 + dist02) / 2.0;
    //                     NekDouble dist1 = (dist11 + dist12) / 2.0;
    //                     if ((dist0 < tolDist || dist1 < tolDist))
    //                     {
    //                         if (dist0 < dist1 &&
    //                             (m_n1_t0 > tmin_max0[0] - tolDist) &&
    //                             (m_n1_t0 < tmin_max0[1] + tolDist) &&
    //                             (m_n2_t0 > tmin_max0[0] - tolDist) &&
    //                             (m_n2_t0 < tmin_max0[1] + tolDist))
    //                         {
    //                             // m_log(VERBOSE)  << "CAD1" << endl ;
    //                             edge->m_parentCAD = Curve_t0;

    //                             edge->m_n1->SetCADCurve(Curve_t0, m_n1_t0);
    //                             edge->m_n2->SetCADCurve(Curve_t0, m_n2_t0);
    //                             // cout << "CADCurve 0 m_n1_t0= " <<
    //                             //
    //                             edge->m_n1->GetCADCurveInfo(Curve_t0->GetId())
    //                             // << endl ; cout << "CADCurve 0 m_n2_t0= "
    //                             <<
    //                             //
    //                             edge->m_n2->GetCADCurveInfo(Curve_t0->GetId())
    //                             // << endl ;
    //                         }
    //                         else if (dist1 < dist0 &&
    //                                  (m_n1_t1 > tmin_max1[0] - tolDist) &&
    //                                  (m_n1_t1 < tmin_max1[1] + tolDist) &&
    //                                  (m_n2_t1 > tmin_max1[0] - tolDist) &&
    //                                  (m_n2_t1 < tmin_max1[1] + tolDist))
    //                         {
    //                             // m_log(VERBOSE) << "CAD2" << endl ;
    //                             edge->m_parentCAD = Curve_t1;
    //                             edge->m_n1->SetCADCurve(Curve_t1, m_n1_t1);
    //                             edge->m_n2->SetCADCurve(Curve_t1, m_n2_t1);
    //                             // cout << "CADCurve 1 m_n1_t1= " <<
    //                             //
    //                             edge->m_n1->GetCADCurveInfo(Curve_t1->GetId())
    //                             // << endl ; cout << "CADCurve 1 m_n2_t1 " <<
    //                             //
    //                             edge->m_n2->GetCADCurveInfo(Curve_t1->GetId())
    //                             // << endl ;
    //                         }
    //                         else
    //                         {
    //                             m_log(WARNING)
    //                                 << "loct circle bug cmn.size()=2 " <<
    //                                 endl;
    //                         }
    //                     }
    //                     else
    //                     {
    //                         m_log(WARNING) << " dist 01 " << dist01
    //                                        << " dist02 " << dist02 << endl;
    //                         m_log(WARNING) << " dist 11 " << dist11
    //                                        << " dist12 " << dist12 << endl;

    //                         cout << "dist of m_n1 to the 2 CAD Curves " <<
    //                         dist0
    //                              << " " << dist1 << endl;
    //                         m_log(WARNING)
    //                             << "Cannot find close CADCurve in CornerCase
    //                             "
    //                                "of Edge with 2 common curves"
    //                             << endl;
    //                     }
    //                     // cout << "end edge commonCADCurves.size()==2 " <<
    //                     endl
    //                     // ;
    //                 }
    //                 else
    //                 {
    //                     // m_log(VERBOSE) << "Edge->m_id =  " << edge->m_id
    //                     << "
    //                     // common CADCurves N= " << commonCADCurves.size()
    //                     // <<endl; m_log(VERBOSE) << edge->m_n1->GetLoc()[0]
    //                     <<
    //                     // " " <<  edge->m_n1->GetLoc()[1] << " " <<
    //                     // edge->m_n1->GetLoc()[2] << endl; m_log(VERBOSE) <<
    //                     // edge->m_n2->GetLoc()[0] << " " <<
    //                     // edge->m_n2->GetLoc()[1] << " " <<
    //                     // edge->m_n2->GetLoc()[2] << endl; m_log(WARNING) <<
    //                     // "Edge has more than 2 CADCurves common, assign
    //                     // automatically one of the surfaces"<<endl;

    //                     // Check all edges of the two surfaces  and select
    //                     one
    //                     // When topologically they are two CAD curves that do
    //                     // not touch
    //                     for (auto surfID : commonSurfacesEdge)
    //                     {
    //                         CADSurfSharedPtr surf =
    //                             m_mesh->m_cad->GetSurf(surfID);
    //                         for (auto edgeloop : surf->GetEdges())
    //                         {
    //                             for (auto curve : edgeloop->edges)
    //                             {
    //                                 NekDouble dist01, dist02, m_n1_t0,
    //                                 m_n2_t0; Array<OneD, NekDouble> tmin_max0
    //                                 =
    //                                     curve->GetBounds();

    //                                 dist01 =
    //                                 curve->loct(edge->m_n1->GetLoc(),
    //                                                      m_n1_t0,
    //                                                      tmin_max0[0],
    //                                                      tmin_max0[1]);
    //                                 dist02 =
    //                                 curve->loct(edge->m_n2->GetLoc(),
    //                                                      m_n2_t0,
    //                                                      tmin_max0[0],
    //                                                      tmin_max0[1]);
    //                                 if (dist01 < 1e-7 && dist02 < 1e-7)
    //                                 {
    //                                     cout << "CADCurve adopted " << endl;
    //                                     edge->m_parentCAD = curve;
    //                                     break;
    //                                 }
    //                             }
    //                         }
    //                         // el->m_parentCAD = NULL;
    //                     }
    //                 }
    //             }
    //             else
    //             {
    //                 m_log(WARNING)
    //                     << "EdgeCommonSurfaces = " <<
    //                     commonSurfacesEdge.size()
    //                     << endl;
    //                 for (auto surfID : commonSurfacesEdge)
    //                 {
    //                     CADSurfSharedPtr surf =
    //                     m_mesh->m_cad->GetSurf(surfID); for (auto edgeloop :
    //                     surf->GetEdges())
    //                     {
    //                         for (auto curve : edgeloop->edges)
    //                         {
    //                             NekDouble dist01, dist02, m_n1_t0, m_n2_t0;
    //                             Array<OneD, NekDouble> tmin_max0 =
    //                                 curve->GetBounds();

    //                             dist01 =
    //                                 curve->loct(edge->m_n1->GetLoc(),
    //                                 m_n1_t0,
    //                                             tmin_max0[0], tmin_max0[1]);
    //                             dist02 =
    //                                 curve->loct(edge->m_n2->GetLoc(),
    //                                 m_n2_t0,
    //                                             tmin_max0[0], tmin_max0[1]);
    //                             if (dist01 < 1e-7 && dist02 < 1e-7)
    //                             {
    //                                 cout << "CADCurve adopted " << endl;
    //                                 edge->m_parentCAD = curve;
    //                                 break;
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     else if (commonSurfacesEL.size() == 0)
    //     {
    //         // CASE 3 (Element Between two NURBS surfaces - use ProjectCAD -
    //         BUG
    //         // #2  )
    //         // m_log(WARNING) << "ProjectCAD legacy" << endl ;
    //         // If CASE 3 (element split between >1 NURBS by StarCCM) use the
    //         // legacy ProjectCAD module Do projectCAD

    //         vector<EdgeSharedPtr> surfEdgesLocal = el->GetEdgeList();
    //         int order                            =
    //         m_config["order"].as<int>();

    //         map<int, vector<int>> eds;

    //         LibUtilities::PointsKey ekey(order + 1,
    //                                      LibUtilities::eGaussLobattoLegendre);
    //         Array<OneD, NekDouble> gll;
    //         LibUtilities::PointsManager()[ekey]->GetPoints(gll);

    //         // make surface edges high-order
    //         for (auto edge : surfEdgesLocal)
    //         {
    //             // m_log(VERBOSE) << " edge->m_id " <<  edge->m_id << endl ;
    //             // m_log(VERBOSE) << "Start edge->m_edgeNodes = " <<
    //             // edge->m_edgeNodes.size() << endl ;
    //             if ((edge->m_edgeNodes.size() != 0) || (edge->m_parentCAD))
    //             {
    //                 // skip projection again Edges that are already projected
    //                 // If not bug on the junction
    //                 // cout << "Skip the edge - already projected or given
    //                 CAD"
    //                 // << endl ;
    //                 continue;
    //             }

    //             if (lockedNodes.count(edge->m_n1) ||
    //                 lockedNodes.count(edge->m_n2))
    //             {
    //                 continue;
    //             }
    //             vector<CADSurfSharedPtr> v1 = edge->m_n1->GetCADSurfs();
    //             vector<CADSurfSharedPtr> v2 = edge->m_n2->GetCADSurfs();

    //             vector<int> vi1, vi2;
    //             for (size_t j = 0; j < v1.size(); ++j)
    //             {
    //                 vi1.push_back(v1[j]->GetId());
    //             }
    //             for (size_t j = 0; j < v2.size(); ++j)
    //             {
    //                 vi2.push_back(v2[j]->GetId());
    //             }

    //             sort(vi1.begin(), vi1.end());
    //             sort(vi2.begin(), vi2.end());

    //             vector<int> cmn;
    //             set_intersection(vi1.begin(), vi1.end(), vi2.begin(),
    //             vi2.end(),
    //                              back_inserter(cmn));
    //             eds[cmn.size()].push_back(0);

    //             edge->m_curveType = LibUtilities::eGaussLobattoLegendre;

    //             if (cmn.size() == 1 || cmn.size() == 2)
    //             {
    //                 m_log(WARNING) << " CASE2 cmn.size() =  " << cmn.size()
    //                                << " " << edge->m_n1->GetLoc()[0] << " "
    //                                << edge->m_n1->GetLoc()[1] << " "
    //                                << edge->m_n1->GetLoc()[2] << endl;

    //                 if (cmn.size() == 1)
    //                 {
    //                     // cout << "asign CADSurf to edge cmn=1 " << endl ;
    //                     edge->m_parentCAD = m_mesh->m_cad->GetSurf(cmn[0]);
    //                     continue;
    //                 }

    //                 for (int j = 0; j < cmn.size(); j++)
    //                 {
    //                     if (m_mesh->m_cad->GetSurf(cmn[j])->IsPlanar())
    //                     {
    //                         // if its planar dont care
    //                         continue;
    //                     }

    //                     Array<OneD, NekDouble> uvb =
    //                         edge->m_n1->GetCADSurfInfo(cmn[j]);
    //                     Array<OneD, NekDouble> uve =
    //                         edge->m_n2->GetCADSurfInfo(cmn[j]);

    //                     // can compare the loction of the projection to the
    //                     // corresponding position of the straight sided edge
    //                     // if the two differ by more than the length of the
    //                     edge
    //                     // something has gone wrong
    //                     NekDouble len = edge->m_n1->Distance(edge->m_n2);

    //                     for (int k = 1; k < order + 1 - 1; k++)
    //                     {
    //                         Array<OneD, NekDouble> uv(2);
    //                         uv[0] = uvb[0] * (1.0 - gll[k]) / 2.0 +
    //                                 uve[0] * (1.0 + gll[k]) / 2.0;
    //                         uv[1] = uvb[1] * (1.0 - gll[k]) / 2.0 +
    //                                 uve[1] * (1.0 + gll[k]) / 2.0;
    //                         Array<OneD, NekDouble> loc;
    //                         loc = m_mesh->m_cad->GetSurf(cmn[j])->P(uv);
    //                         Array<OneD, NekDouble> locT(3);
    //                         locT[0] = edge->m_n1->m_x * (1.0 - gll[k]) / 2.0
    //                         +
    //                                   edge->m_n2->m_x * (1.0 + gll[k]) / 2.0;
    //                         locT[1] = edge->m_n1->m_y * (1.0 - gll[k]) / 2.0
    //                         +
    //                                   edge->m_n2->m_y * (1.0 + gll[k]) / 2.0;
    //                         locT[2] = edge->m_n1->m_z * (1.0 - gll[k]) / 2.0
    //                         +
    //                                   edge->m_n2->m_z * (1.0 + gll[k]) / 2.0;

    //                         NekDouble d =
    //                             sqrt((locT[0] - loc[0]) * (locT[0] - loc[0])
    //                             +
    //                                  (locT[1] - loc[1]) * (locT[1] - loc[1])
    //                                  + (locT[2] - loc[2]) * (locT[2] -
    //                                  loc[2]));

    //                         if (d > len)
    //                         {
    //                             // m_log(WARNING) << "CASE3 cmn1 d > len
    //                             (Legacy
    //                             // ProjectCAD - no problem delete the log KK
    //                             )"
    //                             // << endl ;
    //                             edge->m_edgeNodes.clear();
    //                             break;
    //                         }

    //                         NodeSharedPtr nn = std::shared_ptr<Node>(
    //                             new Node(0, loc[0], loc[1], loc[2]));

    //                         edge->m_edgeNodes.push_back(nn);
    //                     }

    //                     if (edge->m_edgeNodes.size() != 0)
    //                     {
    //                         // it suceeded on this surface so skip the other
    //                         // possibility
    //                         break;
    //                     }
    //                 }
    //             }
    //             else if (cmn.size() == 0)
    //             {
    //                 // m_log(WARNING) << "CASE3 cmn=0 loc= " <<
    //                 // edge->m_n1->GetLoc()[0] << " " <<
    //                 edge->m_n1->GetLoc()[1]
    //                 // << " " << edge->m_n1->GetLoc()[2]   << endl ;
    //                 //  projection, if the projection requires more than two
    //                 //  surfaces including the edge nodes, then,  in theory
    //                 //  projection shouldnt be used

    //                 // Optional FORCING CASE 3 ELEMENTSS to be with a lower p
    //                 to
    //                 // avoid projection problems in C1 discontinuous surfaces
    //                 // (Disabled - orderSafe = order)
    //                 int orderSafe = order;
    //                 LibUtilities::PointsKey ekey(
    //                     orderSafe + 1, LibUtilities::eGaussLobattoLegendre);
    //                 Array<OneD, NekDouble> gll;
    //                 LibUtilities::PointsManager()[ekey]->GetPoints(gll);

    //                 set<int> sused;
    //                 for (int k = 1; k < orderSafe + 1 - 1; k++)
    //                 {
    //                     Array<OneD, NekDouble> locT(3);
    //                     locT[0] = edge->m_n1->m_x * (1.0 - gll[k]) / 2.0 +
    //                               edge->m_n2->m_x * (1.0 + gll[k]) / 2.0;
    //                     locT[1] = edge->m_n1->m_y * (1.0 - gll[k]) / 2.0 +
    //                               edge->m_n2->m_y * (1.0 + gll[k]) / 2.0;
    //                     locT[2] = edge->m_n1->m_z * (1.0 - gll[k]) / 2.0 +
    //                               edge->m_n2->m_z * (1.0 + gll[k]) / 2.0;

    //                     int s;
    //                     // cout << "initial locT " << locT[0] << " " <<
    //                     locT[1]
    //                     // << " " << locT[2] << endl ;
    //                     if (!findAndProject(rtree, locT, s))
    //                     {
    //                         edge->m_edgeNodes.clear();
    //                         break;
    //                     }

    //                     // cout << "end locT " << locT[0] << " " << locT[1]
    //                     << "
    //                     // " << locT[2] << endl ;

    //                     sused.insert(s);

    //                     if (sused.size() > 2)
    //                     {
    //                         edge->m_edgeNodes.clear();
    //                         break;
    //                     }

    //                     NodeSharedPtr nn = std::shared_ptr<Node>(
    //                         new Node(0, locT[0], locT[1], locT[2]));

    //                     edge->m_edgeNodes.push_back(nn);
    //                 }
    //             }

    //             // m_log(VERBOSE) << "Edge  = " << edge->m_id << endl ;
    //             // m_log(WARNING) << "End edge->m_edgeNodes = " <<
    //             // edge->m_edgeNodes.size() << endl ;
    //         }
    //     }
    //     else
    //     {
    //         // TestCase TC09 NACA0012 with very thin TE and a mainplane wing
    //         // constructed by single NURBS CADSurf The TE face has only 1
    //         // element in vertical direction -> hence all vertices are part
    //         of
    //         // both the mainplane NURBS + TE plane The following procedure
    //         will
    //         // work also if the TE surface is curved NURBS. (R of curvature <
    //         TE
    //         // thickness/2.0), otherwise it will pick the NURBS Hence in the
    //         end
    //         // a Jacobian check of the face is constructed, to avoid J=0
    //         // elements

    //         m_log(WARNING)
    //             << "The element[2] vertices are part of too many common "
    //                "surfaces commonSurfacesEL.size()= "
    //             << commonSurfacesEL.size()
    //             << "-> probably Bounding box proximities or CAD Issue " <<
    //             endl;
    //         m_log(WARNING) << commonSurfacesEL.size() << endl;

    //         vector<NekDouble> DistToSurf(commonSurfacesEL.size());
    //         fill(DistToSurf.begin(), DistToSurf.end(), 0.0);

    //         // cout <<  DistToSurf.size() << " initialization of the
    //         vector[0] =
    //         // " << DistToSurf[0] << endl ;

    //         // 0 Calculate the min edge length for the face
    //         NekDouble minEdgeLength = numeric_limits<double>::max();
    //         NekDouble edgeLength;
    //         for (auto edge : el->GetEdgeList())
    //         {
    //             edgeLength = edge->m_n1->Distance(edge->m_n2);
    //             if (edgeLength < minEdgeLength)
    //             {
    //                 minEdgeLength = edgeLength;
    //             }
    //         }

    //         CADSurfSharedPtr Surf;
    //         for (auto edge : el->GetEdgeList())
    //         {
    //             // 1. Calculate edge middle point location
    //             Array<OneD, NekDouble> xyz1 = edge->m_n1->GetLoc();
    //             Array<OneD, NekDouble> xyz2 = edge->m_n2->GetLoc();
    //             Array<OneD, NekDouble> xyz_mid(3);
    //             xyz_mid[0] = (xyz1[0] + xyz2[0]) / 2.0;
    //             xyz_mid[1] = (xyz1[1] + xyz2[1]) / 2.0;
    //             xyz_mid[2] = (xyz1[2] + xyz2[2]) / 2.0;

    //             // 2. Calculate the distance to every of the common surfaces
    //             // Add to the overal DistToSurf[i]
    //             for (int i = 0; i < commonSurfacesEL.size(); i++)
    //             {
    //                 NekDouble dist;
    //                 Surf = m_mesh->m_cad->GetSurf(commonSurfacesEL[i]);
    //                 Surf->locuv(xyz_mid, dist);

    //                 if (dist < minEdgeLength)
    //                 {
    //                     DistToSurf[i] += dist;
    //                     // cout << DistToSurf[i] << endl ;
    //                 }
    //                 else
    //                 {
    //                     m_log(WARNING)
    //                         << "locuv dist= " << dist
    //                         << " minEdgeLength= " << minEdgeLength << endl;
    //                     m_log(WARNING)
    //                         << "The locuv distance on commonSurfacesEL>1 from
    //                         "
    //                            "the edge midpoint to the surface is too big.
    //                            " "The mesh in this region is too coarse for "
    //                            "reconstruction - WARNING instead ? "
    //                         << endl;
    //                 }
    //             }
    //         }

    //         // 3. Sort the distances and pick the minimum surface
    //         auto it = std::minmax_element(DistToSurf.begin(),
    //         DistToSurf.end());

    //         // cout << "The selected surface = " <<
    //         // std::distance(DistToSurf.begin(), it.first) << endl ;

    //         // 4. Assign CADSurf to Element[2] + Edge
    //         Surf = m_mesh->m_cad->GetSurf(
    //             commonSurfacesEL[std::distance(DistToSurf.begin(),
    //             it.first)]);
    //         el->m_parentCAD = Surf;
    //         for (auto edge : el->GetEdgeList())
    //         {
    //             edge->m_parentCAD = Surf;
    //         }
    //     }
    // }
}

} // namespace Nektar::NekMesh
