////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSystem.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <NekMeshUtils/CADSystem/CADSystem.h>
#include <NekMeshUtils/CADSystem/CADObj.h>
#include <NekMeshUtils/CADSystem/CADVert.h>
#include <NekMeshUtils/CADSystem/CADCurve.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

string CADSystem::GetName()
{
    return m_name;
}

void CADSystem::Report()
{
    cout << endl << "CAD report:" << endl;
    cout << "\tCAD has: " << m_curves.size() << " curves." << endl;
    cout << "\tCAD has: " << m_surfs.size() << " surfaces." << endl;
}

Array<OneD, NekDouble> CADSystem::GetBoundingBox()
{
    Array<OneD, NekDouble> bound(6);
    bound[0] = numeric_limits<double>::max(); // xmin
    bound[1] = numeric_limits<double>::min(); // xmax
    bound[2] = numeric_limits<double>::max(); // ymin
    bound[3] = numeric_limits<double>::min(); // ymax
    bound[4] = numeric_limits<double>::max(); // zmin
    bound[5] = numeric_limits<double>::min(); // zmax

    for (int i = 1; i <= m_curves.size(); i++)
    {
        CADCurveSharedPtr c         = GetCurve(i);
        Array<OneD, NekDouble> ends = c->GetMinMax();

        bound[0] = min(bound[0], min(ends[0], ends[3]));
        bound[1] = max(bound[1], max(ends[0], ends[3]));

        bound[2] = min(bound[2], min(ends[1], ends[4]));
        bound[3] = max(bound[3], max(ends[1], ends[4]));

        bound[4] = min(bound[4], min(ends[2], ends[5]));
        bound[5] = max(bound[5], max(ends[2], ends[5]));
    }

    return bound;
}

vector<int> CADSystem::GetBoundarySurfs()
{
    vector<int> ret;

    set<int> surfs;
    vector<CADCurveSharedPtr> cs;

    Array<OneD, NekDouble> bound = GetBoundingBox();

    for (int i = 1; i <= m_curves.size(); i++)
    {
        CADCurveSharedPtr c         = GetCurve(i);
        Array<OneD, NekDouble> ends = c->GetMinMax();

        if((fabs(bound[0] - ends[0]) < 1e-4 ||
            fabs(bound[0] - ends[3]) < 1e-4 ||
            fabs(bound[1] - ends[0]) < 1e-4 ||
            fabs(bound[1] - ends[3]) < 1e-4)
           &&
           (fabs(bound[2] - ends[1]) < 1e-4 ||
            fabs(bound[2] - ends[4]) < 1e-4 ||
            fabs(bound[3] - ends[1]) < 1e-4 ||
            fabs(bound[3] - ends[4]) < 1e-4)
           &&
           (fabs(bound[4] - ends[2]) < 1e-4 ||
            fabs(bound[4] - ends[5]) < 1e-4 ||
            fabs(bound[5] - ends[2]) < 1e-4 ||
            fabs(bound[5] - ends[5]) < 1e-4) )
        {
            //curve touches on bounding box
            cs.push_back(c);
        }

    }

    for(int i = 0; i < cs.size(); i++)
    {
        vector<CADSurfSharedPtr> s = cs[i]->GetAdjSurf();
        for(int j = 0; j < s.size(); j++)
        {
            surfs.insert(s[j]->GetId());
        }
    }

    set<int>::iterator it;
    for(it = surfs.begin(); it != surfs.end(); it++)
    {
        ret.push_back(*it);
    }

    return ret;
}


bool CADSystem::LoadCAD()
{
    if (!boost::filesystem::exists(m_name.c_str()))
    {
        return false;
    }

    string ext;
    size_t pos = m_name.find(".");
    ext        = m_name.substr(pos);

    if (boost::iequals(ext, ".STEP") || boost::iequals(ext, ".STP"))
    {
        // Takes step file and makes OpenCascade shape
        STEPControl_Reader reader;
        reader = STEPControl_Reader();
        reader.ReadFile(m_name.c_str());
        reader.NbRootsForTransfer();
        reader.TransferRoots();
        shape = reader.OneShape();
        if (shape.IsNull())
        {
            return false;
        }
    }
    else if (boost::iequals(ext, ".IGES") || boost::iequals(ext, ".IGS"))
    {
        // Takes IGES file and makes OpenCascade shape
        IGESControl_Reader reader;
        reader = IGESControl_Reader();
        reader.ReadFile(m_name.c_str());
        reader.NbRootsForTransfer();
        reader.TransferRoots();
        shape = reader.OneShape();
        if (shape.IsNull())
        {
            return false;
        }
    }
    else
    {
        return false;
    }

    // faces and verts can be extracted straight from shape
    TopTools_IndexedMapOfShape mapOfVerts, mapOfFaces;
    TopExp::MapShapes(shape, TopAbs_VERTEX, mapOfVerts);
    TopExp::MapShapes(shape, TopAbs_FACE, mapOfFaces);

    // edges need to be built from loops around faces to elimiate degen and
    // hanging edges
    TopTools_IndexedMapOfShape mapOfEdges;

    // build map of verticies
    for (int i = 1; i <= mapOfVerts.Extent(); i++)
    {
        TopoDS_Shape v = mapOfVerts.FindKey(i);
        AddVert(i, v);
    }

    // For each face of the geometry, get the local edges which bound it. If
    // they are valid (their type != 7), then add them to an edge map. This
    // filters out the dummy edges which OCC uses.
    for (int i = 1; i <= mapOfFaces.Extent(); i++)
    {
        TopoDS_Shape face = mapOfFaces.FindKey(i);

        TopTools_IndexedMapOfShape localEdges;
        TopExp::MapShapes(face, TopAbs_EDGE, localEdges);

        for (int j = 1; j <= localEdges.Extent(); j++)
        {
            TopoDS_Shape edge       = localEdges.FindKey(j);
            BRepAdaptor_Curve curve = BRepAdaptor_Curve(TopoDS::Edge(edge));
            if (curve.GetType() != 7)
            {
                if (!(mapOfEdges.Contains(edge)))
                {
                    mapOfEdges.Add(edge);
                }
            }
        }
    }

    map<int, vector<int> > adjsurfmap; // from id of curve to list of ids of
                                       // surfs

    // Adds edges to our type and map
    for (int i = 1; i <= mapOfEdges.Extent(); i++)
    {
        TopoDS_Shape edge = mapOfEdges.FindKey(i);
        TopoDS_Vertex fv =
            TopExp::FirstVertex(TopoDS::Edge(edge), Standard_True);
        TopoDS_Vertex lv =
            TopExp::LastVertex(TopoDS::Edge(edge), Standard_True);

        if (edge.Orientation() == 0)
        {
            AddCurve(i, edge, mapOfVerts.FindIndex(fv),
                     mapOfVerts.FindIndex(lv));
        }
        else
        {
            AddCurve(i, edge, mapOfVerts.FindIndex(lv),
                     mapOfVerts.FindIndex(fv));
        }
    }

    // For each face, examine all the wires (i.e. bounding loops) and
    // investigates the loop. Using this information, connectivity is determined
    // and edges are associated with surfaces.
    for (int i = 1; i <= mapOfFaces.Extent(); i++)
    {
        TopoDS_Shape face = mapOfFaces.FindKey(i);

        TopTools_IndexedMapOfShape mapOfWires;
        TopExp::MapShapes(face, TopAbs_WIRE, mapOfWires);

        // this pice of code does an idiot check on the loops to make sure
        // they dont cross or touch
        if (mapOfWires.Extent() > 1)
        {
            TopoDS_Wire ow = BRepTools::OuterWire(TopoDS::Face(face));

            vector<TopoDS_Shape> wirefacecuts;
            vector<gp_Pnt> centersofcutfaces;

            for (int j = 1; j <= mapOfWires.Extent(); j++)
            {
                TopoDS_Shape wire = mapOfWires.FindKey(j);

                if (wire != ow)
                {
                    BRepBuilderAPI_MakeFace build(
                        BRep_Tool::Surface(TopoDS::Face(face)), 1e-7);
                    build.Add(TopoDS::Wire(wire));
                    TopoDS_Shape newface = build.Shape();
                    wirefacecuts.push_back(newface);
                    BRepAdaptor_Surface b =
                        BRepAdaptor_Surface(TopoDS::Face(newface));
                    NekDouble u, v;
                    u = (b.LastUParameter() - b.FirstUParameter()) / 2.0;
                    v = (b.LastVParameter() - b.FirstVParameter()) / 2.0;
                    centersofcutfaces.push_back(b.Value(u, v));
                }
            }
            for (int j = 0; j < wirefacecuts.size(); j++)
            {
                for (int k = 0; k < wirefacecuts.size(); k++)
                {
                    if (j == k)
                        continue;

                    /// TODO fix this test
                    BRepClass_FaceClassifier fc(TopoDS::Face(wirefacecuts[j]),
                                                centersofcutfaces[k], 1e-7);
                    // ASSERTL0(fc.State() == 1, "Internal face loops make this
                    // cad impossible to mesh");
                }
            }
        }

        vector<EdgeLoop> edgeloops;

        // now we acutally analyse the loops for cad building
        for (int j = 1; j <= mapOfWires.Extent(); j++)
        {
            EdgeLoop edgeloop;

            TopoDS_Shape wire = mapOfWires.FindKey(j);

            ShapeAnalysis_Wire wiretest(TopoDS::Wire(wire), TopoDS::Face(face),
                                        1E-6);

            // calculate the center of the wire
            GProp_GProps massProps;
            BRepGProp::SurfaceProperties(wire, massProps);
            gp_Pnt gPt = massProps.CentreOfMass();
            // alternative locuv methods
            ShapeAnalysis_Surface sas(BRep_Tool::Surface(TopoDS::Face(face)));
            sas.SetDomain(
                BRepAdaptor_Surface(TopoDS::Face(face)).FirstUParameter(),
                BRepAdaptor_Surface(TopoDS::Face(face)).LastUParameter(),
                BRepAdaptor_Surface(TopoDS::Face(face)).FirstVParameter(),
                BRepAdaptor_Surface(TopoDS::Face(face)).LastVParameter());
            gp_Pnt2d p2 = sas.ValueOfUV(gPt, 1e-7);
            Array<OneD, NekDouble> cen(2);
            cen[0]          = p2.X();
            cen[1]          = p2.Y();
            edgeloop.center = cen;

            BRepTools_WireExplorer exp;

            exp.Init(TopoDS::Wire(wire));

            while (exp.More())
            {
                TopoDS_Shape edge = exp.Current();

                if (mapOfEdges.Contains(edge))
                {
                    int e = mapOfEdges.FindIndex(edge);
                    edgeloop.edges.push_back(m_curves[e]);
                    edgeloop.edgeo.push_back(exp.Orientation());
                    adjsurfmap[e].push_back(i);
                }

                exp.Next();
            }

            edgeloops.push_back(edgeloop);
        }

        AddSurf(i, face, edgeloops);
    }

    // attempts to identify properties of the vertex on the degen edge
    for (int i = 1; i <= mapOfFaces.Extent(); i++)
    {
        TopoDS_Shape face = mapOfFaces.FindKey(i);

        TopTools_IndexedMapOfShape localEdges;
        TopExp::MapShapes(face, TopAbs_EDGE, localEdges);

        for (int j = 1; j <= localEdges.Extent(); j++)
        {
            TopoDS_Shape edge = localEdges.FindKey(j);
            if (BRep_Tool::Degenerated(TopoDS::Edge(edge)))
            {
                gp_Pnt2d p1, p2;

                BRep_Tool::UVPoints(TopoDS::Edge(edge), TopoDS::Face(face), p1,
                                    p2);

                m_verts[mapOfVerts.FindIndex(TopExp::FirstVertex(
                            TopoDS::Edge(edge), Standard_True))]
                    ->SetDegen(i, m_surfs[i], (p1.X() + p2.X()) / 2.0,
                               (p1.Y() + p2.Y()) / 2.0);
            }
        }
    }

    // This checks that all edges are bound by two surfaces, sanity check.
    for (map<int, vector<int> >::iterator it = adjsurfmap.begin();
         it != adjsurfmap.end(); it++)
    {
        ASSERTL0(it->second.size() == 2, "no three curve surfaces");
        vector<CADSurfSharedPtr> sfs;
        for (int i = 0; i < it->second.size(); i++)
        {
            sfs.push_back(m_surfs[it->second[i]]);
        }
        m_curves[it->first]->SetAdjSurf(sfs);
    }
    return true;
}

void CADSystem::AddVert(int i, TopoDS_Shape in)
{
    CADVertSharedPtr newVert = MemoryManager<CADVert>::AllocateSharedPtr(i, in);

    m_verts[i] = newVert;
}

void CADSystem::AddCurve(int i, TopoDS_Shape in, int fv, int lv)
{
    CADCurveSharedPtr newCurve =
        MemoryManager<CADCurve>::AllocateSharedPtr(i, in);

    vector<CADVertSharedPtr> vs;
    vs.push_back(m_verts[fv]);
    vs.push_back(m_verts[lv]);
    m_curves[i] = newCurve;
    m_curves[i]->SetVert(vs);
}

void CADSystem::AddSurf(int i, TopoDS_Shape in, vector<EdgeLoop> ein)
{
    CADSurfSharedPtr newSurf =
        MemoryManager<CADSurf>::AllocateSharedPtr(i, in, ein);
    m_surfs[i] = newSurf;

    if (in.Orientation() == 0)
    {
        m_surfs[i]->SetReverseNomral();
    }

    int tote = 0;
    for (int i = 0; i < ein.size(); i++)
    {
        tote += ein[i].edges.size();
    }

    ASSERTL0(tote != 1, "cannot handle periodic curves");

    if (tote == 2)
    {
        m_surfs[i]->SetTwoC();
    }
}

bool CADSystem::InsideShape(Array<OneD, NekDouble> loc)
{
    gp_Pnt p(loc[0] * 1000.0, loc[1] * 1000.0, loc[2] * 1000.0);

    BRepClass3d_SolidClassifier test(shape, p, 1E-7);
    if (test.State() == TopAbs_IN || test.State() == TopAbs_ON)
    {
        return true;
    }
    else
    {
        return false;
    }
}
}
}
