////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMesh.h
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
//  Description: class for indivdual surface meshes
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMeshUtils_SURFACEMESHING_FACEMESH
#define NekMeshUtils_SURFACEMESHING_FACEMESH

#include <NekMeshUtils/MeshElements/Mesh.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>
#include <NekMeshUtils/Octree/Octree.h>
#include <NekMeshUtils/SurfaceMeshing/CurveMesh.h>

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief class for surface meshes on individual surfaces (paramter plane
 * meshes)
 */
class FaceMesh
{
public:
    friend class MemoryManager<FaceMesh>;

    /**
     * @brief Default constructor
     */
    FaceMesh(const int                                id,
             MeshSharedPtr                            m,
             CADSurfSharedPtr                         cad,
             OctreeSharedPtr                          oct,
             const std::map<int, CurveMeshSharedPtr> &cmeshes,
             bool ov)
        : m_mesh(m), m_cadsurf(cad), m_octree(oct), m_curvemeshes(cmeshes),
          m_id(id), over(ov)

    {
        m_edgeloops = m_cadsurf->GetEdges();
    };

    /**
     * @brief mesh exectuation command
     */
    void Mesh(bool remesh = false);

    void QuadRemesh(std::map<NodeSharedPtr, NodeSharedPtr> nmap);

private:
    /**
     * @brief Calculate the paramter plane streching factor
     */
    void Stretching();

    /**
     * @brief performs node smoothing on face
     */
    void Smoothing();

    /**
     * @brief performs diagonal swapping of edges
     */
    void DiagonalSwap();

    /**
     * @brief build a local version of mesh elements
     */
    void BuildLocalMesh();

    /**
     * @brief function which calls the optimisation routines
     */
    void OptimiseLocalMesh();

    /**
     * @brief Validate the surface mesh base on the octree and real
     * dimensions of the edges
     */
    bool Validate();

    /**
     * @brief Get the boundries of the surface and extracts the nodes from
     * the curve meshes in the correct order
     */
    void OrientateCurves();

    /**
     * @brief adds a new stiener point to the triangulation for meshing
     */
    void AddNewPoint(Array<OneD, NekDouble> uv);

    /**
     * @brief adds a quad layer around any interior loops
     */
    void MakeBL();

    /// mesh pointer
    MeshSharedPtr m_mesh;
    /// CAD surface
    CADSurfSharedPtr m_cadsurf;
    /// Octree object
    OctreeSharedPtr m_octree;
    /// Map of the curve meshes which bound the surfaces
    std::map<int, CurveMeshSharedPtr> m_curvemeshes;
    /// data structure containing the edges, their order and oreientation for
    /// the surface
    std::vector<EdgeLoop> m_edgeloops;
    /// id of the surface mesh
    int m_id;
    /// list of boundary nodes in their order loops
    std::vector<std::vector<NodeSharedPtr> > orderedLoops;
    /// list of stiener points in the triangulation
    std::vector<NodeSharedPtr> m_stienerpoints;
    /// pplane stretching
    NekDouble m_str;
    /// triangle connectiviities
    std::vector<std::vector<NodeSharedPtr> > m_connec;
    /// local set of nodes
    NodeSet m_localNodes;
    /// local set of edges
    EdgeSet m_localEdges;
    /// local list of elements
    std::vector<ElementSharedPtr> m_localElements;

    bool over;
};

typedef boost::shared_ptr<FaceMesh> FaceMeshSharedPtr;
}
}

#endif
