///////////////////////////////////////////////////////////////////////////////
//
// File: ShallowWaterSystem.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Generic timestepping for shallow water solvers
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_SHALLOWWATERSYSTEM_H
#define NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_SHALLOWWATERSYSTEM_H

#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

namespace Nektar
{

/// Base class for unsteady solvers.
class ShallowWaterSystem : public SolverUtils::UnsteadySystem
{
public:
    friend class MemoryManager<ShallowWaterSystem>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<ShallowWaterSystem>::AllocateSharedPtr(pSession,
                                                                 pGraph);
        p->InitObject();
        return p;
    }
    /// Name of class
    static std::string className;

    /// Destructor
    ~ShallowWaterSystem() override = default;

protected:
    SolverUtils::RiemannSolverSharedPtr m_riemannSolver;
    SolverUtils::AdvectionSharedPtr m_advection;
    SolverUtils::DiffusionSharedPtr m_diffusion;

    /// Indicates if variables are primitive or conservative
    bool m_primitive;
    /// Indicates if constant depth case
    bool m_constantDepth;
    /// Acceleration of gravity
    NekDouble m_g;
    /// Still water depth
    Array<OneD, NekDouble> m_depth;
    // Bottom slopes
    Array<OneD, Array<OneD, NekDouble>> m_bottomSlope;
    /// Coriolis force
    Array<OneD, NekDouble> m_coriolis;
    // Location of velocity vector.
    Array<OneD, Array<OneD, NekDouble>> m_vecLocs;

    /// Initialises UnsteadySystem class members.
    ShallowWaterSystem(const LibUtilities::SessionReaderSharedPtr &pSession,
                       const SpatialDomains::MeshGraphSharedPtr &pGraph);

    /// Print a summary of time stepping parameters.
    void v_GenerateSummary(SolverUtils::SummaryList &s) override;

    void v_InitObject(bool DeclareFields = true) override;

    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble>> &physarray,
                               NekDouble time);

    void WallBoundary2D(int bcRegion, int cnt,
                        Array<OneD, Array<OneD, NekDouble>> &Fwd);

    void WallBoundary(int bcRegion, int cnt,
                      Array<OneD, Array<OneD, NekDouble>> &Fwd,
                      Array<OneD, Array<OneD, NekDouble>> &physarray);

    void AddCoriolis(const Array<OneD, const Array<OneD, NekDouble>> &physarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray);

    void PrimitiveToConservative();

    void ConservativeToPrimitive();

    NekDouble GetGravity()
    {
        return m_g;
    }

    const Array<OneD, const Array<OneD, NekDouble>> &GetVecLocs()
    {
        return m_vecLocs;
    }

    const Array<OneD, const Array<OneD, NekDouble>> &GetNormals()
    {
        return m_traceNormals;
    }

    const Array<OneD, NekDouble> &GetDepth()
    {
        return m_depth;
    }

    bool IsConstantDepth()
    {
        return m_constantDepth;
    }

    void CopyBoundaryTrace(const Array<OneD, NekDouble> &Fwd,
                           Array<OneD, NekDouble> &Bwd);

private:
    ///
    void EvaluateWaterDepth(void);

    ///
    void EvaluateCoriolis(void);
};
} // namespace Nektar

#endif
