///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingReferenceFrame.h
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
// Description: Allows for a moving frame of reference, through adding c * du/dx
// to the body force, where c is the frame velocity vector
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FORCINGMOVINGREFERENCEFRAME
#define NEKTAR_SOLVERUTILS_FORCINGMOVINGREFERENCEFRAME

#include <string>

#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Filters/FilterAeroForces.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>

namespace Nektar::SolverUtils
{

namespace bn = boost::numeric;

/***
 * Solve the body's motion using Newmark-Beta method
 * M ddx + C dx + K x = F
 * In discrete form
 * CoeffMatrix dx = rhs
 ***/
class Newmark_BetaSolver
{
public:
    Newmark_BetaSolver(){};
    ~Newmark_BetaSolver(){};
    void SetNewmarkBeta(NekDouble beta, NekDouble gamma, NekDouble dt,
                        DNekMatSharedPtr M, DNekMatSharedPtr C,
                        DNekMatSharedPtr K, std::set<int> DirDoFs);
    void Solve(Array<OneD, Array<OneD, NekDouble>> u,
               Array<OneD, NekDouble> force,
               std::map<int, NekDouble> motionPrescribed);
    int m_rows;
    Array<OneD, NekDouble> m_coeffs;
    DNekMatSharedPtr m_M;
    DNekMatSharedPtr m_C;
    DNekMatSharedPtr m_K;
    DNekMatSharedPtr m_coeffMatrix;
    DNekMatSharedPtr m_inverseMatrix;
};

class FrameTransform
{
public:
    FrameTransform();
    ~FrameTransform(){};
    void SetAngle(const Array<OneD, NekDouble> theta);
    void BodyToInerital(const int dim, const Array<OneD, NekDouble> &body,
                        Array<OneD, NekDouble> &inertial);
    void IneritalToBody(const int dim, Array<OneD, NekDouble> &body,
                        const Array<OneD, NekDouble> &inertial);

private:
    Array<OneD, NekDouble> m_matrix;
};

class ForcingMovingReferenceFrame : public Forcing
{

public:
    friend class MemoryManager<ForcingMovingReferenceFrame>;

    /// Creates an instance of this class
    SOLVER_UTILS_EXPORT static ForcingSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int &pNumForcingFields, const TiXmlElement *pForce)
    {
        ForcingSharedPtr p =
            MemoryManager<ForcingMovingReferenceFrame>::AllocateSharedPtr(
                pSession, pEquation);
        p->InitObject(pFields, pNumForcingFields, pForce);
        return p;
    }

    /// Name of the class
    static std::string classNameBody;

protected:
    SOLVER_UTILS_EXPORT void v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int &pNumForcingFields,
        const TiXmlElement *pForce) override;

    SOLVER_UTILS_EXPORT void v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble &time) override;

    SOLVER_UTILS_EXPORT void v_PreApply(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const NekDouble &time) override;

private:
    // name of the function for linear and angular velocities in the session
    // file
    // pivot point
    Array<OneD, NekDouble> m_pivotPoint;
    Array<OneD, NekDouble> m_travelWave;
    FrameTransform m_frame;

    // prescribed functions in the session file
    std::map<int, LibUtilities::EquationSharedPtr> m_frameVelFunction;
    std::map<int, LibUtilities::EquationSharedPtr> m_extForceFunction;
    std::ofstream m_outputStream;
    int m_rank;

    // a boolean switch indicating for which direction the velocities are
    // available. The available velocites could be different from the
    // precscribed one because of the rotation which result in change of basis
    // vector of local frame to the inertial frame.
    Array<OneD, bool> m_hasVel;
    Array<OneD, bool> m_hasOmega;

    // frame linear velocities in inertial frame
    Array<OneD, NekDouble> m_velXYZ;

    // frame linear velocities in local translating-rotating frame
    Array<OneD, NekDouble> m_velxyz;

    // frame angular velocities in inertial frame
    Array<OneD, NekDouble> m_omegaXYZ;

    // frame angular velocities in local translating-rotating frame
    Array<OneD, NekDouble> m_omegaxyz;
    // coordinate vector
    Array<OneD, Array<OneD, NekDouble>> m_coords;
    // externel force
    Array<OneD, NekDouble> m_extForceXYZ;

    // rotation angel
    Array<OneD, NekDouble> m_disp;

    // Projection matrix for transformation of vectors between inertial and
    // moving reference frames
    bn::ublas::matrix<NekDouble> m_ProjMatZ;

    NekDouble m_startTime;
    NekDouble m_timestep;

    bool m_hasRotation;
    bool m_isH1d;
    bool m_hasPlane0;
    bool m_isH2d;
    bool m_hasFreeMotion;
    NekInt m_spacedim;
    NekInt m_expdim;
    unsigned int m_index;
    unsigned int m_outputFrequency;
    Newmark_BetaSolver m_bodySolver;
    Array<OneD, Array<OneD, NekDouble>> m_bodyVel;
    std::set<int> m_DirDoFs;
    FilterAeroForcesSharedPtr m_aeroforceFilter;

    ForcingMovingReferenceFrame(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<EquationSystem> &pEquation);

    ~ForcingMovingReferenceFrame(void) override;

    void Update(const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                const NekDouble &time);
    void UpdateRotMat();
    void CheckForRestartTime(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        NekDouble &time);

    void addRotation(int npoints,
                     const Array<OneD, Array<OneD, NekDouble>> &inarray0,
                     NekDouble angVelScale,
                     const Array<OneD, Array<OneD, NekDouble>> &inarray1,
                     Array<OneD, Array<OneD, NekDouble>> &outarray);
    void InitBodySolver(const TiXmlElement *pForce, const int dim,
                        const int rank, const NekDouble time);
    void SolveBodyMotion(Array<OneD, Array<OneD, NekDouble>> &bodyVel,
                         const Array<OneD, NekDouble> &forcebody,
                         std::map<int, NekDouble> &Dirs);
    void LoadParameters(const TiXmlElement *pForce, const NekDouble time);
    NekDouble EvaluateExpression(std::string expression);
    void InitialiseFilter(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const TiXmlElement *pForce);
    void UpdateBoundaryConditions(NekDouble time);
};

} // namespace Nektar::SolverUtils

#endif
