///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesCFE.cpp
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
// Description: Navier Stokes equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/EquationSystems/NavierStokesCFE.h>

using namespace std;

namespace Nektar
{
    string NavierStokesCFE::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "NavierStokesCFE", NavierStokesCFE::create,
            "NavierStokes equations in conservative variables.");

    NavierStokesCFE::NavierStokesCFE(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          CompressibleFlowSystem(pSession, pGraph)
    {
    }

    NavierStokesCFE::~NavierStokesCFE()
    {

    }

    /**
     * @brief Initialization object for CompressibleFlowSystem class.
     */
    void NavierStokesCFE::v_InitObject()
    {
        CompressibleFlowSystem::v_InitObject();

        // rest of initialisation is in this routine so it can also be called
        // in NavierStokesImplicitCFE initialisation
        InitObject_Explicit();
    }

    void NavierStokesCFE::InitObject_Explicit()
    {
        // Get gas constant from session file and compute Cp
        NekDouble gasConstant;
        m_session->LoadParameter ("GasConstant",   gasConstant,   287.058);
        m_Cp      = m_gamma / (m_gamma - 1.0) * gasConstant;
        m_Cv      = m_Cp / m_gamma;

        m_session->LoadParameter ("Twall", m_Twall, 300.15);

        // Viscosity
        int nPts = m_fields[0]->GetNpoints();
        m_session->LoadSolverInfo("ViscosityType", m_ViscosityType, "Constant");
        m_session->LoadParameter ("mu",            m_muRef,           1.78e-05);
        m_mu = Array<OneD, NekDouble>(nPts, m_muRef);
        if (m_ViscosityType == "Variable")
        {
            m_is_mu_variable = true;
        }

        // Thermal conductivity or Prandtl
        if ( m_session->DefinesParameter("thermalConductivity"))
        {
            ASSERTL0( !m_session->DefinesParameter("Pr"),
                 "Cannot define both Pr and thermalConductivity.");

            m_session->LoadParameter ("thermalConductivity",
                                        m_thermalConductivityRef);
            m_Prandtl = m_Cp * m_muRef / m_thermalConductivityRef;
        }
        else
        {
            m_session->LoadParameter ("Pr", m_Prandtl, 0.72);
            m_thermalConductivityRef = m_Cp * m_muRef / m_Prandtl;
        }
        m_thermalConductivity =
                Array<OneD, NekDouble>(nPts, m_thermalConductivityRef);

        // Artificial viscosity parameter
        m_session->LoadParameter("mu0", m_mu0, 1.0);

        // load smoothing tipe
        m_session->LoadSolverInfo("Smoothing", m_smoothing, "Off");
        if (m_smoothing == "C0")
        {
            m_C0ProjectExp = MemoryManager<MultiRegions::ContField>::
            AllocateSharedPtr(m_session,m_graph,m_session->GetVariable(0));
        }
        // load physical sensor type
        m_session->LoadSolverInfo("PhysicalSensorType", m_physicalSensorType,
            "Off");

        string diffName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDGNS");

        m_diffusion = SolverUtils::GetDiffusionFactory()
                                    .CreateInstance(diffName, diffName);
        if ("InteriorPenalty" == diffName)
        {
            m_is_diffIP = true;
            SetBoundaryConditionsBwdWeight();
        }

        if ("LDGNS" == diffName||
            "LDGNS3DHomogeneous1D" == diffName)
        {
            m_diffusion->SetFluxPenaltyNS(&NavierStokesCFE::
                v_GetFluxPenalty, this);
        }

        if (m_specHP_dealiasing)
        {
            m_diffusion->SetFluxVectorNS(
                &NavierStokesCFE::v_GetViscousFluxVectorDeAlias,
                this);
        }
        else
        {
            m_diffusion->SetFluxVectorNS(&NavierStokesCFE::
                                          v_GetViscousFluxVector, this);
        }

        m_diffusion->SetDiffusionFluxCons(
            &NavierStokesCFE::GetViscousFluxVectorConservVar<false>, this);

        m_diffusion->SetDiffusionFluxConsTrace(
            &NavierStokesCFE::GetViscousFluxVectorConservVar<true>, this);

        m_diffusion->SetSpecialBndTreat(
            &NavierStokesCFE::SpecialBndTreat, this);

        m_diffusion->SetDiffusionSymmFluxCons(
            &NavierStokesCFE::GetViscousSymmtrFluxConservVar, this);

        if (m_shockCaptureType != "Off")
        {
            m_diffusion->SetArtificialDiffusionVector(
                &NavierStokesCFE::GetArtificialViscosity, this);
        }

        m_diffusion->SetCalcViscosity(
                &NavierStokesCFE::CalcViscosity, this);

        // Concluding initialisation of diffusion operator
        m_diffusion->InitObject         (m_session, m_fields);

    }

    void NavierStokesCFE::v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        bool extraFields;
        m_session->MatchSolverInfo("OutputExtraFields","True",
                                   extraFields, true);
        if (extraFields)
        {
            const int nPhys   = m_fields[0]->GetNpoints();
            const int nCoeffs = m_fields[0]->GetNcoeffs();
            Array<OneD, Array<OneD, NekDouble> > tmp(m_fields.size());

            for (int i = 0; i < m_fields.size(); ++i)
            {
                tmp[i] = m_fields[i]->GetPhys();
            }

            Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);
            Array<OneD, Array<OneD, NekDouble> > velFwd  (m_spacedim);
            for (int i = 0; i < m_spacedim; ++i)
            {
                velocity[i] = Array<OneD, NekDouble> (nPhys);
                velFwd[i]   = Array<OneD, NekDouble> (nCoeffs);
            }

            Array<OneD, NekDouble> pressure(nPhys), temperature(nPhys);
            Array<OneD, NekDouble> entropy(nPhys);
            Array<OneD, NekDouble> soundspeed(nPhys), mach(nPhys);
            Array<OneD, NekDouble> sensor(nPhys), SensorKappa(nPhys);

            m_varConv->GetVelocityVector(tmp, velocity);
            m_varConv->GetPressure  (tmp, pressure);
            m_varConv->GetTemperature(tmp, temperature);
            m_varConv->GetEntropy   (tmp, entropy);
            m_varConv->GetSoundSpeed(tmp, soundspeed);
            m_varConv->GetMach      (tmp, soundspeed, mach);

            int sensorOffset;
            m_session->LoadParameter ("SensorOffset", sensorOffset, 1);
            m_varConv->GetSensor (m_fields[0], tmp, sensor, SensorKappa,
                                    sensorOffset);

            Array<OneD, NekDouble> pFwd(nCoeffs), TFwd(nCoeffs);
            Array<OneD, NekDouble> sFwd(nCoeffs);
            Array<OneD, NekDouble> aFwd(nCoeffs), mFwd(nCoeffs);
            Array<OneD, NekDouble> sensFwd(nCoeffs);

            string velNames[3] = {"u", "v", "w"};
            for (int i = 0; i < m_spacedim; ++i)
            {
                m_fields[0]->FwdTrans_IterPerExp(velocity[i], velFwd[i]);
                variables.push_back(velNames[i]);
                fieldcoeffs.push_back(velFwd[i]);
            }

            m_fields[0]->FwdTrans_IterPerExp(pressure,   pFwd);
            m_fields[0]->FwdTrans_IterPerExp(temperature,TFwd);
            m_fields[0]->FwdTrans_IterPerExp(entropy,    sFwd);
            m_fields[0]->FwdTrans_IterPerExp(soundspeed, aFwd);
            m_fields[0]->FwdTrans_IterPerExp(mach,       mFwd);
            m_fields[0]->FwdTrans_IterPerExp(sensor,     sensFwd);

            variables.push_back  ("p");
            variables.push_back  ("T");
            variables.push_back  ("s");
            variables.push_back  ("a");
            variables.push_back  ("Mach");
            variables.push_back  ("Sensor");
            fieldcoeffs.push_back(pFwd);
            fieldcoeffs.push_back(TFwd);
            fieldcoeffs.push_back(sFwd);
            fieldcoeffs.push_back(aFwd);
            fieldcoeffs.push_back(mFwd);
            fieldcoeffs.push_back(sensFwd);

            if (m_artificialDiffusion)
            {
                // reuse pressure
                Array<OneD, NekDouble> sensorFwd(nCoeffs);
                m_artificialDiffusion->GetArtificialViscosity(tmp, pressure);
                m_fields[0]->FwdTrans_IterPerExp(pressure,   sensorFwd);

                variables.push_back  ("ArtificialVisc");
                fieldcoeffs.push_back(sensorFwd);

            }

            if (m_shockCaptureType == "Physical")
            {
                Array<OneD, NekDouble> muavFwd(nCoeffs);
                m_fields[0]->FwdTrans_IterPerExp(m_muav,   muavFwd);
                variables.push_back  ("ArtificialVisc");
                fieldcoeffs.push_back(muavFwd);

                // Debug Ducros
                // div square
                Array<OneD, NekDouble> dv2Fwd(nCoeffs);
                m_fields[0]->FwdTrans_IterPerExp(m_diffusion->m_divVelSquare,
                    dv2Fwd);
                variables.push_back  ("divVelSquare");
                fieldcoeffs.push_back(dv2Fwd);
                // curl square
                Array<OneD, NekDouble> cv2Fwd(nCoeffs);
                m_fields[0]->FwdTrans_IterPerExp(m_diffusion->m_curlVelSquare,
                    cv2Fwd);
                variables.push_back  ("curlVelSquare");
                fieldcoeffs.push_back(cv2Fwd);
                // Ducros
                Array<OneD, NekDouble> duc(nPhys,1.0);
                Ducros(duc);
                Array<OneD, NekDouble> ducFwd(nCoeffs);
                m_fields[0]->FwdTrans_IterPerExp(duc, ducFwd);
                variables.push_back  ("Ducros");
                fieldcoeffs.push_back(ducFwd);

            }
        }
    }

    void NavierStokesCFE::v_DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              Array<OneD,       Array<OneD, NekDouble>> &outarray,
        const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &pBwd)
    {
        size_t nvariables = inarray.size();
        size_t npoints    = GetNpoints();
        size_t nTracePts  = GetTraceTotPoints();

        // this should be preallocated
        Array<OneD, Array<OneD, NekDouble> > outarrayDiff(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(npoints, 0.0);
        }

        // Set artificial viscosity based on NS viscous tensor
        if (m_shockCaptureType == "Physical")
        {
            Array<OneD, NekDouble> div(npoints), curlSquare(npoints);
            GetDivCurlSquared(m_fields, inarray, div, curlSquare,
                pFwd, pBwd);

            // Set volume and trace artificial viscosity
            m_varConv->SetAv(m_fields, inarray, div, curlSquare);
        }

        if (m_is_diffIP)
        {
            if (m_bndEvaluateTime < 0.0)
            {
                NEKERROR(ErrorUtil::efatal, "m_bndEvaluateTime not setup");
            }
            m_diffusion->Diffuse(nvariables, m_fields, inarray, outarrayDiff,
                                m_bndEvaluateTime,
                                pFwd, pBwd);
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarrayDiff[i], 1,
                            outarray[i], 1,
                            outarray[i], 1);
            }
        }
        else
        {
            // Get primitive variables [u,v,w,T]
            Array<OneD, Array<OneD, NekDouble> > inarrayDiff(nvariables-1);
            Array<OneD, Array<OneD, NekDouble> > inFwd(nvariables-1);
            Array<OneD, Array<OneD, NekDouble> > inBwd(nvariables-1);


            for (int i = 0; i < nvariables - 1; ++i)
            {
                inarrayDiff[i] = Array<OneD, NekDouble>{npoints};
                inFwd[i]       = Array<OneD, NekDouble>{nTracePts};
                inBwd[i]       = Array<OneD, NekDouble>{nTracePts};
            }

            // Extract pressure
            // (use inarrayDiff[0] as a temporary storage for the pressure)
            m_varConv->GetPressure(inarray, inarrayDiff[0]);

            // Extract temperature
            m_varConv->GetTemperature(inarray, inarrayDiff[nvariables - 2]);

            // Extract velocities
            m_varConv->GetVelocityVector(inarray, inarrayDiff);

            // Repeat calculation for trace space
            if (pFwd == NullNekDoubleArrayOfArray ||
                pBwd == NullNekDoubleArrayOfArray)
            {
                inFwd = NullNekDoubleArrayOfArray;
                inBwd = NullNekDoubleArrayOfArray;
            }
            else
            {
                m_varConv->GetPressure(pFwd, inFwd[0]);
                m_varConv->GetPressure(pBwd, inBwd[0]);

                m_varConv->GetTemperature(pFwd, inFwd[nvariables - 2]);
                m_varConv->GetTemperature(pBwd, inBwd[nvariables - 2]);

                m_varConv->GetVelocityVector(pFwd, inFwd);
                m_varConv->GetVelocityVector(pBwd, inBwd);
            }

            // Diffusion term in physical rhs form
            m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff,
                                outarrayDiff, inFwd, inBwd);

            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Vadd(npoints,
                            outarrayDiff[i], 1,
                            outarray[i], 1,
                            outarray[i], 1);
            }
        }
    }


    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void NavierStokesCFE::v_GetViscousFluxVector(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
              TensorOfArray3D<NekDouble>                &derivativesO1,
              TensorOfArray3D<NekDouble>                &viscousTensor)
    {
        // Auxiliary variables
        size_t nScalar    = physfield.size();
        size_t nPts       = physfield[0].size();
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Update viscosity and thermal conductivity
        GetViscosityAndThermalCondFromTemp(physfield[nScalar-1], m_mu,
            m_thermalConductivity);

        // Velocity divergence
        for (int j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, divVel, 1, derivativesO1[j][j], 1,
                        divVel, 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, m_mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)
        for (int i = 0; i < m_spacedim; ++i)
        {
            for (int j = i; j < m_spacedim; ++j)
            {
                Vmath::Vadd(nPts, derivativesO1[i][j], 1,
                                  derivativesO1[j][i], 1,
                                  viscousTensor[i][j+1], 1);

                Vmath::Vmul(nPts, m_mu, 1,
                                  viscousTensor[i][j+1], 1,
                                  viscousTensor[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, viscousTensor[i][j+1], 1,
                                  divVel, 1,
                                  viscousTensor[i][j+1], 1);
                }
                else
                {
                    // Copy to make symmetric
                    Vmath::Vcopy(nPts, viscousTensor[i][j+1], 1,
                                       viscousTensor[j][i+1], 1);
                }
            }
        }

        // Terms for the energy equation
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, viscousTensor[i][m_spacedim+1], 1);
            // u_j * tau_ij
            for (int j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vvtvp(nPts, physfield[j], 1,
                               viscousTensor[i][j+1], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
            }
            // Add k*T_i
            Vmath::Vvtvp(nPts, m_thermalConductivity, 1,
                               derivativesO1[i][m_spacedim], 1,
                               viscousTensor[i][m_spacedim+1], 1,
                               viscousTensor[i][m_spacedim+1], 1);
        }
    }

    /**
     * @brief Return the flux vector for the LDG diffusion problem.
     * \todo Complete the viscous flux vector
     */
    void NavierStokesCFE::v_GetViscousFluxVectorDeAlias(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
              TensorOfArray3D<NekDouble>                &derivativesO1,
              TensorOfArray3D<NekDouble>                &viscousTensor)
    {
        // Factor to rescale 1d points in dealiasing.
        NekDouble OneDptscale = 2;
        // Get number of points to dealias a cubic non-linearity
        size_t nScalar   = physfield.size();
        int nPts      = m_fields[0]->Get1DScaledTotPoints(OneDptscale);
        size_t nPts_orig = physfield[0].size();

        // Auxiliary variables
        Array<OneD, NekDouble > divVel             (nPts, 0.0);

        // Stokes hypothesis
        const NekDouble lambda = -2.0/3.0;

        // Update viscosity and thermal conductivity
        GetViscosityAndThermalCondFromTemp(physfield[nScalar-1], m_mu,
            m_thermalConductivity);

        // Interpolate inputs and initialise interpolated output
        Array<OneD, Array<OneD, NekDouble> > vel_interp(m_spacedim);
        TensorOfArray3D<NekDouble>           deriv_interp(m_spacedim);
        TensorOfArray3D<NekDouble>           out_interp(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            // Interpolate velocity
            vel_interp[i]   = Array<OneD, NekDouble> (nPts);
            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, physfield[i], vel_interp[i]);

            // Interpolate derivatives
            deriv_interp[i] = Array<OneD, Array<OneD, NekDouble> >
                             (m_spacedim+1);
            for (int j = 0; j < m_spacedim+1; ++j)
            {
                deriv_interp[i][j] = Array<OneD, NekDouble> (nPts);
                m_fields[0]->PhysInterp1DScaled(
                    OneDptscale, derivativesO1[i][j], deriv_interp[i][j]);
            }

            // Output (start from j=1 since flux is zero for rho)
            out_interp[i] = Array<OneD, Array<OneD, NekDouble> > (m_spacedim+2);
            for (int j = 1; j < m_spacedim+2; ++j)
            {
                out_interp[i][j] = Array<OneD, NekDouble> (nPts);
            }
        }

        // Velocity divergence
        for (int j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vadd(nPts, divVel, 1, deriv_interp[j][j], 1,
                        divVel, 1);
        }

        // Velocity divergence scaled by lambda * mu
        Vmath::Smul(nPts, lambda, divVel, 1, divVel, 1);
        Vmath::Vmul(nPts, m_mu,  1, divVel, 1, divVel, 1);

        // Viscous flux vector for the rho equation = 0 (no need to dealias)
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts_orig, viscousTensor[i][0], 1);
        }

        // Viscous stress tensor (for the momentum equations)
        for (int i = 0; i < m_spacedim; ++i)
        {
            for (int j = i; j < m_spacedim; ++j)
            {
                Vmath::Vadd(nPts, deriv_interp[i][j], 1,
                                  deriv_interp[j][i], 1,
                                  out_interp[i][j+1], 1);

                Vmath::Vmul(nPts, m_mu, 1,
                                  out_interp[i][j+1], 1,
                                  out_interp[i][j+1], 1);

                if (i == j)
                {
                    // Add divergence term to diagonal
                    Vmath::Vadd(nPts, out_interp[i][j+1], 1,
                                  divVel, 1,
                                  out_interp[i][j+1], 1);
                }
                else
                {
                    // Make symmetric
                    out_interp[j][i+1] = out_interp[i][j+1];
                }
            }
        }

        // Terms for the energy equation
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Zero(nPts, out_interp[i][m_spacedim+1], 1);
            // u_j * tau_ij
            for (int j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vvtvp(nPts, vel_interp[j], 1,
                               out_interp[i][j+1], 1,
                               out_interp[i][m_spacedim+1], 1,
                               out_interp[i][m_spacedim+1], 1);
            }
            // Add k*T_i
            Vmath::Vvtvp(nPts, m_thermalConductivity, 1,
                               deriv_interp[i][m_spacedim], 1,
                               out_interp[i][m_spacedim+1], 1,
                               out_interp[i][m_spacedim+1], 1);
        }

        // Project to original space
        for (int i = 0; i < m_spacedim; ++i)
        {
            for (int j = 1; j < m_spacedim+2; ++j)
            {
                m_fields[0]->PhysGalerkinProjection1DScaled(
                    OneDptscale,
                    out_interp[i][j],
                    viscousTensor[i][j]);
            }
        }
    }

<<<<<<< HEAD
      /**
     * @brief Return the flux vector for the IP diffusion problem, based on
     * conservative variables
     */
    void NavierStokesCFE::GetViscousFluxVectorConservVar(
        const int                                              nDim,
        const Array<OneD, Array<OneD, NekDouble> >             &inarray,
        const TensorOfArray3D<NekDouble>                       &qfields,
        TensorOfArray3D<NekDouble>                             &outarray,
        Array< OneD, int >                                     &nonZeroIndex,
        const Array<OneD, Array<OneD, NekDouble> >             &normal)
    {
        size_t nConvectiveFields   = inarray.size();
        size_t nPts=inarray[0].size();
        int n_nonZero   =   nConvectiveFields - 1;
        TensorOfArray3D<NekDouble> fluxVec;
        Array<OneD, Array<OneD, NekDouble>> outtmp{nConvectiveFields};

        for (int i = 0; i < nConvectiveFields; ++i)
        {
            outtmp[i]=Array<OneD, NekDouble>{nPts, 0.0};
        }

        for (int i = 0; i < outarray.size(); ++i)
        {
            for (int j = 0; j < nConvectiveFields; ++j)
            {
                Vmath::Zero(nPts, outarray[i][j], 1);
            }
        }

        if (normal.size())
        {
            for (int nd = 0; nd < nDim; ++nd)
            {
                for (int nderiv = 0; nderiv < nDim; ++nderiv)
                {
                    GetViscousFluxBilinearForm(nDim, nd, nderiv, inarray,
                                                qfields[nderiv], outtmp);

                    for (int j = 0; j < nConvectiveFields; ++j)
                    {
                        Vmath::Vvtvp(nPts, &normal[nd][0], 1,
                                    &outtmp[j][0], 1,
                                    &outarray[0][j][0], 1,
                                    &outarray[0][j][0], 1);
                    }
                }
            }
        }
        else
        {
            fluxVec = outarray;

            for (int nd = 0; nd < nDim; ++nd)
            {
                for (int nderiv = 0; nderiv < nDim; ++nderiv)
                {
                    GetViscousFluxBilinearForm(nDim, nd, nderiv, inarray,
                                                qfields[nderiv], outtmp);

                    for (int j = 0; j < nConvectiveFields; ++j)
                    {
                        Vmath::Vadd(nPts, &outtmp[j][0], 1,
                                    &fluxVec[nd][j][0], 1,
                                    &fluxVec[nd][j][0], 1);
                    }
                }
            }
        }

        nonZeroIndex = Array< OneD, int > {size_t(n_nonZero), 0};

        for (int i = 1; i < n_nonZero + 1; ++i)
        {
            nonZeroIndex[n_nonZero - i] =   nConvectiveFields - i;
        }
    }

    /**
     * @brief For very special treatment. For general boundaries it does nothing
     * But for WallViscous and WallAdiabatic, special treatment is needed
     * because they get the same Bwd value, special treatment is needed for
     * boundary treatment of diffusion flux
     * Note: This special treatment could be removed by seperating
     * WallViscous and WallAdiabatic into two different classes.
     *
     */
    void NavierStokesCFE::SpecialBndTreat(
        Array<OneD, Array<OneD, NekDouble>> &consvar)
    {
        size_t nConvectiveFields = consvar.size();
        int ndens       = 0;
        int nengy       = nConvectiveFields - 1;

        Array<OneD, Array<OneD, NekDouble>> bndCons {nConvectiveFields};
        Array<OneD, NekDouble> bndTotEngy;
        Array<OneD, NekDouble> bndPressure;
        Array<OneD, NekDouble> bndRho;
        Array<OneD, NekDouble> bndIntEndy;
        int nLengthArray = 0;

        int cnt = 0;
        int nBndRegions = m_fields[nengy]->
            GetBndCondExpansions().size();
        for (int j = 0; j < nBndRegions; ++j)
        {
            if (m_fields[nengy]->GetBndConditions()[j]->
                GetBoundaryConditionType() ==
                SpatialDomains::ePeriodic)
            {
                continue;
            }

            size_t nBndEdges = m_fields[nengy]->
            GetBndCondExpansions()[j]->GetExpSize();
            for (int e = 0; e < nBndEdges; ++e)
            {
                size_t nBndEdgePts = m_fields[nengy]->
                GetBndCondExpansions()[j]->GetExp(e)->GetTotPoints();

                int id2 = m_fields[0]->GetTrace()->
                GetPhys_Offset(m_fields[0]->GetTraceMap()->
                            GetBndCondIDToGlobalTraceID(cnt++));

                // Imposing Temperature Twall at the wall
                if (boost::iequals(m_fields[nengy]->GetBndConditions()[j]->
                    GetUserDefined(), "WallViscous"))
                {
                    if (nBndEdgePts != nLengthArray)
                    {
                        for (int i = 0; i < nConvectiveFields; ++i)
                        {
                            bndCons[i] = Array<OneD, NekDouble>
                                        {nBndEdgePts, 0.0};
                        }
                        bndTotEngy  = Array<OneD, NekDouble> {nBndEdgePts, 0.0};
                        bndPressure = Array<OneD, NekDouble> {nBndEdgePts, 0.0};
                        bndRho      = Array<OneD, NekDouble> {nBndEdgePts, 0.0};
                        bndIntEndy  = Array<OneD, NekDouble> {nBndEdgePts, 0.0};
                        nLengthArray = nBndEdgePts;
                    }
                    else
                    {
                        Vmath::Zero(nLengthArray, bndPressure, 1);
                        Vmath::Zero(nLengthArray, bndRho     , 1);
                        Vmath::Zero(nLengthArray, bndIntEndy , 1);
                    }

                    Array<OneD, NekDouble> tmp;

                    for (int k = 0; k < nConvectiveFields; ++k)
                    {
                        Vmath::Vcopy(nBndEdgePts, tmp = consvar[k] + id2, 1,
                                    bndCons[k], 1);
                    }

                    m_varConv->GetPressure(bndCons, bndPressure);
                    Vmath::Fill(nLengthArray, m_Twall, bndTotEngy, 1);
                    m_varConv->GetRhoFromPT(bndPressure, bndTotEngy, bndRho);
                    m_varConv->GetEFromRhoP(bndRho, bndPressure, bndIntEndy);
                    m_varConv->GetDynamicEnergy(bndCons, bndTotEngy);

                    Vmath::Vvtvp(nBndEdgePts, bndIntEndy, 1, bndCons[ndens], 1,
                        bndTotEngy, 1, bndTotEngy, 1);

                    Vmath::Vcopy(nBndEdgePts,
                                bndTotEngy, 1,
                                tmp = consvar[nengy] + id2, 1);
                }
            }
        }
    }

    /**
     *
     * @brief Calculate and return the Symmetric flux in IP method.
     */
    void NavierStokesCFE::GetViscousSymmtrFluxConservVar(
        const int                                           nDim,
        const Array<OneD, Array<OneD, NekDouble> >          &inaverg,
        const Array<OneD, Array<OneD, NekDouble > >         &inarray,
        TensorOfArray3D<NekDouble>                          &outarray,
        Array< OneD, int >                                  &nonZeroIndex,
        const Array<OneD, Array<OneD, NekDouble> >          &normal)
    {
        size_t nConvectiveFields   = inarray.size();
        size_t nPts                = inaverg[nConvectiveFields - 1].size();
        nonZeroIndex = Array<OneD, int>{nConvectiveFields - 1, 0};
        for (int i = 0; i < nConvectiveFields - 1; ++i)
        {
            nonZeroIndex[i] =   i + 1;
        }

        std::vector<NekDouble> inAvgTmp(nConvectiveFields);
        std::vector<NekDouble> inTmp(nConvectiveFields);
        std::vector<NekDouble> outTmp(nConvectiveFields);
        for (int d = 0; d < nDim; ++d)
        {
            for (int nderiv = 0; nderiv < nDim; ++nderiv)
            {
                for (size_t p = 0; p < nPts; ++p)
                {
                    // rearrenge data
                    for (int f = 0; f < nConvectiveFields; ++f)
                    {
                        inAvgTmp[f] = inaverg[f][p];
                        inTmp[f] = inarray[f][p];
                    }

                    // get temp
                    NekDouble temperature = m_varConv->GetTemperature(inTmp.data());
                    // get viscosity
                    NekDouble mu;
                    GetViscosityFromTempKernel(temperature, mu);

                    GetViscousFluxBilinearFormKernel(nDim, d, nderiv,
                        inAvgTmp.data(), inTmp.data(), mu, outTmp.data());

                    for (int f = 0; f < nConvectiveFields; ++f)
                    {
                        outarray[d][f][p] += normal[d][p] * outTmp[f];
                    }
                }
            }
        }
    }


    void NavierStokesCFE::CalcViscosity(
        const Array<OneD, const Array<OneD, NekDouble>> &inaverg,
              Array<OneD, NekDouble>                    &mu)
    {
        int nConvectiveFields = inaverg.size();
        int nPts = inaverg[nConvectiveFields-1].size();

        if (m_ViscosityType == "Variable")
        {
            Array<OneD, NekDouble> tmp(nPts,0.0);
            m_varConv->GetTemperature(inaverg,tmp);
            m_varConv->GetDynamicViscosity(tmp, mu);
        }
        else
        {
            //mu may be on volume or trace
            Vmath::Fill(nPts, m_mu[0], mu, 1);
        }
    }

    /**
     * @brief Return the penalty vector for the LDGNS diffusion problem.
     */
    void NavierStokesCFE::v_GetFluxPenalty(
        const Array<OneD, const Array<OneD, NekDouble>> &uFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &uBwd,
              Array<OneD,       Array<OneD, NekDouble>> &penaltyCoeff)
    {
        unsigned int nTracePts  = uFwd[0].size();

        // Compute average temperature
        unsigned int nVariables = uFwd.size();
        Array<OneD, NekDouble> tAve{nTracePts, 0.0};
        Vmath::Svtsvtp(nTracePts, 0.5, uFwd[nVariables-1], 1,
            0.5, uBwd[nVariables-1], 1, tAve, 1);

        // Get average viscosity and thermal conductivity
        Array<OneD, NekDouble> muAve{nTracePts, 0.0};
        Array<OneD, NekDouble> tcAve{nTracePts, 0.0};

        GetViscosityAndThermalCondFromTemp(tAve, muAve, tcAve);

        // Compute penalty term
        for (int i = 0; i < nVariables; ++i)
        {
            // Get jump of u variables
            Vmath::Vsub(nTracePts, uFwd[i], 1, uBwd[i], 1, penaltyCoeff[i], 1);
            // Multiply by variable coefficient = {coeff} ( u^+ - u^- )
            if ( i < nVariables-1 )
            {
                Vmath::Vmul(nTracePts, muAve, 1, penaltyCoeff[i], 1,
                    penaltyCoeff[i], 1);
            }
            else
            {
                Vmath::Vmul(nTracePts, tcAve, 1, penaltyCoeff[i], 1,
                    penaltyCoeff[i], 1);
            }
        }
    }


    /**
     * @brief Update viscosity
     * todo: add artificial viscosity here
     */
    void NavierStokesCFE::GetViscosityAndThermalCondFromTemp(
        const Array<OneD, NekDouble> &temperature,
              Array<OneD, NekDouble> &mu,
              Array<OneD, NekDouble> &thermalCond)
    {
        auto nPts = temperature.size();

        for (size_t p = 0; p < nPts; ++p)
        {
            GetViscosityAndThermalCondFromTempKernel(temperature[p], mu[p],
                thermalCond[p]);
        }

        // Add artificial viscosity if wanted
        if (m_shockCaptureType == "Physical")
        {
            auto nTracePts = m_fields[0]->GetTrace()->GetTotPoints();
            if (nPts != nTracePts)
            {
                Vmath::Vadd(nPts, mu, 1, m_varConv->GetAv(), 1, mu, 1);
            }
            else
            {
                Vmath::Vadd(nPts, mu, 1, m_varConv->GetAvTrace(), 1, mu, 1);
            }
        }

        // Thermal conductivity
        NekDouble tRa = m_Cp / m_Prandtl;
        Vmath::Smul(nPts, tRa, mu, 1, thermalCond, 1);
    }

    /**
    * @brief Get divergence and curl squared
    *
    * @param input
    *   fields -> expansion list pointer
    *   cnsVar -> conservative variables
    *   cnsVarFwd -> forward trace of conservative variables
    *   cnsVarBwd -> backward trace of conservative variables
    * @paran output
    *   divSquare -> divergence
    *   curlSquare -> curl squared
    *
    */
    void NavierStokesCFE::GetDivCurlSquared(
        const Array<OneD, MultiRegions::ExpListSharedPtr>& fields,
        const Array<OneD, Array<OneD, NekDouble>>& cnsVar,
        Array<OneD, NekDouble>& div,
        Array<OneD, NekDouble>& curlSquare,
        const Array<OneD, Array<OneD, NekDouble>>& cnsVarFwd,
        const Array<OneD, Array<OneD, NekDouble>>& cnsVarBwd)
    {
        auto nDim = fields[0]->GetCoordim(0);
        auto nVar = cnsVar.size();
        auto nPts = cnsVar[0].size();
        auto nPtsTrc = cnsVarFwd[0].size();

        // These should be allocated once
        Array<OneD, Array<OneD, NekDouble>>  primVar(nVar-1),
            primVarFwd(nVar-1), primVarBwd(nVar-1);

        for (unsigned short d = 0; d < nVar-2; ++d)
        {
            primVar[d] = Array<OneD, NekDouble>(nPts, 0.0);
            primVarFwd[d] = Array<OneD, NekDouble>(nPtsTrc, 0.0);
            primVarBwd[d] = Array<OneD, NekDouble>(nPtsTrc, 0.0);
        }
        auto ergLoc = nVar-2;
        primVar[ergLoc] = Array<OneD, NekDouble>(nPts, 0.0);
        primVarFwd[ergLoc] = Array<OneD, NekDouble>(nPtsTrc, 0.0);
        primVarBwd[ergLoc] = Array<OneD, NekDouble>(nPtsTrc, 0.0);

        // Get primitive variables [u,v,w,T=0]
        // Possibly should be changed to [rho,u,v,w,T] to make IP and LDGNS more
        // consistent with each other
        for (unsigned short d = 0; d < nVar-2; ++d)
        {
            // Volume
            for (size_t p = 0; p < nPts; ++p)
            {
                primVar[d][p] = cnsVar[d+1][p] / cnsVar[0][p];
            }
            // Trace
            for (size_t p = 0; p < nPtsTrc; ++p)
            {
                primVarFwd[d][p] = cnsVarFwd[d+1][p] / cnsVarFwd[0][p];
                primVarBwd[d][p] = cnsVarBwd[d+1][p] / cnsVarBwd[0][p];
            }
        }

        // this should be allocated once
        Array<OneD,Array<OneD, Array<OneD, NekDouble>>> primVarDer(nDim);
        for (unsigned short j = 0; j < nDim; ++j)
        {
            primVarDer[j] = Array<OneD, Array<OneD, NekDouble>> (nVar-1);
            for (unsigned short i = 0; i < nVar-1; ++i)
            {
                primVarDer[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
            }
        }

        // Get derivative tensor
        m_diffusion->DiffuseCalculateDerivative(fields, primVar, primVarDer,
            primVarFwd, primVarBwd);

        // Get div curl squared
        GetDivCurlFromDvelT(primVarDer, div, curlSquare);
    }


    /**
     * @brief Get divergence and curl from velocity derivative tensor
     *
     */
    void NavierStokesCFE::GetDivCurlFromDvelT(
        const TensorOfArray3D<NekDouble>& pVarDer,
              Array<OneD, NekDouble>&     div,
              Array<OneD, NekDouble>&     curlSquare)
    {
        auto nDim = pVarDer.size();
        auto nPts = div.size();

        // div velocity
        for (size_t p = 0; p < nPts; ++p)
        {
            NekDouble divTmp = 0;
            for (unsigned short j = 0; j < nDim; ++j)
            {
                divTmp += pVarDer[j][j][p];
            }
            div[p] = divTmp;
        }

        // |curl velocity| ** 2
        if (nDim > 2)
        {
            for (size_t p = 0; p < nPts; ++p)
            {
                // curl[0] 3/2 - 2/3
                NekDouble curl032 = pVarDer[2][1][p]; // load 1x
                NekDouble curl023 = pVarDer[1][2][p]; // load 1x
                NekDouble curl0 = curl032 - curl023;
                // square curl[0]
                NekDouble curl0sqr = curl0 * curl0;

                // curl[1] 3/1 - 1/3
                NekDouble curl131 = pVarDer[2][0][p]; // load 1x
                NekDouble curl113 = pVarDer[0][2][p]; // load 1x
                NekDouble curl1 = curl131 - curl113;
                // square curl[1]
                NekDouble curl1sqr = curl1 * curl1;

                // curl[2] 1/2 - 2/1
                NekDouble curl212 = pVarDer[0][1][p]; // load 1x
                NekDouble curl221 = pVarDer[1][0][p]; // load 1x
                NekDouble curl2 = curl212 - curl221;
                // square curl[2]
                NekDouble curl2sqr = curl2 * curl2;

                // reduce
                curl0sqr += curl1sqr + curl2sqr;
                // store
                curlSquare[p] = curl0sqr; // store 1x

            }
        }
        else if (nDim > 1)
        {
            for (size_t p = 0; p < nPts; ++p)
            {
                // curl[2] 1/2
                NekDouble c212 = pVarDer[0][1][p]; // load 1x
                // curl[2] 2/1
                NekDouble c221 = pVarDer[1][0][p]; // load 1x
                // curl[2] 1/2 - 2/1
                NekDouble curl = c212 - c221;
                // square curl[2]
                curlSquare[p] = curl * curl; // store 1x
            }
        }
        else
        {
            Vmath::Fill(nPts, 0.0, curlSquare, 1);
        }
    }

    void NavierStokesCFE::v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        bool extraFields;
        m_session->MatchSolverInfo("OutputExtraFields","True",
                                   extraFields, true);
        if (extraFields)
        {
            const int nPhys   = m_fields[0]->GetNpoints();
            const int nCoeffs = m_fields[0]->GetNcoeffs();
            Array<OneD, Array<OneD, NekDouble> > cnsVar(m_fields.size());

            for (int i = 0; i < m_fields.size(); ++i)
            {
                cnsVar[i] = m_fields[i]->GetPhys();
            }

            Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);
            Array<OneD, Array<OneD, NekDouble> > velFwd  (m_spacedim);
            for (int i = 0; i < m_spacedim; ++i)
            {
                velocity[i] = Array<OneD, NekDouble> (nPhys);
                velFwd[i]   = Array<OneD, NekDouble> (nCoeffs);
            }

            Array<OneD, NekDouble> pressure(nPhys), temperature(nPhys);
            Array<OneD, NekDouble> entropy(nPhys);
            Array<OneD, NekDouble> soundspeed(nPhys), mach(nPhys);
            Array<OneD, NekDouble> sensor(nPhys), SensorKappa(nPhys);

            m_varConv->GetVelocityVector(cnsVar, velocity);
            m_varConv->GetPressure  (cnsVar, pressure);
            m_varConv->GetTemperature(cnsVar, temperature);
            m_varConv->GetEntropy   (cnsVar, entropy);
            m_varConv->GetSoundSpeed(cnsVar, soundspeed);
            m_varConv->GetMach      (cnsVar, soundspeed, mach);

            int sensorOffset;
            m_session->LoadParameter ("SensorOffset", sensorOffset, 1);
            m_varConv->GetSensor (m_fields[0], cnsVar, sensor, SensorKappa,
                                    sensorOffset);

            Array<OneD, NekDouble> pFwd(nCoeffs), TFwd(nCoeffs);
            Array<OneD, NekDouble> sFwd(nCoeffs);
            Array<OneD, NekDouble> aFwd(nCoeffs), mFwd(nCoeffs);
            Array<OneD, NekDouble> sensFwd(nCoeffs);

            string velNames[3] = {"u", "v", "w"};
            for (int i = 0; i < m_spacedim; ++i)
            {
                m_fields[0]->FwdTrans_IterPerExp(velocity[i], velFwd[i]);
                variables.push_back(velNames[i]);
                fieldcoeffs.push_back(velFwd[i]);
            }

            m_fields[0]->FwdTrans_IterPerExp(pressure,   pFwd);
            m_fields[0]->FwdTrans_IterPerExp(temperature,TFwd);
            m_fields[0]->FwdTrans_IterPerExp(entropy,    sFwd);
            m_fields[0]->FwdTrans_IterPerExp(soundspeed, aFwd);
            m_fields[0]->FwdTrans_IterPerExp(mach,       mFwd);
            m_fields[0]->FwdTrans_IterPerExp(sensor,     sensFwd);

            variables.push_back  ("p");
            variables.push_back  ("T");
            variables.push_back  ("s");
            variables.push_back  ("a");
            variables.push_back  ("Mach");
            variables.push_back  ("Sensor");
            fieldcoeffs.push_back(pFwd);
            fieldcoeffs.push_back(TFwd);
            fieldcoeffs.push_back(sFwd);
            fieldcoeffs.push_back(aFwd);
            fieldcoeffs.push_back(mFwd);
            fieldcoeffs.push_back(sensFwd);

            if (m_artificialDiffusion)
            {
                // reuse pressure
                Array<OneD, NekDouble> sensorFwd(nCoeffs);
                m_artificialDiffusion->GetArtificialViscosity(cnsVar, pressure);
                m_fields[0]->FwdTrans_IterPerExp(pressure,   sensorFwd);

                variables.push_back  ("ArtificialVisc");
                fieldcoeffs.push_back(sensorFwd);
            }

            if (m_shockCaptureType == "Physical")
            {

                Array<OneD, Array<OneD, NekDouble>> cnsVarFwd(m_fields.size()),
                    cnsVarBwd(m_fields.size());

                for (int i = 0; i < m_fields.size(); ++i)
                {
                    cnsVarFwd[i] = Array<OneD, NekDouble>(GetTraceTotPoints());
                    cnsVarBwd[i] = Array<OneD, NekDouble>(GetTraceTotPoints());
                    m_fields[i]->GetFwdBwdTracePhys(cnsVar[i], cnsVarFwd[i], cnsVarBwd[i]);
                }

                Array<OneD, NekDouble> div(nPhys), curlSquare(nPhys);
                GetDivCurlSquared(m_fields, cnsVar, div, curlSquare,
                    cnsVarFwd, cnsVarBwd);

                Array<OneD, NekDouble> divFwd(nCoeffs, 0.0);
                m_fields[0]->FwdTrans_IterPerExp(div, divFwd);
                variables.push_back("div");
                fieldcoeffs.push_back(divFwd);

                Array<OneD, NekDouble> curlFwd(nCoeffs, 0.0);
                m_fields[0]->FwdTrans_IterPerExp(curlSquare, curlFwd);
                variables.push_back("curl^2");
                fieldcoeffs.push_back(curlFwd);

                m_varConv->SetAv(m_fields, cnsVar, div, curlSquare);

                Array<OneD, NekDouble> muavFwd(nCoeffs);
                m_fields[0]->FwdTrans_IterPerExp(m_varConv->GetAv(), muavFwd);
                variables.push_back  ("ArtificialVisc");
                fieldcoeffs.push_back(muavFwd);
            }
        }
    }
}
