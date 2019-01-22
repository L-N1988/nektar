///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionIP.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: IP diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/DiffusionIP.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string DiffusionIP::type = GetDiffusionFactory().
            RegisterCreatorFunction("InteriorPenalty", DiffusionIP::create);

        DiffusionIP::DiffusionIP()
        {
        }

        void DiffusionIP::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            m_session = pSession;

            m_session->LoadSolverInfo("ShockCaptureType",
                                  m_shockCaptureType,    "Off");	

            m_session->LoadParameter("IPSymmFtluxCoeff",
                                  m_IPSymmFtluxCoeff,   0.0);	//0.5

            m_session->LoadParameter("IP2ndDervCoeff",
                                  m_IP2ndDervCoeff,   0.0); // 1.0/12.0	

            // Setting up the normals
            int i;
            int nDim = pFields[0]->GetCoordim(0);
            int nVariable = pFields.num_elements();
            int nTracePts = pFields[0]->GetTrace()->GetTotPoints();
            
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts,0.0);
            }
            m_traceAver = Array<OneD, Array<OneD, NekDouble> >(nVariable);
            m_traceJump = Array<OneD, Array<OneD, NekDouble> >(nVariable);
            for(i = 0; i < nVariable; ++i)
            {
                m_traceAver[i] = Array<OneD, NekDouble> (nTracePts,0.0);
                m_traceJump[i] = Array<OneD, NekDouble> (nTracePts,0.0);
            }

            pFields[0]->GetTrace()->GetNormals(m_traceNormals);
            m_traceNormDirctnElmtLength =   Array<OneD, NekDouble> (nTracePts,0.0);
            pFields[0]->GetTrace()->GetElmtNormalLength(m_traceNormDirctnElmtLength);

            // TODO:: to check parallel case
            Array<OneD, NekDouble> lengthstmp(nTracePts,0.0);
            pFields[0]->PeriodicBwdCopy(m_traceNormDirctnElmtLength,lengthstmp);
            Vmath::Vadd(nTracePts,lengthstmp,1,m_traceNormDirctnElmtLength,1,m_traceNormDirctnElmtLength,1);

            // if(abs(m_IP2ndDervCoeff)>1.0E-14)
            // {
                m_traceNormDirctnElmtLengthRecip =   Array<OneD, NekDouble> (nTracePts,0.0);
                Vmath::Sdiv(nTracePts,1.0,m_traceNormDirctnElmtLength,1,m_traceNormDirctnElmtLengthRecip,1);
            // }

            m_tracBwdWeight  =   Array<OneD, NekDouble> (nTracePts,0.0);
            pFields[0]->GetBwdWeight(m_tracBwdWeight);
            Array<OneD, NekDouble> tmpBwdWeight(nTracePts,0.0);
            for(int i =1; i<nVariable;i++)
            {
                pFields[i]->GetBwdWeight(tmpBwdWeight);
                Vmath::Vsub(nTracePts,tmpBwdWeight,1,m_tracBwdWeight,1,tmpBwdWeight,1);
                Vmath::Vabs(nTracePts,tmpBwdWeight,1,tmpBwdWeight,1);
                NekDouble norm = 0.0;
                for(int j = 0; j<nTracePts; j++)
                {
                    norm += tmpBwdWeight[j];
                }
                ASSERTL0(norm<1.0E-11,"different BWD for different variable not coded yet");
            }

            m_MuVarTrace   =   NullNekDouble1DArray;
            if (m_ArtificialDiffusionVector)
            {
                m_MuVarTrace  =   Array<OneD, NekDouble>(nTracePts, 0.0);
            }
        }
        
        void DiffusionIP::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {

            int nCoeffs   = fields[0]->GetNcoeffs();
            Array<OneD, Array<OneD, NekDouble> > tmp(nConvectiveFields);
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
            }
            DiffusionIP::v_Diffuse_coeff(nConvectiveFields,fields,inarray,tmp,pFwd,pBwd);
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                fields[i]->BwdTrans             (tmp[i], outarray[i]);
            }
        }

        void DiffusionIP::v_Diffuse_coeff(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            int i, j;
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();
            int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();


            Array<OneD, NekDouble>    Fwd(nTracePts,0.0);
            Array<OneD, NekDouble>    Bwd(nTracePts,0.0);

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > elmtFlux(nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > qfield(nDim);
            for (j = 0; j < nDim; ++j)
            {
                qfield[j]       = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                elmtFlux[j]     = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    qfield[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
                }
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    elmtFlux[j][i]   = Array<OneD, NekDouble>(nPts, 0.0);
                }
            }

            Array<OneD, Array<OneD, NekDouble> >    vFwd(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >    vBwd(nConvectiveFields);
            if (pFwd == NullNekDoubleArrayofArray ||
                pBwd == NullNekDoubleArrayofArray)
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    vFwd[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
                    vBwd[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
                }
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    fields[i]->GetFwdBwdTracePhys(inarray[i], vFwd[i], vBwd[i]);
                }
            }
            else
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    vFwd[i]    =   pFwd[i];
                    vBwd[i]    =   pBwd[i];
                }
            }

            DiffuseCalculateDerivative(nConvectiveFields,fields,inarray,qfield,vFwd,vBwd);

            // m_FunctorDerivBndCond(inarray,qfield,m_time,vFwd,tmparray3D);

            Array<OneD, int > nonZeroIndex;
            DiffuseVolumeFlux(nConvectiveFields,fields,inarray,qfield,elmtFlux,nonZeroIndex);

            Array<OneD, Array<OneD, NekDouble> > tmpFluxIprdct(nDim);
            // volume intergration: the nonZeroIndex indicates which flux is nonzero
            for(i = 0; i < nonZeroIndex.num_elements(); ++i)
            {
                int j = nonZeroIndex[i];
                for (int k = 0; k < nDim; ++k)
                {
                    tmpFluxIprdct[k] = elmtFlux[k][j];
                }
                fields[j]->IProductWRTDerivBase(tmpFluxIprdct,outarray[j]);
                Vmath::Neg                      (nCoeffs, outarray[j], 1);
            }
            // release qfield, elmtFlux and muvar;
            for (j = 0; j < nDim; ++j)
            {
                elmtFlux[j]     = NullNekDoubleArrayofArray;
            }

            Array<OneD, Array<OneD, NekDouble > > Traceflux(nConvectiveFields);
            for (int j = 0; j < nConvectiveFields; ++j)
            {
                Traceflux[j]   = Array<OneD, NekDouble>(nTracePts, 0.0);
            }

            DiffuseTraceFlux(nConvectiveFields,fields,inarray,qfield,elmtFlux,Traceflux,vFwd,vBwd,nonZeroIndex);

            for(i = 0; i < nonZeroIndex.num_elements(); ++i)
            {
                int j = nonZeroIndex[i];

                fields[j]->AddTraceIntegral     (Traceflux[j], outarray[j]);
                fields[j]->SetPhysState         (false);
                fields[j]->MultiplyByElmtInvMass(outarray[j], outarray[j]);
            }

            AddDiffusionSymmFluxToCoeff(nConvectiveFields, fields, inarray,qfield,elmtFlux, outarray, vFwd, vBwd);
        }

        void DiffusionIP::v_DiffuseCalculateDerivative(
            const int                                                   nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>           &fields,
            const Array<OneD, Array<OneD, NekDouble> >                  &inarray,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    &qfield,
            const Array<OneD, Array<OneD, NekDouble> >                  &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >                  &pBwd)
        {
            int nDim      = fields[0]->GetCoordim(0);

            Array<OneD, Array<OneD, NekDouble> > qtmp(3);
            for(int nd=0; nd<3; nd++)
            {
                qtmp[nd]    =   NullNekDouble1DArray;
            }
            for(int i = 0; i < nConvectiveFields; ++i)
            {
                for(int nd=0; nd<nDim; nd++)
                {
                    qtmp[nd]    =   qfield[nd][i];
                }
                fields[i]->PhysDeriv(inarray[i], qtmp[0], qtmp[1], qtmp[2]);
            }
        }

        void DiffusionIP::v_DiffuseVolumeFlux(
            const int                                           nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble>>           &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            Array< OneD, int >                                  &nonZeroIndex) 
        {
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();

            Array<OneD, NekDouble> muvar        =   NullNekDouble1DArray;
            if (m_ArtificialDiffusionVector)
            {
                muvar       =   Array<OneD, NekDouble>(nPts, 0.0);
                GetAVmu(fields,inarray,muvar,m_MuVarTrace);
            }

            Array<OneD, Array<OneD, NekDouble> > tmparray2D = NullNekDoubleArrayofArray;

            // TODO: qfield AND elmtFlux share storage????
            m_FunctorDiffusionfluxCons(nConvectiveFields,nDim,inarray,qfield, VolumeFlux,nonZeroIndex,tmparray2D,muvar);
        }
            
        void DiffusionIP::v_DiffuseTraceFlux(
            const int                                           nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble>>           &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            Array<OneD, Array<OneD, NekDouble> >                &TraceFlux,
            const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
            const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
            Array< OneD, int >                                  &nonZeroIndex)
        {
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();
            // int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, Array<OneD, Array<OneD, NekDouble > > > traceflux3D(1);
            traceflux3D[0]  =   TraceFlux;

            CalTraceNumFlux_ReduceComm(
                nConvectiveFields, nDim, nPts, nTracePts, m_IP2ndDervCoeff,
                fields, inarray, qfield, pFwd, pBwd, m_MuVarTrace,
                nonZeroIndex, traceflux3D, m_traceAver, m_traceJump);
        }

        void DiffusionIP::v_DiffuseTraceFlux(
            const int                                                       nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>               &fields,
            const Array<OneD, Array<OneD, NekDouble>>                       &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >  &qfield,
            Array<OneD, Array<OneD, NekDouble> >                            &TraceFlux,
            const Array<OneD, Array<OneD, NekDouble>>                       &pFwd,
            const Array<OneD, Array<OneD, NekDouble>>                       &pBwd,
            const Array<OneD, NekDouble>                                    &MuAVTrace,
            Array< OneD, int >                                              &nonZeroIndex  ,
            const Array<OneD, Array<OneD, NekDouble>>                       &Aver          ,
            const Array<OneD, Array<OneD, NekDouble>>                       &Jump          )
        {
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();
            // int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

            Array<OneD, Array<OneD, Array<OneD, NekDouble > > > traceflux3D(1);
            traceflux3D[0]  =   TraceFlux;

            Array<OneD, Array<OneD, NekDouble> >           pAver;
            Array<OneD, Array<OneD, NekDouble> >           pJump;
            if((Aver.num_elements()&&Jump.num_elements()))
            {
                pAver = Aver;
                pJump = Jump;
            }
            else
            {
                pAver   =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                pJump   =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                for(int i = 0;i<nConvectiveFields;i++)
                {
                    pAver[i]    =   Array<OneD, NekDouble> (nTracePts,0.0);
                    pJump[i]    =   Array<OneD, NekDouble> (nTracePts,0.0);
                }
                // ConsVarAveJump(nConvectiveFields,nTracePts,pFwd,pBwd,pAver,pJump);
            }

            CalTraceNumFlux_ReduceComm(
                nConvectiveFields, nDim, nPts, nTracePts, m_IP2ndDervCoeff,
                fields, inarray, qfield, pFwd, pBwd, m_MuVarTrace,
                nonZeroIndex, traceflux3D, pAver, pJump);
        }

        void DiffusionIP::v_AddDiffusionSymmFluxToCoeff(
            const int                                           nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble> >          &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            Array<OneD, Array<OneD, NekDouble> >                &outarray,
            const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
        {
            if(abs(m_IPSymmFtluxCoeff)>1.0E-12)
            {
                int nDim      = fields[0]->GetCoordim(0);
                int nPts      = fields[0]->GetTotPoints();
                int nTracePts = fields[0]->GetTrace()->GetTotPoints();
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > traceSymflux(nDim);
                for (int nd = 0; nd < nDim; ++nd)
                {
                    traceSymflux[nd]    = Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                    for (int j = 0; j < nConvectiveFields; ++j)
                    {
                        traceSymflux[nd][j]   = Array<OneD, NekDouble>(nTracePts, 0.0);
                    }
                }
                Array< OneD, int >  nonZeroIndex;
                DiffuseTraceSymmFlux(nConvectiveFields,fields,inarray,qfield,VolumeFlux,traceSymflux,pFwd,pBwd,nonZeroIndex);

                AddSymmFluxIntegralToCoeff(nConvectiveFields,nDim,nPts,nTracePts,fields,nonZeroIndex,traceSymflux,outarray);
            }
        }

        void DiffusionIP::v_AddDiffusionSymmFluxToPhys(
            const int                                           nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble> >          &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            Array<OneD, Array<OneD, NekDouble> >                &outarray,
            const Array<OneD, Array<OneD, NekDouble> >          &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >          &pBwd)
        {
            
            if(abs(m_IPSymmFtluxCoeff)>1.0E-12)
            {
                int nDim      = fields[0]->GetCoordim(0);
                int nPts      = fields[0]->GetTotPoints();
                int nTracePts = fields[0]->GetTrace()->GetTotPoints();
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > traceSymflux(nDim);
                for (int nd = 0; nd < nDim; ++nd)
                {
                    traceSymflux[nd]    = Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                    for (int j = 0; j < nConvectiveFields; ++j)
                    {
                        traceSymflux[nd][j]   = Array<OneD, NekDouble>(nTracePts, 0.0);
                    }
                }
                Array< OneD, int >  nonZeroIndex;
                DiffuseTraceSymmFlux(nConvectiveFields,fields,inarray,qfield,VolumeFlux,traceSymflux,pFwd,pBwd,nonZeroIndex);

                AddSymmFluxIntegralToPhys(nConvectiveFields,nDim,nPts,nTracePts,fields,nonZeroIndex,traceSymflux,outarray);
            }
        }
        
        void DiffusionIP::DiffuseTraceSymmFlux(
            const int                                           nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
            const Array<OneD, Array<OneD, NekDouble>>           &inarray,
            Array<OneD,Array<OneD, Array<OneD, NekDouble> > >   &qfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &VolumeFlux,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &SymmFlux,
            const Array<OneD, Array<OneD, NekDouble>>           &pFwd,
            const Array<OneD, Array<OneD, NekDouble>>           &pBwd,
            Array< OneD, int >                                  &nonZeroIndex)
        {
            int nDim      = fields[0]->GetCoordim(0);
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();

            CalTraceSymFlux(nConvectiveFields,nDim,fields,m_traceAver,m_traceJump,
                        nonZeroIndex,SymmFlux);
        }

        void DiffusionIP::CalTraceSymFlux(
            const int                                                           nConvectiveFields,
            const int                                                           nDim,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, Array<OneD, NekDouble> >                          &solution_Aver,
                  Array<OneD, Array<OneD, NekDouble> >                          &solution_jump,
                  Array<OneD, int >                                             &nonZeroIndexsymm,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &traceSymflux)
        {
            int nTracePts = solution_jump[nConvectiveFields-1].num_elements();

            for (int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Smul(nTracePts,m_IPSymmFtluxCoeff,solution_jump[i],1,solution_jump[i],1);
            }
            
            m_FunctorSymmetricfluxCons(nConvectiveFields,nDim,solution_Aver,solution_jump,traceSymflux,nonZeroIndexsymm,m_traceNormals);

            for (int i = 0; i < nConvectiveFields; ++i)
            {
                MultiRegions::ExpListSharedPtr tracelist = fields[i]->GetTrace();
                for(int nd=0;nd<nDim;nd++)
                {
                    tracelist->MultiplyByQuadratureMetric(traceSymflux[nd][i],traceSymflux[nd][i]);
                }
            }
        }

        void DiffusionIP::AddSymmFluxIntegralToCoeff(
            const int                                                           nConvectiveFields,
            const int                                                           nDim,
            const int                                                           nPts,
            const int                                                           nTracePts,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, const int >                                       &nonZeroIndex,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &tracflux,
                  Array<OneD, Array<OneD, NekDouble> >                          &outarray)
        {
            int nCoeffs =   outarray[nConvectiveFields-1].num_elements();
            Array<OneD, NekDouble > tmpCoeff(nCoeffs,0.0);
            Array<OneD, Array<OneD, NekDouble> > tmpfield(nDim);
            for(int i = 0;i<nDim;i++)
            {
                tmpfield[i]    =   Array<OneD, NekDouble>(nPts,0.0);
            }
            int nv = 0;
            for(int j=0;j<nonZeroIndex.num_elements();j++)
            {
                nv  =   nonZeroIndex[j];
                for(int nd=0;nd<nDim;nd++)
                {
                    Vmath::Zero(nPts,tmpfield[nd],1);

                    fields[nv]->AddTraceQuadPhysToField(tracflux[nd][nv],tracflux[nd][nv],tmpfield[nd]);
                    fields[nv]->DividByQuadratureMetric(tmpfield[nd],tmpfield[nd]);
                }
                fields[nv]->IProductWRTDerivBase(tmpfield,tmpCoeff);
                Vmath::Vadd(nCoeffs,tmpCoeff,1,outarray[nv],1,outarray[nv],1);
            }
        }

        void DiffusionIP::AddSymmFluxIntegralToPhys(
            const int                                                           nConvectiveFields,
            const int                                                           nDim,
            const int                                                           nPts,
            const int                                                           nTracePts,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, const int >                                       &nonZeroIndex,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &tracflux,
                  Array<OneD, Array<OneD, NekDouble> >                          &outarray)
        {
            int nCoeffs =   outarray[nConvectiveFields-1].num_elements();
            Array<OneD, NekDouble > tmpCoeff(nCoeffs,0.0);
            Array<OneD, NekDouble > tmpPhysi(nPts,0.0);
            Array<OneD, Array<OneD, NekDouble> > tmpfield(nDim);
            for(int i = 0;i<nDim;i++)
            {
                tmpfield[i]    =   Array<OneD, NekDouble>(nPts,0.0);
            }
            int nv = 0;
            for(int j=0;j<nonZeroIndex.num_elements();j++)
            {
                nv  =   nonZeroIndex[j];
                for(int nd=0;nd<nDim;nd++)
                {
                    Vmath::Zero(nPts,tmpfield[nd],1);

                    fields[nv]->AddTraceQuadPhysToField(tracflux[nd][nv],tracflux[nd][nv],tmpfield[nd]);
                    fields[nv]->DividByQuadratureMetric(tmpfield[nd],tmpfield[nd]);
                }
                fields[nv]->IProductWRTDerivBase(tmpfield,tmpCoeff);
                fields[nv]->BwdTrans            (tmpCoeff,tmpPhysi);
                Vmath::Vadd(nPts,tmpPhysi,1,outarray[nv],1,outarray[nv],1);
            }
        }
        
        void DiffusionIP::GetPenaltyFactor(
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                  Array<OneD, NekDouble >                       &factor)
        {
            MultiRegions::ExpListSharedPtr tracelist = fields[0]->GetTrace();
            std::shared_ptr<LocalRegions::ExpansionVector> traceExp= tracelist->GetExp();
            int ntotTrac            = (*traceExp).size();
            int nTracPnt,noffset;
            
            const MultiRegions::LocTraceToTraceMapSharedPtr locTraceToTraceMap = fields[0]->GetlocTraceToTraceMap();
            
            const Array<OneD, const Array<OneD, int >> LRAdjExpid  =   locTraceToTraceMap->GetLeftRightAdjacentExpId();
            const Array<OneD, const Array<OneD, bool>> LRAdjflag   =   locTraceToTraceMap->GetLeftRightAdjacentExpFlag();

            std::shared_ptr<LocalRegions::ExpansionVector> fieldExp= fields[0]->GetExp();

            Array<OneD, NekDouble > factorFwdBwd(2,0.0);

            NekDouble spaceDim    =   NekDouble( fields[0]->GetCoordim(0) );

            for(int ntrace = 0; ntrace < ntotTrac; ++ntrace)
            {
                noffset     = tracelist->GetPhys_Offset(ntrace);
                nTracPnt    = tracelist->GetTotPoints(ntrace);

                factorFwdBwd[0] =   0.0;
                factorFwdBwd[1] =   0.0;
                
                for(int  nlr = 0; nlr < 2; nlr++)
                {
                    if(LRAdjflag[nlr][ntrace])
                    {
                        int numModes        = fields[0]->GetNcoeffs(LRAdjExpid[nlr][ntrace]);  
                        NekDouble numModesdir     = pow(NekDouble(numModes),(1.0/spaceDim));
                        factorFwdBwd[nlr]   =   1.0 * numModesdir * (numModesdir + 1.0);
                    }
                }

                for(int np = 0; np < nTracPnt; ++np)
                {
                    factor[noffset+np]    =   max(factorFwdBwd[0],factorFwdBwd[1]);
                }
            }
        }

        void DiffusionIP::GetPenaltyFactor_const(
            const Array<OneD, MultiRegions::ExpListSharedPtr>   &fields,
                  Array<OneD, NekDouble >                       &factor)
        {
            Vmath::Fill(factor.num_elements(),4.0,factor,1);
        }

        void DiffusionIP::v_ConsVarAveJump(
            const int                                           nConvectiveFields,
            const int                                           npnts,
            const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
            const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                  Array<OneD,       Array<OneD, NekDouble> >    &aver,
                  Array<OneD,       Array<OneD, NekDouble> >    &jump)
        {
            ConsVarAve(nConvectiveFields,npnts,vFwd,vBwd,aver);

            m_SpecialBndTreat(nConvectiveFields,aver);

            // note: here the jump is 2.0*(aver-vFwd) 
            //       because Viscous wall use a symmetry value as the Bwd, not the target one   
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Vsub(npnts,aver[i],1,vFwd[i],1,jump[i],1);
                Vmath::Smul(npnts,2.0,jump[i],1,jump[i],1);
            }
        }
        
        void DiffusionIP::ConsVarAve(
            const int                                           nConvectiveFields,
            const int                                           npnts,
            const Array<OneD, const Array<OneD, NekDouble> >    &vFwd,
            const Array<OneD, const Array<OneD, NekDouble> >    &vBwd,
                  Array<OneD,       Array<OneD, NekDouble> >    &aver)
        {
            NekDouble LinternalEngy =0.0;
            NekDouble RinternalEngy =0.0;
            NekDouble AinternalEngy =0.0;

            Array<OneD, NekDouble> Fweight (npnts,1.0);
            Array<OneD, NekDouble> Bweight;

            Bweight = m_tracBwdWeight;

            Vmath::Vsub(npnts,Fweight,1,Bweight,1,Fweight,1);

            for (int i = 0; i < nConvectiveFields-1; ++i)
            {
                Vmath::Vmul (npnts,Fweight,1,vFwd[i],1,aver[i],1);
                Vmath::Vvtvp(npnts,Bweight,1,vBwd[i],1,aver[i],1,aver[i],1);
            }
            
            int nengy = nConvectiveFields-1;
            int nvelst    = 1;
            int nveled    = nengy;
            for (int nt = 0; nt < npnts; ++nt)
            {
                LinternalEngy =0.0;
                for(int j=nvelst;j<nveled;j++)
                {
                    LinternalEngy += vFwd[j][nt]*vFwd[j][nt];
                }
                LinternalEngy *= -0.5/vFwd[0][nt];
                LinternalEngy += vFwd[nengy][nt];

                RinternalEngy =0.0;
                for(int j=nvelst;j<nveled;j++)
                {
                    RinternalEngy += vBwd[j][nt]*vBwd[j][nt];
                }
                RinternalEngy *= -0.5/vBwd[0][nt];
                RinternalEngy += vBwd[nengy][nt];

                AinternalEngy =0.0;
                aver[nengy][nt] = Fweight[nt]*LinternalEngy + Bweight[nt]*RinternalEngy;
                for(int j=nvelst;j<nveled;j++)
                {
                    AinternalEngy += aver[j][nt]*aver[j][nt];
                }
                aver[nengy][nt] += AinternalEngy*(0.5/aver[0][nt]);
            }
        }

        void DiffusionIP::CalTraceNumFlux_ReduceComm(
            const int                                                           nConvectiveFields,
            const int                                                           nDim,
            const int                                                           nPts,
            const int                                                           nTracePts,
            const NekDouble                                                     PenaltyFactor2,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, Array<OneD, NekDouble> >                          &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qfield,
            const Array<OneD, Array<OneD, NekDouble> >                          &vFwd,
            const Array<OneD, Array<OneD, NekDouble> >                          &vBwd,
            const Array<OneD, NekDouble >                                       &MuVarTrace,
                  Array<OneD, int >                                             &nonZeroIndexflux,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &traceflux,
                  Array<OneD, Array<OneD, NekDouble> >                          &solution_Aver,
                  Array<OneD, Array<OneD, NekDouble> >                          &solution_jump)
        {
            const MultiRegions::AssemblyMapDGSharedPtr                      TraceMap=fields[0]->GetTraceMap();

            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    numDerivBwd(nDim);
            //Fwd is also used for final numerical results
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    numDerivFwd(nDim);
            for (int nd = 0; nd < nDim; ++nd)
            {
                numDerivBwd[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                numDerivFwd[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    numDerivBwd[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
                    numDerivFwd[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
                }
            }

            // Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  tmparray3D = NullNekDoubleArrayofArrayofArray;
            // m_FunctorDerivBndCond(inarray,qfield,m_time,vFwd,tmparray3D);

            Array<OneD, NekDouble> Fwd(nTracePts,0.0);
            Array<OneD, NekDouble> Bwd(nTracePts,0.0);

            if(abs(PenaltyFactor2)>1.0E-12)
            {
                AddSecondDerivTOTrace_ReduceComm(nConvectiveFields,nDim,nPts,nTracePts,PenaltyFactor2,fields,qfield,numDerivFwd,numDerivBwd);
            }

            for (int nd = 0; nd < nDim; ++nd)
            {
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    Vmath::Zero(nTracePts, Bwd,1);
                    Vmath::Zero(nTracePts, Fwd,1);
                    fields[i]->GetFwdBwdTracePhysDeriv_serial(nd,qfield[nd][i], Fwd, Bwd);
                    Vmath::Svtvp(nTracePts,0.5,Bwd,1,numDerivBwd[nd][i],1,numDerivBwd[nd][i],1);
                    Vmath::Svtvp(nTracePts,0.5,Fwd,1,numDerivFwd[nd][i],1,numDerivFwd[nd][i],1);
                    TraceMap->UniversalTraceAssemble(numDerivBwd[nd][i]);
                    TraceMap->UniversalTraceAssemble(numDerivFwd[nd][i]);
                    Vmath::Vadd(nTracePts,numDerivFwd[nd][i],1,numDerivBwd[nd][i],1,numDerivFwd[nd][i],1);
                }
            }
        
            ConsVarAveJump(nConvectiveFields,nTracePts,vFwd,vBwd,solution_Aver,solution_jump);

            Array<OneD, NekDouble>  jumpTmp         =   Fwd;
            Array<OneD, NekDouble>  PenaltyFactor   =   Bwd;
            GetPenaltyFactor_const(fields,PenaltyFactor);

            Vmath::Vmul(nTracePts,PenaltyFactor,1, m_traceNormDirctnElmtLengthRecip,1,PenaltyFactor,1);
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Vmul(nTracePts,solution_jump[i],1, PenaltyFactor,1,jumpTmp,1);
                for (int nd = 0; nd < nDim; ++nd)
                {
                    Vmath::Vvtvp(nTracePts, m_traceNormals[nd],1,jumpTmp,1, numDerivFwd[nd][i],1, numDerivFwd[nd][i],1);
                }
            }
            jumpTmp         =   NullNekDouble1DArray;
            PenaltyFactor   =   NullNekDouble1DArray;

            // Calculate normal viscous flux
            m_FunctorDiffusionfluxCons(nConvectiveFields,nDim,solution_Aver,numDerivFwd,traceflux,nonZeroIndexflux,m_traceNormals,MuVarTrace);
        }
        
        void DiffusionIP::AddSecondDerivTOTrace_ReduceComm(
            const int                                                           nConvectiveFields,
            const int                                                           nDim,
            const int                                                           nPts,
            const int                                                           nTracePts,
            const NekDouble                                                     PenaltyFactor2,
            const Array<OneD, MultiRegions::ExpListSharedPtr>                   &fields,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > >      &qfield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &numDerivFwd,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >            &numDerivBwd)
        {
            Array<OneD, NekDouble> Fwd(nTracePts,0.0);
            Array<OneD, NekDouble> Bwd(nTracePts,0.0);
            Array<OneD,NekDouble>  tmp(nTracePts,0.0);

            Array<OneD, Array<OneD, NekDouble> > elmt2ndDerv(nDim);
            for(int nd1=0; nd1<nDim; nd1++)
            {
                elmt2ndDerv[nd1]    =   Array<OneD, NekDouble>(nPts,0.0);
            }

            Array<OneD, Array<OneD, NekDouble> > qtmp(3);
            for(int nd=0; nd<3; nd++)
            {
                qtmp[nd]    =   NullNekDouble1DArray;
            }
            for(int nd2=0; nd2<nDim; nd2++)
            {
                qtmp[nd2]    =   elmt2ndDerv[nd2];
            }

            Vmath::Smul(nTracePts,PenaltyFactor2,m_traceNormDirctnElmtLength,1,tmp,1);
            // the derivatives are assumed to be exchangable  
            for(int nd1=0; nd1<nDim; nd1++)
            {
                for(int i = 0; i < nConvectiveFields; ++i)
                {
                    fields[i]->PhysDeriv(qfield[nd1][i], qtmp[0], qtmp[1], qtmp[2]);

                    for(int nd2=nd1; nd2<nDim; nd2++)
                    {
                        Vmath::Zero(nTracePts,Bwd,1);
                        fields[i]->GetFwdBwdTracePhysDeriv_serial(nd2,elmt2ndDerv[nd2], Fwd, Bwd);
                        Vmath::Vmul(nTracePts,tmp,1,Bwd,1,Bwd,1);
                        Vmath::Vvtvp(nTracePts,m_traceNormals[nd2],1,Bwd,1,numDerivBwd[nd1][i],1,numDerivBwd[nd1][i],1);
                        Vmath::Vmul(nTracePts,tmp,1,Fwd,1,Fwd,1);
                        Vmath::Vvtvm(nTracePts,m_traceNormals[nd2],1,Fwd,1,numDerivFwd[nd1][i],1,numDerivFwd[nd1][i],1);
                        Vmath::Neg(nTracePts,numDerivFwd[nd1][i],1);

                        if(nd2!=nd1)
                        {
                            Vmath::Vvtvp(nTracePts,m_traceNormals[nd1],1,Bwd,1,numDerivBwd[nd2][i],1,numDerivBwd[nd2][i],1);
                            Vmath::Vvtvm(nTracePts,m_traceNormals[nd1],1,Fwd,1,numDerivFwd[nd2][i],1,numDerivFwd[nd2][i],1);
                            Vmath::Neg(nTracePts,numDerivFwd[nd2][i],1);
                        }
                    }
                }
            }
        }

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF
        void DiffusionIP::v_MinusVolumDerivJacToMat( 
            const int                                               nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr>       &pFields,
            const Array<OneD, const Array<OneD, DNekMatSharedPtr> > &ElmtJac,
            const int                                               nfluxDir, 
            const int                                               nDervDir, 
                  Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >    &gmtxarray)
        {
            MultiRegions::ExpListSharedPtr explist = pFields[0];
                std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
            int ntotElmt            = (*pexp).size();
            int nElmtPnt,nElmtCoef;

            NekDouble tmp;
            DNekMatSharedPtr        tmpGmtx,ElmtMat;

            Array<OneD, DNekMatSharedPtr>  mtxPerVar(ntotElmt);
            Array<OneD, Array<OneD, NekDouble> > JacArray(ntotElmt);
            Array<OneD, int > elmtpnts(ntotElmt);
            Array<OneD, int > elmtcoef(ntotElmt);
            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                nElmtCoef           = (*pexp)[nelmt]->GetNcoeffs();
                nElmtPnt            = (*pexp)[nelmt]->GetTotPoints();
                elmtpnts[nelmt]     =   nElmtPnt;
                elmtcoef[nelmt]     =   nElmtCoef;
                mtxPerVar[nelmt]    =MemoryManager<DNekMat>
                                    ::AllocateSharedPtr(nElmtCoef, nElmtPnt);
                JacArray[nelmt]    =Array<OneD, NekDouble>(nElmtPnt,0.0);
            }

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtPnt            = elmtpnts[nelmt];
                        for(int npnt = 0; npnt < nElmtPnt; npnt++)
                        {
                            JacArray[nelmt][npnt]   =   (*(ElmtJac[nelmt][npnt]))(m,n);
                        }
                    }

                    // explist->GetMatIpwrtdbWeightBwd(JacArray,nDirctn,mtxPerVar);
                    explist->GetMatIpwrtDeriveBase(JacArray,nfluxDir,mtxPerVar);
                    explist->AddRightIPTPhysDerivBase(nDervDir,mtxPerVar,mtxPerVar);

                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtCoef       = elmtcoef[nelmt];
                        nElmtPnt        = elmtpnts[nelmt];

                        tmpGmtx         = gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                        ElmtMat         = mtxPerVar[nelmt];

                        for(int ncl = 0; ncl < nElmtCoef; ncl++)
                        {
                            for(int nrw = 0; nrw < nElmtCoef; nrw++)
                            {
                                tmp   =   (*tmpGmtx)(nrw,ncl) - (*ElmtMat)(nrw,ncl);
                                tmpGmtx->SetValue(nrw,ncl,tmp);
                            }
                        }
                    }
                }
            }
        }
#endif

    }
}
