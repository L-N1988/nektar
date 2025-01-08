///////////////////////////////////////////////////////////////////////////////
//
// File: SpatialDomains.h
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
// Description: NekPy configuration for SpatialDomains to defined opaque types.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBRARY_SPATIALDOMAINS_PYTHON_SPATIALDOMAINS_H
#define NEKTAR_LIBRARY_SPATIALDOMAINS_PYTHON_SPATIALDOMAINS_H

#include <LibUtilities/Python/NekPyConfig.hpp>

#include <SpatialDomains/MeshGraph.h>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

// Define common opaque types
PYBIND11_MAKE_OPAQUE(LibUtilities::FieldMetaDataMap);
PYBIND11_MAKE_OPAQUE(std::vector<GeometrySharedPtr>);
PYBIND11_MAKE_OPAQUE(PointGeomMap);
PYBIND11_MAKE_OPAQUE(SegGeomMap);
PYBIND11_MAKE_OPAQUE(QuadGeomMap);
PYBIND11_MAKE_OPAQUE(TriGeomMap);
PYBIND11_MAKE_OPAQUE(TetGeomMap);
PYBIND11_MAKE_OPAQUE(PrismGeomMap);
PYBIND11_MAKE_OPAQUE(PyrGeomMap);
PYBIND11_MAKE_OPAQUE(HexGeomMap);
PYBIND11_MAKE_OPAQUE(CurveMap);
PYBIND11_MAKE_OPAQUE(CompositeMap);
PYBIND11_MAKE_OPAQUE(std::map<int, CompositeMap>);

#endif
