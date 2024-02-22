////////////////////////////////////////////////////////////////////////////////
//
//  File: Interfaces.cpp
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
//  Description: Python wrapper for Interface.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Python/NekPyConfig.hpp>
#include <SpatialDomains/Movement/InterfaceInterpolation.h>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace Nektar;
using namespace Nektar::SpatialDomains;

std::shared_ptr<Interface> Interface_Init(int indx, const CompositeMap &edge,
                                          bool skipCoordCheck)
{
    return std::make_shared<Interface>(indx, edge, skipCoordCheck);
}

std::shared_ptr<InterfacePair> InterfacePair_Init(
    const InterfaceShPtr &leftInterface, const InterfaceShPtr &rightInterface)
{
    return std::make_shared<InterfacePair>(leftInterface, rightInterface);
}

void export_Interfaces()
{
    py::class_<std::vector<unsigned int>>("UIntList")
        .def(py::vector_indexing_suite<std::vector<unsigned int>, true>());

    py::class_<Interface, std::shared_ptr<Interface>>("Interface", py::no_init)
        .def("__init__", py::make_constructor(&Interface_Init))
        .def<const std::map<int, GeometrySharedPtr> &(Interface::*)() const>(
            "GetEdge", &Interface::GetEdge, py::return_internal_reference<>())
        .def("IsEmpty", &Interface::IsEmpty)
        .def("GetId", &Interface::GetId,
             py::return_value_policy<py::copy_non_const_reference>())
        .def("GetOppInterace", &Interface::GetOppInterface,
             py::return_internal_reference<>())
        .def("GetCompositeIDs", &Interface::GetCompositeIDs,
             py::return_value_policy<py::copy_const_reference>());

    py::class_<InterfacePair, std::shared_ptr<InterfacePair>>("InterfacePair",
                                                              py::no_init)
        .def("__init__", py::make_constructor(&InterfacePair_Init))
        .def("GetLeftInterface", &InterfacePair::GetLeftInterface,
             py::return_value_policy<py::copy_const_reference>())
        .def("GetRightInterface", &InterfacePair::GetRightInterface,
             py::return_value_policy<py::copy_const_reference>());
}
