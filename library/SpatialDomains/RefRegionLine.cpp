////////////////////////////////////////////////////////////////////////////////
//
//  File: RefRegionLine.cpp
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
//  Description: Linear surface for the refinement region.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/RefRegionLine.h>

using namespace std;

namespace Nektar::SpatialDomains
{

RefRegionLine::RefRegionLine(const unsigned int coordim, NekDouble radius,
                             std::vector<NekDouble> coord1,
                             std::vector<NekDouble> coord2,
                             std::vector<unsigned int> numModes,
                             std::vector<unsigned int> numPoints)
    : RefRegion(coordim, radius, coord1, coord2, numModes, numPoints)
{
}

RefRegionLine::~RefRegionLine()
{
}

/**
 * @brief Check if vertex is inside the line.
 *
 * @param coords    coordinates of the vertex
 * @return          true or false depending on if the vertex is inside
 *                  or not of the surface defined by the user.
 */
bool RefRegionLine::v_Contains(const Array<OneD, NekDouble> &coords)
{
    // Is P between the region below?
    // xa, xb plus the radius
    if (m_coord1[0] < m_coord2[0])
    {
        if ((m_coord1[0] - m_radius < coords[0]) &&
            (m_coord2[0] + m_radius > coords[0]))
        {
            return true;
        }
    }
    else
    {
        if ((m_coord1[0] + m_radius > coords[0]) &&
            (m_coord2[0] - m_radius < coords[0]))
        {
            return true;
        }
    }

    return false;
}

} // namespace Nektar::SpatialDomains
