/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#define INDEXVALUEMAPPER_CPP_

#include "Indices2ValuesMapper.inl"
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace engine
{

using namespace sofa;
using namespace sofa::defaulttype;

SOFA_DECL_CLASS(Indices2ValuesMapper)

int Indices2ValuesMapperClass = core::RegisterObject("?")
#ifndef SOFA_FLOAT
        .add< Indices2ValuesMapper<Vec3dTypes> >(true)
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
        .add< Indices2ValuesMapper<Vec3fTypes> >()
#endif //SOFA_DOUBLE
        ;

#ifndef SOFA_FLOAT
template class SOFA_ENGINE_API Indices2ValuesMapper<Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_ENGINE_API Indices2ValuesMapper<Vec3fTypes>;
#endif //SOFA_DOUBLE


} // namespace engine

} // namespace component

} // namespace sofa