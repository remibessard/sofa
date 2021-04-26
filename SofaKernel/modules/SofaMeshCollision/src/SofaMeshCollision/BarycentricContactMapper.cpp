/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#define SOFA_COMPONENT_COLLISION_BARYCENTRICCONTACTMAPPER_CPP
#include <SofaMeshCollision/BarycentricContactMapper.inl>
#include <sofa/helper/Factory.inl>

namespace sofa::component::collision
{

using namespace defaulttype;

ContactMapperCreator< ContactMapper<LineCollisionModel<sofa::defaulttype::Vec3Types>> > LineContactMapperClass("default",true);
ContactMapperCreator< ContactMapper<TriangleCollisionModel<sofa::defaulttype::Vec3Types>> > TriangleContactMapperClass("default",true);
ContactMapperCreator< ContactMapper<CapsuleCollisionModel<sofa::defaulttype::Vec3Types>> > CapsuleContactMapperClass("default",true);
ContactMapperCreator< ContactMapper<CylinderCollisionModel<sofa::defaulttype::Vec3Types>> > CylinderContactMapperClass("default", true);

template class SOFA_SOFAMESHCOLLISION_API ContactMapper<LineCollisionModel<sofa::defaulttype::Vec3Types>, sofa::defaulttype::Vec3Types>;
template class SOFA_SOFAMESHCOLLISION_API ContactMapper<TriangleCollisionModel<sofa::defaulttype::Vec3Types>, sofa::defaulttype::Vec3Types>;
template class SOFA_SOFAMESHCOLLISION_API ContactMapper<CapsuleCollisionModel<sofa::defaulttype::Vec3Types>, sofa::defaulttype::Vec3Types>;
template class SOFA_SOFAMESHCOLLISION_API ContactMapper<CylinderCollisionModel<sofa::defaulttype::Vec3Types>, sofa::defaulttype::Vec3Types>;

} //namespace sofa::component::collision
