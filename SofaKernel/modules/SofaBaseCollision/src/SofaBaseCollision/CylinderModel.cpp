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
#define SOFA_COMPONENT_COLLISION_CYLINDERCOLLISIONMODEL_CPP
#include <SofaBaseCollision/CylinderModel.inl>

namespace sofa::component::collision
{

using namespace sofa::defaulttype;
using namespace sofa::core::collision;
using namespace helper;

int RigidCylinderCollisionModelClass = core::RegisterObject("Collision model which represents a set of rigid cylinders")
        .add<  CylinderCollisionModel<defaulttype::Vec3dTypes> >()
        //.add<  CylinderCollisionModel<defaulttype::Rigid3Types> >()

        .addAlias("Cylinder")
        .addAlias("CylinderModel")
        ;

/////
///// RIGID3 IMPLEMENTATION
/////
const sofa::defaulttype::Quaternion CylinderCollisionModel<defaulttype::Rigid3Types >::orientation(Index index)const {
    return m_mstate->read(core::ConstVecCoordId::position())->getValue()[index].getOrientation();
}
// 
///--------------------------------------------------------------------------------------------------------------------///
///
/// VEC3 IMPLEMENTATION
///

void CylinderCollisionModel<defaulttype::Vec3dTypes>::init()
{
    this->CollisionModel::init();
    m_mstate = dynamic_cast<core::behavior::MechanicalState<DataTypes>*> (getContext()->getMechanicalState());
    //this->getContext()->get(mpoints);

    if (m_mstate == nullptr)
    {
        msg_error() << "LineModel requires a Vec3 Mechanical Model";
        d_componentState.setValue(ComponentState::Invalid);
        return;
    }

    //simulation::Node* node = dynamic_cast<simulation::Node*>(this->getContext());
    //if (node != 0)
    //{
    //    m_lmdFilter = node->getNodeObject< LineLocalMinDistanceFilter >();
    //}

    if (l_topology.empty())
    {
        msg_info() << "link to Topology container should be set to ensure right behavior. First Topology found in current context will be used.";
        l_topology.set(this->getContext()->getMeshTopologyLink());
    }

    m_topology = l_topology.get();
    msg_info() << "Topology path used: '" << l_topology.getLinkedPath() << "'";

    if (!m_topology)
    {
        msg_error() << "No topology component found at path: " << l_topology.getLinkedPath() << ", nor in current context: " << this->getContext()->name << ". LineCollisionModel<sofa::defaulttype::Vec3Types> requires a MeshTopology";
        sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }
    resize(m_topology->getNbEdges());
}


//void CylinderCollisionModel<defaulttype::Vec3dTypes>::resize(Size size)
//{
//    //this->core::CollisionModel::resize(size);
//
//    VecReal & Cylinder_radii = *d_cylinder_radii.beginEdit();
//    //VecReal & Cylinder_heights = *d_cylinder_heights.beginEdit();
//    //VecAxisCoord & Cylinder_local_axes = *d_cylinder_local_axes.beginEdit();
//
//    if (Cylinder_radii.size() < size)
//    {
//        while (Cylinder_radii.size() < size)
//            Cylinder_radii.push_back(d_default_radius.getValue());
//    }
//    else
//    {
//        Cylinder_radii.reserve(size);
//    }
//
//    //if (Cylinder_heights.size() < size)
//    //{
//    //    while (Cylinder_heights.size() < size)
//    //        Cylinder_heights.push_back(d_default_height.getValue());
//    //}
//    //else
//    //{
//    //    Cylinder_heights.reserve(size);
//    //}
//
//    //if (Cylinder_local_axes.size() < size)
//    //{
//    //    while (Cylinder_local_axes.size() < size)
//    //        Cylinder_local_axes.push_back(d_default_local_axis.getValue());
//    //}
//    //else
//    //{
//    //    Cylinder_local_axes.reserve(size);
//    //}
//
//    d_cylinder_radii.endEdit();
//    //d_cylinder_heights.endEdit();
//    //d_cylinder_local_axes.endEdit();
//}


typename CylinderCollisionModel<defaulttype::Vec3dTypes>::Coord CylinderCollisionModel< defaulttype::Vec3dTypes >::point1(Index i) const
{
    return this->m_mstate->read(core::ConstVecCoordId::position())->getValue()[m_topology->getEdge(i)[0]];
    //return  center(i) - axis(i) * height(i) / 2.0;
}

typename CylinderCollisionModel<defaulttype::Vec3dTypes>::Coord CylinderCollisionModel< defaulttype::Vec3dTypes >::point2(Index i) const
{
    return this->m_mstate->read(core::ConstVecCoordId::position())->getValue()[m_topology->getEdge(i)[1]];
    //return  center(i) + axis(i) * height(i) / 2.0;
}

typename CylinderCollisionModel<defaulttype::Vec3dTypes>::Coord CylinderCollisionModel<defaulttype::Vec3dTypes >::axis(Index index) const {
    //Coord ax = d_cylinder_local_axes.getValue()[index];
    //const sofa::defaulttype::Quaternion & ori = orientation(index);
    //return ori.rotate(ax);
    return (this->m_mstate->read(core::ConstVecCoordId::position())->getValue()[m_topology->getEdge(index)[1]]
        - this->m_mstate->read(core::ConstVecCoordId::position())->getValue()[m_topology->getEdge(index)[0]]).normalized();
}

typename CylinderCollisionModel<defaulttype::Vec3dTypes>::Coord CylinderCollisionModel< defaulttype::Vec3dTypes >::center(Index i)const {
    return
        (this->m_mstate->read(core::ConstVecCoordId::position())->getValue()[m_topology->getEdge(i)[0]]
            + this->m_mstate->read(core::ConstVecCoordId::position())->getValue()[m_topology->getEdge(i)[1]])
        / 2.0;
    //return DataTypes::getCPos((m_mstate->read(core::ConstVecCoordId::position())->getValue())[i]);
}

typename TCylinder<defaulttype::Vec3dTypes>::Coord TCylinder< defaulttype::Vec3dTypes >::point1() const
{
    return this->model->point1(this->index);
}

typename TCylinder<defaulttype::Vec3dTypes>::Coord TCylinder<defaulttype::Vec3dTypes >::point2() const
{
    return this->model->point2(this->index);
}

typename TCylinder<defaulttype::Vec3dTypes>::Coord TCylinder<defaulttype::Vec3dTypes >::axis() const {
    return this->model->axis(this->index);
}

typename TCylinder<defaulttype::Vec3dTypes>::Coord TCylinder< defaulttype::Vec3dTypes >::center() const
{
    return this->model->center(this->index);
}



const typename CylinderCollisionModel<defaulttype::Vec3dTypes>::Coord CylinderCollisionModel<defaulttype::Vec3dTypes >::velocity(Index i) const {
    return (this->m_mstate->read(core::ConstVecDerivId::velocity())->getValue()[m_topology->getEdge(i)[0]]
        + this->m_mstate->read(core::ConstVecDerivId::velocity())->getValue()[m_topology->getEdge(i)[1]])
        / ((Real)(2.0));
    //return DataTypes::getDPos(((m_mstate->read(core::ConstVecDerivId::velocity())->getValue()))[index]);
}

//const typename TCylinder<defaulttype::Vec3dTypes>::Coord& TCylinder<defaulttype::Vec3dTypes >::v() const { return this->model->velocity(this->index); }


///TODO compute correct orientation for Vec3d
const sofa::defaulttype::Quaternion CylinderCollisionModel<defaulttype::Vec3dTypes >::orientation(Index index)const {
    return sofa::defaulttype::Quaternion::createFromRotationVector(0, 1, 0);// test;// new Quat(0, 0, 0, 1);//  m_mstate->read(core::ConstVecCoordId::position())->getValue()[index];
}


typename CylinderCollisionModel<defaulttype::Vec3dTypes>::Real CylinderCollisionModel<defaulttype::Vec3dTypes>::height(Index index) const {
    //if (this->m_mstate == nullptr)
    //    return 0.0;
    m_topology->getEdge(index);
    //std::cout << "index: " << index << "   and m_topology->getEdge(index):" << m_topology->getEdge(index) 
    //    << "    and mstatesize: " << this->m_mstate->getSize()<< std::endl;
    return (this->m_mstate->read(core::ConstVecCoordId::position())->getValue()[m_topology->getEdge(index)[1]]
        - this->m_mstate->read(core::ConstVecCoordId::position())->getValue()[m_topology->getEdge(index)[0]]).norm();
    //return ((d_cylinder_heights.getValue()))[index];
}

typename TCylinder<defaulttype::Vec3dTypes>::Real TCylinder<defaulttype::Vec3dTypes >::radius() const
{
    return this->model->radius(this->index);
}

template class SOFA_SOFABASECOLLISION_API TCylinder<defaulttype::Vec3dTypes>;
template class SOFA_SOFABASECOLLISION_API CylinderCollisionModel<defaulttype::Vec3dTypes>;
//template class SOFA_SOFABASECOLLISION_API TCylinder<defaulttype::Rigid3Types>;
//template class SOFA_SOFABASECOLLISION_API CylinderCollisionModel<defaulttype::Rigid3Types>;

} // namespace sofa::component::collision
