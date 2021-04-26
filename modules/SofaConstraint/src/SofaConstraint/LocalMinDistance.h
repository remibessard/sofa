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
#pragma once
#include <SofaConstraint/config.h>

#include <SofaBaseCollision/BaseProximityIntersection.h>
#include <sofa/helper/FnDispatcher.h>

#include <SofaBaseCollision/SphereModel.h>
#include <SofaMeshCollision/TriangleModel.h>
#include <SofaMeshCollision/LineModel.h>
#include <SofaMeshCollision/PointModel.h>
#include <SofaBaseCollision/CubeModel.h>
#include <SofaUserInteraction/RayModel.h>
#include <SofaBaseCollision/CapsuleModel.h>
#include <SofaBaseCollision/CylinderModel.h>

namespace sofa::component::collision
{

class SOFA_SOFACONSTRAINT_API LocalMinDistance : public BaseProximityIntersection
{
public:
    SOFA_CLASS(LocalMinDistance,BaseProximityIntersection);

    typedef core::collision::IntersectorFactory<LocalMinDistance> IntersectorFactory;

    Data<bool> filterIntersection; ///< Activate LMD filter
    Data<double> angleCone; ///< Filtering cone extension angle
    Data<double> coneFactor; ///< Factor for filtering cone angle computation
    Data<bool> useLMDFilters; ///< Use external cone computation (Work in Progress)

    int m_nb;
    defaulttype::Vector3 m_A;
    defaulttype::Vector3 m_B;
    defaulttype::Vector3 m_C;
    defaulttype::Vector3 m_T1,m_T2,m_T3;
    defaulttype::Vector3 m_P;
    defaulttype::Vector3 m_H;
    defaulttype::Vector3 m_X;
    defaulttype::Vector3 m_N;

protected:
    LocalMinDistance();
public:
    void init() override;

    bool testIntersection(Cube& ,Cube&);

    bool testIntersection(Point&, Point&);
    bool testIntersection(Sphere&, Point&);
    bool testIntersection(Sphere&, Sphere&);
    bool testIntersection(Line&, Point&);
    bool testIntersection(Line&, Sphere&);
    bool testIntersection(Line&, Line&);
    bool testIntersection(Triangle&, Point&);
    bool testIntersection(Triangle&, Sphere&);
    bool testIntersection(Ray&, Sphere&);
    bool testIntersection(Ray&, Triangle&);

    bool testIntersection(Cylinder&, Point&);
    //bool testIntersection(Cylinder&, Sphere&);
    //bool testIntersection(Cylinder&, Line&);
    //bool testIntersection(Cylinder&, Triangle&);
    //bool testIntersection(Cylinder&, Cylinder&);

    bool testIntersection(Capsule&, Point&);
    //bool testIntersection(Capsule&, Sphere&);
    bool testIntersection(Capsule&, Line&);
    bool testIntersection(Capsule&, Triangle&);
    //bool testIntersection(Capsule&, Capsule&);


    int computeIntersection(Cube&, Cube&, OutputVector*);
    int computeIntersection(Point&, Point&, OutputVector*);
    int computeIntersection(Sphere&, Point&, OutputVector*);
    int computeIntersection(Sphere&, Sphere&, OutputVector*);
    int computeIntersection(Line&, Point&, OutputVector*);
    int computeIntersection(Line&, Sphere&, OutputVector*);
    int computeIntersection(Line&, Line&, OutputVector*);
    int computeIntersection(Triangle&, Point&, OutputVector*);
    int doIntersectionTrianglePoint(SReal alarmDist, Triangle& e2, core::CollisionElementIterator& e1, const defaulttype::Vector3 p, OutputVector* contacts);
    int computeIntersection(Triangle&, Sphere&, OutputVector*);
    int computeIntersection(Ray&, Sphere&, OutputVector*);
    int computeIntersection(Ray&, Triangle&, OutputVector*);

    int computeIntersection(Cylinder&, Point&, OutputVector*);
    //int computeIntersection(Cylinder&, Sphere&, OutputVector*);
    //int computeIntersection(Cylinder&, Line&, OutputVector*);
    //int computeIntersection(Cylinder&, Triangle&, OutputVector*);
    //int computeIntersection(Cylinder&, Cylinder&, OutputVector*);

    int computeIntersection(Capsule&, Point&, OutputVector*);
    //int computeIntersection(Capsule&, Sphere&, OutputVector*);
    int computeIntersection(Capsule&, Line&, OutputVector*);
    int doIntersectionCapsuleLine(Capsule& e2, core::CollisionElementIterator& e1, const defaulttype::Vector3 l1, const defaulttype::Vector3 l2, OutputVector* contacts);
    int computeIntersection(Capsule&, Triangle&, OutputVector*);
    //int computeIntersection(Capsule&, Capsule&, OutputVector*);

    /// These methods check the validity of a found intersection.
    /// According to the local configuration around the found intersected primitive,
    /// we build a "Region Of Interest" geometric cone.
    /// Pertinent intersections have to belong to this cone, others are not taking into account anymore.
    bool testValidity(Sphere&, const defaulttype::Vector3&) { return true; }
    bool testValidity(Point&, const defaulttype::Vector3&);
    bool testValidity(Line&, const defaulttype::Vector3&);
    bool testValidity(Triangle&, const defaulttype::Vector3&);

    void draw(const core::visual::VisualParams* vparams) override;

    /// Actions to accomplish when the broadPhase is started. By default do nothing.
    void beginBroadPhase() override {}

    int beginIntersection(sofa::core::CollisionModel* /*model1*/, sofa::core::CollisionModel* /*model2*/, OutputVector* /*contacts*/)
    {
        return 0;
    }

private:
    double mainAlarmDistance;
    double mainContactDistance;
};

} // namespace sofa::component::collision
