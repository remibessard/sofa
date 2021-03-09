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
#include <SofaMeshCollision/config.h>

#include <SofaBaseCollision/NewProximityIntersection.h>
#include <SofaMeshCollision/TriangleModel.h>
#include <SofaMeshCollision/LineModel.h>
#include <SofaMeshCollision/MeshIntTool.h>
#include <SofaBaseCollision/IntrUtility3.h>
#include <SofaBaseCollision/BaseIntTool.h>

namespace sofa::component::collision
{

class SOFA_SOFAMESHCOLLISION_API MeshNewProximityIntersection : public core::collision::BaseIntersector
{
    typedef NewProximityIntersection::OutputVector OutputVector;

public:
    MeshNewProximityIntersection(NewProximityIntersection* object, bool addSelf=true);


    template <class T1,class T2>
    bool testIntersection(T1 & e1,T2 & e2){
        return BaseIntTool::testIntersection(e1,e2,intersection->getAlarmDistance());
    }


    int computeIntersection(Point&, Point&, OutputVector*);
    template <class T> int computeIntersection(TSphere<T>&, Point&, OutputVector*);
    int computeIntersection(Line&, Point&, OutputVector*);
    template <class T> int computeIntersection(Line&, TSphere<T>&, OutputVector*);
    int computeIntersection(Line&, Line&, OutputVector*);
    int computeIntersection(Triangle&, Point&, OutputVector*);

    template <class T> int computeIntersection(Triangle&, TSphere<T>&, OutputVector*);
    int computeIntersection(Triangle&, Line&, OutputVector*);

    int computeIntersection(Triangle&, Triangle&, OutputVector*);

    template <class T1,class T2>
    int computeIntersection(T1 & e1,T2 & e2,OutputVector* contacts){
        return MeshIntTool::computeIntersection(e1,e2,e1.getProximity() + e2.getProximity() + intersection->getAlarmDistance(),e1.getProximity() + e2.getProximity() + intersection->getContactDistance(),contacts);
    }

    static inline int doIntersectionLineLine(SReal dist2, const defaulttype::Vector3& p1, const defaulttype::Vector3& p2, const defaulttype::Vector3& q1, const defaulttype::Vector3& q2, OutputVector* contacts, int id, const defaulttype::Vector3& n=defaulttype::Vector3(), bool useNormal=false);

	static inline int doBarycentricIntersectionLinePoint(SReal dist2, const defaulttype::Vector3& p1, const defaulttype::Vector3& p2, const defaulttype::Vector3& barycentre, const defaulttype::Vector3& q, OutputVector* contacts, int id, bool swapElems = false);

    static inline int doIntersectionLinePoint(SReal dist2, const defaulttype::Vector3& p1, const defaulttype::Vector3& p2, const defaulttype::Vector3& q, OutputVector* contacts, int id, bool swapElems = false);

    static inline int doIntersectionTrianglePoint(SReal dist2, int flags, const defaulttype::Vector3& p1, const defaulttype::Vector3& p2, const defaulttype::Vector3& p3, const defaulttype::Vector3& n, const defaulttype::Vector3& q, OutputVector* contacts, int id, bool swapElems = false, bool useNormal=false);

    static inline int doIntersectionTrianglePoint2(SReal dist2, int flags, const defaulttype::Vector3& p1, const defaulttype::Vector3& p2, const defaulttype::Vector3& p3, const defaulttype::Vector3& n, const defaulttype::Vector3& q, OutputVector* contacts, int id, bool swapElems = false);

protected:

    NewProximityIntersection* intersection;
};

} //namespace sofa::component::collision
