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
#include <SofaMeshCollision/MeshNewProximityIntersection.h>
#include <SofaBaseCollision/NewProximityIntersection.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/core/collision/Intersection.inl>

namespace sofa::component::collision
{

inline int MeshNewProximityIntersection::doIntersectionLineLine(SReal dist2, const defaulttype::Vector3& p1, const defaulttype::Vector3& p2, const defaulttype::Vector3& q1, const defaulttype::Vector3& q2, OutputVector* contacts, int id, const defaulttype::Vector3& /*n*/, bool /*useNormal*/)
{  
    defaulttype::Vector3 p,q;
    IntrUtil<SReal>::segNearestPoints(p1,p2,q1,q2,p,q);
	SReal alpha1 = (p - p1).norm() / (p2 - p1).norm();
	SReal alpha2 = (q - q1).norm() / (q2 - q1).norm();
    defaulttype::Vector3 pq = q-p;
    SReal norm2 = pq.norm2();

    if (norm2 >= dist2)
        return 0;

    contacts->resize(contacts->size()+1);
    core::collision::DetectionOutput *detection = &*(contacts->end()-1);
    detection->id = id;
    detection->point[0]=p;
    detection->point[1]=q;
    detection->value = helper::rsqrt(norm2);
    detection->normal = pq / detection->value;

#ifdef DETECTIONOUTPUT_BARYCENTRICINFO
	detection->baryCoords[0] = defaulttype::Vector3(alpha1, 0, 0);
	detection->baryCoords[1] = defaulttype::Vector3(alpha2, 0, 0);
#endif
	std::cout << "ll ";
    //detection->value -= contactDist;
    return 1;
}

inline int MeshNewProximityIntersection::doBarycentricIntersectionLinePoint(SReal dist2, const defaulttype::Vector3& p1, const defaulttype::Vector3& p2, const defaulttype::Vector3& barycentre, const defaulttype::Vector3& q, OutputVector* contacts, int id, bool swapElems)
{
	defaulttype::Vector3 p;
	p = IntrUtil<SReal>::nearestPointOnSeg(p1, p2, q);
	SReal alpha = (p - p1).norm() / (p2 - p1).norm();
	defaulttype::Vector3 pq = q - p;
	SReal norm2 = pq.norm2();
	if (norm2 >= dist2)
		return 0;
	if (alpha < -1e-7 || 1 - alpha < -1e-7)
		return 0;

	//const SReal contactDist = getContactDistance() + e1.getProximity() + e2.getProximity();
	contacts->resize(contacts->size() + 1);
	core::collision::DetectionOutput *detection = &*(contacts->end() - 1);

	//detection.elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e2, e1);
	detection->id = id;
	detection->value = helper::rsqrt(norm2);
	if (swapElems)
	{
		detection->point[0] = q;
		detection->point[1] = p;
		detection->normal = -pq / detection->value;
	}
	else
	{
		detection->point[0] = p;
		detection->point[1] = q;
		detection->normal = pq / detection->value;
	}
#ifdef DETECTIONOUTPUT_BARYCENTRICINFO
	detection->baryCoords[1] = defaulttype::Vector3(alpha, 0, 0);
	detection->baryCoords[0] = barycentre;
#endif

	//detection.value -= contactDist;
	return 1;
}

inline int MeshNewProximityIntersection::doIntersectionLinePoint(SReal dist2, const defaulttype::Vector3& p1, const defaulttype::Vector3& p2, const defaulttype::Vector3& q, OutputVector* contacts, int id, bool swapElems)
{
    defaulttype::Vector3 p;
    p = IntrUtil<SReal>::nearestPointOnSeg(p1,p2,q);

    defaulttype::Vector3 pq = q-p;
    SReal norm2 = pq.norm2();
    if (norm2 >= dist2)
        return 0;

    contacts->resize(contacts->size()+1);
    core::collision::DetectionOutput *detection = &*(contacts->end()-1);

    detection->id = id;
    detection->value = helper::rsqrt(norm2);
    if (swapElems)
    {
        detection->point[0]=q;
        detection->point[1]=p;
        detection->normal = -pq / detection->value;
    }
    else
    {
        detection->point[0]=p;
        detection->point[1]=q;
        detection->normal = pq / detection->value;
    }

    return 1;
}

inline int MeshNewProximityIntersection::doIntersectionTrianglePoint2(SReal dist2, int flags, const defaulttype::Vector3& p1, const defaulttype::Vector3& p2, const defaulttype::Vector3& p3, const defaulttype::Vector3& /*n*/, const defaulttype::Vector3& q, OutputVector* contacts, int id, bool swapElems)
{
    const defaulttype::Vector3 AB = p2-p1;
    const defaulttype::Vector3 AC = p3-p1;
    const defaulttype::Vector3 AQ = q -p1;
    defaulttype::Matrix2 A;
    defaulttype::Vector2 b;
    A[0][0] = AB*AB;
    A[1][1] = AC*AC;
    A[0][1] = A[1][0] = AB*AC;
    b[0] = AQ*AB;
    b[1] = AQ*AC;
    const SReal det = defaulttype::determinant(A);

    SReal alpha = 0.5;
    SReal beta = 0.5;

    alpha = (b[0]*A[1][1] - b[1]*A[0][1])/det;
    beta  = (b[1]*A[0][0] - b[0]*A[1][0])/det;
    if (alpha < 0.000001 || beta < 0.000001 || alpha + beta > 0.999999)
    {
        // nearest point is on an edge or corner
        // barycentric coordinate on AB
        SReal pAB = b[0] / A[0][0]; // AQ*AB / AB*AB
        // barycentric coordinate on AC
        SReal pAC = b[1] / A[1][1]; // AQ*AB / AB*AB
        if (pAB < 0.000001 && pAC < 0.0000001)
        {
            // closest point is A
            if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P1)) return 0; // this corner is not considered
            alpha = 0.0;
            beta = 0.0;
        }
        else if (pAB < 0.999999 && beta < 0.000001)
        {
            // closest point is on AB
            if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_E12)) return 0; // this edge is not considered
            alpha = pAB;
            beta = 0.0;
        }
        else if (pAC < 0.999999 && alpha < 0.000001)
        {
            // closest point is on AC
            if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_E12)) return 0; // this edge is not considered
            alpha = 0.0;
            beta = pAC;
        }
        else
        {
            // barycentric coordinate on BC
            // BQ*BC / BC*BC = (AQ-AB)*(AC-AB) / (AC-AB)*(AC-AB) = (AQ*AC-AQ*AB + AB*AB-AB*AC) / (AB*AB+AC*AC-2AB*AC)
            SReal pBC = (b[1] - b[0] + A[0][0] - A[0][1]) / (A[0][0] + A[1][1] - 2*A[0][1]); // BQ*BC / BC*BC
            if (pBC < 0.000001)
            {
                // closest point is B
                if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P2)) return 0; // this edge is not considered
                alpha = 1.0;
                beta = 0.0;
            }
            else if (pBC > 0.999999)
            {
                // closest point is C
                if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P3)) return 0; // this edge is not considered
                alpha = 0.0;
                beta = 1.0;
            }
            else
            {
                // closest point is on BC
                if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_E31)) return 0; // this edge is not considered
                alpha = 1.0-pBC;
                beta = pBC;
            }
        }
    }

    defaulttype::Vector3 p, pq;
    p = p1 + AB * alpha + AC * beta;
    pq = q-p;
    SReal norm2 = pq.norm2();
    if (pq.norm2() >= dist2)
        return 0;

    contacts->resize(contacts->size()+1);
    core::collision::DetectionOutput *detection = &*(contacts->end()-1);
    detection->id = id;
    detection->value = helper::rsqrt(norm2);
    if (swapElems)
    {
        detection->point[1]=p;
        detection->normal = -pq / detection->value;
    }
    else
    {
        detection->point[0]=p;
        detection->normal = pq / detection->value;
    }
    return 1;
}



inline int MeshNewProximityIntersection::doIntersectionTrianglePoint(SReal dist2, int flags, const defaulttype::Vector3& p1, const defaulttype::Vector3& p2, const defaulttype::Vector3& p3, const defaulttype::Vector3& barycoord, const defaulttype::Vector3& q, OutputVector* contacts, int id, bool swapElems, bool /*useNormal*/)
{
    const defaulttype::Vector3 AB = p2-p1;
    const defaulttype::Vector3 AC = p3-p1;
    const defaulttype::Vector3 AQ = q -p1;
    defaulttype::Matrix2 A;
    defaulttype::Vector2 b;
    A[0][0] = AB*AB;
    A[1][1] = AC*AC;
    A[0][1] = A[1][0] = AB*AC;
    b[0] = AQ*AB;
    b[1] = AQ*AC;
    const SReal det = defaulttype::determinant(A);

    SReal alpha = 0.5;
    SReal beta = 0.5;
    const SReal epsilon=std::numeric_limits<SReal>::epsilon();

    alpha = (b[0]*A[1][1] - b[1]*A[0][1])/det;
    beta  = (b[1]*A[0][0] - b[0]*A[1][0])/det;
    if (alpha < epsilon || beta < epsilon || alpha + beta > 1 - epsilon)
    {
            //return 0;
        // nearest point is on an edge or corner
        // barycentric coordinate on AB
        SReal pAB = b[0] / A[0][0]; // AQ*AB / AB*AB
        // barycentric coordinate on AC
        SReal pAC = b[1] / A[1][1]; // AQ*AB / AB*AB
        if (pAB < epsilon && pAC < epsilon)
        {
            // closest point is A
            if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P1)) return 0; // this corner is not considered
            alpha = 0.0;
            beta = 0.0;
        }
        else if (pAB < 1 - epsilon && pAB >= epsilon && beta < epsilon)
        {
            // closest point is on AB
            if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_E12)) return 0; // this edge is not considered
            alpha = pAB;
            beta = 0.0;
        }
        else if (pAC < 1 - epsilon && pAC >= epsilon && alpha < epsilon)
        {
            // closest point is on AC
            if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_E31)) return 0; // this edge is not considered
            alpha = 0.0;
            beta = pAC;
        }
        else
        {
            // barycentric coordinate on BC
            // BQ*BC / BC*BC = (AQ-AB)*(AC-AB) / (AC-AB)*(AC-AB) = (AQ*AC-AQ*AB + AB*AB-AB*AC) / (AB*AB+AC*AC-2AB*AC)
            SReal pBC = (b[1] - b[0] + A[0][0] - A[0][1]) / (A[0][0] + A[1][1] - 2*A[0][1]); // BQ*BC / BC*BC
            if (pBC < epsilon)
            {
                // closest point is B
                if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P2)) return 0; // this edge is not considered
                alpha = 1.0;
                beta = 0.0;
            }
            else if (pBC > 1 - epsilon)
            {
                // closest point is C
                if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P3)) return 0; // this edge is not considered
                alpha = 0.0;
                beta = 1.0;
            }
            else
            {
                // closest point is on BC
                if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_E23)) return 0; // this edge is not considered
                alpha = 1.0-pBC;
                beta = pBC;
            }
        }
    }

    defaulttype::Vector3 p, pq;
    p = p1 + AB * alpha + AC * beta;
    pq = q-p;
    SReal norm2 = pq.norm2();
    if (pq.norm2() >= dist2)
        return 0;

    contacts->resize(contacts->size()+1);
    core::collision::DetectionOutput *detection = &*(contacts->end()-1);
    detection->id = id;
    detection->value = helper::rsqrt(norm2);

    if (swapElems)
    {
        detection->point[0]=q;
        detection->point[1]=p;
        detection->normal = -pq / detection->value;
#ifdef DETECTIONOUTPUT_BARYCENTRICINFO
		detection->baryCoords[1] = defaulttype::Vector3(alpha, beta, 0);
		detection->baryCoords[0] = barycoord;
#endif
    }
    else
    {
        detection->point[0]=p;
        detection->point[1]=q;
        detection->normal = pq / detection->value;
#ifdef DETECTIONOUTPUT_BARYCENTRICINFO
		detection->baryCoords[0] = defaulttype::Vector3(alpha, beta, 0);
		detection->baryCoords[1] = barycoord;
#endif
    }
    return 1;
}

template<class T>
int MeshNewProximityIntersection::computeIntersection(TSphere<T>& e1, Point& e2, OutputVector* contacts)
{
    const SReal alarmDist = intersection->getAlarmDistance() + e1.getProximity() + e2.getProximity() + e1.r();
    int n = intersection->doIntersectionPointPoint(alarmDist*alarmDist, e1.center(), e2.p(), contacts, (e1.getCollisionModel()->getSize() > e2.getCollisionModel()->getSize()) ? e1.getIndex() : e2.getIndex());
    if (n>0)
    {
        const SReal contactDist = intersection->getContactDistance() + e1.getProximity() + e2.getProximity() + e1.r();
        for (OutputVector::iterator detection = contacts->end()-n; detection != contacts->end(); ++detection)
        {
            detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
            detection->value -= contactDist;
        }
    }
    return n;
}

template<class T>
int MeshNewProximityIntersection::computeIntersection(Line& e1, TSphere<T>& e2, OutputVector* contacts)
{
    const SReal alarmDist = intersection->getAlarmDistance() + e1.getProximity() + e2.getProximity() + e2.r();
    int n = doIntersectionLinePoint(alarmDist*alarmDist, e1.p1(),e1.p2(), e2.center(), contacts, e2.getIndex());
    if (n>0)
    {
        const SReal contactDist = intersection->getContactDistance() + e1.getProximity() + e2.getProximity() + e2.r();
        for (OutputVector::iterator detection = contacts->end()-n; detection != contacts->end(); ++detection)
        {
            detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
            detection->value -= contactDist;
        }
    }
    return n;
}

template<class T>
int MeshNewProximityIntersection::computeIntersection(Triangle& e1, TSphere<T>& e2, OutputVector* contacts)
{
    int flags = e1.flags();
    const SReal alarmDist = intersection->getAlarmDistance() + e1.getProximity() + e2.getProximity() + e2.r();
    const SReal dist2 = alarmDist*alarmDist;

    const defaulttype::Vector3 AB = e1.p2() - e1.p1();
    const defaulttype::Vector3 AC = e1.p3() - e1.p1();
    const defaulttype::Vector3 AQ = e2.center() - e1.p1();
    defaulttype::Matrix2 A;
    defaulttype::Vector2 b;
    A[0][0] = AB*AB;
    A[1][1] = AC*AC;
    A[0][1] = A[1][0] = AB*AC;
    b[0] = AQ*AB;
    b[1] = AQ*AC;
    const SReal det = defaulttype::determinant(A);

    SReal alpha = 0.5;
    SReal beta = 0.5;

    alpha = (b[0]*A[1][1] - b[1]*A[0][1])/det;
    beta  = (b[1]*A[0][0] - b[0]*A[1][0])/det;
    if (alpha < 0.000001 || beta < 0.000001 || alpha + beta > 0.999999)
    {
        // nearest point is on an edge or corner
        // barycentric coordinate on AB
        SReal pAB = b[0] / A[0][0]; // AQ*AB / AB*AB
        // barycentric coordinate on AC
        SReal pAC = b[1] / A[1][1]; // AQ*AB / AB*AB
        if (pAB < 0.000001 && pAC < 0.0000001)
        {
            // closest point is A
            if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P1)) return 0; // this corner is not considered
            alpha = 0.0;
            beta = 0.0;
        }
        else if (pAB < 0.999999 && beta < 0.000001)
        {
            // closest point is on AB
            if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_E12)) return 0; // this edge is not considered
            alpha = pAB;
            beta = 0.0;
        }
        else if (pAC < 0.999999 && alpha < 0.000001)
        {
            // closest point is on AC
            if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_E12)) return 0; // this edge is not considered
            alpha = 0.0;
            beta = pAC;
        }
        else
        {
            // barycentric coordinate on BC
            // BQ*BC / BC*BC = (AQ-AB)*(AC-AB) / (AC-AB)*(AC-AB) = (AQ*AC-AQ*AB + AB*AB-AB*AC) / (AB*AB+AC*AC-2AB*AC)
            SReal pBC = (b[1] - b[0] + A[0][0] - A[0][1]) / (A[0][0] + A[1][1] - 2*A[0][1]); // BQ*BC / BC*BC
            if (pBC < 0.000001)
            {
                // closest point is B
                if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P2)) return 0; // this edge is not considered
                alpha = 1.0;
                beta = 0.0;
            }
            else if (pBC > 0.999999)
            {
                // closest point is C
                if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P3)) return 0; // this edge is not considered
                alpha = 0.0;
                beta = 1.0;
            }
            else
            {
                // closest point is on BC
                if (!(flags&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_E31)) return 0; // this edge is not considered
                alpha = 1.0-pBC;
                beta = pBC;
            }
        }
    }

    defaulttype::Vector3 p, pq;
    p = e1.p1() + AB * alpha + AC * beta;
    pq = e2.center() - p;
    SReal norm2 = pq.norm2();
    if (pq.norm2() >= dist2)
        return 0;

    contacts->resize(contacts->size()+1);
    core::collision::DetectionOutput *detection = &*(contacts->end()-1);
    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    detection->id = e2.getIndex();
    detection->value = helper::rsqrt(norm2) ;

    if(detection->value>1e-15)
    {
        detection->normal = pq / detection->value;
    }
    else
    {
        msg_warning(intersection) <<"Null distance between contact detected";
        detection->normal= defaulttype::Vector3(1,0,0);
    }

    detection->value -= (intersection->getContactDistance() + e1.getProximity() + e2.getProximity() + e2.r());
    detection->point[0]=p;
    detection->point[1]= e2.getContactPointByNormal(detection->normal);

    return 1;
}

} //namespace sofa::component::collision
