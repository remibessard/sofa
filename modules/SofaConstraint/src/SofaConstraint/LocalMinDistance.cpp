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
#define SOFA_COMPONENT_COLLISION_LOCALMINDISTANCE_CPP
#include <SofaConstraint/LocalMinDistance.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/proximity.h>
#include <sofa/core/collision/Intersection.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/Node.h>

#include <SofaBaseCollision/CylinderModel.cpp>
#include <SofaBaseCollision/CapsuleModel.h>

#include <sofa/helper/types/RGBAColor.h>


#define DYNAMIC_CONE_ANGLE_COMPUTATION
#define EMIT_EXTRA_DEBUG_MESSAGE false

namespace sofa::core::collision
{
    template class SOFA_SOFACONSTRAINT_API IntersectorFactory<component::collision::LocalMinDistance>;

} // namespace sofa::core::collision

namespace sofa::component::collision
{

using namespace sofa::core::collision;
using namespace helper;
using namespace sofa::defaulttype;

using core::topology::BaseMeshTopology;


int LocalMinDistanceClass = core::RegisterObject("A set of methods to compute (for constraint methods) if two primitives are close enough to consider they collide")
        .add< LocalMinDistance >()
        ;

LocalMinDistance::LocalMinDistance()
    : BaseProximityIntersection()
    , filterIntersection(initData(&filterIntersection, true, "filterIntersection","Activate LMD filter"))
    , angleCone(initData(&angleCone, 0.0, "angleCone","Filtering cone extension angle"))
    , coneFactor(initData(&coneFactor, 0.5, "coneFactor", "Factor for filtering cone angle computation"))
    , useLMDFilters(initData(&useLMDFilters, false, "useLMDFilters", "Use external cone computation (Work in Progress)"))
{
    m_A = Vector3(0, 0, 0);
    m_B = Vector3(0, 0, 0);
    m_C = Vector3(0, 0, 0);
    m_T1 = Vector3(0, 0, 0);
    m_T2 = Vector3(0, 0, 0);
    m_T3 = Vector3(0, 0, 0);
    m_P = Vector3(0, 0, 0);
    m_H = Vector3(0, 0, 0);
    m_X = Vector3(0, 0, 0);
    m_N = Vector3(0, 0, 0);
    m_nb = 0;
}

void LocalMinDistance::init()
{
    intersectors.add<CubeCollisionModel, CubeCollisionModel, LocalMinDistance>(this);
    intersectors.add<SphereCollisionModel<sofa::defaulttype::Vec3Types>, SphereCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this); // sphere-sphere is always activated
    intersectors.add<SphereCollisionModel<sofa::defaulttype::Vec3Types>, PointCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this); // sphere-point is always activated

    intersectors.add<PointCollisionModel<sofa::defaulttype::Vec3Types>, PointCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this); // point-point is always activated
    intersectors.add<LineCollisionModel<sofa::defaulttype::Vec3Types>, LineCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this);
    intersectors.add<LineCollisionModel<sofa::defaulttype::Vec3Types>, PointCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this);
    intersectors.add<LineCollisionModel<sofa::defaulttype::Vec3Types>, SphereCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this);
    intersectors.add<TriangleCollisionModel<sofa::defaulttype::Vec3Types>, PointCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this);
    intersectors.add<TriangleCollisionModel<sofa::defaulttype::Vec3Types>, SphereCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this);

    intersectors.ignore<TriangleCollisionModel<sofa::defaulttype::Vec3Types>, LineCollisionModel<sofa::defaulttype::Vec3Types>>();			// never the case with LMD
    intersectors.ignore<TriangleCollisionModel<sofa::defaulttype::Vec3Types>, TriangleCollisionModel<sofa::defaulttype::Vec3Types>>();		// never the case with LMD

    intersectors.ignore<RayCollisionModel, PointCollisionModel<sofa::defaulttype::Vec3Types>>();
    intersectors.ignore<RayCollisionModel, LineCollisionModel<sofa::defaulttype::Vec3Types>>();
    intersectors.add<RayCollisionModel, TriangleCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this);
    intersectors.add<RayCollisionModel, SphereCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this);

    intersectors.add<CylinderCollisionModel<sofa::defaulttype::Vec3Types>, PointCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this);

    intersectors.add<CapsuleCollisionModel<sofa::defaulttype::Vec3Types>, PointCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this);
    intersectors.add<CapsuleCollisionModel<sofa::defaulttype::Vec3Types>, LineCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this);
    intersectors.add<CapsuleCollisionModel<sofa::defaulttype::Vec3Types>, TriangleCollisionModel<sofa::defaulttype::Vec3Types>, LocalMinDistance>(this);


    IntersectorFactory::getInstance()->addIntersectors(this);

    BaseProximityIntersection::init();
}

bool LocalMinDistance::testIntersection(Cube &cube1, Cube &cube2)
{
    const Vector3& minVect1 = cube1.minVect();
    const Vector3& minVect2 = cube2.minVect();
    const Vector3& maxVect1 = cube1.maxVect();
    const Vector3& maxVect2 = cube2.maxVect();

    const double alarmDist = getAlarmDistance() + cube1.getProximity() + cube2.getProximity();

    for (int i=0; i<3; i++)
    {
        if ( minVect1[i] > maxVect2[i] + alarmDist || minVect2[i]> maxVect1[i] + alarmDist )
            return false;
    }

    return true;
}

int LocalMinDistance::computeIntersection(Cube&, Cube&, OutputVector* /*contacts*/)
{
    return 0; /// \todo
}

bool LocalMinDistance::testIntersection(Line& e1, Line& e2)
{
    if(!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return false;

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();

    const Vector3 AB = e1.p2()-e1.p1();
    const Vector3 CD = e2.p2()-e2.p1();
    const Vector3 AC = e2.p1()-e1.p1();
    Matrix2 A;
    Vector2 b;

    A[0][0] = AB*AB;
    A[1][1] = CD*CD;
    A[0][1] = A[1][0] = -CD*AB;
    b[0] = AB*AC;
    b[1] = -CD*AC;

    const double det = defaulttype::determinant(A);

    double alpha = 0.5;
    double beta = 0.5;

    if (det < -1.0e-30 || det > 1.0e-30)
    {
        alpha = (b[0]*A[1][1] - b[1]*A[0][1])/det;
        beta  = (b[1]*A[0][0] - b[0]*A[1][0])/det;
        if (alpha < 1e-15 || alpha > (1.0-1e-15) ||
            beta  < 1e-15  || beta  > (1.0-1e-15) )
            return false;
    }

    Vector3 PQ = AC + CD * beta - AB * alpha;

    if (PQ.norm2() < alarmDist*alarmDist)
    {
        // filter for LMD

        if (!useLMDFilters.getValue())
        {
            if (!testValidity(e1, PQ))
                return false;

            Vector3 QP = -PQ;
            return testValidity(e2, QP);
        }
        else
        {
            return true;
        }
    }
    else
        return false;
}

int LocalMinDistance::computeIntersection(Line& e1, Line& e2, OutputVector* contacts)
{

    if(!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
    {
        dmsg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
            <<" not activated" ;
        return 0;
    }

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();

    // E1 => A-->B
    // E2 => C-->D
    const Vector3 AB = e1.p2()-e1.p1();
    const Vector3 CD = e2.p2()-e2.p1();
    const Vector3 AC = e2.p1()-e1.p1();
    Matrix2 A;
    Vector2 b;

    A[0][0] = AB*AB;
    A[1][1] = CD*CD;
    A[0][1] = A[1][0] = -CD*AB;
    b[0] = AB*AC;
    b[1] = -CD*AC;
    const double det = defaulttype::determinant(A);

    double alpha = 0.5;
    double beta = 0.5;

    if (det < -1.0e-30 || det > 1.0e-30)
    {
        alpha = (b[0]*A[1][1] - b[1]*A[0][1])/det;
        beta  = (b[1]*A[0][0] - b[0]*A[1][0])/det;
        if (alpha < 1e-15 || alpha > (1.0-1e-15) ||
            beta  < 1e-15  || beta  > (1.0-1e-15) )
            return 0;
    }
    else
    {
        // several possibilities :
        // -one point in common (auto-collision) => return false !
        // -no point in common but line are // => we can continue to test
        msg_warning() <<"Determinant is null";
    }

    Vector3 P,Q,PQ;
    P = e1.p1() + AB * alpha;
    Q = e2.p1() + CD * beta;
    PQ = Q-P;

    if (PQ.norm2() >= alarmDist*alarmDist)
        return 0;

    // filter for LMD //

    if (!useLMDFilters.getValue())
    {
        if (!testValidity(e1, PQ))
        {
            dmsg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                 <<" testValidity rejected for the first segment" ;
            return 0;
        }

        Vector3 QP = -PQ;

        if (!testValidity(e2, QP))
        {
            dmsg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                <<" testValidity rejected for the second segment";
            return 0;
        }
    }

    contacts->resize(contacts->size() + 1);
    DetectionOutput *detection = &*(contacts->end() - 1);

#ifdef SOFA_DETECTIONOUTPUT_FREEMOTION

    if (e1.hasFreePosition() && e2.hasFreePosition())
    {
        Vector3 Pfree, Qfree, ABfree, CDfree;
        ABfree = e1.p2Free()-e1.p1Free();
        CDfree = e2.p2Free()-e2.p1Free();
        Pfree = e1.p1Free() + ABfree * alpha;
        Qfree = e2.p1Free() + CDfree * beta;
        detection->freePoint[0] = Pfree;
        detection->freePoint[1] = Qfree;
    }

#endif

    const double contactDist = getContactDistance() + e1.getProximity() + e2.getProximity();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    detection->id = (e1.getCollisionModel()->getSize() > e2.getCollisionModel()->getSize()) ? e1.getIndex() : e2.getIndex();
    detection->point[0] = P;
    detection->point[1] = Q;
    detection->normal = PQ;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= contactDist;

    return 1;
}

bool LocalMinDistance::testIntersection(Triangle& e2, Point& e1)
{
    if(!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return false;

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();

    const Vector3 AB = e2.p2()-e2.p1();
    const Vector3 AC = e2.p3()-e2.p1();
    const Vector3 AP = e1.p() -e2.p1();
    Matrix2 A;
    Vector2 b;

    // We want to find alpha,beta so that:
    // AQ = AB*alpha+AC*beta
    // PQ.AB = 0 and PQ.AC = 0
    // (AQ-AP).AB = 0 and (AQ-AP).AC = 0
    // AQ.AB = AP.AB and AQ.AC = AP.AC
    //
    // (AB*alpha+AC*beta).AB = AP.AB and
    // (AB*alpha+AC*beta).AC = AP.AC
    //
    // AB.AB*alpha + AC.AB*beta = AP.AB and
    // AB.AC*alpha + AC.AC*beta = AP.AC
    //
    // A . [alpha beta] = b
    A[0][0] = AB*AB;
    A[1][1] = AC*AC;
    A[0][1] = A[1][0] = AB*AC;
    b[0] = AP*AB;
    b[1] = AP*AC;
    const double det = defaulttype::determinant(A);

    double alpha = 0.5;
    double beta = 0.5;


    alpha = (b[0]*A[1][1] - b[1]*A[0][1])/det;
    beta  = (b[1]*A[0][0] - b[0]*A[1][0])/det;
    if (alpha < 0.000001 ||
            beta  < 0.000001 ||
            alpha + beta  > 0.999999)
        return false;

    const Vector3 PQ = AB * alpha + AC * beta - AP;

    if (PQ.norm2() < alarmDist*alarmDist)
    {
        //filter for LMD
        if (!useLMDFilters.getValue())
        {
            if (!testValidity(e1, PQ))
                return false;

            Vector3 QP = -PQ;
            return testValidity(e2, QP);
        }
        else
        {
            return true;
        }
        // end filter
    }
    else
        return false;
}

int LocalMinDistance::doIntersectionTrianglePoint(SReal alarmDist, Triangle& e2, core::CollisionElementIterator& e1, const Vector3 p, OutputVector* contacts)
{
    const Vector3 AB = e2.p2() - e2.p1();
    const Vector3 AC = e2.p3() - e2.p1();
    const Vector3 AP = p - e2.p1();
    Matrix2 A;
    Vector2 b;

    A[0][0] = AB * AB;
    A[1][1] = AC * AC;
    A[0][1] = A[1][0] = AB * AC;
    b[0] = AP * AB;
    b[1] = AP * AC;

    const double det = defaulttype::determinant(A);

    double alpha = 0.5;
    double beta = 0.5;


    alpha = (b[0] * A[1][1] - b[1] * A[0][1]) / det;
    beta = (b[1] * A[0][0] - b[0] * A[1][0]) / det;
    if (alpha < 0.000001 ||
        beta  < 0.000001 ||
        alpha + beta  > 0.999999)
        return 0;


    Vector3 P, Q; //PQ
    P = p;
    Q = e2.p1() + AB * alpha + AC * beta;
    Vector3 PQ = Q - P;
    Vector3 QP = -PQ;

    if (PQ.norm2() >= alarmDist * alarmDist)
        return 0;


    // filter for LMD

    if (!useLMDFilters.getValue())
    {
        Point* point = dynamic_cast<TPoint<sofa::defaulttype::Vec3Types>*>(e1.getCollisionModel());
        if (point != NULL)
        {
            if (!testValidity(*point, PQ))
                return 0;
        }

        if (!testValidity(e2, QP))
            return 0;
    }

    //end filter

    contacts->resize(contacts->size() + 1);
    DetectionOutput *detection = &*(contacts->end() - 1);

    const double contactDist = getContactDistance() + e1.getProximity() + e2.getProximity();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e2, e1); //Tri , Point
    Capsule* cap = dynamic_cast<Capsule*>(e1.getCollisionModel());
    if(cap)
        std::cout << "e1.getIndex: " << e1.getIndex() << "  and pos:" << p << std::endl;
    detection->id = e1.getIndex(); // Point
    detection->point[0] = Q; // Bary in Tri
    detection->point[1] = P; //Point
    detection->normal = QP;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= contactDist;
    return 1;
}

int LocalMinDistance::computeIntersection(Triangle& e2, Point& e1, OutputVector* contacts)
{
    if(!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();

    doIntersectionTrianglePoint(alarmDist, e2, static_cast<core::CollisionElementIterator>(e1), e1.p(), contacts);
}


bool LocalMinDistance::testIntersection(Triangle& e2, Sphere& e1)
{
    if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return false;

    const double alarmDist = getAlarmDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    const Vector3 AB = e2.p2()-e2.p1();
    const Vector3 AC = e2.p3()-e2.p1();
    const Vector3 AP = e1.p() -e2.p1();
    Matrix2 A;
    Vector2 b;

    // We want to find alpha,beta so that:
    // AQ = AB*alpha+AC*beta
    // PQ.AB = 0 and PQ.AC = 0
    // (AQ-AP).AB = 0 and (AQ-AP).AC = 0
    // AQ.AB = AP.AB and AQ.AC = AP.AC
    //
    // (AB*alpha+AC*beta).AB = AP.AB and
    // (AB*alpha+AC*beta).AC = AP.AC
    //
    // AB.AB*alpha + AC.AB*beta = AP.AB and
    // AB.AC*alpha + AC.AC*beta = AP.AC
    //
    // A . [alpha beta] = b
    A[0][0] = AB*AB;
    A[1][1] = AC*AC;
    A[0][1] = A[1][0] = AB*AC;
    b[0] = AP*AB;
    b[1] = AP*AC;
    const double det = defaulttype::determinant(A);

    double alpha = 0.5;
    double beta = 0.5;


    alpha = (b[0]*A[1][1] - b[1]*A[0][1])/det;
    beta  = (b[1]*A[0][0] - b[0]*A[1][0])/det;
    if (alpha < 0.000001 ||
            beta  < 0.000001 ||
            alpha + beta  > 0.999999)
        return false;


    const Vector3 PQ = AB * alpha + AC * beta - AP;

    if (PQ.norm2() < alarmDist*alarmDist)
    {

        //filter for LMD

        if (!useLMDFilters.getValue())
        {
            if (!testValidity(e1, PQ))
                return false;

            Vector3 QP = -PQ;
            return testValidity(e2, QP);
        }
        else
        {
            return true;
        }

        // end filter
    }
    else
        return false;
}

int LocalMinDistance::computeIntersection(Triangle& e2, Sphere& e1, OutputVector* contacts)
{
	if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
		return false;

    const double alarmDist = getAlarmDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    const Vector3 AB = e2.p2()-e2.p1();
    const Vector3 AC = e2.p3()-e2.p1();
    const Vector3 AP = e1.p() -e2.p1();
    Matrix2 A;
    Vector2 b;

    A[0][0] = AB*AB;
    A[1][1] = AC*AC;
    A[0][1] = A[1][0] = AB*AC;
    b[0] = AP*AB;
    b[1] = AP*AC;



    const double det = defaulttype::determinant(A);

    double alpha = 0.5;
    double beta = 0.5;

    alpha = (b[0]*A[1][1] - b[1]*A[0][1])/det;
    beta  = (b[1]*A[0][0] - b[0]*A[1][0])/det;
    if (alpha < 0.000001 ||
            beta  < 0.000001 ||
            alpha + beta  > 0.999999)
        return 0;

    Vector3 P,Q; //PQ
    P = e1.p();
    Q = e2.p1() + AB * alpha + AC * beta;
    Vector3 PQ = Q-P;
    Vector3 QP = -PQ;

    if (PQ.norm2() >= alarmDist*alarmDist)
        return 0;


    // filter for LMD

    if (!useLMDFilters.getValue())
    {
        if (!testValidity(e1, PQ))
            return 0;

        if (!testValidity(e2, QP))
            return 0;
    }

    //end filter

    contacts->resize(contacts->size()+1);
    DetectionOutput *detection = &*(contacts->end()-1);

#ifdef SOFA_DETECTIONOUTPUT_FREEMOTION
    if (e1.hasFreePosition() && e2.hasFreePosition())
    {
        Vector3 Pfree,Qfree,ABfree,ACfree;
        ABfree = e2.p2Free()-e2.p1Free();
        ACfree = e2.p3Free()-e2.p1Free();
        Pfree = e1.pFree();
        Qfree = e2.p1Free() + ABfree * alpha + ACfree * beta;

        detection->freePoint[0] = Qfree;
        detection->freePoint[1] = Pfree;
    }
#endif

    const double contactDist = getContactDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e2, e1);
    detection->id = e1.getIndex();
    detection->point[0] = Q;
    detection->point[1] = P;
    detection->normal = QP;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= contactDist;
    return 1;
}

bool LocalMinDistance::testIntersection(Line& e2, Point& e1)
{

    if(!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return false;

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();
    const Vector3 AB = e2.p2()-e2.p1();
    const Vector3 AP = e1.p()-e2.p1();

    double A;
    double b;
    A = AB*AB;
    b = AP*AB;

    double alpha = 0.5;


    alpha = b/A;
    if (alpha < 0.000001 || alpha > 0.999999)
        return false;


    Vector3 P,Q,PQ;
    P = e1.p();
    Q = e2.p1() + AB * alpha;
    PQ = Q-P;

    if (PQ.norm2() < alarmDist*alarmDist)
    {
        // filter for LMD

        if (!useLMDFilters.getValue())
        {
            if (!testValidity(e1, PQ))
                return false;

            Vector3 QP = -PQ;
            return testValidity(e2, QP);
        }
        else
        {
            return true;
        }

        // end filter
    }
    else
        return false;
}

int LocalMinDistance::computeIntersection(Line& e2, Point& e1, OutputVector* contacts)
{
    if(!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();
    const Vector3 AB = e2.p2()-e2.p1();
    const Vector3 AP = e1.p()-e2.p1();

    if (AB.norm()<0.000000000001*AP.norm())
    {
        return 0;
    }


    double A;
    double b;
    A = AB*AB;
    b = AP*AB;

    double alpha = 0.5;

    alpha = b/A;
    if (alpha < 0.000001 || alpha > 0.999999)
        return 0;

    Vector3 P,Q;
    P = e1.p();
    Q = e2.p1() + AB * alpha;
    Vector3 PQ = Q - P;
    Vector3 QP = -PQ;

    if (PQ.norm2() >= alarmDist*alarmDist)
        return 0;

    // filter for LMD
    if (!useLMDFilters.getValue())
    {
        if (!testValidity(e1, PQ))
            return 0;

        if (!testValidity(e2, QP))
            return 0;
    }

    // end filter

    contacts->resize(contacts->size()+1);
    DetectionOutput *detection = &*(contacts->end()-1);

#ifdef SOFA_DETECTIONOUTPUT_FREEMOTION
    if (e1.hasFreePosition() && e2.hasFreePosition())
    {
        Vector3 ABfree = e2.p2Free() - e2.p1Free();
        Vector3 Pfree = e1.pFree();
        Vector3 Qfree = e2.p1Free() + ABfree * alpha;

        detection->freePoint[0] = Qfree;
        detection->freePoint[1] = Pfree;
    }
#endif

    const double contactDist = getContactDistance() + e1.getProximity() + e2.getProximity();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e2, e1);
    detection->id = e1.getIndex();
    detection->point[0]=Q;
    detection->point[1]=P;
    detection->normal=QP;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= contactDist;

    return 1;
}


bool LocalMinDistance::testIntersection(Line& e2, Sphere& e1)
{
	if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
		return false;

    const double alarmDist = getAlarmDistance() + e1.r() + e1.getProximity() + e2.getProximity();
    const Vector3 AB = e2.p2()-e2.p1();
    const Vector3 AP = e1.p()-e2.p1();

    double A;
    double b;
    A = AB*AB;
    b = AP*AB;

    double alpha = 0.5;


    alpha = b/A;
    if (alpha < 0.000001 || alpha > 0.999999)
        return false;


    Vector3 P,Q,PQ;
    P = e1.p();
    Q = e2.p1() + AB * alpha;
    PQ = Q-P;

    if (PQ.norm2() < alarmDist*alarmDist)
    {
        // filter for LMD

        if (!useLMDFilters.getValue())
        {
            if (!testValidity(e1, PQ))
                return false;

            Vector3 QP = -PQ;
            return testValidity(e2, QP);
        }
        else
        {
            return true;
        }

        // end filter
    }
    else
        return false;
}

int LocalMinDistance::computeIntersection(Line& e2, Sphere& e1, OutputVector* contacts)
{
	if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
		return false;

    const double alarmDist = getAlarmDistance() + e1.r() + e1.getProximity() + e2.getProximity();
    const Vector3 AB = e2.p2()-e2.p1();
    const Vector3 AP = e1.p()-e2.p1();

    if (AB.norm()<0.000000000001*AP.norm())
    {
        return 0;
    }


    double A;
    double b;
    A = AB*AB;
    b = AP*AB;

    double alpha = 0.5;


    alpha = b/A;
    if (alpha < 0.000001 || alpha > 0.999999)
        return 0;


    Vector3 P,Q;
    P = e1.p();
    Q = e2.p1() + AB * alpha;
    Vector3 PQ = Q - P;
    Vector3 QP = -PQ;

    if (PQ.norm2() >= alarmDist*alarmDist)
        return 0;

    // filter for LMD
    if (!useLMDFilters.getValue())
    {
        if (!testValidity(e1, PQ))
            return 0;

        if (!testValidity(e2, QP))
            return 0;
    }


    // end filter

    contacts->resize(contacts->size()+1);
    DetectionOutput *detection = &*(contacts->end()-1);

#ifdef SOFA_DETECTIONOUTPUT_FREEMOTION
    if (e1.hasFreePosition() && e2.hasFreePosition())
    {
        Vector3 ABfree = e2.p2Free() - e2.p1Free();
        Vector3 Pfree = e1.pFree();
        Vector3 Qfree = e2.p1Free() + ABfree * alpha;

        detection->freePoint[0] = Qfree;
        detection->freePoint[1] = Pfree;
    }
#endif

    const double contactDist = getContactDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e2, e1);
    detection->id = e1.getIndex();
    detection->point[0]=Q;
    detection->point[1]=P;
    detection->normal=QP;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= contactDist;

    return 1;
}

bool LocalMinDistance::testIntersection(Point& e1, Point& e2)
{
    if(!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();

    Vector3 PQ = e2.p()-e1.p();

    if (PQ.norm2() < alarmDist*alarmDist)
    {
        // filter for LMD

        if (!useLMDFilters.getValue())
        {
            if (!testValidity(e1, PQ))
                return false;

            Vector3 QP = -PQ;
            return testValidity(e2, QP);
        }
        else
        {
            return true;
        }

        // end filter
    }
    else
        return false;
}

int LocalMinDistance::computeIntersection(Point& e1, Point& e2, OutputVector* contacts)
{
    if(!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();

    Vector3 P,Q,PQ;
    P = e1.p();
    Q = e2.p();
    PQ = Q-P;


    if (PQ.norm2() >= alarmDist*alarmDist)
        return 0;

    // filter for LMD

    if (!useLMDFilters.getValue())
    {
        if (!testValidity(e1, PQ))
            return 0;

        Vector3 QP = -PQ;

        if (!testValidity(e2, QP))
            return 0;
    }

    // end filter

    contacts->resize(contacts->size()+1);
    DetectionOutput *detection = &*(contacts->end()-1);

#ifdef SOFA_DETECTIONOUTPUT_FREEMOTION
    if (e1.hasFreePosition() && e2.hasFreePosition())
    {
        Vector3 Pfree,Qfree;
        Pfree = e1.pFree();
        Qfree = e2.pFree();

        detection->freePoint[0] = Pfree;
        detection->freePoint[1] = Qfree;
    }
#endif

    const double contactDist = getContactDistance() + e1.getProximity() + e2.getProximity();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    detection->id = (e1.getCollisionModel()->getSize() > e2.getCollisionModel()->getSize()) ? e1.getIndex() : e2.getIndex();
    detection->point[0]=P;
    detection->point[1]=Q;
    detection->normal=PQ;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= contactDist;
    return 1;
}

bool LocalMinDistance::testIntersection(Sphere& e1, Point& e2)
{
    if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;

    const double alarmDist = getAlarmDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    Vector3 PQ = e2.p()-e1.p();

    if (PQ.norm2() < alarmDist*alarmDist)
    {
        // filter for LMD

        if (!useLMDFilters.getValue())
        {
            if (!testValidity(e1, PQ))
                return false;

            Vector3 QP = -PQ;
            return testValidity(e2, QP);
        }
        else
        {
            return true;
        }

        // end filter
    }
    else
        return false;
}

int LocalMinDistance::computeIntersection(Sphere& e1, Point& e2, OutputVector* contacts)
{
    if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;

    const double alarmDist = getAlarmDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    Vector3 P,Q,PQ;
    P = e1.p();
    Q = e2.p();
    PQ = Q-P;


    if (PQ.norm2() >= alarmDist*alarmDist)
        return 0;

    // filter for LMD

    if (!useLMDFilters.getValue())
    {
        if (!testValidity(e1, PQ))
            return 0;

        Vector3 QP = -PQ;

        if (!testValidity(e2, QP))
            return 0;
    }

    // end filter

    contacts->resize(contacts->size()+1);
    DetectionOutput *detection = &*(contacts->end()-1);

#ifdef SOFA_DETECTIONOUTPUT_FREEMOTION
    if (e1.hasFreePosition() && e2.hasFreePosition())
    {
        Vector3 Pfree,Qfree;
        Pfree = e1.pFree();
        Qfree = e2.pFree();

        detection->freePoint[0] = Pfree;
        detection->freePoint[1] = Qfree;
    }
#endif

    const double contactDist = getContactDistance() + e1.r() + e1.getProximity() + e2.getProximity();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    detection->id = (e1.getCollisionModel()->getSize() > e2.getCollisionModel()->getSize()) ? e1.getIndex() : e2.getIndex();
    detection->point[0]=P;
    detection->point[1]=Q;
    detection->normal=PQ;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= contactDist;
    return 1;
}

bool LocalMinDistance::testIntersection(Sphere& e1, Sphere& e2)
{
    if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;

    const double alarmDist = getAlarmDistance() + e1.r() + e1.getProximity() + e2.r() + e2.getProximity();

    Vector3 PQ = e2.p()-e1.p();

    if (PQ.norm2() < alarmDist*alarmDist)
    {
        // filter for LMD

        if (!useLMDFilters.getValue())
        {
            if (!testValidity(e1, PQ))
                return false;

            Vector3 QP = -PQ;
            return testValidity(e2, QP);
        }
        else
        {
            return true;
        }

        // end filter
    }
    else
        return false;
}

int LocalMinDistance::computeIntersection(Sphere& e1, Sphere& e2, OutputVector* contacts)
{
    if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;

    const double alarmDist = getAlarmDistance() + e1.r() + e1.getProximity() + e2.r() + e2.getProximity();

    Vector3 P,Q,PQ;
    P = e1.p();
    Q = e2.p();
    PQ = Q-P;


    if (PQ.norm2() >= alarmDist*alarmDist)
        return 0;

    // filter for LMD

    if (!useLMDFilters.getValue())
    {
        if (!testValidity(e1, PQ))
            return 0;

        Vector3 QP = -PQ;

        if (!testValidity(e2, QP))
            return 0;
    }

    // end filter

    contacts->resize(contacts->size()+1);
    DetectionOutput *detection = &*(contacts->end()-1);

#ifdef SOFA_DETECTIONOUTPUT_FREEMOTION
    if (e1.hasFreePosition() && e2.hasFreePosition())
    {
        Vector3 Pfree,Qfree;
        Pfree = e1.pFree();
        Qfree = e2.pFree();

        detection->freePoint[0] = Pfree;
        detection->freePoint[1] = Qfree;
    }
#endif

    const double contactDist = getContactDistance() + e1.r() + e1.getProximity() + e2.r() + e2.getProximity();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    detection->id = (e1.getCollisionModel()->getSize() > e2.getCollisionModel()->getSize()) ? e1.getIndex() : e2.getIndex();
    detection->point[0]=P;
    detection->point[1]=Q;
    detection->normal=PQ;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= contactDist;
    return 1;
}


bool LocalMinDistance::testIntersection(Ray &t1,Triangle &t2)
{
    Vector3 P,Q,PQ;
    static DistanceSegTri proximitySolver;

    const double alarmDist = getAlarmDistance() + t1.getProximity() + t2.getProximity();

    if (fabs(t2.n() * t1.direction()) < 0.000001)
        return false; // no intersection for edges parallel to the triangle

    Vector3 A = t1.origin();
    Vector3 B = A + t1.direction() * t1.l();

    proximitySolver.NewComputation( t2.p1(), t2.p2(), t2.p3(), A, B,P,Q);
    PQ = Q-P;

    if (PQ.norm2() < alarmDist*alarmDist)
    {
        return true;
    }
    else
        return false;
}

int LocalMinDistance::computeIntersection(Ray &t1, Triangle &t2, OutputVector* contacts)
{
    if (!t1.isActive(t2.getCollisionModel()) || !t2.isActive(t1.getCollisionModel()))
        return 0;

    const double alarmDist = getAlarmDistance() + t1.getProximity() + t2.getProximity();


    if (fabs(t2.n() * t1.direction()) < 0.000001)
        return false; // no intersection for edges parallel to the triangle

    Vector3 A = t1.origin();
    Vector3 B = A + t1.direction() * t1.l();

    Vector3 P,Q,PQ;
    static DistanceSegTri proximitySolver;

    proximitySolver.NewComputation( t2.p1(), t2.p2(), t2.p3(), A,B,P,Q);
    PQ = Q-P;

    if (PQ.norm2() >= alarmDist*alarmDist)
        return 0;

    const double contactDist = alarmDist;
    contacts->resize(contacts->size()+1);
    DetectionOutput *detection = &*(contacts->end()-1);

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(t1, t2);
    detection->id = t1.getIndex();
    detection->point[1]=P;
    detection->point[0]=Q;
#ifdef SOFA_DETECTIONOUTPUT_FREEMOTION
    detection->freePoint[1] = P;
    detection->freePoint[0] = Q;
#endif
    detection->normal=-t2.n();
    detection->value = PQ.norm();
    detection->value -= contactDist;
    return 1;
}


bool LocalMinDistance::testIntersection(Ray &ray1,Sphere &sph2)
{
    // Center of the sphere
    const Vector3 sph2Pos(sph2.center());
    // Radius of the sphere
    const double radius1 = sph2.r();

    const Vector3 ray1Origin(ray1.origin());
    const Vector3 ray1Direction(ray1.direction());
    const double length2 = ray1.l();
    const Vector3 tmp = sph2Pos - ray1Origin;
    const double rayPos = tmp*ray1Direction;
    const double rayPosInside = std::max(std::min(rayPos,length2),0.0);
    const double dist2 = tmp.norm2() - (rayPosInside*rayPosInside);
    return (dist2 < (radius1*radius1));
}

int LocalMinDistance::computeIntersection(Ray &ray1, Sphere &sph2, OutputVector* contacts)
{
    if (!ray1.isActive(sph2.getCollisionModel()) || !sph2.isActive(ray1.getCollisionModel()))
        return 0;

    // Center of the sphere
    const Vector3 sph2Pos(sph2.center());
    // Radius of the sphere
    const double radius1 = sph2.r();

    const Vector3 ray1Origin(ray1.origin());
    const Vector3 ray1Direction(ray1.direction());
    const double length2 = ray1.l();
    const Vector3 tmp = sph2Pos - ray1Origin;
    const double rayPos = tmp*ray1Direction;
    const double rayPosInside = std::max(std::min(rayPos,length2),0.0);
    const double dist2 = tmp.norm2() - (rayPosInside*rayPosInside);
    if (dist2 >= (radius1*radius1))
        return 0;

    const double dist = sqrt(dist2);

    contacts->resize(contacts->size()+1);
    DetectionOutput *detection = &*(contacts->end()-1);

    detection->point[0] = ray1Origin + ray1Direction*rayPosInside;
    detection->normal = sph2Pos - detection->point[0];
    detection->normal /= dist;
    detection->point[1] = sph2Pos - detection->normal * radius1;
    detection->value = dist - radius1;
    detection->elem.first = ray1;
    detection->elem.second = sph2;
    detection->id = ray1.getIndex();
    return 1;
}

bool LocalMinDistance::testIntersection(Cylinder& e1, Point& e2)
{
    if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return false;

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();

    const Vector3 AB = e1.point2() - e1.point1();
    const Vector3 AP = e2.p() - e1.point1();

    std::cout << "on test Cylinder point !" << std::endl;
    return (cross(AB, AP).norm() / AB.norm()) <= e1.radius();
}

int LocalMinDistance::computeIntersection(Cylinder& e1, Point& e2, OutputVector* contacts)
{
    if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();

    //                  P   
    // _________________|______
    // |                |X    |
    // |                |     |
    // |                |     |
    // |__________'_____'_____|
    // A          C     H     B

    const Vector3 AB = e1.point2() - e1.point1();
    const Vector3 PC = e1.center() - e2.p();
    const double HCnorm = dot(PC,AB.normalized());
    const Vector3 H = e1.center() - HCnorm * e1.axis().normalized();
    const double PCnorm2 = PC.norm2();
    const Vector3 HP = e2.p() - H;
    const double HPnorm2 = PCnorm2 - HCnorm*HCnorm;

    Vector3 X;
    if (fabs(HCnorm) < e1.height() / 2) // H is on the line segment AB
    {
        X = H + e1.radius() * ((e2.p() - H).normalized());
    }
    else // we need to compute distance to disc
    {
        if (HPnorm2 < e1.radius()*e1.radius()) // P projects onto the disc
        {
            if (HCnorm > 0)
            {
                X = e1.point1() + HP;
            }
            else
            {
                X = e1.point2() + HP;
            }
        }
        else // P projects onto a circle 
        {
            if (HCnorm > 0)
            {
                X = e1.point1() + e1.radius() * HP.normalized();
            }
            else
            {
                X = e1.point2() + e1.radius() * HP.normalized();
            }
        }
    }
    m_A = e1.point1();
    m_B = e1.point2();
    m_C = e1.center();
    m_P = e2.p();
    m_H = H;
    m_X = X;

    if ((e2.p()-X).norm2() >= alarmDist * alarmDist)
        return 0;


    //// filter for LMD

    //if (!useLMDFilters.getValue())
    //{
    //    if (!testValidity(e1, PQ))
    //        return 0;

    //    if (!testValidity(e2, QP))
    //        return 0;
    //}

    ////end filter

    contacts->resize(contacts->size() + 1);
    DetectionOutput *detection = &*(contacts->end() - 1);

    const double contactDist = getContactDistance() + e1.getProximity() + e2.getProximity();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    detection->id = e2.getIndex();
    detection->point[0] = X;
    detection->point[1] = e2.p();
    detection->normal = e2.p() - X;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= contactDist;
    return 1;
}

bool LocalMinDistance::testIntersection(Capsule& e1, Point& e2)
{
    std::cout << "on test Capsule point !" << std::endl;

    const defaulttype::Vector3 p1 = e1.point1();
    const defaulttype::Vector3 p2 = e1.point2();
    const defaulttype::Vector3 q = e2.p();
    const defaulttype::Vector3 AB = p2 - p1;
    const defaulttype::Vector3 AQ = q - p1;
    const SReal cap_rad = e1.radius();
    SReal A;
    SReal b;
    A = AB * AB;
    b = AQ * AB;

    SReal alpha = 0.5;

    alpha = b / A;//projection of the point on the capsule segment such as the projected point P = p1 + AB * alpha
    if (alpha < 0.0) alpha = 0.0;//if the projection is out the segment, we associate it to a segment apex
    else if (alpha > 1.0) alpha = 1.0;

    defaulttype::Vector3 p, pq;
    p = p1 + AB * alpha;
    pq = q - p;

    SReal enough_to_touch = getAlarmDistance() + cap_rad;
    if (pq.norm2() >= enough_to_touch * enough_to_touch)
        return false;
    return true;
}

int LocalMinDistance::computeIntersection(Capsule& e1, Point& e2, OutputVector* contacts)
{
    if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();

    const defaulttype::Vector3 p1 = e1.point1();
    const defaulttype::Vector3 p2 = e1.point2();
    const defaulttype::Vector3 q  = e2.p();
    const defaulttype::Vector3 AB = p2 - p1;
    const defaulttype::Vector3 AQ = q - p1;
    SReal A;
    SReal b;
    A = AB * AB;
    b = AQ * AB;
    SReal cap_rad = e1.radius();

    SReal alpha = 0.5;

    alpha = b / A;//projection of the point on the capsule segment such as the projected point P = p1 + AB * alpha
    if (alpha < 0.0) alpha = 0.0;//if the projection is out the segment, we associate it to a segment apex
    else if (alpha > 1.0) alpha = 1.0;

    defaulttype::Vector3 p, pq;
    p = p1 + AB * alpha;
    pq = q - p;

    m_A = p1;
    m_B = p2;
    m_C = p1+(p2-p1)/2;
    m_P = q;
    //m_H = H;
    m_X = p + cap_rad * pq.normalized();

    SReal enough_to_touch = alarmDist + cap_rad;
    if (pq.norm2() >= enough_to_touch * enough_to_touch)
        return 0;

    //const SReal contactDist = getContactDistance() + e1.getProximity() + e2.getProximity();
    contacts->resize(contacts->size() + 1);
    DetectionOutput *detection = &*(contacts->end() - 1);
    
    const double contactDist = getContactDistance() + e1.getProximity() + e2.getProximity();
    
    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    detection->id = e2.getIndex();
    detection->point[0] = p;
    detection->point[1] = q;
    detection->normal = pq;

    detection->value = detection->normal.norm();
    detection->normal /= detection->value;

    detection->value -= (contactDist + cap_rad);

    return 1;
}

bool LocalMinDistance::testIntersection(Capsule& e1, Line& e2)
{
    std::cout << "on test Capsule Line !" << std::endl;
    return true;
}

int LocalMinDistance::doIntersectionCapsuleLine(Capsule& e1, core::CollisionElementIterator& e2, const defaulttype::Vector3 l1, const defaulttype::Vector3 l2, OutputVector* contacts)
{
    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();

    const defaulttype::Vector3 p1 = e1.point1();
    const defaulttype::Vector3 p2 = e1.point2();
    const defaulttype::Vector3 q1 = l1;
    const defaulttype::Vector3 q2 = l2;
    const Vector3 AB = p2 - p1; // capsule segment
    const Vector3 CD = q2 - q1; // line segment
    const Vector3 AC = q1 - p1;
    Matrix2 A;
    Vector2 b;
    A[0][0] = AB * AB;
    A[1][1] = CD * CD;
    A[0][1] = A[1][0] = -CD * AB;
    b[0] = AB * AC;
    b[1] = -CD * AC;
    const SReal det = defaulttype::determinant(A);
    const SReal cap_rad = e1.radius();

    SReal alpha = 0.5;
    SReal beta = 0.5;

    if (det < -0.000000000001 || det > 0.000000000001)//AB and CD are not on the same plane
    {
        alpha = (b[0] * A[1][1] - b[1] * A[0][1]) / det;
        beta = (b[1] * A[0][0] - b[0] * A[1][0]) / det;

        if (alpha < 0)
            alpha = 0;
        else if (alpha > 1)
            alpha = 1;

        if (beta < 0)
            beta = 0;
        else if (beta > 1)
            beta = 1;
    }
    else {//Segments on a same plane. Here the idea to find the nearest points
        //is to project segment apexes on the other segment.
        //Visual example with semgents AB and CD :
        //            A----------------B
        //                     C----------------D
        //After projection :
        //            A--------c-------B
        //                     C-------b--------D
        //So the nearest points are p and q which are respecively in the middle of cB and Cb:
        //            A--------c---p---B
        //                     C---q---b--------D

        Vector3 AD = q2 - p1;
        Vector3 CB = p2 - q1;

        SReal AB_norm2 = AB.norm2();
        SReal CD_norm2 = CD.norm2();
        SReal c_proj = b[0] / AB_norm2;//alpha = (AB * AC)/AB_norm2
        SReal d_proj = (AB * AD) / AB_norm2;
        SReal a_proj = b[1];//beta = (-CD*AC)/CD_norm2
        SReal b_proj = (CD*CB) / CD_norm2;

        if (c_proj >= 0 && c_proj <= 1) {//projection of C on AB is lying on AB
            if (d_proj > 1) {//case :
                           //             A----------------B
                           //                      C---------------D
                alpha = (1.0 + c_proj) / 2.0;
                beta = b_proj / 2.0;
            }
            else if (d_proj < 0) {//case :
                                //             A----------------B
                                //     D----------------C
                alpha = c_proj / 2.0;
                beta = (1 + a_proj) / 2.0;
            }
            else {//case :
                //             A----------------B
                //                 C------D
                alpha = (c_proj + d_proj) / 2.0;
                beta = 0.5;
            }
        }
        else if (d_proj >= 0 && d_proj <= 1) {
            if (c_proj < 0) {//case :
                           //             A----------------B
                           //     C----------------D
                alpha = d_proj / 2.0;
                beta = (1 + a_proj) / 2.0;
            }
            else {//case :
                 //          A---------------B
                 //                 D-------------C
                alpha = (1 + d_proj) / 2.0;
                beta = b_proj / 2.0;
            }
        }
        else {
            if (c_proj * d_proj < 0) {//case :
                                    //           A--------B
                                    //       D-----------------C
                alpha = 0.5;
                beta = (a_proj + b_proj) / 2.0;
            }
            else {
                if (c_proj < 0) {//case :
                               //                    A---------------B
                               // C-------------D

                    alpha = 0;
                }
                else {
                    alpha = 1;
                }

                if (a_proj < 0) {//case :
                               // A---------------B
                               //                     C-------------D
                }
                else {//case :
                     //                     A---------------B
                     //   C-------------D
                    beta = 1;
                }
            }
        }
    }

    bool ignore_p1 = true;
    bool ignore_p2 = true;
    if (ignore_p1 && beta == 0)
        return 0;
    if (ignore_p2 && beta == 1)
        return 0;

    SReal enough_to_touch = alarmDist + cap_rad;
    Vector3 p, q, pq;
    p = p1 + AB * alpha;
    q = q1 + CD * beta;
    pq = q - p;
    if (pq.norm2() >= enough_to_touch * enough_to_touch)
        return 0;

    contacts->resize(contacts->size() + 1);
    DetectionOutput *detection = &*(contacts->end() - 1);

    const double contactDist = getContactDistance() + e1.getProximity() + e2.getProximity();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    detection->id = e2.getIndex();
    detection->point[0] = p;
    detection->point[1] = q;
    detection->normal = pq;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= (contactDist + cap_rad);
    return 1;
}

int LocalMinDistance::computeIntersection(Capsule& e1, Line& e2, OutputVector* contacts)
{
    if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();

    const defaulttype::Vector3 p1 = e1.point1();
    const defaulttype::Vector3 p2 = e1.point2();
    const defaulttype::Vector3 q1 = e2.p1();
    const defaulttype::Vector3 q2 = e2.p2();
    const Vector3 AB = p2 - p1; // capsule segment
    const Vector3 CD = q2 - q1; // line segment
    const Vector3 AC = q1 - p1;
    Matrix2 A;
    Vector2 b;
    A[0][0] = AB * AB;
    A[1][1] = CD * CD;
    A[0][1] = A[1][0] = -CD * AB;
    b[0] = AB * AC;
    b[1] = -CD * AC;
    const SReal det = defaulttype::determinant(A);
    const SReal cap_rad = e1.radius();

    SReal alpha = 0.5;
    SReal beta = 0.5;

    if (det < -0.000000000001 || det > 0.000000000001)//AB and CD are not on the same plane
    {
        alpha = (b[0] * A[1][1] - b[1] * A[0][1]) / det;
        beta = (b[1] * A[0][0] - b[0] * A[1][0]) / det;

        if (alpha < 0)
            alpha = 0;
        else if (alpha > 1)
            alpha = 1;

        if (beta < 0)
            beta = 0;
        else if (beta > 1)
            beta = 1;
    }
    else {//Segments on a same plane. Here the idea to find the nearest points
        //is to project segment apexes on the other segment.
        //Visual example with semgents AB and CD :
        //            A----------------B
        //                     C----------------D
        //After projection :
        //            A--------c-------B
        //                     C-------b--------D
        //So the nearest points are p and q which are respecively in the middle of cB and Cb:
        //            A--------c---p---B
        //                     C---q---b--------D

        Vector3 AD = q2 - p1;
        Vector3 CB = p2 - q1;

        SReal AB_norm2 = AB.norm2();
        SReal CD_norm2 = CD.norm2();
        SReal c_proj = b[0] / AB_norm2;//alpha = (AB * AC)/AB_norm2
        SReal d_proj = (AB * AD) / AB_norm2;
        SReal a_proj = b[1];//beta = (-CD*AC)/CD_norm2
        SReal b_proj = (CD*CB) / CD_norm2;

        if (c_proj >= 0 && c_proj <= 1) {//projection of C on AB is lying on AB
            if (d_proj > 1) {//case :
                           //             A----------------B
                           //                      C---------------D
                alpha = (1.0 + c_proj) / 2.0;
                beta = b_proj / 2.0;
            }
            else if (d_proj < 0) {//case :
                                //             A----------------B
                                //     D----------------C
                alpha = c_proj / 2.0;
                beta = (1 + a_proj) / 2.0;
            }
            else {//case :
                //             A----------------B
                //                 C------D
                alpha = (c_proj + d_proj) / 2.0;
                beta = 0.5;
            }
        }
        else if (d_proj >= 0 && d_proj <= 1) {
            if (c_proj < 0) {//case :
                           //             A----------------B
                           //     C----------------D
                alpha = d_proj / 2.0;
                beta = (1 + a_proj) / 2.0;
            }
            else {//case :
                 //          A---------------B
                 //                 D-------------C
                alpha = (1 + d_proj) / 2.0;
                beta = b_proj / 2.0;
            }
        }
        else {
            if (c_proj * d_proj < 0) {//case :
                                    //           A--------B
                                    //       D-----------------C
                alpha = 0.5;
                beta = (a_proj + b_proj) / 2.0;
            }
            else {
                if (c_proj < 0) {//case :
                               //                    A---------------B
                               // C-------------D

                    alpha = 0;
                }
                else {
                    alpha = 1;
                }

                if (a_proj < 0) {//case :
                               // A---------------B
                               //                     C-------------D
                }
                else {//case :
                     //                     A---------------B
                     //   C-------------D
                    beta = 1;
                }
            }
        }
    }

    bool ignore_p1 = true;
    bool ignore_p2 = true;
    if (ignore_p1 && beta == 0)
        return 0;
    if (ignore_p2 && beta == 1)
        return 0;

    SReal enough_to_touch = alarmDist + cap_rad;
    Vector3 p, q, pq;
    p = p1 + AB * alpha;
    q = q1 + CD * beta;
    pq = q - p;
    if (pq.norm2() >= enough_to_touch * enough_to_touch)
        return 0;

    // filter for LMD //

    if (!useLMDFilters.getValue())
    {
        //if (!testValidity(e1, pq))
        //{
        //    dmsg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
        //        << " testValidity rejected for the first segment";
        //    return 0;
        //}
        Vector3 qp = -pq;
        if (!testValidity(e2, qp))
        {
            dmsg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                << " testValidity rejected for the second segment";
            return 0;
        }
    }


    contacts->resize(contacts->size() + 1);
    DetectionOutput *detection = &*(contacts->end() - 1);

    const double contactDist = getContactDistance() + e1.getProximity() + e2.getProximity();

    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    detection->id = e2.getIndex();
    detection->point[0] = p;
    detection->point[1] = q;
    detection->normal = pq;
    detection->value = detection->normal.norm();
    detection->normal /= detection->value;
    detection->value -= (contactDist + cap_rad);
    return 1;
}


bool LocalMinDistance::testIntersection(Capsule& e1, Triangle& e2)
{
    std::cout << "on test Capsule Triangle !" << std::endl;
    return true;
}

int LocalMinDistance::computeIntersection(Capsule& e1, Triangle& e2, OutputVector* contacts)
{
    if (!e1.isActive(e2.getCollisionModel()) || !e2.isActive(e1.getCollisionModel()))
        return 0;
    //std::cout << "seems active: " << e1.getIndex() << std::endl;

    const double alarmDist = getAlarmDistance() + e1.getProximity() + e2.getProximity();

    const defaulttype::Vector3 cap_p1 = e1.point1();
    const defaulttype::Vector3 cap_p2 = e1.point2();
    const defaulttype::Vector3 tri_p1 = e2.p1();
    const defaulttype::Vector3 tri_p2 = e2.p2();
    const defaulttype::Vector3 tri_p3 = e2.p3();
    const defaulttype::Vector3 AB = cap_p2 - cap_p1;
    const SReal cap_rad = e1.radius();

    const int tri_flg = e2.flags();

    int id = e1.getIndex();
    int n = 0;


    m_A = e1.point1();
    m_B = e1.point2();
    m_T1 = tri_p1;
    m_T2 = tri_p2;
    m_T3 = tri_p3;


    SReal dist2 = (alarmDist + cap_rad);// *(alarmDist + cap_rad);

    
    n += doIntersectionTrianglePoint(dist2, e2, static_cast<core::CollisionElementIterator>(e1), e1.point1(), contacts);
    n += doIntersectionTrianglePoint(dist2, e2, static_cast<core::CollisionElementIterator>(e1), e1.point2(), contacts);
    //n += doIntersectionTrianglePoint(dist2, tri_flg, tri_p1, tri_p2, tri_p3, cap_p1, contacts, true);
    //n += doIntersectionTrianglePoint(dist2, tri_flg, tri_p1, tri_p2, tri_p3, cap_p2, contacts, true);

    const double contactDist = getContactDistance() + e1.getProximity() + e2.getProximity();
    SReal substract_dist = contactDist + cap_rad;

    //if (n != 0)
    //    std::cout << "booya! " << n << std::endl;
    m_nb = n;
    if (n == 2) {
        OutputVector::iterator detection1 = contacts->end() - 2;
        OutputVector::iterator detection2 = contacts->end() - 1;

        if (detection1->value > detection2->value - 1e-15 && detection1->value < detection2->value + 1e-15) // capsule colliding alongside the triangle
        {                         
            // contact created at the middle of the capsule
            detection1->point[0] = (detection1->point[0] + detection2->point[0]) / 2.0;
            detection1->point[1] = (detection1->point[1] + detection2->point[1]) / 2.0;
            detection1->normal = (detection1->normal + detection2->normal) / 2.0;
            detection1->value = (detection1->value + detection2->value) / 2.0 - substract_dist;
            detection1->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
            m_X = detection1->point[0];
            m_N = m_X + (detection1->normal).normalized();

            contacts->pop_back();
            n = 1;
        }
        else 
        {
            for (OutputVector::iterator detection = contacts->end() - n; detection != contacts->end(); ++detection) 
            {
                detection->normal = -detection->normal;
                detection->value -= substract_dist;
                detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
                detection->id = id;

                m_X = detection->point[0];
                m_N = m_X + (detection1->normal).normalized();
            }
        }
    }
    else 
    {
        for (OutputVector::iterator detection = contacts->end() - n; detection != contacts->end(); ++detection) 
        {
            detection->normal = -detection->normal;
            detection->value -= substract_dist;
            detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
            detection->id = id;

            m_X = detection->point[0];
            m_N = m_X + (detection->normal).normalized();
        }
    }

    return n;
    //int old_n = n;
    //n = 0;


    //if (tri_flg&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_E12)
    //    n += doIntersectionCapsuleLine(e1, static_cast<core::CollisionElementIterator>(e2), tri_p1, tri_p2, contacts);
    //    //n += doCapLineInt(cap_p1, cap_p2, cap_rad, tri_p1, tri_p2, alarmDist, contactDist, contacts, !(tri_flg&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P1), !(tri_flg&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P2));
    //if (tri_flg&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_E23)
    //    n += doIntersectionCapsuleLine(e1, static_cast<core::CollisionElementIterator>(e2), tri_p2, tri_p3, contacts);
    //    //n += doCapLineInt(cap_p1, cap_p2, cap_rad, tri_p2, tri_p3, alarmDist, contactDist, contacts, !(tri_flg&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P2), !(tri_flg&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P3));
    //if (tri_flg&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_E31)
    //    n += doIntersectionCapsuleLine(e1, static_cast<core::CollisionElementIterator>(e2), tri_p3, tri_p1, contacts);
    //    //n += doCapLineInt(cap_p1, cap_p2, cap_rad, tri_p3, tri_p1, alarmDist, contactDist, contacts, !(tri_flg&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P3), !(tri_flg&TriangleCollisionModel<sofa::defaulttype::Vec3Types>::FLAG_P1));

    //for (OutputVector::iterator detection = contacts->end() - n; detection != contacts->end(); ++detection) {
    //    detection->elem = std::pair<core::CollisionElementIterator, core::CollisionElementIterator>(e1, e2);
    //    detection->id = id;
    //}

    //return n + old_n;
}


bool LocalMinDistance::testValidity(Point &p, const Vector3 &PQ)
{
    if (!filterIntersection.getValue())
        return true;

    Vector3 pt = p.p();

    sofa::simulation::Node* node = dynamic_cast<sofa::simulation::Node*>(p.getCollisionModel()->getContext());
    if ( !(node->get< LineCollisionModel<sofa::defaulttype::Vec3Types> >()) )
        return true;

    BaseMeshTopology* topology = p.getCollisionModel()->getCollisionTopology();
    const helper::vector<Vector3>& x =(p.getCollisionModel()->getMechanicalState()->read(core::ConstVecCoordId::position())->getValue());

    const auto& trianglesAroundVertex = topology->getTrianglesAroundVertex(p.getIndex());
    const auto& edgesAroundVertex = topology->getEdgesAroundVertex(p.getIndex());
    Vector3 nMean;

    for (unsigned int i=0; i<trianglesAroundVertex.size(); i++)
    {
        unsigned int t = trianglesAroundVertex[i];
        const auto& ptr = topology->getTriangle(t);
        Vector3 nCur = (x[ptr[1]]-x[ptr[0]]).cross(x[ptr[2]]-x[ptr[0]]);
        nCur.normalize();
        nMean += nCur;
    }

    if (trianglesAroundVertex.size()==0)
    {
        for (unsigned int i=0; i<edgesAroundVertex.size(); i++)
        {
            unsigned int e = edgesAroundVertex[i];
            const auto& ped = topology->getEdge(e);
            Vector3 l = (pt - x[ped[0]]) + (pt - x[ped[1]]);
            l.normalize();
            nMean += l;
        }
    }



    if (nMean.norm()> 0.0000000001)
    {
        /// validity test with nMean, except if bothSide
        PointCollisionModel<sofa::defaulttype::Vec3Types> *pM = p.getCollisionModel();
        bool bothSide_computation = pM->bothSide.getValue();
        nMean.normalize();
        if (dot(nMean, PQ) < -angleCone.getValue()*PQ.norm() && !bothSide_computation)
        {
            return false;
        }
    }

    for (unsigned int i=0; i<edgesAroundVertex.size(); i++)
    {
        unsigned int e = edgesAroundVertex[i];
        const auto& ped = topology->getEdge(e);
        Vector3 l = (pt - x[ped[0]]) + (pt - x[ped[1]]);
        l.normalize();
        double computedAngleCone = dot(nMean , l) * coneFactor.getValue();
        if (computedAngleCone<0)
            computedAngleCone=0.0;
        computedAngleCone+=angleCone.getValue();
        if (dot(l , PQ) < -computedAngleCone*PQ.norm())
        {
            return false;
        }
    }

    return true;
}

bool LocalMinDistance::testValidity(Line &l, const Vector3 &PQ)
{
    if (!filterIntersection.getValue())
        return true;

    LineCollisionModel<sofa::defaulttype::Vec3Types> *lM = l.getCollisionModel();
    bool bothSide_computation = lM->bothSide.getValue();

    Vector3 nMean;
    Vector3 n1, n2;
    Vector3 t1, t2;

    const Vector3 &pt1 = l.p1();
    const Vector3 &pt2 = l.p2();

    Vector3 AB = pt2 - pt1;
    AB.normalize();

    BaseMeshTopology* topology = l.getCollisionModel()->getCollisionTopology();
    const helper::vector<Vector3>& x =(l.getCollisionModel()->getMechanicalState()->read(core::ConstVecCoordId::position())->getValue());
    const auto& trianglesAroundEdge = topology->getTrianglesAroundEdge(l.getIndex());

    if ( trianglesAroundEdge.size() == 2)
    {

        // which triangle is left ?
        const BaseMeshTopology::Triangle& triangle0 = topology->getTriangle(trianglesAroundEdge[0]);
        bool triangle0_is_left=false;
        if ( (l.i1()==triangle0[0]&&l.i2()==triangle0[1]) || (l.i1()==triangle0[1]&&l.i2()==triangle0[2]) || (l.i1()==triangle0[2]&&l.i2()==triangle0[0]) )
        {
            triangle0_is_left=true;
        }


        // compute the normal of the triangle situated on the right
        const BaseMeshTopology::Triangle& triangleRight = triangle0_is_left ? topology->getTriangle(trianglesAroundEdge[1]): topology->getTriangle(trianglesAroundEdge[0]);
        n1 = cross(x[triangleRight[1]]-x[triangleRight[0]], x[triangleRight[2]]-x[triangleRight[0]]);
        n1.normalize();
        nMean = n1;
        t1 = cross(n1, AB);
        t1.normalize(); // necessary ?

        // compute the normal of the triangle situated on the left
        const BaseMeshTopology::Triangle& triangleLeft = triangle0_is_left ? topology->getTriangle(trianglesAroundEdge[0]): topology->getTriangle(trianglesAroundEdge[1]);
        n2 = cross(x[triangleLeft[1]]-x[triangleLeft[0]], x[triangleLeft[2]]-x[triangleLeft[0]]);
        n2.normalize();
        nMean += n2;
        t2 = cross(AB, n2);
        t2.normalize(); // necessary ?

        nMean.normalize();

        if ((nMean*PQ) < 0  && !bothSide_computation) // test
        {
            msg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                    <<" rejected because of nMean: "<<nMean ;
            return false;
        }

        // compute the angle for the cone to filter contacts using the normal of the triangle situated on the right
        double computedAngleCone = (nMean * t1) * coneFactor.getValue();
        if (computedAngleCone<0)
            computedAngleCone=0.0;
        computedAngleCone+=angleCone.getValue();

        if (t1*PQ < -computedAngleCone*PQ.norm())
        {
            msg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                    <<" rejected because of right triangle normal: "<<n1<<" tang "<< t1 ;
            return false;
        }

        // compute the angle for the cone to filter contacts using the normal of the triangle situated on the left
        computedAngleCone = (nMean * t2) * coneFactor.getValue();
        if (computedAngleCone<0)
            computedAngleCone=0.0;
        computedAngleCone+=angleCone.getValue();

        if (t2*PQ < -computedAngleCone*PQ.norm())
        {
            msg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                <<" rejected because of left triangle normal: "<<n2 ;

            return false;
        }
    }
    else
    {
        n1 = PQ;
        n1.normalize();
        if (fabs(dot(AB,n1)) > angleCone.getValue() + 0.0001 )		// dot(AB,n1) should be equal to 0
        {
            // means that proximity was detected with a null determinant
            // in function computeIntersection
            msg_info_when(EMIT_EXTRA_DEBUG_MESSAGE)
                <<"bad case detected  -  abs(dot(AB,n1)) ="<<fabs(dot(AB,n1)) ;
            return false;
        }
    }
    return true;
}

bool LocalMinDistance::testValidity(Triangle &t, const Vector3 &PQ)
{
    TriangleCollisionModel<sofa::defaulttype::Vec3Types> *tM = t.getCollisionModel();
    bool bothSide_computation = tM->d_bothSide.getValue();

    if (!filterIntersection.getValue()  || bothSide_computation)
        return true;

    const Vector3& pt1 = t.p1();
    const Vector3& pt2 = t.p2();
    const Vector3& pt3 = t.p3();

    Vector3 n = cross(pt2-pt1,pt3-pt1);

    return ( (n*PQ) >= 0.0);
}


void LocalMinDistance::draw(const core::visual::VisualParams* vparams)
{
    //if (!vparams->displayFlags().getShowCollisionModels())
    //    return;

    //vparams->drawTool()->drawLine(m_A, m_B, Vector4(1, 0, 0, 1));
    //vparams->drawTool()->drawSphere(m_A, 0.001, Vector4(1, 0, 0, 1));
    //vparams->drawTool()->drawSphere(m_B, 0.001, Vector4(1, 0, 0, 1));
    //vparams->drawTool()->drawSphere(m_C, 0.001, Vector4(1, 0, 0, 1));


 /*   Vector3 n = cross(m_T2 - m_T1, m_T3 - m_T1);
    vparams->drawTool()->drawTriangle(m_T1, m_T2, m_T3, n, helper::types::RGBAColor(1, 0, 0, 1));
    std::ostringstream oss;
    oss << m_nb;
    vparams->drawTool()->draw3DText(m_A, 0.01, helper::types::RGBAColor(0.2, 0.2, 0, 1), oss.str().c_str());
    vparams->drawTool()->drawSphere(m_X, 0.001, Vector4(0, 1, 0, 1));
    vparams->drawTool()->drawArrow(m_X, m_N, 0.001, Vector4(0, 1, 0, 1));*/


    /*vparams->drawTool()->drawSphere(m_P, 0.001, Vector4(0, 1, 0, 1));
    vparams->drawTool()->drawLine(m_P, m_H, Vector4(0, 1, 0, 1));
    vparams->drawTool()->drawSphere(m_H, 0.001, Vector4(1, 0, 1, 1));
    vparams->drawTool()->drawSphere(m_X, 0.001,Vector4(0, 0, 1, 1));*/
}

} //namespace sofa::component::collision
