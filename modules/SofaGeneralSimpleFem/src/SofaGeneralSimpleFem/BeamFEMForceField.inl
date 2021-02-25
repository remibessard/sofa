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
#include <SofaBaseTopology/TopologyData.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <SofaBaseTopology/GridTopology.h>
#include <sofa/helper/rmath.h>
#include <cassert>
#include <iostream>
#include <set>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/defaulttype/RigidTypes.h>

#include "StiffnessContainer.h"
#include "PoissonContainer.h"
#include "BeamFEMForceField.h"


namespace sofa::component::forcefield::_beamfemforcefield_
{

using core::objectmodel::BaseContext;
using defaulttype::Quat;

template<class DataTypes>
BeamFEMForceField<DataTypes>::BeamFEMForceField()
    : m_beamsData(initData(&m_beamsData, "beamsData", "Internal element data"))
    , m_indexedElements(nullptr)
    , d_poissonRatio(initData(&d_poissonRatio,(Real)0.49f,"poissonRatio","Potion Ratio"))
    , d_youngModulus(initData(&d_youngModulus,(Real)5000,"youngModulus","Young Modulus"))
    , d_radius(initData(&d_radius,(Real)0.1,"radius","radius of the section"))
    , d_radiusInner(initData(&d_radiusInner,(Real)0.0,"radiusInner","inner radius of the section for hollow beams"))
    , d_listSegment(initData(&d_listSegment,"listSegment", "apply the forcefield to a subset list of beam segments. If no segment defined, forcefield applies to the whole topology"))
    , d_useSymmetricAssembly(initData(&d_useSymmetricAssembly,false,"useSymmetricAssembly","use symmetric assembly of the matrix K"))
	, l_topology(initLink("topology", "link to the topology container"))
	, d_alternativeBeamDescription(initData(&d_alternativeBeamDescription, false, "alternativeBeamDescription", "another beam description (base on interpolated center of the beam instead of 1st extremity of the beam)"))
	, m_partialListSegment(false)
    , m_updateStiffnessMatrix(true)
    , m_assembling(false)
    , m_edgeHandler(nullptr)
{
    m_edgeHandler = new BeamFFEdgeHandler(this, &m_beamsData);

    d_poissonRatio.setRequired(true);
    d_youngModulus.setReadOnly(true);
}

template<class DataTypes>
BeamFEMForceField<DataTypes>::BeamFEMForceField(Real poissonRatio, Real youngModulus, Real radius, Real radiusInner)
    : m_beamsData(initData(&m_beamsData, "beamsData", "Internal element data"))
    , m_indexedElements(nullptr)
    , d_poissonRatio(initData(&d_poissonRatio,(Real)poissonRatio,"poissonRatio","Potion Ratio"))
    , d_youngModulus(initData(&d_youngModulus,(Real)youngModulus,"youngModulus","Young Modulus"))
    , d_radius(initData(&d_radius,(Real)radius,"radius","radius of the section"))
    , d_radiusInner(initData(&d_radiusInner,(Real)radiusInner,"radiusInner","inner radius of the section for hollow beams"))
    , d_listSegment(initData(&d_listSegment,"listSegment", "apply the forcefield to a subset list of beam segments. If no segment defined, forcefield applies to the whole topology"))
    , d_useSymmetricAssembly(initData(&d_useSymmetricAssembly,false,"useSymmetricAssembly","use symmetric assembly of the matrix K"))
    , l_topology(initLink("topology", "link to the topology container"))
	, d_alternativeBeamDescription(initData(&d_alternativeBeamDescription, false, "alternativeBeamDescription", "another beam description (base on interpolated center of the beam instead of 1st extremity of the beam)"))
    , m_partialListSegment(false)
    , m_updateStiffnessMatrix(true)
    , m_assembling(false)
    , m_edgeHandler(nullptr)
{
    m_edgeHandler = new BeamFFEdgeHandler(this, &m_beamsData);

    d_poissonRatio.setRequired(true);
    d_youngModulus.setReadOnly(true);
}

template<class DataTypes>
BeamFEMForceField<DataTypes>::~BeamFEMForceField()
{
    if(m_edgeHandler) delete m_edgeHandler;
}

template <class DataTypes>
void BeamFEMForceField<DataTypes>::bwdInit()
{
    core::behavior::BaseMechanicalState* state = this->getContext()->getMechanicalState();
    if(!state)
        msg_warning() << "Missing mechanical state";
    m_lastUpdatedStep=-1.0;
}



template <class DataTypes>
void BeamFEMForceField<DataTypes>::init()
{
    Inherit1::init();
    
    if (l_topology.empty())
    {
        msg_info() << "link to Topology container should be set to ensure right behavior. First Topology found in current context will be used.";
        l_topology.set(this->getContext()->getMeshTopologyLink());
    }

    m_topology = l_topology.get();
    msg_info() << "Topology path used: '" << l_topology.getLinkedPath() << "'";

    if (m_topology == nullptr)
    {
        msg_error() << "No topology component found at path: " << l_topology.getLinkedPath() << ", nor in current context: " << this->getContext()->name << ". Object must have a BaseMeshTopology (i.e. EdgeSetTopology or MeshTopology)";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    BaseContext* context = this->getContext();
    m_stiffnessContainer = context->BaseContext::get<container::StiffnessContainer>();
    m_poissonContainer = context->BaseContext::get<container::PoissonContainer>();

    if(m_topology->getNbEdges()==0)
    {
        msg_error() << "Topology is empty.";
        return;
    }
    m_indexedElements = &m_topology->getEdges();
    if (d_listSegment.getValue().size() == 0)
    {
        msg_info() <<"Forcefield named "<<this->getName()<<" applies to the wholo topo.";
        m_partialListSegment = false;
    }
    else
    {
        msg_info() <<"Forcefield named "<<this->getName()<<" applies to a subset of edges.";
        m_partialListSegment = true;

        for (unsigned int j=0; j<d_listSegment.getValue().size(); j++)
        {
            unsigned int i = d_listSegment.getValue()[j];
            if (i>=m_indexedElements->size())
            {
                msg_warning() <<"Defined listSegment is not compatible with topology";
                m_partialListSegment = false;
            }
        }
    }

    m_beamsData.createTopologicalEngine(m_topology,m_edgeHandler);
    m_beamsData.registerTopologicalData();

    reinit();
}

template <class DataTypes>
void BeamFEMForceField<DataTypes>::reinit()
{
    unsigned int n = m_indexedElements->size();
    m_forces.resize( this->mstate->getSize() );

    initBeams( n );
    for (unsigned int i=0; i<n; ++i)
        reinitBeam(i);
    msg_info() << "Reinit OK, "<<n<<" elements." ;
}

template <class DataTypes>
void BeamFEMForceField<DataTypes>::reinitBeam(Index i)
{
    double stiffness, length, radius, poisson, radiusInner;
    Index a = (*m_indexedElements)[i][0];
    Index b = (*m_indexedElements)[i][1];

    const VecCoord& x0 = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();
    if (m_stiffnessContainer)
        stiffness = m_stiffnessContainer->getStiffness(i) ;
    else
        stiffness =  d_youngModulus.getValue() ;

    length = (x0[a].getCenter()-x0[b].getCenter()).norm() ;

    radius = d_radius.getValue() ;
    radiusInner = d_radiusInner.getValue();
    poisson = d_poissonRatio.getValue() ;


    setBeam(i, stiffness, length, poisson, radius, radiusInner);

    computeStiffness(i,a,b);

    initLarge(i,a,b);
}

template< class DataTypes>
void BeamFEMForceField<DataTypes>::BeamFFEdgeHandler::applyCreateFunction(Index edgeIndex, BeamInfo &ei,
                                                                          const core::topology::BaseMeshTopology::Edge &,
                                                                          const sofa::helper::vector<Index> &,
                                                                          const sofa::helper::vector<double> &)
{
    if(ff)
    {
        ff->reinitBeam(edgeIndex);
        ei = ff->m_beamsData.getValue()[edgeIndex];
    }
}

template<class DataTypes>
Quat& BeamFEMForceField<DataTypes>::beamQuat(int i)
{
    helper::vector<BeamInfo>& bd = *(m_beamsData.beginEdit());
    return bd[i].quat;
}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams,
                                            DataVecDeriv &  dataF,
                                            const DataVecCoord &  dataX ,
                                            const DataVecDeriv &dataV )
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(dataV);

    VecDeriv& f = *(dataF.beginEdit());
    const VecCoord& p=dataX.getValue();
    f.resize(p.size());

    //// First compute each node rotation
    typename VecElement::const_iterator it;

    if (m_partialListSegment)
    {

        for (unsigned int j=0; j<d_listSegment.getValue().size(); j++)
        {
            unsigned int i = d_listSegment.getValue()[j];
            Element edge= (*m_indexedElements)[i];
            Index a = edge[0];
            Index b = edge[1];
            //initLarge(i,a,b);
            accumulateForceLarge( f, p, i, a, b );
        }
    }
    else
    {
        unsigned int i;
        for(it=m_indexedElements->begin(),i=0; it!=m_indexedElements->end(); ++it,++i)
        {

            Index a = (*it)[0];
            Index b = (*it)[1];

            //initLarge(i,a,b);
            accumulateForceLarge( f, p, i, a, b );
        }
    }

    dataF.endEdit();
}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams *mparams, DataVecDeriv& datadF , const DataVecDeriv& datadX)
{
    VecDeriv& df = *(datadF.beginEdit());
    const VecDeriv& dx=datadX.getValue();
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    df.resize(dx.size());

    if (m_partialListSegment)
    {
        for (unsigned int j=0; j<d_listSegment.getValue().size(); j++)
        {
            unsigned int i = d_listSegment.getValue()[j];
            Element edge= (*m_indexedElements)[i];
            Index a = edge[0];
            Index b = edge[1];

            applyStiffnessLarge(df, dx, i, a, b, kFactor);
        }
    }
    else
    {
        typename VecElement::const_iterator it;
        unsigned int i = 0;
        for(it = m_indexedElements->begin() ; it != m_indexedElements->end() ; ++it, ++i)
        {
            Index a = (*it)[0];
            Index b = (*it)[1];

            applyStiffnessLarge(df, dx, i, a, b, kFactor);
        }
    }

    datadF.endEdit();
}

template<class DataTypes>
typename BeamFEMForceField<DataTypes>::Real BeamFEMForceField<DataTypes>::pseudoDeterminantForCoef ( const defaulttype::Mat<2, 3, Real>&  M )
{
    return  M[0][1]*M[1][2] - M[1][1]*M[0][2] -  M[0][0]*M[1][2] + M[1][0]*M[0][2] + M[0][0]*M[1][1] - M[1][0]*M[0][1];
}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::computeStiffness(int i, Index , Index )
{
    Real   phiy, phiz;
    Real _L = (Real)m_beamsData.getValue()[i]._L;
    Real _A = (Real)m_beamsData.getValue()[i]._A;
    Real _nu = (Real)m_beamsData.getValue()[i]._nu;
    Real _E = (Real)m_beamsData.getValue()[i]._E;
    Real _Iy = (Real)m_beamsData.getValue()[i]._Iy;
    Real _Iz = (Real)m_beamsData.getValue()[i]._Iz;
    Real _Asy = (Real)m_beamsData.getValue()[i]._Asy;
    Real _Asz = (Real)m_beamsData.getValue()[i]._Asz;
    Real _G = (Real)m_beamsData.getValue()[i]._G;
    Real _J = (Real)m_beamsData.getValue()[i]._J;
    Real L2 = (Real) (_L * _L);
    Real L3 = (Real) (L2 * _L);
    Real EIy = (Real)(_E * _Iy);
    Real EIz = (Real)(_E * _Iz);

    // Find shear-deformation parameters
    if (_Asy == 0)
        phiy = 0.0;
    else
        phiy = (Real)(24.0*(1.0+_nu)*_Iz/(_Asy*L2));

    if (_Asz == 0)
        phiz = 0.0;
    else
        phiz = (Real)(24.0*(1.0+_nu)*_Iy/(_Asz*L2));

    helper::vector<BeamInfo>& bd = *(m_beamsData.beginEdit());
    StiffnessMatrix& k_loc = bd[i]._k_loc;

    // Define stiffness matrix 'k' in local coordinates
    k_loc.clear();
    k_loc[6][6]   = k_loc[0][0]   = _E*_A/_L;
    k_loc[7][7]   = k_loc[1][1]   = (Real)(12.0*EIz/(L3*(1.0+phiy)));
    k_loc[8][8]   = k_loc[2][2]   = (Real)(12.0*EIy/(L3*(1.0+phiz)));
    k_loc[9][9]   = k_loc[3][3]   = _G*_J/_L;
    k_loc[10][10] = k_loc[4][4]   = (Real)((4.0+phiz)*EIy/(_L*(1.0+phiz)));
    k_loc[11][11] = k_loc[5][5]   = (Real)((4.0+phiy)*EIz/(_L*(1.0+phiy)));

    k_loc[4][2]   = (Real)(-6.0*EIy/(L2*(1.0+phiz)));
    k_loc[5][1]   = (Real)( 6.0*EIz/(L2*(1.0+phiy)));
    k_loc[6][0]   = -k_loc[0][0];
    k_loc[7][1]   = -k_loc[1][1];
    k_loc[7][5]   = -k_loc[5][1];
    k_loc[8][2]   = -k_loc[2][2];
    k_loc[8][4]   = -k_loc[4][2];
    k_loc[9][3]   = -k_loc[3][3];
    k_loc[10][2]  = k_loc[4][2];
    k_loc[10][4]  = (Real)((2.0-phiz)*EIy/(_L*(1.0+phiz)));
    k_loc[10][8]  = -k_loc[4][2];
    k_loc[11][1]  = k_loc[5][1];
    k_loc[11][5]  = (Real)((2.0-phiy)*EIz/(_L*(1.0+phiy)));
    k_loc[11][7]  = -k_loc[5][1];

    for (int i=0; i<=10; i++)
        for (int j=i+1; j<12; j++)
            k_loc[i][j] = k_loc[j][i];

    m_beamsData.endEdit();
}

inline defaulttype::Quat qDiff(defaulttype::Quat a, const defaulttype::Quat& b)
{
    if (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]<0)
    {
        a[0] = -a[0];
        a[1] = -a[1];
        a[2] = -a[2];
        a[3] = -a[3];
    }
    defaulttype::Quat q = b.inverse() * a;
    return q;
}

////////////// large displacements method
template<class DataTypes>
void BeamFEMForceField<DataTypes>::initLarge(int i, Index a, Index b)
{
    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();

    defaulttype::Quat quatA, quatB, dQ;
    Vec3 dW;

    quatA = x[a].getOrientation();
    quatB = x[b].getOrientation();

    quatA.normalize();
    quatB.normalize();

    dQ = qDiff(quatB, quatA);
    dQ.normalize();

    dW = dQ.quatToRotationVector();     // TODO(e.coevoet) remove before v20:
                                        // Use of quatToRotationVector instead of toEulerVector:	    dW = dQ.quatToRotationVector();
                                        // this is done to keep the old behavior (before the
                                        // correction of the toEulerVector  function). If the
                                        // purpose was to obtain the Eulerian vector and not the
                                        // rotation vector please use the following line instead
                                        // dW = dQ.toEulerVector();

    SReal Theta = dW.norm();

    if(Theta>(SReal)0.0000001)
    {
        dW.normalize();

        beamQuat(i) = quatA*dQ.axisToQuat(dW, Theta/2);
        beamQuat(i).normalize();
    }
    else
        beamQuat(i)= quatA;

    m_beamsData.endEdit();
}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::accumulateForceLarge( VecDeriv& f, const VecCoord & x, int i, Index a, Index b )
{
    const VecCoord& x0 = this->mstate->read(core::ConstVecCoordId::restPosition())->getValue();

	// local displacement
	Displacement depl;

	if (d_alternativeBeamDescription.getValue())
	{
		Vec3 Pmiddle, Pmiddle_rest;
		Vec3 P0, P1, P2, P3;
		Real L = (x[b].getCenter() - x[a].getCenter()).norm();
		P0 = x[a].getCenter();
		P3 = x[b].getCenter();
		P1 = P0 + x[a].getOrientation().rotate(Vec3(1.0, 0, 0)*(L / 3.0));
		P2 = P3 + x[b].getOrientation().rotate(Vec3(-1.0, 0, 0)*(L / 3.0));
		Pmiddle = 0.125 * (P0 + P3) + 0.375 * (P1 + P2);  //bezier curve interpolation

		Pmiddle_rest = x0[a].getCenter() + (x0[b].getCenter() - x0[a].getCenter()) / 2.0; //assuming straight beam at rest position

		beamQuat(i).slerp(x[a].getOrientation(), x[b].getOrientation(), (float)0.5, true);
		beamQuat(i).normalize();

		Vec3 u0, u3;

		// translations //
		Vec3 PmidP0, PmidP0_rest, PmidP3, PmidP3_rest;
		PmidP0_rest = x0[a].getCenter() - Pmiddle_rest;
		PmidP0_rest = x0[a].getOrientation().inverseRotate(PmidP0_rest);
		PmidP0 = P0 - Pmiddle;
		PmidP0 = beamQuat(i).inverseRotate(PmidP0);
		u0 = PmidP0 - PmidP0_rest;

		PmidP3_rest = x0[b].getCenter() - Pmiddle_rest;
		PmidP3_rest = x0[a].getOrientation().inverseRotate(PmidP3_rest);
		PmidP3 = P3 - Pmiddle;
		PmidP3 = beamQuat(i).inverseRotate(PmidP3);
		u3 = PmidP3 - PmidP3_rest;

		depl[0] = u0[0]; depl[1] = u0[1]; depl[2] = u0[2];
		depl[6] = u3[0]; depl[7] = u3[1]; depl[8] = u3[2];


		// rotations //
		Quat Qmiddle, Qmiddle_rest;
		Qmiddle = beamQuat(i);
		Qmiddle_rest = x0[a].getOrientation();  //assuming straight beam at rest position

		Quat QmidQ0_rest, QmidQ0, QmidQ3_rest, QmidQ3;
		QmidQ0_rest = qDiff(x0[a].getOrientation(), Qmiddle_rest);
		QmidQ0_rest.normalize();
		QmidQ0 = qDiff(x[a].getOrientation(), Qmiddle);
		QmidQ0.normalize();
		Quat dQQ_rest = qDiff(QmidQ0, QmidQ0_rest);
		dQQ_rest.normalize();
		u0 = dQQ_rest.quatToRotationVector();

		QmidQ3_rest = qDiff(x0[b].getOrientation(), Qmiddle_rest);
		QmidQ3_rest.normalize();
		QmidQ3 = qDiff(Qmiddle, x[b].getOrientation());
		QmidQ3.normalize();
		dQQ_rest = qDiff(QmidQ3, QmidQ3_rest);
		dQQ_rest.normalize();
		dQQ_rest = dQQ_rest.inverse();
		u3 = dQQ_rest.quatToRotationVector();

		depl[3] = u0[0]; depl[4] = u0[1]; depl[5] = u0[2];
		depl[9] = u3[0]; depl[10] = u3[1]; depl[11] = u3[2];

		double Rot[4][4];
		beamQuat(i).inverse().buildRotationMatrix(Rot);
		defaulttype::Mat<6, 6, Real> J;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				J[i][j] = Rot[i][j];
				J[i + 3][j + 3] = Rot[i][j];
				J[i][j + 3] = 0.0;// R_Origin[i][j];
				J[i + 3][j] = 0.0;
			}
		}

		defaulttype::VecNoInit<6, Real> JtD;
		JtD[0] = J[0][0] * depl[0] +/*J[ 1][0]*depl[ 1]+  J[ 2][0]*depl[ 2]+*/
			J[3][0] * depl[3] +/*J[ 4][0]*depl[ 4]+  J[ 5][0]*depl[ 5]+*/
			J[0][0] * depl[6] +/*J[ 7][0]*depl[ 7]+  J[ 8][0]*depl[ 8]+*/
			J[3][0] * depl[9] /*J[10][0]*depl[10]+  J[11][0]*depl[11]*/;
		JtD[1] = /*J[ 0][1]*depl[ 0]+*/J[1][1] * depl[1] +/*J[ 2][1]*depl[ 2]+*/
			/*J[ 3][1]*depl[ 3]+*/J[4][1] * depl[4] +/*J[ 5][1]*depl[ 5]+*/
			/*J[ 6][1]*depl[ 6]+*/J[1][1] * depl[7] +/*J[ 8][1]*depl[ 8]+*/
			/*J[ 9][1]*depl[ 9]+*/J[4][1] * depl[10] /*J[11][1]*depl[11]*/;
		JtD[2] = /*J[ 0][2]*depl[ 0]+  J[ 1][2]*depl[ 1]+*/J[2][2] * depl[2] +
			/*J[ 3][2]*depl[ 3]+  J[ 4][2]*depl[ 4]+*/J[5][2] * depl[5] +
			/*J[ 6][2]*depl[ 6]+  J[ 7][2]*depl[ 7]+*/J[2][2] * depl[8] +
			/*J[ 9][2]*depl[ 9]+  J[10][2]*depl[10]+*/J[5][2] * depl[11];
		JtD[3] = J[0][3] * depl[0] + J[1][3] * depl[1] +/*J[ 2][3]*depl[ 2]+*/
			J[3][3] * depl[3] + J[4][3] * depl[4] +/*J[ 5][3]*depl[ 5]+*/
			J[0][3] * depl[6] + J[1][3] * depl[7] +/*J[ 8][3]*depl[ 8]+*/
			J[3][3] * depl[9] + J[4][3] * depl[10] /*J[11][3]*depl[11]*/;
		JtD[4] = /*J[ 0][4]*depl[ 0]+*/J[1][4] * depl[1] + J[2][4] * depl[2] +
			/*J[ 3][4]*depl[ 3]+*/J[4][4] * depl[4] + J[5][4] * depl[5] +
			/*J[ 6][4]*depl[ 6]+*/J[1][4] * depl[7] + J[2][4] * depl[8] +
			/*J[ 9][4]*depl[ 9]+*/J[4][4] * depl[10] + J[5][4] * depl[11];
		JtD[5] = J[0][5] * depl[0] +/*J[ 1][5]*depl[ 1]*/ J[2][5] * depl[2] +
			J[3][5] * depl[3] +/*J[ 4][5]*depl[ 4]*/ J[5][5] * depl[5] +
			J[0][5] * depl[6] +/*J[ 7][5]*depl[ 7]*/ J[2][5] * depl[8] +
			J[3][5] * depl[9] +/*J[10][5]*depl[10]*/ J[5][5] * depl[11];

		//std::cout << "norm jtd: " << JtD.norm() << std::endl;
		// eventually remove a part of the strain to simulate plasticity
		if (d_plasticMaxThreshold.getValue() > 0)
		{
			defaulttype::VecNoInit<6, Real> elasticStrain = JtD; // JtD is the total strain
			elasticStrain -= m_plasticStrains[i]; // totalStrain = elasticStrain + plasticStrain

			if (elasticStrain.norm2() > d_plasticYieldThreshold.getValue()*d_plasticYieldThreshold.getValue())
				m_plasticStrains[i] += d_plasticCreep.getValue() * elasticStrain;

			Real plasticStrainNorm2 = m_plasticStrains[i].norm2();
			if (plasticStrainNorm2 > d_plasticMaxThreshold.getValue()*d_plasticMaxThreshold.getValue())
				m_plasticStrains[i] *= d_plasticMaxThreshold.getValue() / helper::rsqrt(plasticStrainNorm2);

			// remaining elasticStrain = totatStrain - plasticStrain
			JtD -= m_plasticStrains[i];
		}
	}
	else
	{
		beamQuat(i) = x[a].getOrientation();
		beamQuat(i).normalize();

		defaulttype::Vec<3, Real> u, P1P2, P1P2_0;

		// translations //
		P1P2_0 = x0[b].getCenter() - x0[a].getCenter();
		P1P2_0 = x0[a].getOrientation().inverseRotate(P1P2_0);
		P1P2 = x[b].getCenter() - x[a].getCenter();
		P1P2 = beamQuat(i).inverseRotate(P1P2);
		u = P1P2 - P1P2_0;

		depl[0] = 0.0; 	depl[1] = 0.0; 	depl[2] = 0.0;
		depl[6] = u[0]; depl[7] = u[1]; depl[8] = u[2];

		// rotations //
		defaulttype::Quat dQ0, dQ;

		dQ0 = qDiff(x0[b].getOrientation(), x0[a].getOrientation());
		dQ = qDiff(x[b].getOrientation(), x[a].getOrientation());

		dQ0.normalize();
		dQ.normalize();

		defaulttype::Quat tmpQ = qDiff(dQ, dQ0);
		tmpQ.normalize();

		u = tmpQ.quatToRotationVector();	// Use of quatToRotationVector instead of toEulerVector:	    u = tmpQ.quatToRotationVector();
	  // this is done to keep the old behavior (before the
	  // correction of the toEulerVector  function). If the
	  // purpose was to obtain the Eulerian vector and not the
	  // rotation vector please use the following line instead
	  // u = tmpQ.toEulerVector();


		depl[3] = 0.0; 	depl[4] = 0.0; 	depl[5] = 0.0;
		depl[9] = u[0]; depl[10] = u[1]; depl[11] = u[2];
	}

	m_beamsData.endEdit();

	// this computation can be optimised: (we know that half of "depl" is null)
	Displacement force = m_beamsData.getValue()[i]._k_loc * depl;


	// Apply lambda transpose (we use the rotation value of point a for the beam)
	Vec3 fa1 = beamQuat(i).rotate(defaulttype::Vec3d(force[0], force[1], force[2]));
	Vec3 fa2 = beamQuat(i).rotate(defaulttype::Vec3d(force[3], force[4], force[5]));

	Vec3 fb1 = beamQuat(i).rotate(defaulttype::Vec3d(force[6], force[7], force[8]));
	Vec3 fb2 = beamQuat(i).rotate(defaulttype::Vec3d(force[9], force[10], force[11]));

	f[a] += Deriv(-fa1, -fa2);
	f[b] += Deriv(-fb1, -fb2);

}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::applyStiffnessLarge(VecDeriv& df, const VecDeriv& dx, int i, Index a, Index b, double fact)
{
    Displacement local_depl;
    defaulttype::Vec<3,Real> u;
    defaulttype::Quat& q = beamQuat(i);
    q.normalize();

    u = q.inverseRotate(getVCenter(dx[a]));
    local_depl[0] = u[0];
    local_depl[1] = u[1];
    local_depl[2] = u[2];

    u = q.inverseRotate(getVOrientation(dx[a]));
    local_depl[3] = u[0];
    local_depl[4] = u[1];
    local_depl[5] = u[2];

    u = q.inverseRotate(getVCenter(dx[b]));
    local_depl[6] = u[0];
    local_depl[7] = u[1];
    local_depl[8] = u[2];

    u = q.inverseRotate(getVOrientation(dx[b]));
    local_depl[9] = u[0];
    local_depl[10] = u[1];
    local_depl[11] = u[2];

    Displacement local_force = m_beamsData.getValue()[i]._k_loc * local_depl;

    Vec3 fa1 = q.rotate(defaulttype::Vec3d(local_force[0],local_force[1] ,local_force[2] ));
    Vec3 fa2 = q.rotate(defaulttype::Vec3d(local_force[3],local_force[4] ,local_force[5] ));
    Vec3 fb1 = q.rotate(defaulttype::Vec3d(local_force[6],local_force[7] ,local_force[8] ));
    Vec3 fb2 = q.rotate(defaulttype::Vec3d(local_force[9],local_force[10],local_force[11]));


    df[a] += Deriv(-fa1,-fa2) * fact;
    df[b] += Deriv(-fb1,-fb2) * fact;
}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::addKToMatrix(const sofa::core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix )
{
    sofa::core::behavior::MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    Real k = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    defaulttype::BaseMatrix* mat = r.matrix;

    if (r)
    {
        unsigned int i=0;

        unsigned int &offset = r.offset;

        if (m_partialListSegment)
        {

            for (unsigned int j=0; j<d_listSegment.getValue().size(); j++)
            {

                i = d_listSegment.getValue()[j];
                Element edge= (*m_indexedElements)[i];
                Index a = edge[0];
                Index b = edge[1];

                defaulttype::Quat& q = beamQuat(i);
                q.normalize();
                Transformation R,Rt;
                q.toMatrix(R);
                Rt.transpose(R);
                const StiffnessMatrix& K0 = m_beamsData.getValue()[i]._k_loc;
                StiffnessMatrix K;
                for (int x1=0; x1<12; x1+=3)
                    for (int y1=0; y1<12; y1+=3)
                    {
                        defaulttype::Mat<3,3,Real> m;
                        K0.getsub(x1,y1, m);
                        m = R*m*Rt;
                        K.setsub(x1,y1, m);
                    }
                int index[12];
                for (int x1=0; x1<6; x1++)
                    index[x1] = offset+a*6+x1;
                for (int x1=0; x1<6; x1++)
                    index[6+x1] = offset+b*6+x1;
                for (int x1=0; x1<12; ++x1)
                    for (int y1=0; y1<12; ++y1)
                        mat->add(index[x1], index[y1], - K(x1,y1)*k);

            }

        }
        else
        {
            typename VecElement::const_iterator it;
            for(it = m_indexedElements->begin() ; it != m_indexedElements->end() ; ++it, ++i)
            {
                Index a = (*it)[0];
                Index b = (*it)[1];

                defaulttype::Quat& q = beamQuat(i);
                q.normalize();
                Transformation R,Rt;
                q.toMatrix(R);
                Rt.transpose(R);
                const StiffnessMatrix& K0 = m_beamsData.getValue()[i]._k_loc;
                StiffnessMatrix K;
                bool exploitSymmetry = d_useSymmetricAssembly.getValue();

                if (exploitSymmetry) {
                    for (int x1=0; x1<12; x1+=3) {
                        for (int y1=x1; y1<12; y1+=3)
                        {
                            defaulttype::Mat<3,3,Real> m;
                            K0.getsub(x1,y1, m);
                            m = R*m*Rt;

                            for (int i=0; i<3; i++)
                                for (int j=0; j<3; j++) {
                                    K.elems[i+x1][j+y1] += m[i][j];
                                    K.elems[j+y1][i+x1] += m[i][j];
                                }
                            if (x1 == y1)
                                for (int i=0; i<3; i++)
                                    for (int j=0; j<3; j++)
                                        K.elems[i+x1][j+y1] *= double(0.5);

                        }
                    }
                } else  {
                    for (int x1=0; x1<12; x1+=3) {
                        for (int y1=0; y1<12; y1+=3)
                        {
                            defaulttype::Mat<3,3,Real> m;
                            K0.getsub(x1,y1, m);
                            m = R*m*Rt;
                            K.setsub(x1,y1, m);
                        }
                    }
                }

                int index[12];
                for (int x1=0; x1<6; x1++)
                    index[x1] = offset+a*6+x1;
                for (int x1=0; x1<6; x1++)
                    index[6+x1] = offset+b*6+x1;
                for (int x1=0; x1<12; ++x1)
                    for (int y1=0; y1<12; ++y1)
                        mat->add(index[x1], index[y1], - K(x1,y1)*k);

            }
        }

    }

}

template<class DataTypes>
SReal BeamFEMForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams* mparams, const DataVecCoord&  x) const
{
    SOFA_UNUSED(x);
    SOFA_UNUSED(mparams);
    msg_warning() << "Method getPotentialEnergy not implemented yet.";
    return 0.0;
}


template<class DataTypes>
void BeamFEMForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowForceFields()) return;
    if (!this->mstate) return;

    vparams->drawTool()->saveLastState();

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();

    std::vector< defaulttype::Vector3 > points[3];

    if (m_partialListSegment)
    {
        for (unsigned int j=0; j<d_listSegment.getValue().size(); j++)
            drawElement(d_listSegment.getValue()[j], points, x);
    }
    else
    {
        for (unsigned int i=0; i<m_indexedElements->size(); ++i)
            drawElement(i, points, x);
    }
    vparams->drawTool()->drawLines(points[0], 1, sofa::helper::types::RGBAColor::red());
    vparams->drawTool()->drawLines(points[1], 1, sofa::helper::types::RGBAColor::green());
    vparams->drawTool()->drawLines(points[2], 1, sofa::helper::types::RGBAColor::blue());

    vparams->drawTool()->restoreLastState();
}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::computeBBox(const core::ExecParams* params, bool onlyVisible)
{
    SOFA_UNUSED(params);

    if( !onlyVisible ) return;


    static const Real max_real = std::numeric_limits<Real>::max();
    static const Real min_real = std::numeric_limits<Real>::lowest();
    Real maxBBox[3] = {min_real,min_real,min_real};
    Real minBBox[3] = {max_real,max_real,max_real};


    const size_t npoints = this->mstate->getSize();
    const VecCoord& p = this->mstate->read(core::ConstVecCoordId::position())->getValue();

    for (size_t i=0; i<npoints; i++)
    {
        const defaulttype::Vector3 &pt = p[i].getCenter();

        for (int c=0; c<3; c++)
        {
            if (pt[c] > maxBBox[c]) maxBBox[c] = pt[c];
            else if (pt[c] < minBBox[c]) minBBox[c] = pt[c];
        }
    }

    this->f_bbox.setValue(sofa::defaulttype::TBoundingBox<Real>(minBBox,maxBBox));

}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::drawElement(int i, std::vector< defaulttype::Vector3 >* points, const VecCoord& x)
{
    Index a = (*m_indexedElements)[i][0];
    Index b = (*m_indexedElements)[i][1];
    defaulttype::Vec3d p; p = (x[a].getCenter()+x[b].getCenter())*0.5;
    defaulttype::Vec3d beamVec;
    beamVec[0]=m_beamsData.getValue()[i]._L*0.5; beamVec[1] = 0.0; beamVec[2] = 0.0;

    const defaulttype::Quat& q = beamQuat(i);

    // axis X
    points[0].push_back(p - q.rotate(beamVec) );
    points[0].push_back(p + q.rotate(beamVec) );

    // axis Y
    beamVec[0]=0.0; beamVec[1] = m_beamsData.getValue()[i]._r*0.5;
    points[1].push_back(p );
    points[1].push_back(p + q.rotate(beamVec) );

    // axis Z
    beamVec[1]=0.0; beamVec[2] = m_beamsData.getValue()[i]._r*0.5;
    points[2].push_back(p);
    points[2].push_back(p + q.rotate(beamVec) );
}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::initBeams(std::size_t size)
{
    helper::vector<BeamInfo>& bd = *(m_beamsData.beginEdit());
    bd.resize(size);
    m_beamsData.endEdit();
}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::setUpdateStiffnessMatrix(bool val)
{
    this->m_updateStiffnessMatrix = val;
}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::setComputeGlobalMatrix(bool val)
{
    this->m_assembling= val;
}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::setBeam(Index i, double E, double L, double nu, double r, double rInner)
{
    helper::vector<BeamInfo>& bd = *(m_beamsData.beginEdit());
    bd[i].init(E,L,nu,r,rInner);
    m_beamsData.endEdit();
}

template<class DataTypes>
void BeamFEMForceField<DataTypes>::BeamInfo::init(double E, double L, double nu, double r, double rInner)
{
    _E = E;
    _E0 = E;
    _nu = nu;
    _L = L;
    _r = r;
    _rInner = rInner;

    _G=_E/(2.0*(1.0+_nu));
    _Iz = M_PI*(r*r*r*r - rInner*rInner*rInner*rInner)/4.0;

    _Iy = _Iz ;
    _J = _Iz+_Iy;
    _A = M_PI*(r*r - rInner*rInner);


    _Asy = 0.0;
    _Asz = 0.0;
}

} // namespace  sofa::component::forcefield::_beamfemforcefield_
