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
#include <SofaBaseMechanics/BarycentricMappers/TopologyBarycentricMapper.h>

#include <SofaBaseTopology/RegularGridTopology.h>

namespace sofa::component::mapping
{

using core::visual::VisualParams;
using sofa::defaulttype::BaseMatrix;
using sofa::defaulttype::Vec3dTypes;
using sofa::defaulttype::Vec3fTypes;
using topology::RegularGridTopology;
using topology::PointSetTopologyContainer;

/// Class allowing barycentric mapping computation on a RegularGridTopology
template<class In, class Out>
class BarycentricMapperRegularGridTopology : public TopologyBarycentricMapper<In,Out>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(BarycentricMapperRegularGridTopology,In,Out),SOFA_TEMPLATE2(TopologyBarycentricMapper,In,Out));
    typedef typename Inherit1::Real Real;
    typedef typename Inherit1::OutReal OutReal;
    typedef typename Inherit1::OutDeriv  OutDeriv;
    typedef typename Inherit1::InDeriv  InDeriv;
    typedef typename Inherit1::CubeData CubeData;
    typedef typename Inherit1::MBloc MBloc;
    typedef typename Inherit1::MatrixType MatrixType;
    typedef typename MatrixType::Index MatrixTypeIndex;
    typedef typename Inherit1::ForceMask ForceMask;

    using Index = sofa::Index;

    enum { NIn = Inherit1::NIn };
    enum { NOut = Inherit1::NOut };

public:
    ~BarycentricMapperRegularGridTopology() override;
    void clear(std::size_t reserve=0) override;
    void resize( core::State<Out>* toModel ) override;
    virtual bool isEmpty() {return this->m_map.size() == 0;}
    virtual void setTopology(topology::RegularGridTopology* _topology) {this->m_fromTopology = _topology;}
    RegularGridTopology *getTopology() {return dynamic_cast<topology::RegularGridTopology *>(this->m_fromTopology);}
    Index addPointInCube(const Index cubeIndex, const SReal* baryCoords) override;

    void init(const typename Out::VecCoord& out, const typename In::VecCoord& in) override;
    void draw(const core::visual::VisualParams*,const typename Out::VecCoord& out, const typename In::VecCoord& in) override;
    void apply( typename Out::VecCoord& out, const typename In::VecCoord& in ) override;
    void applyJ( typename Out::VecDeriv& out, const typename In::VecDeriv& in ) override;
    void applyJT( typename In::VecDeriv& out, const typename Out::VecDeriv& in ) override;
    void applyJT( typename In::MatrixDeriv& out, const typename Out::MatrixDeriv& in ) override;
    const BaseMatrix* getJ(int outSize, int inSize) override;

    template<class I, class O>
    friend std::istream& operator >> ( std::istream& in, BarycentricMapperRegularGridTopology<I, O> &b );
    template<class I, class O>
    friend std::ostream& operator << ( std::ostream& out, const BarycentricMapperRegularGridTopology<I, O> & b );

protected:
    BarycentricMapperRegularGridTopology(RegularGridTopology* fromTopology,
                                         PointSetTopologyContainer* toTopology);
    //void addMatrixContrib(MatrixType* m, int row, int col, Real value);

    helper::vector<CubeData> m_map;
    RegularGridTopology* m_fromTopology   {nullptr};
    MatrixType* m_matrixJ                 {nullptr};
    bool m_updateJ                        {false};
};

#if !defined(SOFA_COMPONENT_MAPPING_BARYCENTRICMAPPERREGULARGRIDTOPOLOGY_CPP)
extern template class SOFA_SOFABASEMECHANICS_API BarycentricMapperRegularGridTopology< Vec3dTypes, Vec3dTypes >;


#endif

} // namespace sofa::component::mapping
