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
#include <SofaBaseMechanics/BarycentricMappers/BarycentricMapperTopologyContainer.h>
#include <sofa/core/State.h>
#include <sofa/core/visual/VisualParams.h>

namespace sofa::component::mapping::_barycentricmappertopologycontainer_
{

using defaulttype::Vec3d;
using defaulttype::Vec3i;
typedef typename core::topology::BaseMeshTopology::SeqEdges SeqEdges;

template <class In, class Out, class MappingDataType, class Element>
BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::BarycentricMapperTopologyContainer(core::topology::BaseMeshTopology* fromTopology,
                                                                                                       topology::PointSetTopologyContainer* toTopology)
     : Inherit1(fromTopology, toTopology),
       d_map(initData(&d_map,"map", "mapper data")),
       m_matrixJ(nullptr),
       m_updateJ(true)
 {}



template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::clear(std::size_t size)
{
    helper::vector<MappingDataType>& vectorData = *(d_map.beginEdit());
    vectorData.clear();
    if ( size>0 ) vectorData.reserve ( size );
    d_map.endEdit();
}


template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::initHashing(const typename In::VecCoord& in)
{
    computeHashingCellSize(in);
    computeHashTable(in);
}


template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::computeHashingCellSize(const typename In::VecCoord& in )
{
    // The grid cell size is set to the average edge length of all elements
    const SeqEdges& edges = m_fromTopology->getEdges();
    Real averageLength=0.;

    if(edges.size()>0)
    {
        for(unsigned int i=0; i<edges.size(); i++)
        {
            Edge edge = edges[i];
            averageLength += (in[edge[0]]-in[edge[1]]).norm();
        }
        averageLength/=Real(edges.size());
    }
    else
    {
        const helper::vector<Element>& elements = getElements();

        for(unsigned int i=0; i<elements.size(); i++)
        {
            Element element = elements[i];
            averageLength += (in[element[0]]-in[element[1]]).norm();
        }
        averageLength/=Real(elements.size());
    }

    m_gridCellSize = averageLength;
    m_convFactor = 1./Real(m_gridCellSize);
}


template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::computeHashTable( const typename In::VecCoord& in )
{
    const helper::vector<Element>& elements = getElements();
    m_hashTableSize = elements.size()*2; // Next prime number would be better
    m_hashTable.clear();
    if(m_hashTable.size()<m_hashTableSize)
        m_hashTable.reserve(m_hashTableSize);

    for(unsigned int i=0; i<elements.size(); i++)
    {
        Element element = elements[i];
        Vector3 min=in[element[0]], max=in[element[0]];

        for(unsigned int j=0; j<element.size(); j++)
        {
            unsigned int pointId = element[j];
            for(int k=0; k<3; k++)
            {
                if(in[pointId][k]<min[k]) min[k]=in[pointId][k];
                if(in[pointId][k]>max[k]) max[k]=in[pointId][k];
            }
        }

        Vec3i i_min=getGridIndices(min);
        Vec3i i_max=getGridIndices(max);

        for(int j=i_min[0]; j<=i_max[0]; j++)
            for(int k=i_min[1]; k<=i_max[1]; k++)
                for(int l=i_min[2]; l<=i_max[2]; l++)
                    m_hashTable[Key(j,k,l)].push_back(i);
    }
}


template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::init ( const typename Out::VecCoord& out, const typename In::VecCoord& in )
{
    initHashing(in);
    this->clear ( int(out.size()) );
    computeBasesAndCenters(in);

    // Compute distances to get nearest element and corresponding bary coef
    const helper::vector<Element>& elements = getElements();
    for ( unsigned int i=0; i<out.size(); i++ )
    {
        Vector3 outPos = Out::getCPos(out[i]);
        NearestParams nearestParams;

        // Search nearest element in grid cell
        Vec3i gridIds = getGridIndices(outPos);
        Key key(gridIds[0],gridIds[1],gridIds[2]);

        auto it_entries = m_hashTable.find(key);
        if( it_entries != m_hashTable.end() )
        {
            for(auto entry : it_entries->second)
            {
                Vector3 inPos = in[elements[entry][0]];
                checkDistanceFromElement(entry, outPos, inPos, nearestParams);
            }
        }

        if(nearestParams.elementId==std::numeric_limits<unsigned int>::max()) // No element in grid cell, perform exhaustive search
        {
            for ( unsigned int e = 0; e < elements.size(); e++ )
            {
                Vector3 inPos = in[elements[e][0]];
                checkDistanceFromElement(e, outPos, inPos, nearestParams);
            }
            addPointInElement(nearestParams.elementId, nearestParams.baryCoords.ptr());
        }
        else if(fabs(nearestParams.distance)>m_gridCellSize/2.) // Nearest element in grid cell may not be optimal, check neighbors
        {
            Vec3i centerGridIds = gridIds;
            for(int xId=-1; xId<=1; xId++)
                for(int yId=-1; yId<=1; yId++)
                    for(int zId=-1; zId<=1; zId++)
                    {
                        gridIds = Vec3i(centerGridIds[0]+xId,centerGridIds[1]+yId,centerGridIds[2]+zId);
                        Key key(gridIds[0],gridIds[1],gridIds[2]);

                        auto it_entries = m_hashTable.find(key);
                        if( it_entries != m_hashTable.end() )
                        {
                            for(auto entry : it_entries->second)
                            {
                                Vector3 inPos = in[elements[entry][0]];
                                checkDistanceFromElement(entry, outPos, inPos, nearestParams);
                            }
                        }
                    }
            addPointInElement(nearestParams.elementId, nearestParams.baryCoords.ptr());
        }
        else
            addPointInElement(nearestParams.elementId, nearestParams.baryCoords.ptr());
    }
}


template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::computeBasesAndCenters( const typename In::VecCoord& in )
{
    const helper::vector<Element>& elements = getElements();
    m_bases.resize ( elements.size() );
    m_centers.resize ( elements.size() );

    for ( unsigned int e = 0; e < elements.size(); e++ )
    {
        Element element = elements[e];

        Mat3x3d base;
        computeBase(base,in,element);
        m_bases[e] = base;

        Vector3 center;
        computeCenter(center,in,element);
        m_centers[e] = center;
    }
}


template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::checkDistanceFromElement(unsigned int e,
                                                                                                  const Vector3& outPos,
                                                                                                  const Vector3& inPos,
                                                                                                  NearestParams& nearestParams)
{
    Vector3 bary = m_bases[e] * ( outPos - inPos);
    double dist;
    computeDistance(dist, bary);
    if ( dist>0 )
        dist = ( outPos-m_centers[e] ).norm2();
    if ( dist<nearestParams.distance )
    {
        nearestParams.baryCoords = bary;
        nearestParams.distance = dist;
        nearestParams.elementId = e;
    }
};


template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::applyJT ( typename In::MatrixDeriv& out, const typename Out::MatrixDeriv& in )
{
    typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();
    const helper::vector< Element >& elements = getElements();

    for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename Out::MatrixDeriv::ColConstIterator colItEnd = rowIt.end();
        typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin();

        if (colIt != colItEnd)
        {
            typename In::MatrixDeriv::RowIterator o = out.writeLine(rowIt.index());

            for ( ; colIt != colItEnd; ++colIt)
            {
                unsigned indexIn = colIt.index();
                InDeriv data = InDeriv(Out::getDPos(colIt.val()));

                const Element& element = elements[d_map.getValue()[indexIn].in_index];

                helper::vector<SReal> baryCoef = getBaryCoef(d_map.getValue()[indexIn].baryCoords);
                for (unsigned int j=0; j<element.size(); j++)
                    o.addCol(element[j], data*baryCoef[j]);
            }
        }
    }
}



template <class In, class Out, class MappingDataType, class Element>
const defaulttype::BaseMatrix* BarycentricMapperTopologyContainer<In,Out,MappingDataType, Element>::getJ(int outSize, int inSize)
{
	return nullptr;
#if 0
    if (m_matrixJ && !m_updateJ)
        return m_matrixJ;

    if (!m_matrixJ) m_matrixJ = new MatrixType;
    if (m_matrixJ->rowBSize() != MatrixTypeIndex(outSize) || m_matrixJ->colBSize() != MatrixTypeIndex(inSize))
        m_matrixJ->resize(outSize*NOut, inSize*NIn);
    else
        m_matrixJ->clear();

    const helper::vector<Element>& elements = getElements();

    for( size_t outId=0 ; outId<this->maskTo->size() ; ++outId)
    {
        if( !this->maskTo->getEntry(outId) ) continue;

        const Element& element = elements[d_map.getValue()[outId].in_index];

        helper::vector<SReal> baryCoef = getBaryCoef(d_map.getValue()[outId].baryCoords);
        for (unsigned int j=0; j<element.size(); j++)
            this->addMatrixContrib(m_matrixJ, int(outId), element[j], baryCoef[j]);
    }

    m_matrixJ->compress();
    m_updateJ = false;
    return m_matrixJ;
#endif
}


template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::applyJT ( typename In::VecDeriv& out, const typename Out::VecDeriv& in )
{
    const helper::vector<Element>& elements = getElements();

    ForceMask& mask = *this->maskFrom;
    for( size_t i=0 ; i<this->maskTo->size() ; ++i)
    {
        if( !this->maskTo->getEntry(i) ) continue;

        Index index = d_map.getValue()[i].in_index;
        const Element& element = elements[index];

        const typename Out::DPos inPos = Out::getDPos(in[i]);
        helper::vector<SReal> baryCoef = getBaryCoef(d_map.getValue()[i].baryCoords);
        for (unsigned int j=0; j<element.size(); j++)
        {
            out[element[j]] += inPos * baryCoef[j];
            mask.insertEntry(element[j]);
        }
    }
}

template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::applyJ ( typename Out::VecDeriv& out, const typename In::VecDeriv& in )
{
    out.resize( d_map.getValue().size() );

    const helper::vector<Element>& elements = getElements();

    for( size_t i=0 ; i<this->maskTo->size() ; ++i)
    {
        if( this->maskTo->isActivated() && !this->maskTo->getEntry(i) ) continue;

        Index index = d_map.getValue()[i].in_index;
        const Element& element = elements[index];

        helper::vector<SReal> baryCoef = getBaryCoef(d_map.getValue()[i].baryCoords);
        InDeriv inPos{0.,0.,0.};
        for (unsigned int j=0; j<element.size(); j++)
            inPos += in[element[j]] * baryCoef[j];

        Out::setDPos(out[i] , inPos);
    }
}


template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::resize( core::State<Out>* toModel )
{
    toModel->resize(Size(d_map.getValue().size()));
}

template<class In, class Out, class MappingDataType, class Element>
bool BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::isEmpty()
{
    return d_map.getValue().empty();
}



template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::apply ( typename Out::VecCoord& out, const typename In::VecCoord& in )
{
    out.resize( d_map.getValue().size() );

    const helper::vector<Element>& elements = getElements();
    for ( unsigned int i=0; i<d_map.getValue().size(); i++ )
    {
        Index index = d_map.getValue()[i].in_index;
        const Element& element = elements[index];

        helper::vector<SReal> baryCoef = getBaryCoef(d_map.getValue()[i].baryCoords);
        InDeriv inPos{0.,0.,0.};
        for (unsigned int j=0; j<element.size(); j++)
            inPos += in[element[j]] * baryCoef[j];

        Out::setCPos(out[i] , inPos);
    }
}


template <class In, class Out, class MappingDataType, class Element>
void BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::draw  (const core::visual::VisualParams* vparams,
                                                                                const typename Out::VecCoord& out,
                                                                                const typename In::VecCoord& in )
{
    // Draw line between mapped node (out) and nodes of nearest element (in)
    const helper::vector<Element>& elements = getElements();

    std::vector< Vector3 > points;
    {
        for ( unsigned int i=0; i<d_map.getValue().size(); i++ )
        {
            Index index = d_map.getValue()[i].in_index;
            const Element& element = elements[index];
            helper::vector<SReal> baryCoef = getBaryCoef(d_map.getValue()[i].baryCoords);
            for ( unsigned int j=0; j<element.size(); j++ )
            {
                if ( baryCoef[j]<=-0.0001 || baryCoef[j]>=0.0001 )
                {
                    points.push_back ( Out::getCPos(out[i]) );
                    points.push_back ( in[element[j]] );
                }
            }
        }
    }
    vparams->drawTool()->drawLines ( points, 1, sofa::helper::types::RGBAColor::green());
}


template <class In, class Out, class MappingDataType, class Element>
Vec3i BarycentricMapperTopologyContainer<In,Out,MappingDataType,Element>::getGridIndices(const Vector3& pos)
{
    Vec3i i_x;
    for(int i=0; i<3; i++)
        i_x[i]=int(std::floor(Real(pos[i]*m_convFactor)));

    return i_x;
}


template<class In, class Out, class MappingData, class Element>
std::istream& operator >> ( std::istream& in, BarycentricMapperTopologyContainer<In, Out, MappingData, Element> &b )
{
    unsigned int size_vec;

    in >> size_vec;
    helper::vector<MappingData>& m = *(b.d_map.beginEdit());
    m.clear();

    MappingData value;
    for (unsigned int i=0; i<size_vec; i++)
    {
        in >> value;
        m.push_back(value);
    }
    b.d_map.endEdit();
    return in;
}

template<class In, class Out, class MappingData, class Element>
std::ostream& operator << ( std::ostream& out, const BarycentricMapperTopologyContainer<In, Out, MappingData, Element> & b )
{
    out << b.d_map.getValue().size();
    out << " " ;
    out << b.d_map;

    return out;
}

} // namespace sofa::component::mapping::_barycentricmappertopologycontainer_
