/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2020 MGH, INRIA, USTL, UJF, CNRS, InSimo            *
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
#ifndef SOFA_DEFAULTTYPE_COMPRESSEDROWSPARSEMATRIXMECHANICAL_INL
#define SOFA_DEFAULTTYPE_COMPRESSEDROWSPARSEMATRIXMECHANICAL_INL

#include <SofaBaseLinearSolver/CompressedRowSparseMatrixMechanical.h>

namespace sofa::component::linearsolver
{

template <> template <>
inline void CompressedRowSparseMatrixMechanical<double>::filterValues(CompressedRowSparseMatrixMechanical<defaulttype::Mat<3,3,double> >& M, filter_fn* filter, const Real ref, bool keepEmptyRows)
{
    M.compress();
    nRow = M.rowSize();
    nCol = M.colSize();
    nBlocRow = 1;
    nBlocCol = 1;
    rowIndex.clear();
    rowBegin.clear();
    colsIndex.clear();
    colsValue.clear();
    btemp.clear();
    skipCompressZero = true;
    rowIndex.reserve(M.rowIndex.size()*3);
    rowBegin.reserve(M.rowBegin.size()*3);
    colsIndex.reserve(M.colsIndex.size()*9);
    colsValue.reserve(M.colsValue.size()*9);

    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.clear();

    Index vid = 0;
    for (std::size_t rowId = 0; rowId < M.rowIndex.size(); ++rowId)
    {
        Index i = M.rowIndex[rowId] * 3;

        Range rowRange(M.rowBegin[rowId], M.rowBegin[rowId+1]);

        for (Index lb = 0; lb<3 ; lb++)
        {
            rowIndex.push_back(i+lb);
            rowBegin.push_back(vid);

            for (std::size_t xj = static_cast<std::size_t>(rowRange.begin()); xj < static_cast<std::size_t>(rowRange.end()); ++xj)
            {
                Index j = M.colsIndex[xj] * 3;
                defaulttype::Mat<3,3,double> b = M.colsValue[xj];
                if ((*filter)(i+lb,j+0,b[lb][0],ref))
                {
                    colsIndex.push_back(j+0);
                    colsValue.push_back(b[lb][0]);
                    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.push_back(true);
                    ++vid;
                }
                if ((*filter)(i+lb,j+1,b[lb][1],ref))
                {
                    colsIndex.push_back(j+1);
                    colsValue.push_back(b[lb][1]);
                    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.push_back(true);
                    ++vid;
                }
                if ((*filter)(i+lb,j+2,b[lb][2],ref))
                {
                    colsIndex.push_back(j+2);
                    colsValue.push_back(b[lb][2]);
                    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.push_back(true);
                    ++vid;
                }
            }

            if (!keepEmptyRows && rowBegin.back() == vid)   // row was empty
            {
                rowIndex.pop_back();
                rowBegin.pop_back();
            }
        }
    }
    rowBegin.push_back(vid); // end of last row
}

template <> template <>
inline void CompressedRowSparseMatrixMechanical<double>::filterValues(CompressedRowSparseMatrixMechanical<defaulttype::Mat<3,3,float> >& M, filter_fn* filter, const Real ref, bool keepEmptyRows)
{
    M.compress();
    nRow = M.rowSize();
    nCol = M.colSize();
    nBlocRow = 1;
    nBlocCol = 1;
    rowIndex.clear();
    rowBegin.clear();
    colsIndex.clear();
    colsValue.clear();
    skipCompressZero = true;
    btemp.clear();
    rowIndex.reserve(M.rowIndex.size()*3);
    rowBegin.reserve(M.rowBegin.size()*3);
    colsIndex.reserve(M.colsIndex.size()*9);
    colsValue.reserve(M.colsValue.size()*9);

    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.clear();

    Index vid = 0;
    for (std::size_t rowId = 0; rowId < M.rowIndex.size(); ++rowId)
    {
        Index i = M.rowIndex[rowId] * 3;

        Range rowRange(M.rowBegin[rowId], M.rowBegin[rowId+1]);

        for (Index lb = 0; lb<3 ; lb++)
        {
            rowIndex.push_back(i+lb);
            rowBegin.push_back(vid);

            for (std::size_t xj = static_cast<std::size_t>(rowRange.begin()); xj < static_cast<std::size_t>(rowRange.end()); ++xj)
            {
                Index j = M.colsIndex[xj] * 3;
                defaulttype::Mat<3,3,double> b = M.colsValue[xj];
                if ((*filter)(i+lb,j+0,b[lb][0],ref))
                {
                    colsIndex.push_back(j+0);
                    colsValue.push_back(b[lb][0]);
                    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.push_back(true);
                    ++vid;
                }
                if ((*filter)(i+lb,j+1,b[lb][1],ref))
                {
                    colsIndex.push_back(j+1);
                    colsValue.push_back(b[lb][1]);
                    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.push_back(true);
                    ++vid;
                }
                if ((*filter)(i+lb,j+2,b[lb][2],ref))
                {
                    colsIndex.push_back(j+2);
                    colsValue.push_back(b[lb][2]);
                    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.push_back(true);
                    ++vid;
                }
            }

            if (!keepEmptyRows && rowBegin.back() == vid)   // row was empty
            {
                rowIndex.pop_back();
                rowBegin.pop_back();
            }
        }
    }
    rowBegin.push_back(vid); // end of last row
}

template <> template <>
inline void CompressedRowSparseMatrixMechanical<float>::filterValues(CompressedRowSparseMatrixMechanical<defaulttype::Mat<3,3,float> >& M, filter_fn* filter, const Real ref, bool keepEmptyRows)
{
    M.compress();
    nRow = M.rowSize();
    nCol = M.colSize();
    nBlocRow = 1;
    nBlocCol = 1;
    rowIndex.clear();
    rowBegin.clear();
    colsIndex.clear();
    colsValue.clear();
    skipCompressZero = true;
    btemp.clear();
    rowIndex.reserve(M.rowIndex.size()*3);
    rowBegin.reserve(M.rowBegin.size()*3);
    colsIndex.reserve(M.colsIndex.size()*9);
    colsValue.reserve(M.colsValue.size()*9);

    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.clear();

    Index vid = 0;
    for (std::size_t rowId = 0; rowId < M.rowIndex.size(); ++rowId)
    {
        Index i = M.rowIndex[rowId] * 3;

        Range rowRange(M.rowBegin[rowId], M.rowBegin[rowId+1]);

        for (Index lb = 0; lb<3 ; lb++)
        {
            rowIndex.push_back(i+lb);
            rowBegin.push_back(vid);

            for (std::size_t xj = static_cast<std::size_t>(rowRange.begin()); xj < static_cast<std::size_t>(rowRange.end()); ++xj)
            {
                Index j = M.colsIndex[xj] * 3;
                defaulttype::Mat<3,3,float> b = M.colsValue[xj];
                if ((*filter)(i+lb,j+0,b[lb][0],ref))
                {
                    colsIndex.push_back(j+0);
                    colsValue.push_back(b[lb][0]);
                    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.push_back(true);
                    ++vid;
                }
                if ((*filter)(i+lb,j+1,b[lb][1],ref))
                {
                    colsIndex.push_back(j+1);
                    colsValue.push_back(b[lb][1]);
                    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.push_back(true);
                    ++vid;
                }
                if ((*filter)(i+lb,j+2,b[lb][2],ref))
                {
                    colsIndex.push_back(j+2);
                    colsValue.push_back(b[lb][2]);
                    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.push_back(true);
                    ++vid;
                }
            }

            if (!keepEmptyRows && rowBegin.back() == vid)   // row was empty
            {
                rowIndex.pop_back();
                rowBegin.pop_back();
            }
        }
    }
    rowBegin.push_back(vid); // end of last row
}

template <> template <>
inline void CompressedRowSparseMatrixMechanical<float>::filterValues(CompressedRowSparseMatrixMechanical<defaulttype::Mat<3,3,double> >& M, filter_fn* filter, const Real ref, bool keepEmptyRows)
{
    M.compress();
    nRow = M.rowSize();
    nCol = M.colSize();
    nBlocRow = 1;
    nBlocCol = 1;
    rowIndex.clear();
    rowBegin.clear();
    colsIndex.clear();
    colsValue.clear();
    skipCompressZero = true;
    btemp.clear();
    rowIndex.reserve(M.rowIndex.size()*3);
    rowBegin.reserve(M.rowBegin.size()*3);
    colsIndex.reserve(M.colsIndex.size()*9);
    colsValue.reserve(M.colsValue.size()*9);

    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.clear();

    Index vid = 0;
    for (std::size_t rowId = 0; rowId < M.rowIndex.size(); ++rowId)
    {
        Index i = M.rowIndex[rowId] * 3;

        Range rowRange(M.rowBegin[rowId], M.rowBegin[rowId+1]);

        for (Index lb = 0; lb<3 ; lb++)
        {
            rowIndex.push_back(i+lb);
            rowBegin.push_back(vid);

            for (std::size_t xj = static_cast<std::size_t>(rowRange.begin()); xj < static_cast<std::size_t>(rowRange.end()); ++xj)
            {
                Index j = M.colsIndex[xj] * 3;
                defaulttype::Mat<3,3,float> b = M.colsValue[xj];
                if ((*filter)(i+lb,j+0,b[lb][0],ref))
                {
                    colsIndex.push_back(j+0);
                    colsValue.push_back(b[lb][0]);
                    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.push_back(true);
                    ++vid;
                }
                if ((*filter)(i+lb,j+1,b[lb][1],ref))
                {
                    colsIndex.push_back(j+1);
                    colsValue.push_back(b[lb][1]);
                    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.push_back(true);
                    ++vid;
                }
                if ((*filter)(i+lb,j+2,b[lb][2],ref))
                {
                    colsIndex.push_back(j+2);
                    colsValue.push_back(b[lb][2]);
                    //SOFA_IF_CONSTEXPR (Policy::StoreTouchFlags) touchedBloc.push_back(true);
                    ++vid;
                }
            }

            if (!keepEmptyRows && rowBegin.back() == vid)   // row was empty
            {
                rowIndex.pop_back();
                rowBegin.pop_back();
            }
        }
    }
    rowBegin.push_back(vid); // end of last row
}


} // namespace defaulttype

#endif
