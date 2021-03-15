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
#ifndef SOFA_DEFAULTTYPE_COMPRESSEDROWSPARSEMATRIXMECHANICAL_H
#define SOFA_DEFAULTTYPE_COMPRESSEDROWSPARSEMATRIXMECHANICAL_H

//#include <sofa/SofaFramework.h>
#include <SofaBaseLinearSolver/config.h>

#include <SofaBaseLinearSolver/matrix_bloc_traits.h>
#include <SofaBaseLinearSolver/CompressedRowSparseMatrix.h>
#include <sofa/defaulttype/BaseMatrix.h>
#include <SofaBaseLinearSolver/MatrixExpr.h>
#include <SofaBaseLinearSolver/FullVector.h>
#include <sofa/defaulttype/Mat.h>

namespace sofa::component::linearsolver
{

/// Mechanical policy type, showing the types and flags to give to CompressedRowSparseMatrixMechanical
/// for its second template type.
class CRSMechanicalPolicy : public CRSDefaultPolicy
{
public:
    static constexpr bool CompressZeros = true;
    static constexpr bool IsAlwaysSquare = true;
    static constexpr bool IsAlwaysSymmetric = true;
    static constexpr bool OrderedInsertion = true;
    static constexpr bool StoreLowerTriangularBloc = true;

    static constexpr int matrixType = 1;
};

class CRSMechanicalStoreTouchedPolicy : public CRSDefaultPolicy
{
public:
    static constexpr bool CompressZeros            = true;
    static constexpr bool IsAlwaysSquare           = true;
    static constexpr bool IsAlwaysSymmetric        = true;
    static constexpr bool OrderedInsertion         = true;
    static constexpr bool StoreLowerTriangularBloc = true;
    static constexpr bool StoreTouchFlags          = true;
    static constexpr int  matrixType               = 1;
};

template<typename TBloc, typename TPolicy = CRSMechanicalPolicy >
class CompressedRowSparseMatrixMechanical final // final is used to allow the compiler to inline virtual methods
    : public CompressedRowSparseMatrix<TBloc, TPolicy>, public sofa::defaulttype::BaseMatrix
{
public:
    typedef CompressedRowSparseMatrixMechanical<TBloc, TPolicy> Matrix;
    //SOFA_MATRIX_CLASS_EXTERNAL((Matrix),((defaulttype::BaseMatrix)));

    typedef CompressedRowSparseMatrix<TBloc, TPolicy> CRSMatrix;
    typedef typename CRSMatrix::Policy Policy;

    using Bloc     = TBloc;
    using VecBloc  = typename CRSBlocTraits<Bloc>::VecBloc;
    using VecIndex = typename CRSBlocTraits<Bloc>::VecIndex;
    using VecFlag  = typename CRSBlocTraits<Bloc>::VecFlag;
    using Index    = typename VecIndex::value_type;

    typedef typename CRSMatrix::Bloc Data;
    typedef typename CRSMatrix::Range Range;
    typedef typename CRSMatrix::traits traits;
    typedef typename CRSMatrix::Real Real;
    typedef typename CRSMatrix::Index KeyType;
    typedef typename CRSMatrix::IndexedBloc IndexedBloc;
    typedef typename CRSMatrix::VecIndexedBloc VecIndexedBloc;

    typedef Matrix Expr;
    typedef CompressedRowSparseMatrixMechanical<double> matrix_type;
    enum { category = MATRIX_SPARSE };
    enum { operand = 1 };

    enum { NL = CRSMatrix::NL };  ///< Number of rows of a block
    enum { NC = CRSMatrix::NC };  ///< Number of columns of a block

    /// Size
    Index nRow,nCol;         ///< Mathematical size of the matrix, in scalars
    static_assert(!Policy::AutoSize,
        "CompressedRowSparseMatrixMechanical cannot use AutoSize policy to make sure bloc-based and scalar-based sizes match");

    CompressedRowSparseMatrixMechanical()
        : CRSMatrix()
        , nRow(0), nCol(0)
    {
		//nRow = 0;
		//nCol = 0;
    }

    CompressedRowSparseMatrixMechanical(Index nbRow, Index nbCol)
        : CRSMatrix((nbRow + NL-1) / NL, (nbCol + NC-1) / NC)
        ,nRow(nbRow), nCol(nbCol)
    {
		//nRow = nbRow;
		//nCol = nbCol;
    }

    static void split_row_index(Index& index, Index& modulo) { bloc_index_func<NL>::split(index, modulo); }
    static void split_col_index(Index& index, Index& modulo) { bloc_index_func<NC>::split(index, modulo); }

    void compress() override
    {
        CRSMatrix::compress();
    }

    void swap(Matrix& m)
    {
        Index t;
        t = nRow; nRow = m.nRow; m.nRow = t;
        t = nCol; nCol = m.nCol; m.nCol = t;
        CRSMatrix::swap(m);
    }

    /// Make sure all diagonal entries are present even if they are zero
    template< typename = typename std::enable_if< Policy::IsAlwaysSquare> >
    void fullDiagonal()
    {
        if constexpr (Policy::LogTrace) this->logCall(FnEnum::fullDiagonal);
        compress();
        Index ndiag = 0;
        for (Index r = 0; r < static_cast<Index>(this->rowIndex.size()); ++r)
        {
            Index i = this->rowIndex[r];
            Index b = this->rowBegin[r];
            Index e = this->rowBegin[r+1];
            Index t = b;
            while (b < e && this->colsIndex[t] != i)
            {
                if (this->colsIndex[t] < i)
                    b = t+1;
                else
                    e = t;
                t = (b+e)>>1;
            }
            if (b<e) ++ndiag;
        }
        if (ndiag == this->nBlocRow) return;

        this->oldRowIndex.swap(this->rowIndex);
        this->oldRowBegin.swap(this->rowBegin);
        this->oldColsIndex.swap(this->colsIndex);
        this->oldColsValue.swap(this->colsValue);
        this->rowIndex.resize(this->nBlocRow);
        this->rowBegin.resize(this->nBlocRow + 1);
        this->colsIndex.resize(this->oldColsIndex.size() + this->nBlocRow-ndiag);
        this->colsValue.resize(this->oldColsValue.size() + this->nBlocRow-ndiag);
        if constexpr(Policy::StoreTouchFlags)
        {
            this->touchedBloc.resize(this->oldColsValue.size() + this->nBlocRow - ndiag);
            std::fill(this->touchedBloc.begin(), this->touchedBloc.end(), false);
        }
        Index nv = 0;
        for (Index i = 0; i < this->nBlocRow; ++i) this->rowIndex[i] = i;
        Index j = 0;
        for (Index i = 0; i < static_cast<Index>(this->oldRowIndex.size()); ++i)
        {
            for (; j < this->oldRowIndex[i]; ++j)
            {
                this->rowBegin[j] = nv;
                this->colsIndex[nv] = j;
                traits::clear(this->colsValue[nv]);
                ++nv;
            }
            this->rowBegin[j] = nv;
            Index b = this->oldRowBegin[i];
            Index e = this->oldRowBegin[i+1];
            for (; b < e && this->oldColsIndex[b] < j; ++b)
            {
                this->colsIndex[nv] = this->oldColsIndex[b];
                this->colsValue[nv] = this->oldColsValue[b];
                ++nv;
            }
            if (b >= e || this->oldColsIndex[b] > j)
            {
                this->colsIndex[nv] = j;
                traits::clear(this->colsValue[nv]);
                ++nv;
            }
            for (; b < e; ++b)
            {
                this->colsIndex[nv] = this->oldColsIndex[b];
                this->colsValue[nv] = this->oldColsValue[b];
                ++nv;
            }
            ++j;
        }
        for (; j < this->nBlocRow; ++j)
        {
            this->rowBegin[j] = nv;
            this->colsIndex[nv] = j;
            traits::clear(this->colsValue[nv]);
            ++nv;
        }
        this->rowBegin[j] = nv;
    }

    ///< Mathematical size of the matrix
    Index rowSize() const
    {
        return nRow;
    }

    ///< Mathematical size of the matrix
    Index colSize() const
    {
        return nCol;
    }

    /// This override classic resizeBloc to fill nRow and nCol values.
    void resizeBloc(Index nbBRow, Index nbBCol)
    {
        CRSMatrix::resizeBloc(nbBRow, nbBCol);
        nRow = NL * nbBRow;
        nCol = NC * nbBCol;
    }

    void resize(Index nbRow, Index nbCol)
    {
        if constexpr (Policy::Verbose)
        {
            if (nbRow != rowSize() || nbCol != colSize()) std::cout << this->Name()  << ": resize("<<nbRow<<","<<nbCol<<")"<<std::endl;
        }
        this->resizeBloc((nbRow + NL-1) / NL, (nbCol + NC-1) / NC);
        nRow = nbRow;
        nCol = nbCol;
    }

    void extend(Index nbRow, Index nbCol)
    {
        if constexpr(Policy::Check)
        {
            if (nbRow < rowSize() || nbCol < colSize() )
            {
                std::cerr << "ERROR: extend("<<nbRow<<","<<nbCol<<") reduces the size of the matrix " << this->Name() << " : size ("<<rowSize()<<","<<colSize()<<")" << std::endl;
            }
        }
        nRow = nbRow;
        nCol = nbCol;
        this->nBlocRow = (nbRow + NL-1) / NL;
        this->nBlocCol = (nbCol + NC-1) / NC;
    }

    /**
    * \brief get scalar element i, j of matrix
    **/
    SReal element(Index i, Index j) const override
    {
        if constexpr (Policy::AutoCompress) const_cast<Matrix*>(this)->compress(); /// \warning this violates the const-ness of the method !
        if constexpr (Policy::Check)
        {
            if (i >= rowSize() || j >= colSize())
            {
                std::cerr << "ERROR: invalid read access to element ("<<i<<","<<j<<") in "<< this->Name() << " of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
                return 0.0;
            }
        }
        if constexpr (!Policy::StoreLowerTriangularBloc) if ((i / NL) > (j / NC)) std::swap(i,j);

        Index bi=0, bj=0; split_row_index(i, bi); split_col_index(j, bj);
        return traits::v(this->bloc(i, j), bi, bj);
    }

    /**
    * \brief set scalar element i, j of matrix
    **/
    void set(Index i, Index j, double v) override
    {
        if constexpr (Policy::LogTrace) this->logCall(FnEnum::setVal, i, j, v);
        if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") = "<<v<<std::endl;
        if constexpr (Policy::Check)
        {
            if (i >= rowSize() || j >= colSize())
            {
                std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<< this->Name() << " of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
                return;
            }
        }
        if constexpr (!Policy::StoreLowerTriangularBloc) if ((i / NL) > (j / NC)) return;

        Index bi=0, bj=0; split_row_index(i, bi); split_col_index(j, bj);
        traits::vset(*this->wbloc(i, j, true), bi, bj, static_cast<Real>(v) );
        if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<this->rowBSize()<<"*"<<NL<<","<<this->colBSize()<<"*"<<NC<<"): bloc("<<i<<","<<j<<")["<<bi<<","<<bj<<"] = "<<v<<std::endl;
    }

    /**
    * \brief add scalar v at element i, j of matrix
    **/
    void add(Index i, Index j, double v) override
    {
        if constexpr (Policy::LogTrace) this->logCall(FnEnum::addVal, i, j, v);
        if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") += "<<v<<std::endl;
        if constexpr (Policy::Check)
        {
            if (i >= rowSize() || j >= colSize())
            {
                std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<< this->Name() << " of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
                return;
            }
        }
        if constexpr (!Policy::StoreLowerTriangularBloc) if ((i / NL) > (j / NC)) return;

        Index bi=0, bj=0; split_row_index(i, bi); split_col_index(j, bj);
        if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<this->rowBSize()<<"*"<<NL<<","<<this->colBSize()<<"*"<<NC<<"): bloc("<<i<<","<<j<<")["<<bi<<","<<bj<<"] += "<<v<<std::endl;
        traits::vadd(*this->wbloc(i,j,true), bi, bj, static_cast<Real>(v) );
    }

    /**
    * \brief set scalar element i, j of matrix when rowId and colId are known
    **/
    void set(Index i, Index j, int& rowId, int& colId, double v)
    {
        if constexpr (Policy::LogTrace) this->logCall(FnEnum::setValId, i, j, rowId, colId, v);
        if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") = "<<v<<std::endl;
        if constexpr (Policy::Check)
        {
            if (i >= rowSize() || j >= colSize())
            {
                std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<< this->Name() << " of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
                return;
            }
        }
        if constexpr (!Policy::StoreLowerTriangularBloc) if ((i / NL) > (j / NC)) return;

        Index bi=0, bj=0; split_row_index(i, bi); split_col_index(j, bj);
        traits::vset(*this->wbloc(i,j,rowId,colId,true), bi, bj, static_cast<Real>(v) );
        if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<this->rowBSize()<<"*"<<NL<<","<<this->colBSize()<<"*"<<NC<<"): bloc("<<i<<","<<j<<")["<<bi<<","<<bj<<"] = "<<v<<std::endl;
    }

    /**
    * \brief add scalar v at element i, j when rowId and colId are known
    **/
    void add(Index i, Index j, int& rowId, int& colId, double v)
    {
        if constexpr (Policy::LogTrace) this->logCall(FnEnum::addValId, i, j, rowId, colId, v);
        if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") += "<<v<<std::endl;
        if constexpr (Policy::Check)
        {
            if (i >= rowSize() || j >= colSize())
            {
                std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<< this->Name() << " of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
                return;
            }
        }
        if constexpr (!Policy::StoreLowerTriangularBloc) if ((i / NL) > (j / NC)) return;

        Index bi=0, bj=0; split_row_index(i, bi); split_col_index(j, bj);
        if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<this->rowBSize()<<"*"<<NL<<","<<this->colBSize()<<"*"<<NC<<"): bloc("<<i<<","<<j<<")["<<bi<<","<<bj<<"] += "<<v<<std::endl;
        traits::vadd(*this->wbloc(i,j,rowId,colId,true), bi, bj, static_cast<Real>(v) );
    }

    /**
    * \brief clear scalar at element i, j of matrix
    **/
    void clear(Index i, Index j)
    {
        if constexpr (Policy::LogTrace) this->logCall(FnEnum::clearIndex, i, j);
        if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<rowSize()<<","<<colSize()<<"): element("<<i<<","<<j<<") = 0"<<std::endl;
        if constexpr (Policy::Check)
        {
            if (i >= rowSize() || j >= colSize())
            {
                std::cerr << "ERROR: invalid write access to element ("<<i<<","<<j<<") in "<< this->Name() << " of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
                return;
            }
        }
        if constexpr (Policy::AutoCompress) this->compress();

        Index bi=0, bj=0; split_row_index(i, bi); split_col_index(j, bj);
        Bloc* b = this->wbloc(i,j,false);
        if (b) traits::vset(*b, bi, bj, 0);
    }

    /**
    * \brief Clear row scalar method. Clear all col of this line.
    * @param i : Line index considering size of matrix in scalar.
    * \warning If you want to clear all value of a bloc, it is better to call clearRowBloc
    **/
    void clearRow(Index i)
    {
        if constexpr (Policy::LogTrace) this->logCall(FnEnum::clearRow, i);
        if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<rowSize()<<","<<colSize()<<"): row("<<i<<") = 0"<<std::endl;
        if constexpr (Policy::Check)
        {
            if (i >= rowSize())
            {
                std::cerr << "ERROR: invalid write access to row "<<i<<" in "<< this->Name() << " of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
                return;
            }
        }
        if constexpr (Policy::AutoCompress) this->compress(); /// If AutoCompress policy is activated, we neeed to be sure not missing btemp registered value.

        Index bi=0; split_row_index(i, bi);
        Index rowId = i * this->rowIndex.size() / this->nBlocRow;
        if (this->sortedFind(this->rowIndex, i, rowId))
        {
            Range rowRange(this->rowBegin[rowId], this->rowBegin[rowId+1]);
            for (Index xj = rowRange.begin(); xj < rowRange.end(); ++xj)
            {
                Bloc& b = this->colsValue[xj];
                for (Index bj = 0; bj < NC; ++bj)
                    traits::vset(b, bi, bj, 0);
                if constexpr (Policy::StoreTouchFlags) this->touchedBloc[xj] = true;
            }
        }
    }

    /**
    * \brief Clear col scalar method. Clear this col in all row of matrix.
    * @param j : Col index considering size of matrix in scalar.
    * \warning If you want to clear all value of a bloc, it is better to call clearColBloc
    **/
    void clearCol(Index j)
    {
        if constexpr (Policy::LogTrace) this->logCall(FnEnum::clearCol, j);
        if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<rowSize()<<","<<colSize()<<"): col("<<j<<") = 0"<<std::endl;
        if constexpr (Policy::Check)
        {
            if (j >= colSize())
            {
                std::cerr << "ERROR: invalid write access to col "<<j<<" in "<< this->Name() << " of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
                return;
            }
        }
        /// If AutoCompress policy is activated, we neeed to be sure not missing btemp registered value.
        if constexpr (Policy::AutoCompress) this->compress();

        Index bj = 0; split_col_index(j, bj);
        for (Index i = 0; i < this->nBlocRow; ++i)
        {
            Bloc* b = this->wbloc(i,j,false);
            if (b)
            {
                for (Index bi = 0; bi < NL; ++bi)
                    traits::vset(*b, bi, bj, 0);
            }
        }
    }

    /**
    * \brief Clear both row i and column i in a square matrix
    * @param i : Row and Col index considering size of matrix in scalar.
    **/
    void clearRowCol(Index i)
    {
        if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<rowSize()<<","<<colSize()<<"): row("<<i<<") = 0 and col("<<i<<") = 0"<<std::endl;

        if constexpr (!Policy::IsAlwaysSquare || !Policy::StoreLowerTriangularBloc)
        {
            clearRow(i);
            clearCol(i);
        }
        else
        {
            if constexpr (Policy::Check)
            {
                if (i >= rowSize() || i >= colSize())
                {
                    std::cerr << "ERROR: invalid write access to row and column "<<i<<" in "<< this->Name() << " of size ("<<rowSize()<<","<<colSize()<<")"<<std::endl;
                    return;
                }
            }

            /// If AutoCompress policy is activated, we need to be sure that we are not missing btemp registered value.
            if constexpr (Policy::AutoCompress) this->compress();

            if constexpr (Policy::LogTrace) this->logCall(FnEnum::clearRowCol, i);

            Index bi=0; split_row_index(i, bi);
            Index rowId = i * this->rowIndex.size() / this->nBlocRow;
            if (this->sortedFind(this->rowIndex, i, rowId))
            {
                Range rowRange(this->rowBegin[rowId], this->rowBegin[rowId+1]);
                for (Index xj = rowRange.begin(); xj < rowRange.end(); ++xj)
                {
                    Bloc* b = &this->colsValue[xj];
                    // first clear (i,j)
                    for (Index bj = 0; bj < NC; ++bj)
                        traits::vset(*b, bi, bj, 0);

                    if constexpr(Policy::StoreTouchFlags) this->touchedBloc[xj] = true;

                    // then clear (j,i) 
                    Index j = this->colsIndex[xj];
                    
                    if (j != i)
                    {
                        Range jrowRange(this->rowBegin[j], this->rowBegin[j + 1]);
                        Index colId = 0;

                        // look for column i
                        if (this->sortedFind(this->colsIndex, jrowRange, i, colId))
                        {
                            if constexpr(Policy::StoreTouchFlags) this->touchedBloc[colId] = true;
                            b = &this->colsValue[colId];
                        }
                    }

                    for (Index bj = 0; bj < NL; ++bj)
                        traits::vset(*b, bj, bi, 0);
             
                }
            }
        }
    }

    /**
    * \brief Completely clear the matrix
    **/
    void clear() override /// Need implement clear to override BaseMatrix one.
    {
        CRSMatrix::clear();
    }

/// @name BlocMatrixWriter operators
/// @{
    /// Override CRSMatrix add method to avoid mis-understanding by compilator with other add method overriding BaseMatrix.
    void add(unsigned int bi, unsigned int bj, const Bloc& b)
    {
        CRSMatrix::add(bi, bj, b);
    }
    void add(unsigned int bi, unsigned int bj, int& rowId, int& colId, const Bloc& b)
    {
        CRSMatrix::add(bi, bj, rowId, colId, b);
    }
/// @}

/// @name Get information about the content and structure of this matrix (diagonal, band, sparse, full, block size, ...)
/// @{

    /// @return type of elements stored in this matrix
    virtual ElementType getElementType() const override { return traits::getElementType(); }

    /// @return size of elements stored in this matrix
    virtual std::size_t getElementSize() const override { return sizeof(Real); }

    /// @return the category of this matrix
    virtual MatrixCategory getCategory() const override { return MATRIX_SPARSE; }

    /// @return the number of rows in each block, or 1 of there are no fixed block size
    virtual Index getBlockRows() const override { return NL; }

    /// @return the number of columns in each block, or 1 of there are no fixed block size
    virtual Index getBlockCols() const override { return NC; }

    /// @return the number of rows of blocks
    virtual Index bRowSize() const override{ return this->rowBSize(); }

    /// @return the number of columns of blocks
    virtual Index bColSize() const override { return this->colBSize(); }

    /// @return the width of the band on each side of the diagonal (only for band matrices)
    virtual Index getBandWidth() const override { return NC-1; }

/// @}

/// @name Filtering-out part of a matrix
/// @{

    typedef bool filter_fn     (Index   i  , Index   j  , Bloc& val, const Real   ref  );
    static bool nonzeros(Index /*i*/, Index /*j*/, Bloc& val, const Real /*ref*/) { return (!traits::empty(val)); }
    static bool nonsmall(Index /*i*/, Index /*j*/, Bloc& val, const Real   ref  )
    {
        for (Index bi = 0; bi < NL; ++bi)
            for (Index bj = 0; bj < NC; ++bj)
                if (helper::rabs(traits::v(val, bi, bj)) >= ref) return true;
        return false;
    }
    static bool upper(Index   i  , Index   j  , Bloc& val, const Real /*ref*/)
    {
        if (NL>1 && i*NL == j*NC)
        {
            for (Index bi = 1; bi < NL; ++bi)
                for (Index bj = 0; bj < bi; ++bj)
                    traits::vset(val, bi, bj, 0);
        }
        return i*NL <= j*NC;
    }
    static bool lower(Index   i  , Index   j  , Bloc& val, const Real /*ref*/)
    {
        if (NL>1 && i*NL == j*NC)
        {
            for (Index bi = 0; bi < NL-1; ++bi)
                for (Index bj = bi+1; bj < NC; ++bj)
                    traits::vset(val, bi, bj, 0);
        }
        return i*NL >= j*NC;
    }
    static bool upper_nonzeros(Index   i  , Index   j  , Bloc& val, const Real   ref  ) { return upper(i,j,val,ref) && nonzeros(i,j,val,ref); }
    static bool lower_nonzeros(Index   i  , Index   j  , Bloc& val, const Real   ref  ) { return lower(i,j,val,ref) && nonzeros(i,j,val,ref); }
    static bool upper_nonsmall(Index   i  , Index   j  , Bloc& val, const Real   ref  ) { return upper(i,j,val,ref) && nonsmall(i,j,val,ref); }
    static bool lower_nonsmall(Index   i  , Index   j  , Bloc& val, const Real   ref  ) { return lower(i,j,val,ref) && nonsmall(i,j,val,ref); }

    template<class TMatrix>
    void filterValues(TMatrix& M, filter_fn* filter = &nonzeros, const Real ref = Real(), bool keepEmptyRows=false)
    {
        M.compress();
        this->nBlocRow = M.rowBSize();
        this->nBlocCol = M.colBSize();
        this->rowIndex.clear();
        this->rowBegin.clear();
        this->colsIndex.clear();
        this->colsValue.clear();
        this->skipCompressZero = true;
        this->btemp.clear();
        this->rowIndex.reserve(M.rowIndex.size());
        this->rowBegin.reserve(M.rowBegin.size());
        this->colsIndex.reserve(M.colsIndex.size());
        this->colsValue.reserve(M.colsValue.size());

        Index vid = 0;
        for (Index rowId = 0; rowId < static_cast<Index>(M.rowIndex.size()); ++rowId)
        {
            Index i = M.rowIndex[rowId];
            this->rowIndex.push_back(i);
            this->rowBegin.push_back(vid);
            Range rowRange(M.rowBegin[rowId], M.rowBegin[rowId+1]);
            for (Index xj = rowRange.begin(); xj < rowRange.end(); ++xj)
            {
                Index j = M.colsIndex[xj];
                Bloc b = M.colsValue[xj];
                if ((*filter)(i,j,b,ref))
                {
                    this->colsIndex.push_back(j);
                    this->colsValue.push_back(b);
                    ++vid;
                }
            }
            if (!keepEmptyRows && this->rowBegin.back() == vid) // row was empty
            {
                this->rowIndex.pop_back();
                this->rowBegin.pop_back();
            }
        }
        this->rowBegin.push_back(vid); // end of last row
    }

    template <class TMatrix>
    void copyNonZeros(TMatrix& M, bool keepEmptyRows=false)
    {
        filterValues(M, nonzeros, Real(), keepEmptyRows);
    }

    template <class TMatrix>
    void copyNonSmall(TMatrix& M, const Real ref, bool keepEmptyRows=false)
    {
        filterValues(M, nonsmall, ref, keepEmptyRows);
    }

    void copyUpper(Matrix& M, bool keepEmptyRows=false)
    {
        filterValues(M, upper, Real(), keepEmptyRows);
    }

    void copyLower(Matrix& M, bool keepEmptyRows=false)
    {
        filterValues(M, lower, Real(), keepEmptyRows);
    }

    template <class TMatrix>
    void copyUpperNonZeros(TMatrix& M, bool keepEmptyRows=false)
    {
        filterValues(M, upper_nonzeros, Real(), keepEmptyRows);
    }

    template <class TMatrix>
    void copyLowerNonZeros(TMatrix& M, bool keepEmptyRows=false)
    {
        filterValues(M, lower_nonzeros, Real(), keepEmptyRows);
    }

    void copyUpperNonSmall(Matrix& M, const Real ref, bool keepEmptyRows=false)
    {
        filterValues(M, upper_nonsmall, ref, keepEmptyRows);
    }

    void copyLowerNonSmall(Matrix& M, const Real ref, bool keepEmptyRows=false)
    {
        filterValues(M, lower_nonsmall, ref, keepEmptyRows);
    }

/// @}

/// @name Virtual iterator classes and methods
/// @{

protected:
    virtual void bAccessorDelete(const InternalBlockAccessor* /*b*/) const {}
    virtual void bAccessorCopy(InternalBlockAccessor* /*b*/) const {}
    virtual SReal bAccessorElement(const InternalBlockAccessor* b, Index i, Index j) const
    {
        //return element(b->row * getBlockRows() + i, b->col * getBlockCols() + j);
        Index index = b->data;
        const Bloc& data = (index >= 0) ? this->colsValue[index] : this->btemp[-index-1].value;
        return static_cast<SReal>(traits::v(data, i, j));
    }
    virtual void bAccessorSet(InternalBlockAccessor* b, Index i, Index j, double v)
    {
        //set(b->row * getBlockRows() + i, b->col * getBlockCols() + j, v);
        Index index = b->data;
        Bloc& data = (index >= 0) ? this->colsValue[index] : this->btemp[-index-1].value;
        traits::vset(data, i, j, static_cast<Real>(v) );
    }
    virtual void bAccessorAdd(InternalBlockAccessor* b, Index i, Index j, double v)
    {
        //add(b->row * getBlockRows() + i, b->col * getBlockCols() + j, v);
        Index index = b->data;
        Bloc& data = (index >= 0) ? this->colsValue[index] : this->btemp[-index-1].value;
        traits::vadd(data, i, j, static_cast<Real>(v) );
    }

    template<class T>
    const T* bAccessorElementsCSRImpl(const InternalBlockAccessor* b, T* buffer) const
    {
        Index index = b->data;
        const Bloc& data = (index >= 0) ? this->colsValue[index] : this->btemp[-index-1].value;
        for (Index l=0; l<NL; ++l)
            for (Index c=0; c<NC; ++c)
                buffer[l*NC+c] = static_cast<T>(traits::v(data, l, c));
        return buffer;
    }
    virtual const float* bAccessorElements(const InternalBlockAccessor* b, float* buffer) const
    {
        return bAccessorElementsCSRImpl<float>(b, buffer);
    }
    virtual const double* bAccessorElements(const InternalBlockAccessor* b, double* buffer) const
    {
        return bAccessorElementsCSRImpl<double>(b, buffer);
    }
    virtual const int* bAccessorElements(const InternalBlockAccessor* b, int* buffer) const
    {
        return bAccessorElementsCSRImpl<int>(b, buffer);
    }

    template<class T>
    void bAccessorSetCSRImpl(InternalBlockAccessor* b, const T* buffer)
    {
        Index index = b->data;
        Bloc& data = (index >= 0) ? this->colsValue[index] : this->btemp[-index-1].value;
        for (Index l=0; l<NL; ++l)
            for (Index c=0; c<NC; ++c)
                traits::vset(data, l, c, static_cast<Real>(buffer[l*NC+c]) );
    }
    virtual void bAccessorSet(InternalBlockAccessor* b, const float* buffer)
    {
        bAccessorSetCSRImpl<float>(b, buffer);
    }
    virtual void bAccessorSet(InternalBlockAccessor* b, const double* buffer)
    {
        bAccessorSetCSRImpl<double>(b, buffer);
    }
    virtual void bAccessorSet(InternalBlockAccessor* b, const int* buffer)
    {
        bAccessorSetCSRImpl<int>(b, buffer);
    }

    template<class T>
    void bAccessorAddCSRImpl(InternalBlockAccessor* b, const T* buffer)
    {
        Index index = b->data;
        Bloc& data = (index >= 0) ? this->colsValue[index] : this->btemp[-index-1].value;
        for (Index l=0; l<NL; ++l)
            for (Index c=0; c<NC; ++c)
                traits::vadd(data, l, c,static_cast<Real>(buffer[l*NC+c]) );
    }
    virtual void bAccessorAdd(InternalBlockAccessor* b, const float* buffer)
    {
        bAccessorAddCSRImpl<float>(b, buffer);
    }
    virtual void bAccessorAdd(InternalBlockAccessor* b, const double* buffer)
    {
        bAccessorAddCSRImpl<double>(b, buffer);
    }
    virtual void bAccessorAdd(InternalBlockAccessor* b, const int* buffer)
    {
        bAccessorAddCSRImpl<int>(b, buffer);
    }

public:

    /// Get read access to a bloc
    virtual BlockConstAccessor blocGet(Index i, Index j) const
    {
        if constexpr (Policy::AutoCompress) const_cast<Matrix*>(this)->compress(); /// \warning this violates the const-ness of the method !

        Index rowId = i * this->rowIndex.size() / this->nBlocRow;
        if (this->sortedFind(this->rowIndex, i, rowId))
        {
            Range rowRange(this->rowBegin[rowId], this->rowBegin[rowId+1]);
            Index colId = rowRange.begin() + j * rowRange.size() / this->nBlocCol;
            if (this->sortedFind(this->colsIndex, rowRange, j, colId))
            {
                return createBlockConstAccessor(i, j, colId);
            }
        }
        return createBlockConstAccessor(-1-i, -1-j, static_cast<Index>(0));
    }

    /// Get write access to a bloc
    virtual BlockAccessor blocGetW(Index i, Index j)
    {
        if constexpr (Policy::AutoCompress) compress();

        Index rowId = i * this->rowIndex.size() / this->nBlocRow;
        if (this->sortedFind(this->rowIndex, i, rowId))
        {
            Range rowRange(this->rowBegin[rowId], this->rowBegin[rowId+1]);
            Index colId = rowRange.begin() + j * rowRange.size() / this->nBlocCol;
            if (this->sortedFind(this->colsIndex, rowRange, j, colId))
            {
                return createBlockAccessor(i, j, colId);
            }
        }
        return createBlockAccessor(-1-i, -1-j, static_cast<Index>(0));
    }

    /// Get write access to a bloc, possibly creating it
    virtual BlockAccessor blocCreate(Index i, Index j)
    {
        Index rowId = i * this->rowIndex.size() / this->nBlocRow;
        if (this->sortedFind(this->rowIndex, i, rowId))
        {
            Range rowRange(this->rowBegin[rowId], this->rowBegin[rowId+1]);
            Index colId = rowRange.begin() + j * rowRange.size() / this->nBlocCol;
            if (this->sortedFind(this->colsIndex, rowRange, j, colId))
            {
                if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<this->rowBSize()<<"*"<<NL<<","<<this->colBSize()<<"*"<<NC<<"): bloc("<<i<<","<<j<<") found at "<<colId<<" (line "<<rowId<<")."<<std::endl;
                return createBlockAccessor(i, j, colId);
            }
        }
        if (this->btemp.empty() || this->btemp.back().l != i || this->btemp.back().c != j)
        {
            if constexpr (Policy::Verbose) std::cout << this->Name()  << "("<<rowSize()<<","<<colSize()<<"): new temp bloc ("<<i<<","<<j<<")"<<std::endl;
            this->btemp.push_back(IndexedBloc(i,j));
            traits::clear(this->btemp.back().value);
        }
        return createBlockAccessor(i, j, -static_cast<Index>(this->btemp.size()));
    }

protected:
    virtual void itCopyColBlock(InternalColBlockIterator* /*it*/) const {}
    virtual void itDeleteColBlock(const InternalColBlockIterator* /*it*/) const {}
    virtual void itAccessColBlock(InternalColBlockIterator* it, BlockConstAccessor* b) const
    {
        Index index = it->data;
        setMatrix(b);
        getInternal(b)->row = it->row;
        getInternal(b)->data = index;
        getInternal(b)->col = this->colsIndex[index];
    }
    virtual void itIncColBlock(InternalColBlockIterator* it) const
    {
        Index index = it->data;
        ++index;
        it->data = index;
    }
    virtual void itDecColBlock(InternalColBlockIterator* it) const
    {
        Index index = it->data;
        --index;
        it->data = index;
    }
    virtual bool itEqColBlock(const InternalColBlockIterator* it, const InternalColBlockIterator* it2) const
    {
        Index index = it->data;
        Index index2 = it2->data;
        return index == index2;
    }
    virtual bool itLessColBlock(const InternalColBlockIterator* it, const InternalColBlockIterator* it2) const
    {
        Index index = it->data;
        Index index2 = it2->data;
        return index < index2;
    }

public:
    /// Get the iterator corresponding to the beginning of the given row of blocks
    virtual ColBlockConstIterator bRowBegin(Index ib) const
    {
        if constexpr (Policy::AutoCompress) const_cast<Matrix*>(this)->compress(); /// \warning this violates the const-ness of the method !
        Index rowId = ib * this->rowIndex.size() / this->nBlocRow;
        Index index = 0;
        if (this->sortedFind(this->rowIndex, ib, rowId))
        {
            index = this->rowBegin[rowId];
        }
        return createColBlockConstIterator(ib, index);
    }

    /// Get the iterator corresponding to the end of the given row of blocks
    virtual ColBlockConstIterator bRowEnd(Index ib) const
    {
        if constexpr (Policy::AutoCompress) const_cast<Matrix*>(this)->compress(); /// \warning this violates the const-ness of the method !
        Index rowId = ib * this->rowIndex.size() / this->nBlocRow;
        Index index2 = 0;
        if (this->sortedFind(this->rowIndex, ib, rowId))
        {
            index2 = this->rowBegin[rowId+1];
        }
        return createColBlockConstIterator(ib, index2);
    }

    /// Get the iterators corresponding to the beginning and end of the given row of blocks
    virtual std::pair<ColBlockConstIterator, ColBlockConstIterator> bRowRange(Index ib) const
    {
        if constexpr (Policy::AutoCompress) const_cast<Matrix*>(this)->compress(); /// \warning this violates the const-ness of the method !
        Index rowId = ib * this->rowIndex.size() / this->nBlocRow;
        Index index = 0, index2 = 0;
        if (this->sortedFind(this->rowIndex, ib, rowId))
        {
            index = this->rowBegin[rowId];
            index2 = this->rowBegin[rowId+1];
        }
        return std::make_pair(createColBlockConstIterator(ib, index ),
                createColBlockConstIterator(ib, index2));
    }


protected:
    virtual void itCopyRowBlock(InternalRowBlockIterator* /*it*/) const {}
    virtual void itDeleteRowBlock(const InternalRowBlockIterator* /*it*/) const {}
    virtual Index itAccessRowBlock(InternalRowBlockIterator* it) const
    {
        Index rowId = it->data[0];
        return this->rowIndex[rowId];
    }
    virtual ColBlockConstIterator itBeginRowBlock(InternalRowBlockIterator* it) const
    {
        Index rowId = it->data[0];
        Index row = this->rowIndex[rowId];
        Index index = this->rowBegin[rowId];
        return createColBlockConstIterator(row, index);
    }
    virtual ColBlockConstIterator itEndRowBlock(InternalRowBlockIterator* it) const
    {
        Index rowId = it->data[0];
        Index row = this->rowIndex[rowId];
        Index index2 = this->rowBegin[rowId+1];
        return createColBlockConstIterator(row, index2);
    }
    virtual std::pair<ColBlockConstIterator, ColBlockConstIterator> itRangeRowBlock(InternalRowBlockIterator* it) const
    {
        Index rowId = it->data[0];
        Index row = this->rowIndex[rowId];
        Index index = this->rowBegin[rowId];
        Index index2 = this->rowBegin[rowId+1];
        return std::make_pair(createColBlockConstIterator(row, index ),
                createColBlockConstIterator(row, index2));
    }

    virtual void itIncRowBlock(InternalRowBlockIterator* it) const
    {
        Index rowId = it->data[0];
        ++rowId;
        it->data[0] = rowId;
    }
    virtual void itDecRowBlock(InternalRowBlockIterator* it) const
    {
        Index rowId = it->data[0];
        --rowId;
        it->data[0] = rowId;
    }
    virtual bool itEqRowBlock(const InternalRowBlockIterator* it, const InternalRowBlockIterator* it2) const
    {
        Index rowId = it->data[0];
        Index rowId2 = it2->data[0];
        return rowId == rowId2;
    }
    virtual bool itLessRowBlock(const InternalRowBlockIterator* it, const InternalRowBlockIterator* it2) const
    {
        Index rowId = it->data[0];
        Index rowId2 = it2->data[0];
        return rowId < rowId2;
    }

public:
    /// Get the iterator corresponding to the beginning of the rows of blocks
    virtual RowBlockConstIterator bRowsBegin() const
    {
        if constexpr (Policy::AutoCompress) const_cast<Matrix*>(this)->compress(); /// \warning this violates the const-ness of the method !
        return createRowBlockConstIterator(0, 0);
    }

    /// Get the iterator corresponding to the end of the rows of blocks
    virtual RowBlockConstIterator bRowsEnd() const
    {
        if constexpr (Policy::AutoCompress) const_cast<Matrix*>(this)->compress(); /// \warning this violates the const-ness of the method !
        return createRowBlockConstIterator(this->rowIndex.size(), 0);
    }

    /// Get the iterators corresponding to the beginning and end of the given row of blocks
    virtual std::pair<RowBlockConstIterator, RowBlockConstIterator> bRowsRange() const
    {
        if constexpr (Policy::AutoCompress) const_cast<Matrix*>(this)->compress(); /// \warning this violates the const-ness of the method !
        return std::make_pair(createRowBlockConstIterator(0, 0),
                createRowBlockConstIterator(this->rowIndex.size(), 0));
    }

/// @}

protected:

/// @name setter/getter & product methods on template vector types
/// @{

    template<class Vec> static Real vget(const Vec& vec, Index i, Index j, Index k) { return vget( vec, i*j+k ); }
    template<class Vec> static Real vget(const helper::vector<Vec>&vec, Index i, Index /*j*/, Index k) { return vec[i][k]; }

                          static Real  vget(const defaulttype::BaseVector& vec, Index i) { return vec.element(i); }
    template<class Real2> static Real2 vget(const FullVector<Real2>& vec, Index i) { return vec[i]; }


    template<class Vec> static void vset(Vec& vec, Index i, Index j, Index k, Real v) { vset( vec, i*j+k, v ); }
    template<class Vec> static void vset(helper::vector<Vec>&vec, Index i, Index /*j*/, Index k, Real v) { vec[i][k] = v; }

                          static void vset(defaulttype::BaseVector& vec, Index i, Real v) { vec.set(i, v); }
    template<class Real2> static void vset(FullVector<Real2>& vec, Index i, Real2 v) { vec[i] = v; }


    template<class Vec> static void vadd(Vec& vec, Index i, Index j, Index k, Real v) { vadd( vec, i*j+k, v ); }
    template<class Vec> static void vadd(helper::vector<Vec>&vec, Index i, Index /*j*/, Index k, Real v) { vec[i][k] += v; }

                          static void vadd(defaulttype::BaseVector& vec, Index i, Real v) { vec.add(i, v); }
    template<class Real2> static void vadd(FullVector<Real2>& vec, Index i, Real2 v) { vec[i] += v; }

    template<class Vec> static void vresize(Vec& vec, Index /*blockSize*/, Index totalSize) { vec.resize( totalSize ); }
    template<class Vec> static void vresize(helper::vector<Vec>&vec, Index blockSize, Index /*totalSize*/) { vec.resize( blockSize ); }


    /** Product of the matrix with a templated vector res = this * vec*/
    template<class Real2, class V1, class V2>
    void tmul(V1& res, const V2& vec) const
    {
        assert( vec.size() % bColSize() == 0 ); // vec.size() must be a multiple of block size.

        if constexpr (Policy::AutoCompress) const_cast<Matrix*>(this)->compress(); /// \warning this violates the const-ness of the method !

        vresize( res, this->rowBSize(), rowSize() );

        for (Index xi = 0; xi < static_cast<Index>(this->rowIndex.size()); ++xi)  // for each non-empty block row
        {
            Range rowRange(this->rowBegin[xi], this->rowBegin[xi+1]);
            for (Index xj = rowRange.begin(); xj < rowRange.end(); ++xj)
            {
                sofa::defaulttype::Vec<NL,Real2> vi;
                const Bloc& b = this->colsValue[xj];
                Index rowIndex = this->rowIndex[xi] * NL;
                Index colIndex = this->colsIndex[xj] * NC;
                std::copy(vec.begin() + colIndex, vec.begin() + colIndex + NC, vi.begin());
                for (Index bi = 0; bi < NL; ++bi)
                    for (Index bj = 0; bj < NC; ++bj)
                        res[rowIndex + bi] += traits::v(b, bi, bj) * vi[bj];

                if constexpr (!Policy::StoreLowerTriangularBloc)
                {
                    if (colIndex != rowIndex)
                    {
                        sofa::defaulttype::Vec<NL,Real2> vj;
                        std::copy(vec.begin() + rowIndex, vec.begin() + rowIndex + NL, vj.begin());
                        for (Index bi = 0; bi < NL; ++bi)
                            for (Index bj = 0; bj < NC; ++bj)
                                res[colIndex + bi] += traits::v(b, bj, bi) * vj[bj];
                    }
                }
            }
        }
    }


    /** Product of the matrix with a templated vector res += this * vec*/
    template<class Real2, class V1, class V2>
    void taddMul(V1& res, const V2& vec) const
    {
        assert( vec.size()%bColSize() == 0 ); // vec.size() must be a multiple of block size.

        if constexpr (Policy::AutoCompress) const_cast<Matrix*>(this)->compress(); /// \warning this violates the const-ness of the method !
        vresize( res, this->rowBSize(), rowSize() );

        for (Index xi = 0; xi < static_cast<Index>(this->rowIndex.size()); ++xi)
        {
            Range rowRange(this->rowBegin[xi], this->rowBegin[xi+1]);
            for (Index xj = rowRange.begin(); xj < rowRange.end(); ++xj)
            {
                sofa::defaulttype::Vec<NL,Real2> vi;
                const Bloc& b = this->colsValue[xj];
                Index rowIndex = this->rowIndex[xi] * NL;
                Index colIndex = this->colsIndex[xj] * NC;
                std::copy(vec.begin() + colIndex, vec.begin() + colIndex + NC, vi.begin());
                for (Index bi = 0; bi < NL; ++bi)
                    for (Index bj = 0; bj < NC; ++bj)
                        res[rowIndex + bi] += traits::v(b, bi, bj) * vi[bj];

                if constexpr (!Policy::StoreLowerTriangularBloc)
                {
                    if (colIndex != rowIndex)
                    {
                        sofa::defaulttype::Vec<NL,Real2> vj;
                        std::copy(vec.begin() + rowIndex, vec.begin() + rowIndex + NL, vj.begin());
                        for (Index bi = 0; bi < NL; ++bi)
                            for (Index bj = 0; bj < NC; ++bj)
                                res[colIndex + bi] += traits::v(b, bj, bi) * vj[bj];
                    }
                }
            }
        }
    }


    /** Product of the matrix with a templated vector that have the size of the bloc res += this * [vec,...,vec]^T */
    template<class Real2, class V1, class V2>
    void taddMul_by_line(V1& res, const V2& vec) const
    {
        assert( vec.size() == NC ); // vec.size() must have the block size.

        if constexpr (Policy::AutoCompress) const_cast<Matrix*>(this)->compress(); /// \warning this violates the const-ness of the method !
        vresize( res, this->rowBSize(), rowSize() );

        for (Index xi = 0; xi < static_cast<Index>(this->rowIndex.size()); ++xi)
        {
            Range rowRange(this->rowBegin[xi], this->rowBegin[xi+1]);
            for (Index xj = rowRange.begin(); xj < rowRange.end(); ++xj)
            {
                const Bloc& b = this->colsValue[xj];
                Index rowIndex = this->rowIndex[xi] * NL;
                Index colIndex = this->colsIndex[xj] * NC;
                for (Index bi = 0; bi < NL; ++bi)
                    for (Index bj = 0; bj < NC; ++bj)
                        res[rowIndex + bi] += traits::v(b, bi, bj) * vec[bj];

                if constexpr (!Policy::StoreLowerTriangularBloc)
                {
                    if (colIndex != rowIndex)
                    {
                        for (Index bi = 0; bi < NL; ++bi)
                            for (Index bj = 0; bj < NC; ++bj)
                                res[colIndex + bi] += traits::v(b, bj, bi) * vec[bj];
                    }
                }
            }
        }
    }

    /** Product of the transpose with a templated vector and add it to res   res += this^T * vec */
    template<class Real2, class V1, class V2>
    void taddMulTranspose(V1& res, const V2& vec) const
    {
        assert( vec.size()%bRowSize() == 0 ); // vec.size() must be a multiple of block size.

        if constexpr (Policy::AutoCompress) const_cast<Matrix*>(this)->compress(); /// \warning this violates the const-ness of the method !
        vresize( res, this->colBSize(), colSize() );

        if constexpr (Policy::IsAlwaysSymetric) /// In symetric case this^T = this
        {
            taddMul(res, vec);
            return;
        }

        for (Index xi = 0; xi < this->rowIndex.size(); ++xi) // for each non-empty block row (i.e. column of the transpose)
        {
            // copy the corresponding chunk of the input to a local vector
            defaulttype::Vec<NL,Real2> v;
            //Index iN = rowIndex[xi] * NL;    // index of the row in the vector
            for (Index bi = 0; bi < NL; ++bi)
                v[bi] = vget(vec, this->rowIndex[xi], NL, bi);

            // accumulate the product of the column with the local vector
            Range rowRange(this->rowBegin[xi], this->rowBegin[xi+1]);
            for (Index xj = rowRange.begin(); xj < rowRange.end(); ++xj) // for each non-empty block in the row
            {
                const Bloc& b = this->colsValue[xj]; // non-empty block

                defaulttype::Vec<NC,Real2> r;  // local vector to store the product
                //Index jN = colsIndex[xj] * NC;

                // columnwise bloc-vector product
                for (Index bj = 0; bj < NC; ++bj)
                    r[bj] = traits::v(b, 0, bj) * v[0];
                for (Index bi = 1; bi < NL; ++bi)
                    for (Index bj = 0; bj < NC; ++bj)
                        r[bj] += traits::v(b, bi, bj) * v[bi];

                // accumulate the product to the result
                for (Index bj = 0; bj < NC; ++bj)
                    vadd(res, this->colsIndex[xj], NC, bj, r[bj]);
            }
        }
    }

/// @}


public:

/// @name Matrix operators
/// @{

    /// equal result = this * v
    /// @warning The block sizes must be compatible ie v.size() must be a multiple of block size.
    template< typename V1, typename V2 >
    void mul( V2& result, const V1& v ) const
    {
        this-> template tmul< Real, V2, V1 >(result, v);
    }


    /// equal result += this^T * v
    /// @warning The block sizes must be compatible ie v.size() must be a multiple of block size.
    template< typename V1, typename V2 >
    void addMultTranspose( V1& result, const V2& v ) const
    {
        this-> template taddMulTranspose< Real, V1, V2 >(result, v);
    }

    /// @returns this * v
    /// @warning The block sizes must be compatible ie v.size() must be a multiple of block size.
    template<class Vec>
    Vec operator*(const Vec& v) const
    {
        Vec res;
        mul( res, v );
        return res;
    }

    /// result += this * (v,...,v)^T
    /// v has the size of one bloc
    template< typename V, typename Real2 >
    void addMul_by_line( V& res, const defaulttype::Vec<NC,Real2>& v ) const
    {
        this-> template taddMul_by_line< Real2,V,defaulttype::Vec<NC,Real2> >( res, v );
    }
    template< typename Real, typename V, typename V2 >
    void addMul_by_line( V& res, const V2& v ) const
    {
        this-> template taddMul_by_line< Real,V,V2 >( res, v );
    }

    /// result += this * v
    template< typename V1, typename V2 >
    void addMul( V1& res, const V2& v ) const
    {
        taddMul< Real,V1,V2 >( res, v );
    }

    /// @}

    // methods for MatrixExpr support

    template<class M2>
    bool hasRef(const M2* m) const
    {
        return const_cast<void*>(this) == const_cast<void*>(m);
    }

    std::string expr() const
    {
        return std::string(Name());
    }

    bool valid() const
    {
        return true;
    }


    /// dest += this
    /// different bloc types possible
    /// @todo how to optimize when same bloc types
    template<class Dest>
    void addTo(Dest* dest) const
    {
        for (Index xi = 0; xi < static_cast<Index>(this->rowIndex.size()); ++xi)
        {
            Index iN = this->rowIndex[xi] * NL;
            Range rowRange(this->rowBegin[xi], this->rowBegin[xi+1]);
            for (Index xj = rowRange.begin(); xj < rowRange.end(); ++xj)
            {
                Index jN = this->colsIndex[xj] * NC;
                const Bloc& b = this->colsValue[xj];
                for (Index bi = 0; bi < NL; ++bi)
                    for (Index bj = 0; bj < NC; ++bj)
                        dest->add(iN+bi, jN+bj, traits::v(b, bi, bj));
            }
        }
        if (!this->btemp.empty())
        {
            for (typename VecIndexedBloc::const_iterator it = this->btemp.begin(), itend = this->btemp.end(); it != itend; ++it)
            {
                Index iN = it->l * NL;
                Index jN = it->c * NC;
                const Bloc& b = it->value;
                for (Index bi = 0; bi < NL; ++bi)
                    for (Index bj = 0; bj < NC; ++bj)
                        dest->add(iN+bi, jN+bj, traits::v(b, bi, bj));
            }
        }
    }

protected:

    /// add ? this += m : this = m
    /// m can be the same as this
    template<class M>
    void equal( const M& m, bool add = false )
    {
        if (m.hasRef(this))
        {
            Matrix tmp;
            tmp.resize(m.rowSize(), m.colSize());
            m.addTo(&tmp);
            if (add)
                tmp.addTo(this);
            else
                swap(tmp);
        }
        else
        {
            if (!add)
                resize(m.rowSize(), m.colSize());
            m.addTo(this);
        }
    }

    /// this += m
    template<class M>
    inline void addEqual( const M& m )
    {
        equal( m, true );
    }



public:

    template<class TBloc2, class TPolicy2>
    void operator=(const CompressedRowSparseMatrixMechanical<TBloc2, TPolicy2>& m)
    {
        if (&m == this) return;
        resize(m.rowSize(), m.colSize());
        m.addTo(this);
    }

    template<class TBloc2, class TPolicy2>
    void operator+=(const CompressedRowSparseMatrixMechanical<TBloc2, TPolicy2>& m)
    {
        addEqual(m);
    }

    template<class TBloc2, class TPolicy2>
    void operator-=(const CompressedRowSparseMatrixMechanical<TBloc2, TPolicy2>& m)
    {
        equal(MatrixExpr< MatrixNegative< CompressedRowSparseMatrixMechanical<TBloc2, TPolicy2> > >(MatrixNegative< CompressedRowSparseMatrixMechanical<TBloc2, TPolicy2> >(m)), true);
    }

    template<class Expr2>
    void operator=(const MatrixExpr< Expr2 >& m)
    {
        equal(m, false);
    }

    template<class Expr2>
    void operator+=(const MatrixExpr< Expr2 >& m)
    {
        addEqual(m);
    }

    template<class Expr2>
    void operator-=(const MatrixExpr< Expr2 >& m)
    {
        addEqual(MatrixExpr< MatrixNegative< Expr2 > >(MatrixNegative< Expr2 >(m)));
    }

    MatrixExpr< MatrixTranspose< Matrix > > t() const
    {
        return MatrixExpr< MatrixTranspose< Matrix > >(MatrixTranspose< Matrix >(*this));
    }


    MatrixExpr< MatrixNegative< Matrix > > operator-() const
    {
        return MatrixExpr< MatrixNegative< Matrix > >(MatrixNegative< Matrix >(*this));
    }

    MatrixExpr< MatrixScale< Matrix, double > > operator*(const double& r) const
    {
        return MatrixExpr< MatrixScale< Matrix, double > >(MatrixScale< Matrix, double >(*this, r));
    }


    static const char* Name()
    {
        // Note: to preserve backward compatibility, CompressedRowSparseMatrixMechanical keeps the same
        // name as CompressedRowSparseMatrix. We could change it later but it requires either being
        // sure all old code/scenes are updated, or add an alias mechanism in template names.
        return CRSMatrix::Name();
        // static std::string name = std::string("CompressedRowSparseMatrixMechanical") + std::string(traits::Name());
        // return name.c_str();
    }
};


//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<float>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat1x1f>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat2x2f>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat3x3f>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat4x4f>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat<6, 6, float> >), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat<8, 8, float> >), SOFA_DEFAULTTYPE_API);
//
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<double>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat1x1d>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat2x2d>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat3x3d>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat4x4d>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat<6, 6, double> >), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat<8, 8, double> >), SOFA_DEFAULTTYPE_API);
//
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<float, CRSMechanicalStoreTouchedPolicy>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat1x1f, CRSMechanicalStoreTouchedPolicy>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat2x2f, CRSMechanicalStoreTouchedPolicy>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat3x3f, CRSMechanicalStoreTouchedPolicy>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat4x4f, CRSMechanicalStoreTouchedPolicy>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat<6, 6, float>, CRSMechanicalStoreTouchedPolicy >), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat<8, 8, float>, CRSMechanicalStoreTouchedPolicy >), SOFA_DEFAULTTYPE_API);
//
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<double, CRSMechanicalStoreTouchedPolicy>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat1x1d, CRSMechanicalStoreTouchedPolicy>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat2x2d, CRSMechanicalStoreTouchedPolicy>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat3x3d, CRSMechanicalStoreTouchedPolicy>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat4x4d, CRSMechanicalStoreTouchedPolicy>), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat<6, 6, double>, CRSMechanicalStoreTouchedPolicy >), SOFA_DEFAULTTYPE_API);
//SOFA_TEMPLATE_MATRIX_CLASS_DECL_EXPORT((CompressedRowSparseMatrixMechanical<Mat<8, 8, double>, CRSMechanicalStoreTouchedPolicy >), SOFA_DEFAULTTYPE_API);


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_BASE_LINEAR_SOLVER) 
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<float>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat1x1f>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat2x2f>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat3x3f>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat4x4f>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat<6, 6, float> >;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat<8, 8, float> >;

extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<double>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat1x1d>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat2x2d>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat3x3d>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat4x4d>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat<6, 6, double> >;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat<8, 8, double> >;

extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<float, CRSMechanicalStoreTouchedPolicy>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat1x1f, CRSMechanicalStoreTouchedPolicy>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat2x2f, CRSMechanicalStoreTouchedPolicy>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat3x3f, CRSMechanicalStoreTouchedPolicy>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat4x4f, CRSMechanicalStoreTouchedPolicy>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat<6, 6, float>, CRSMechanicalStoreTouchedPolicy >;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat<8, 8, float>, CRSMechanicalStoreTouchedPolicy >;

extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<double, CRSMechanicalStoreTouchedPolicy>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat1x1d, CRSMechanicalStoreTouchedPolicy>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat2x2d, CRSMechanicalStoreTouchedPolicy>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat3x3d, CRSMechanicalStoreTouchedPolicy>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat4x4d, CRSMechanicalStoreTouchedPolicy>;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat<6, 6, double>, CRSMechanicalStoreTouchedPolicy >;
extern template class SOFA_SOFABASELINEARSOLVER_API CompressedRowSparseMatrixMechanical<defaulttype::Mat<8, 8, double>, CRSMechanicalStoreTouchedPolicy >;

#endif

} // namespace defaulttype

#endif
