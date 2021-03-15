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
#ifndef SOFA_COMPONENT_LINEARSOLVER_SPARSECHOLESKYSOLVER_H
#define SOFA_COMPONENT_LINEARSOLVER_SPARSECHOLESKYSOLVER_H
#include <SofaSparseSolver/config.h>

#include <sofa/core/behavior/LinearSolver.h>
#include <SofaBaseLinearSolver/MatrixLinearSolver.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <SofaBaseLinearSolver/SparseMatrix.h>
#include <SofaBaseLinearSolver/CompressedRowSparseMatrixMechanical.h>
#include <sofa/helper/map.h>
#include <cmath>
#include <csparse.h>

namespace sofa
{

namespace component
{

namespace linearsolver
{

/// Direct linear solver based on Sparse Cholesky factorization, implemented with the CSPARSE library
template<class TMatrix, class TVector>
class SparseCholeskySolver : public sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(SparseCholeskySolver,TMatrix,TVector),SOFA_TEMPLATE2(sofa::component::linearsolver::MatrixLinearSolver,TMatrix,TVector));

    typedef TMatrix Matrix;
    typedef TVector Vector;
    typedef sofa::component::linearsolver::MatrixLinearSolver<TMatrix,TVector> Inherit;

    Data<bool> f_verbose; ///< Dump system state at each iteration

    SparseCholeskySolver();
    ~SparseCholeskySolver();
    void solve (Matrix& M, Vector& x, Vector& b) override;
    void invert(Matrix& M) override;

public :
    cs A;
    css *S;
    csn *N;
    int * A_i;
    int * A_p;
    helper::vector<double> A_x,z_tmp,r_tmp,tmp;

    void solveT(double * z, double * r);
    void solveT(float * z, float * r);
};

#if  !defined(SOFA_COMPONENT_LINEARSOLVER_SPARSECHOLESKYSOLVER_CPP)
extern template class SOFA_SOFASPARSESOLVER_API SparseCholeskySolver< CompressedRowSparseMatrixMechanical<double>,FullVector<double> >;
extern template class SOFA_SOFASPARSESOLVER_API SparseCholeskySolver< CompressedRowSparseMatrixMechanical<float>,FullVector<float> >;
#endif

} // namespace linearsolver

} // namespace component

} // namespace sofa


#endif
