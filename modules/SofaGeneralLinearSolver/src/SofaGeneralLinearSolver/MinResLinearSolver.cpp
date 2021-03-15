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
#define SOFA_COMPONENT_LINEARSOLVER_MINRESLINEARSOLVER_CPP
#include <SofaGeneralLinearSolver/MinResLinearSolver.inl>

#include <sofa/core/ObjectFactory.h>
#include <iostream>

namespace sofa::component::linearsolver
{

using namespace sofa::defaulttype;


int MinResLinearSolverClass = core::RegisterObject("Linear system solver using the MINRES iterative algorithm")
        .add< MinResLinearSolver< GraphScatteredMatrix, GraphScatteredVector > >(true)
        .add< MinResLinearSolver< FullMatrix<double>, FullVector<double> > >()
        .add< MinResLinearSolver< SparseMatrix<double>, FullVector<double> > >()
        .add< MinResLinearSolver< CompressedRowSparseMatrixMechanical<double>, FullVector<double> > >()
        .add< MinResLinearSolver< CompressedRowSparseMatrixMechanical<Mat<2,2,double> >, FullVector<double> > >()
        .add< MinResLinearSolver< CompressedRowSparseMatrixMechanical<Mat<3,3,double> >, FullVector<double> > >()
        .add< MinResLinearSolver< CompressedRowSparseMatrixMechanical<Mat<4,4,double> >, FullVector<double> > >()
        .add< MinResLinearSolver< CompressedRowSparseMatrixMechanical<Mat<6,6,double> >, FullVector<double> > >()
        .add< MinResLinearSolver< CompressedRowSparseMatrixMechanical<Mat<8,8,double> >, FullVector<double> > >()

        .addAlias("MINRESSolver")
        .addAlias("MinResSolver")
        ;

template class SOFA_SOFAGENERALLINEARSOLVER_API MinResLinearSolver< GraphScatteredMatrix, GraphScatteredVector >;
template class SOFA_SOFAGENERALLINEARSOLVER_API MinResLinearSolver< FullMatrix<double>, FullVector<double> >;
template class SOFA_SOFAGENERALLINEARSOLVER_API MinResLinearSolver< SparseMatrix<double>, FullVector<double> >;
template class SOFA_SOFAGENERALLINEARSOLVER_API MinResLinearSolver< CompressedRowSparseMatrixMechanical<double>, FullVector<double> >;
template class SOFA_SOFAGENERALLINEARSOLVER_API MinResLinearSolver< CompressedRowSparseMatrixMechanical<Mat<2,2,double> >, FullVector<double> >;
template class SOFA_SOFAGENERALLINEARSOLVER_API MinResLinearSolver< CompressedRowSparseMatrixMechanical<Mat<3,3,double> >, FullVector<double> >;
template class SOFA_SOFAGENERALLINEARSOLVER_API MinResLinearSolver< CompressedRowSparseMatrixMechanical<Mat<4,4,double> >, FullVector<double> >;
template class SOFA_SOFAGENERALLINEARSOLVER_API MinResLinearSolver< CompressedRowSparseMatrixMechanical<Mat<6,6,double> >, FullVector<double> >;
template class SOFA_SOFAGENERALLINEARSOLVER_API MinResLinearSolver< CompressedRowSparseMatrixMechanical<Mat<8,8,double> >, FullVector<double> >;


} //namespace sofa::component::linearsolver
