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
#ifndef SOFA_CORE_BEHAVIOR_LINEARSOLVER_H
#define SOFA_CORE_BEHAVIOR_LINEARSOLVER_H

#include <sofa/core/objectmodel/BaseObject.h>

#include <sofa/core/MechanicalParams.h>
#include <sofa/defaulttype/BaseMatrix.h>

#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/core/ConstraintParams.h>


namespace sofa
{

namespace core
{

namespace behavior
{


/**
 *  \brief Abstract base class (as type identifier) for linear system solvers without any API
 *
 */
class SOFA_CORE_API BaseLinearSolver : public objectmodel::BaseObject
{
public:
    SOFA_ABSTRACT_CLASS(BaseLinearSolver, objectmodel::BaseObject);
    SOFA_BASE_CAST_IMPLEMENTATION(BaseLinearSolver)

    bool insertInNode( objectmodel::BaseNode* node ) override;
    bool removeInNode( objectmodel::BaseNode* node ) override;

    virtual bool isMultiGroup() const
    {
        return false;
    }

}; // class BaseLinearSolver




/**
 *  \brief Abstract interface for linear system solvers
 *
 */
class SOFA_CORE_API LinearSolver : public BaseLinearSolver
{
public:
    SOFA_ABSTRACT_CLASS(LinearSolver, BaseLinearSolver);
    SOFA_BASE_CAST_IMPLEMENTATION(LinearSolver)
protected:
    LinearSolver();

    ~LinearSolver() override;
public:
    /// Reset the current linear system.
    virtual void resetSystem() = 0;

    /// Set the linear system matrix, combining the mechanical M,B,K matrices using the given coefficients
    ///
    /// @todo Should we put this method in a specialized class for mechanical systems, or express it using more general terms (i.e. coefficients of the second order ODE to solve)
    virtual void setSystemMBKMatrix(const MechanicalParams* mparams) = 0;

    /// Rebuild the system using a mass and force factor
    /// Experimental API used to investigate convergence issues.
    virtual void rebuildSystem(SReal /*massFactor*/, SReal /*forceFactor*/){}

    /// Indicate if the solver update the system in parallel
    virtual bool isAsyncSolver() { return false; }

    /// Indicate if the solver updated the system after the last call of setSystemMBKMatrix (should return true if isParallelSolver return false)
    virtual bool hasUpdatedMatrix() { return true; }

    /// This function is use for the preconditioner it must be called at each time step event if setSystemMBKMatrix is not called
    virtual void updateSystemMatrix() {}

    /// Set the linear system right-hand term vector, from the values contained in the (Mechanical/Physical)State objects
    virtual void setSystemRHVector(core::MultiVecDerivId v) = 0;

    /// Set the initial estimate of the linear system left-hand term vector, from the values contained in the (Mechanical/Physical)State objects
    /// This vector will be replaced by the solution of the system once solveSystem is called
    virtual void setSystemLHVector(core::MultiVecDerivId v) = 0;

    /// Solve the system as constructed using the previous methods
    virtual void solveSystem() = 0;

    /// call of solveSystem during ode integrated step when graph propagation are forbiden; Default behavior
    virtual void directSolveSystem() {
        if (!isIterativeSolver)
        {
            solveSystem();
        }
    }

    /// call of solveSystem during ode assemble step when graph propagation is allowed; need to be override
    virtual void iterativeSolveSystem() {
        if (isIterativeSolver)
        {
            solveSystem();
        }
    }
	/// write system resolution in mstates
    virtual void writeSolution() = 0;

    ///
    virtual void init_partial_solve() { msg_warning() << "partial_solve is not implemented yet."; }

    ///
    virtual void partial_solve(std::list<defaulttype::BaseMatrix::Index>& /*I_last_Disp*/, std::list<defaulttype::BaseMatrix::Index>& /*I_last_Dforce*/, bool /*NewIn*/) { msg_warning() << "partial_solve is not implemented yet"; }

    /// Invert the system, this method is optional because it's called when solveSystem() is called for the first time
    virtual void invertSystem() {}

    /// Multiply the inverse of the system matrix by the transpose of the given matrix J
    ///
    /// @param result the variable where the result will be added
    /// @param J the matrix J to use
    /// @param fact integrator parameter
    /// @return false if the solver does not support this operation, of it the system matrix is not invertible
    virtual bool addMInvJt(defaulttype::BaseMatrix* result, defaulttype::BaseMatrix* J, SReal fact)
    {
        SOFA_UNUSED(result);
        SOFA_UNUSED(J);
        SOFA_UNUSED(fact);
        return false;
    }

    /// Build the jacobian of the constraints using a visitor
    ///
    /// @param cparams contains the MultiMatrixDerivId  which allows to retrieve the constraint jacobian to use for 
    ///        each mechanical object. 
    /// @param result the variable where the result will be added
    /// @param fact integrator parameter
    /// @return false if the solver does not support this operation, of it the system matrix is not invertible
    virtual bool buildComplianceMatrix(const sofa::core::ConstraintParams* cparams, defaulttype::BaseMatrix* result, SReal fact)
    {
        SOFA_UNUSED(cparams);
        SOFA_UNUSED(result);
        SOFA_UNUSED(fact);
        msg_error() << "buildComplianceMatrix has not been implemented.";
        return false;
    }

    /// Apply the contactforce dx = Minv * J^t * f and store the resut in dx VecId
    virtual void applyConstraintForce(const sofa::core::ConstraintParams* /*cparams*/,sofa::core::MultiVecDerivId /*dx*/, const defaulttype::BaseVector* /*f*/) {
        msg_error() << "applyConstraintForce has not been implemented.";
    }

    /// Compute the residual in the newton iterations due to the constraints forces
    /// i.e. compute mparams->dF() = J^t lambda
    /// the result is written in mparams->dF()
    virtual void computeResidual(const core::ExecParams* /*params*/, defaulttype::BaseVector* /*f*/) {
        msg_error() << "computeResidual has not been implemented.";
    }


    /// Multiply the inverse of the system matrix by the transpose of the given matrix, and multiply the result with the given matrix J
    ///
    /// @param result the variable where the result will be added
    /// @param J the matrix J to use
    /// @param fact integrator parameter
    /// @return false if the solver does not support this operation, of it the system matrix is not invertible
    virtual bool addJMInvJt(defaulttype::BaseMatrix* result, defaulttype::BaseMatrix* J, SReal fact)
    {
        SOFA_UNUSED(result);
        SOFA_UNUSED(J);
        SOFA_UNUSED(fact);
        return false;
    }

    /// Get the linear system matrix, or nullptr if this solver does not build it
    virtual defaulttype::BaseMatrix* getSystemBaseMatrix() { return nullptr; }

    /// Get the MultiMatrix view of the linear system, or nullptr if this solved does not build it
    virtual const behavior::MultiMatrixAccessor* getSystemMultiMatrixAccessor() const { return nullptr; }

    /// Get the linear system right-hand term vector, or nullptr if this solver does not build it
    virtual defaulttype::BaseVector* getSystemRHBaseVector() { return nullptr; }

    /// Get the linear system left-hand term vector, or nullptr if this solver does not build it
    virtual defaulttype::BaseVector* getSystemLHBaseVector() { return nullptr; }

    /// Get the linear system inverse matrix, or nullptr if this solver does not build it
    virtual defaulttype::BaseMatrix* getSystemInverseBaseMatrix() { return nullptr; }

    /// Read the Matrix solver from a file
    virtual bool readFile(std::istream& /*in*/) { return false;}

    /// Read the Matrix solver from a file
    virtual bool writeFile(std::ostream& /*out*/) {return false;}

    /// Ask the solver to no longer update the system matrix
    virtual void freezeSystemMatrix() { frozen = true; }



protected:

    bool frozen;
    bool isIterativeSolver = false;
};

} // namespace behavior

} // namespace core

} // namespace sofa

#endif
