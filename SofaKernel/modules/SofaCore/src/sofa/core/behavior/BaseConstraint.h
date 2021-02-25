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
#ifndef SOFA_CORE_BEHAVIOR_BASECONSTRAINT_H
#define SOFA_CORE_BEHAVIOR_BASECONSTRAINT_H

#include <sofa/core/behavior/BaseConstraintSet.h>



namespace sofa
{

namespace core
{

namespace behavior
{

/**
 *  \brief Object computing a constraint resolution within a Gauss-Seidel algorithm
 */

class SOFA_CORE_API ConstraintResolution
{
public:
    ConstraintResolution(unsigned int nbLines = 3, double tolerance = 0.0);

    virtual ~ConstraintResolution() {}

    /// The resolution object can do precomputation with the compliance matrix, and give an initial guess.
    virtual void init(int /*line*/, double** /*w*/, double* /*force*/);

    /// The resolution object can provide an initial guess
    virtual void initForce(int /*line*/, double* /*force*/);

    /// Resolution of the constraint for one Gauss-Seidel iteration
    virtual void resolution(int line, double** w, double* d, double* force, double * dFree)
    {
        SOFA_UNUSED(line);
        SOFA_UNUSED(w);
        SOFA_UNUSED(d);
        SOFA_UNUSED(force);
        SOFA_UNUSED(dFree);
        dmsg_error("ConstraintResolution")
                << "resolution(int , double** , double* , double* , double * ) not implemented." ;
    }

    /// Called after Gauss-Seidel last iteration, in order to store last computed forces for the inital guess
    virtual void store(int /*line*/, double* /*force*/, bool /*convergence*/) {}

    inline void setNbLines(unsigned int nbLines)
    {
        m_nbLines = nbLines;
    }

    inline unsigned int getNbLines() const
    {
        return m_nbLines;
    }

    inline void setTolerance(double tolerance)
    {
        m_tolerance = tolerance;
    }

    inline double getTolerance() const
    {
        return m_tolerance;
    }

private:
    /// Number of dof used by this particular constraint. To be modified in the object's constructor.
    unsigned int m_nbLines;

    /// Custom tolerance, used for the convergence of this particular constraint instead of the global tolerance
    double m_tolerance;
};

/**
 *  \brief Component computing constraints within a simulated body.
 *
 *  This class define the abstract API common to all constraints.
 *  A BaseConstraint computes constraints applied to one or more simulated body
 *  given its current position and velocity.
 *
 *  Constraints can be internal to a given body (attached to one MechanicalState,
 *  see the Constraint class), or link several bodies together (such as contacts,
 *  see the InteractionConstraint class).
 *
 */
class SOFA_CORE_API BaseConstraint : public BaseConstraintSet
{
public:
    SOFA_ABSTRACT_CLASS(BaseConstraint, BaseConstraintSet);
    SOFA_BASE_CAST_IMPLEMENTATION(BaseConstraint)

protected:
    BaseConstraint() {}
    ~BaseConstraint() override {}

private:
    BaseConstraint(const BaseConstraint& n) ;
    BaseConstraint& operator=(const BaseConstraint& n) ;

public:
    /// Get the ID of the group containing this constraint. This ID is used to specify which constraints are solved by which solver, by specifying in each solver which groups of constraints it should handle.
    int getGroup() const { return group.getValue(); }

    /// Set the ID of the group containing this constraint. This ID is used to specify which constraints are solved by which solver, by specifying in each solver which groups of constraints it should handle.
    void setGroup(int g) { group.setValue(g); }

    typedef long long PersistentID;
    typedef helper::vector<PersistentID> VecPersistentID;
    typedef defaulttype::Vec<3,int> ConstCoord;
    typedef helper::vector<ConstCoord> VecConstCoord;
    typedef defaulttype::Vec<3,double> ConstDeriv;
    typedef helper::vector<ConstDeriv> VecConstDeriv;
    typedef double ConstArea;
    typedef helper::vector<ConstArea> VecConstArea;

    class ConstraintBlockInfo
    {
    public:
        BaseConstraint* parent;
        int const0; ///< index of first constraint
        int nbLines; ///< how many dofs (i.e. lines in the matrix) are used by each constraint
        int nbGroups; ///< how many groups of constraints are active
        bool hasId; ///< true if this constraint has persistent ID information
        bool hasPosition; ///< true if this constraint has coordinates information
        bool hasDirection; ///< true if this constraint has direction information
        bool hasArea; ///< true if this constraint has area information
        int offsetId; ///< index of first constraint group info in vector of persistent ids and coordinates
        int offsetPosition; ///< index of first constraint group info in vector of coordinates
        int offsetDirection; ///< index of first constraint info in vector of directions
        int offsetArea; ///< index of first constraint group info in vector of areas
        ConstraintBlockInfo() : parent(nullptr), const0(0), nbLines(1), nbGroups(0), hasId(false), hasPosition(false), hasDirection(false), hasArea(false), offsetId(0), offsetPosition(0), offsetDirection(0), offsetArea(0)
        {}
    };
    typedef helper::vector<ConstraintBlockInfo> VecConstraintBlockInfo;

    /// Get information for each constraint: pointer to parent BaseConstraint, unique persistent ID, 3D position
    /// \param cParams defines the state vectors to use for positions and velocities. Also defines the order of the constraint (POS, VEL, ACC) and resolution parameters (smoothness, ...)
    virtual void getConstraintInfo(const ConstraintParams* cParams, VecConstraintBlockInfo& blocks, VecPersistentID& ids, VecConstCoord& positions, VecConstDeriv& directions, VecConstArea& areas)
    {
        SOFA_UNUSED(cParams);
        SOFA_UNUSED(blocks);
        SOFA_UNUSED(ids);
        SOFA_UNUSED(positions);
        SOFA_UNUSED(directions);
        SOFA_UNUSED(areas);

    }

    /// Add the corresponding ConstraintResolution using the offset parameter
    /// \param cParams defines the state vectors to use for positions and velocities. Also defines the order of the constraint (POS, VEL, ACC) and resolution parameters (smoothness, ...)
    /// \param resTab is the result vector that contains the contraint resolution algorithms
    virtual void getConstraintResolution(const ConstraintParams* cParams, std::vector<ConstraintResolution*> &resTab, unsigned int &offset)
    {

        getConstraintResolution(resTab, offset);
        SOFA_UNUSED(cParams);

    }

    virtual void getConstraintResolution(std::vector<ConstraintResolution*> &resTab, unsigned int &offset)
    {
        SOFA_UNUSED(resTab);
        SOFA_UNUSED(offset);
    }



    /// Store the constraint lambda at the constraint dofs at the given VecDerivId location. 
    /// res = J^t * lambda. 
    /// J is the sparse matrix containing the constraint jacobian that was used to build the constraint matrix ( see BaseConstraintSet::buildConstraintMatrix ).
    /// \param cParams stores the id of the state vectors used during the constraint solving step. Mostly it helps retrieving the MatrixDerivId where
    ///        the constraint jacobian J is stored.
    /// \param res is the state vector Id where to store the result.
    /// \param lambda is the vector of scalar constraint impulses. The direction are stored in the MatrixDerivId stored in the cParams.
    virtual void storeLambda(const ConstraintParams* cParams, MultiVecDerivId res, const sofa::defaulttype::BaseVector* lambda) = 0;
};

} // namespace behavior

} // namespace core

} // namespace sofa

#endif
