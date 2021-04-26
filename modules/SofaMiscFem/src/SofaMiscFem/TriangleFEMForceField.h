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

#include <SofaMiscFem/config.h>



#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>


// corotational triangle from
// @InProceedings{NPF05,
//   author       = "Nesme, Matthieu and Payan, Yohan and Faure, Fran\c{c}ois",
//   title        = "Efficient, Physically Plausible Finite Elements",
//   booktitle    = "Eurographics (short papers)",
//   month        = "august",
//   year         = "2005",
//   editor       = "J. Dingliana and F. Ganovelli",
//   keywords     = "animation, physical model, elasticity, finite elements",
//   url          = "http://www-evasion.imag.fr/Publications/2005/NPF05"
// }



namespace sofa::component::forcefield
{

/** Triangle FEM force field using the QR decomposition of the deformation gradient, inspired from http://www-evasion.imag.fr/Publications/2005/NPF05 , to handle large displacements.
  The material properties are uniform across the domain.
  Two methods are proposed, one for small displacements and one for large displacements.
  The method for small displacements has not been validated and we suspect that it is broke. Use it very carefully, and compare with the method for large displacements.
  */
template<class DataTypes>
class TriangleFEMForceField : public core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TriangleFEMForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherited;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;

    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;

    typedef sofa::core::topology::BaseMeshTopology::Index Index;
    typedef sofa::core::topology::BaseMeshTopology::Triangle Element;
    typedef sofa::core::topology::BaseMeshTopology::SeqTriangles VecElement;

    static const int SMALL = 1;										///< Symbol of small displacements triangle solver
    static const int LARGE = 0;										///< Symbol of large displacements triangle solver

protected:
    typedef defaulttype::Vec<6, Real> Displacement;								///< the displacement vector

    typedef defaulttype::Mat<3, 3, Real> MaterialStiffness;						///< the matrix of material stiffness
    typedef sofa::helper::vector<MaterialStiffness> VecMaterialStiffness;    ///< a vector of material stiffness matrices
    VecMaterialStiffness _materialsStiffnesses;						///< the material stiffness matrices vector

    typedef defaulttype::Mat<6, 3, Real> StrainDisplacement;						///< the strain-displacement matrix (the transpose, actually)
    typedef sofa::helper::vector<StrainDisplacement> VecStrainDisplacement;	///< a vector of strain-displacement matrices
    VecStrainDisplacement _strainDisplacements;						///< the strain-displacement matrices vector

    typedef defaulttype::Mat<3, 3, Real > Transformation;						///< matrix for rigid transformations like rotations

    /// Stiffness matrix ( = RJKJtRt  with K the Material stiffness matrix, J the strain-displacement matrix, and R the transformation matrix if any )
    typedef defaulttype::Mat<9, 9, Real> StiffnessMatrix;
    
    const VecElement *_indexedElements;
    Data< VecCoord > _initialPoints; ///< the intial positions of the points


    TriangleFEMForceField();
    virtual ~TriangleFEMForceField();

    sofa::core::topology::BaseMeshTopology* m_topology;

public:

    void init() override;
    void reinit() override;
    void addForce(const core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override;
    void addDForce(const core::MechanicalParams* mparams, DataVecDeriv& df, const DataVecDeriv& dx) override;
    SReal getPotentialEnergy(const core::MechanicalParams* /*mparams*/, const DataVecCoord&  /* x */) const override
    {
        msg_warning() << "Method getPotentialEnergy not implemented yet.";
        return 0.0;
    }

    void draw(const core::visual::VisualParams* vparams) override;

    int method;
    Data<std::string> f_method; ///< Choice of method: 0 for small, 1 for large displacements
    Data<Real> f_poisson;       ///< Poisson ratio of the material
    Data<Real> f_young;         ///< Young modulus of the material
    Data<Real> f_thickness;     ///< Thickness of the elements
//    Data<Real> f_damping;       ///< Damping coefficient of the material, currently unused
    Data<bool> f_planeStrain; ///< compute material stiffness corresponding to the plane strain assumption, or to the plane stress otherwise.

    Real getPoisson() { return f_poisson.getValue(); }
    void setPoisson(Real val) { f_poisson.setValue(val); }
    Real getYoung() { return f_young.getValue(); }
    void setYoung(Real val) { f_young.setValue(val); }
    Real getThickness() { return f_thickness.getValue(); }
    void setThickness(Real val) { f_thickness.setValue(val); }

//    Real getDamping() { return f_damping.getValue(); }
//    void setDamping(Real val) { f_damping.setValue(val); }
    int  getMethod() { return method; }
    void setMethod(int val) { method = val; }

    /// Link to be set to the topology container in the component graph.
    SingleLink<TriangleFEMForceField<DataTypes>, sofa::core::topology::BaseMeshTopology, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_topology;

protected :

    /// f += Kx where K is the stiffness matrix and x a displacement
    virtual void applyStiffness( VecCoord& f, Real h, const VecCoord& x, const SReal &kFactor );
    void computeStrainDisplacement( StrainDisplacement &J, Coord a, Coord b, Coord c);
    void computeMaterialStiffnesses();
    void computeForce( Displacement &F, const Displacement &Depl, const MaterialStiffness &K, const StrainDisplacement &J );

    ////////////// small displacements method
    void initSmall();
    void accumulateForceSmall( VecCoord& f, const VecCoord & p, Index elementIndex, bool implicit = false );
//    void accumulateDampingSmall( VecCoord& f, Index elementIndex );
    void applyStiffnessSmall( VecCoord& f, Real h, const VecCoord& x, const SReal &kFactor );

    ////////////// large displacements method
    sofa::helper::vector< helper::fixed_array <Coord, 3> > _rotatedInitialElements;   ///< The initials positions in its frame
    sofa::helper::vector< Transformation > _rotations;
    void initLarge();
    void computeRotationLarge( Transformation &r, const VecCoord &p, const Index &a, const Index &b, const Index &c);
    void accumulateForceLarge( VecCoord& f, const VecCoord & p, Index elementIndex, bool implicit=false );
//    void accumulateDampingLarge( VecCoord& f, Index elementIndex );
    void applyStiffnessLarge( VecCoord& f, Real h, const VecCoord& x, const SReal &kFactor );

    //// stiffness matrix assembly
    void computeElementStiffnessMatrix( StiffnessMatrix& S, StiffnessMatrix& SR, const MaterialStiffness &K, const StrainDisplacement &J, const Transformation& Rot );
    void addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal k, unsigned int &offset) override; // compute and add all the element stiffnesses to the global stiffness matrix

    /** Accumulate an element matrix to a global assembly matrix. This is a helper for addKToMatrix, to accumulate each (square) element matrix in the (square) assembled matrix.
      \param bm the global assembly matrix
      \param offset start index of the local DOFs within the global matrix
      \param nodeIndex indices of the nodes of the element within the local nodes, as stored in the topology
      \param em element matrix, typically a stiffness, damping, mass, or weighted sum thereof
      \param scale weight applied to the matrix, typically ±params->kfactor() for a stiffness matrix
      */
    template<class IndexArray, class ElementMat>
    void addToMatrix(sofa::defaulttype::BaseMatrix* bm, unsigned offset, const IndexArray& nodeIndex, const ElementMat& em, SReal scale )
    {
        const unsigned  S = DataTypes::deriv_total_size; // size of node blocks
        for (unsigned n1=0; n1<nodeIndex.size(); n1++)
        {
            for(unsigned i=0; i<S; i++)
            {
                unsigned ROW = offset + S*nodeIndex[n1] + i;  // i-th row associated with node n1 in BaseMatrix
                unsigned row = S*n1+i;                        // i-th row associated with node n1 in the element matrix

                for (unsigned n2=0; n2<nodeIndex.size(); n2++)
                {
                    for (unsigned j=0; j<S; j++)
                    {
                        unsigned COLUMN = offset + S*nodeIndex[n2] +j; // j-th column associated with node n2 in BaseMatrix
                        unsigned column = 3*n2+j;                      // j-th column associated with node n2 in the element matrix
                        bm->add( ROW,COLUMN, em[row][column]* scale );
                    }
                }
            }
        }
    }
};


#if  !defined(SOFA_COMPONENT_FORCEFIELD_TRIANGLEFEMFORCEFIELD_CPP)

extern template class SOFA_SOFAMISCFEM_API TriangleFEMForceField<sofa::defaulttype::Vec3Types>;


#endif //  !defined(SOFA_COMPONENT_FORCEFIELD_TRIANGLEFEMFORCEFIELD_CPP)

} // namespace sofa::component::forcefield
