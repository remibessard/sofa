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
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>

#include <sofa/core/visual/VisualParams.h>
#include <SofaBaseLinearSolver/CompressedRowSparseMatrix.h>

using sofa::core::behavior::BaseMechanicalState;

/// This line registers the DefaultMultiMatrixAccessor to the messaging system
/// allowing to use msg_info() instead of msg_info("DefaultMultiMatrixAccessor")
MSG_REGISTER_CLASS(sofa::component::linearsolver::DefaultMultiMatrixAccessor, "DefaultMultiMatrixAccessor")

namespace sofa::component::linearsolver
{


DefaultMultiMatrixAccessor::DefaultMultiMatrixAccessor()
    : globalMatrix(nullptr)
    , globalDim(0)
{
}


DefaultMultiMatrixAccessor::~DefaultMultiMatrixAccessor()
{
    this->clear();
}

void DefaultMultiMatrixAccessor::clear()
{
    globalDim = 0;
    for (auto it = realStateOffsets.begin(), itend = realStateOffsets.end(); it != itend; ++it)
        it->second = -1;

    for (std::map< const sofa::core::behavior::BaseMechanicalState*, defaulttype::BaseMatrix* >::iterator it = mappedMatrices.begin(), itend = mappedMatrices.end(); it != itend; ++it)
        if (it->second != nullptr) delete it->second;
    mappedMatrices.clear();
    diagonalStiffnessBloc.clear();

    for (std::map< std::pair<const BaseMechanicalState*, const BaseMechanicalState*>, InteractionMatrixRef >::iterator it = interactionStiffnessBloc.begin(), itend = interactionStiffnessBloc.end(); it != itend; ++it)
        if (it->second.matrix != nullptr && it->second.matrix != globalMatrix) delete it->second.matrix;

    interactionStiffnessBloc.clear();
    mappingList.clear();

}


void DefaultMultiMatrixAccessor::setGlobalMatrix(defaulttype::BaseMatrix* matrix)
{
    this->globalMatrix = matrix;
}

void DefaultMultiMatrixAccessor::addMechanicalState(const sofa::core::behavior::BaseMechanicalState* mstate)
{
    auto dim = mstate->getMatrixSize();
    realStateOffsets[mstate] = globalDim;
    globalDim += dim;

    if(m_doPrintInfo)/////////////////////////////////////////////////////////
    {
        msg_info() << "Adding '" << mstate->getPathName()
                << "' in the global matrix ["<<dim<<"."<<dim
                <<"] at offset (" << realStateOffsets[mstate] <<","<< realStateOffsets[mstate]<<")";
    }
}


void DefaultMultiMatrixAccessor::addMechanicalMapping(sofa::core::BaseMapping* mapping)
{
    const sofa::defaulttype::BaseMatrix* jmatrix = nullptr;
    if (mapping->isMechanical() && mapping->areMatricesMapped())
        jmatrix = mapping->getJ();

    if (jmatrix)
    {

        const BaseMechanicalState* mappedState  = const_cast<const BaseMechanicalState*>(mapping->getMechTo()[0]);
        defaulttype::BaseMatrix* mappedstiffness;
        mappedstiffness = mapping->createMappedMatrix(mappedState,mappedState,&DefaultMultiMatrixAccessor::createMatrix);
        mappedMatrices[mappedState]=mappedstiffness;

        mappingList.push_back(mapping);

        if(m_doPrintInfo)/////////////////////////////////////////////////////////
        {
            msg_info() << "Adding validated MechanicalMapping '" << mapping->getPathName()
                       << "' with J["<< jmatrix->rowSize()<<"."<<jmatrix->colSize()<<"]" ;
        }
    }
}


void DefaultMultiMatrixAccessor::addMappedMechanicalState(const sofa::core::behavior::BaseMechanicalState* /*mstate*/)
{
    // do not add the mapped mechanical state here because
    // a mapped mechanical state is added if and only if it has its own stiffness matrix and its mapping must satistify several conditions
    // so if and only if there are a forcefield or other component call getMatrix(mappedstate)
}

void DefaultMultiMatrixAccessor::setupMatrices()
{
    auto it = realStateOffsets.begin(), itend = realStateOffsets.end();
    while (it != itend)
    {
        if (globalMatrix)
        {
            MatrixRef& r = diagonalStiffnessBloc[it->first];
            r.matrix = globalMatrix;
            r.offset = it->second;
        }
        ++it;
    }

    if(m_doPrintInfo)
    {
        msg_info() << "Setting up the Global Matrix [" << globalDim << "." << globalDim << "] for " << realStateOffsets.size() << " real mechanical state(s)." ;
    }
}

DefaultMultiMatrixAccessor::Index DefaultMultiMatrixAccessor::getGlobalDimension() const
{
    return globalDim;
}

int DefaultMultiMatrixAccessor::getGlobalOffset(const sofa::core::behavior::BaseMechanicalState* mstate) const
{
    auto it = realStateOffsets.find(mstate);
    if (it != realStateOffsets.end())
        return it->second;
    return -1;
}

DefaultMultiMatrixAccessor::MatrixRef DefaultMultiMatrixAccessor::getMatrix(const sofa::core::behavior::BaseMechanicalState* mstate) const
{
#if 0
    MatrixRef r;
    auto itRealState = realStateOffsets.find(mstate);

    if (itRealState != realStateOffsets.end()) //case where mechanical state is a non mapped state
    {
        if (globalMatrix)
        {
            r.matrix = globalMatrix;
            r.offset = itRealState->second;
        }
    }
    else //case where mechanical state is a mapped state
    {
        std::map< const sofa::core::behavior::BaseMechanicalState*, defaulttype::BaseMatrix*>::iterator itmapped = mappedMatrices.find(mstate);
        if (itmapped != mappedMatrices.end()) // this mapped state and its matrix has been already added and created
        {
            r.matrix = itmapped->second;
            r.offset = 0;
        }
        else // this mapped state and its matrix hasnt been created we creat it and its matrix by "createMatrix"
        {

            defaulttype::BaseMatrix* m = createMatrixImpl(mstate,mstate, m_doPrintInfo);
            r.matrix = m;
            r.offset = 0;
            //when creating an matrix, it dont have to be added before
            assert(diagonalStiffnessBloc.find(mstate) == diagonalStiffnessBloc.end());
            mappedMatrices[mstate]=r.matrix;
        }
    }

    diagonalStiffnessBloc[mstate] = r;

    if(m_doPrintInfo)
    {
        if (r.matrix != nullptr)
        {
            msg_info() << "Giving Stiffness Matrix [" << r.matrix->rowSize() << "." << r.matrix->colSize() << "] for state '" << mstate->getPathName()
                    << "' at offset (" << r.offset  <<","<< r.offset <<")" ;
        }
        else
            msg_warning() << "nullptr matrix found for state " << mstate->getName() ;
    }
    return r;
#endif

	MatrixRef r;

	/*std::map< const sofa::core::behavior::BaseMechanicalState*, int >::const_iterator*/ auto itRealState = realStateOffsets.find(mstate);

	if (itRealState != realStateOffsets.end()) //case where mechanical state is a non mapped state
	{
		if (globalMatrix)
		{
			r.matrix = globalMatrix;
			r.offset = itRealState->second;
		}
	}
	else //case where mechanical state is a mapped state
	{

		std::map< const sofa::core::behavior::BaseMechanicalState*, defaulttype::BaseMatrix*>::iterator itmapped = mappedMatrices.find(mstate);
		if (itmapped != mappedMatrices.end()) // this mapped state and its matrix has been already added and created
		{
			r.matrix = itmapped->second;
			r.offset = 0;
		}
		else // this mapped state and its matrix hasnt been created we creat it and its matrix by "createMatrix"
		{
			defaulttype::BaseMatrix* m = createMatrix(mstate, mstate);
			r.matrix = m;
			r.offset = 0;
			//when creating an matrix, it dont have to be added before
			assert(diagonalStiffnessBloc.find(mstate) == diagonalStiffnessBloc.end());
			mappedMatrices[mstate] = r.matrix;
		}
	}

	diagonalStiffnessBloc[mstate] = r;

	return r;
}

DefaultMultiMatrixAccessor::InteractionMatrixRef DefaultMultiMatrixAccessor::getMatrix(const sofa::core::behavior::BaseMechanicalState* mstate1, const sofa::core::behavior::BaseMechanicalState* mstate2) const
{
    InteractionMatrixRef r2;
    if (mstate1 == mstate2)// case where state1 == state2, interaction matrix is on the diagonal stiffness bloc
    {
        MatrixRef r = diagonalStiffnessBloc.find(mstate1)->second;
        r2.matrix = r.matrix;
        r2.offRow = r.offset;
        r2.offCol = r.offset;

        if(m_doPrintInfo)///////////////////////////////////////////
        {
            if (r2.matrix != nullptr)
            {
                msg_info() << "Giving Interaction Stiffness Matrix ["
                        << r2.matrix->rowSize() << "." << r2.matrix->colSize()
                        <<"] at offset ("<<r2.offRow <<","<< r2.offCol<<") for self-interaction : "
                        <<mstate1->getName()<<"["<<mstate1->getMatrixSize()<<"] --- "<<mstate2->getName()<<"["<<mstate2->getMatrixSize()<<"]" ;
            }else
            {
                msg_warning() << "Giving nullptr matrix for self-interaction "<<
                        mstate1->getName()<<"["<<mstate1->getMatrixSize()<<"] --- "<<mstate2->getName()<<"["<<mstate2->getMatrixSize()<<"]" ;
            }
        }
    }
    else// case where state1 # state2
    {
        std::pair<const BaseMechanicalState*,const BaseMechanicalState*> pairMS = std::make_pair(mstate1,mstate2);

        std::map< std::pair<const BaseMechanicalState*,const BaseMechanicalState*>, InteractionMatrixRef >::iterator it = interactionStiffnessBloc.find(pairMS);
        if (it != interactionStiffnessBloc.end())// the interaction is already added
        {
            if(it->second.matrix != nullptr)
            {
                r2 = it->second;
            }

            if(m_doPrintInfo)///////////////////////////////////////////
            {
                if(r2.matrix != nullptr)
                {
                    msg_info() << "Giving Interaction Stiffness Matrix ["
                            << r2.matrix->rowSize() << "." << r2.matrix->colSize()
                            <<"] at offset ("<<r2.offRow <<","<< r2.offCol<<")  for interaction : "
                            <<mstate1->getName()<<"["<<mstate1->getMatrixSize()<<"] --- "<<mstate2->getName()<<"["<<mstate2->getMatrixSize()<<"]" ;
                }
                else
                {
                    msg_warning() << "Giving nullptr matrix  for interaction "
                            <<mstate1->getName()<<"["<<mstate1->getMatrixSize()<<"] --- "<<mstate2->getName()<<"["<<mstate2->getMatrixSize()<<"]" ;
                }
            }
        }
        else// the interaction is not added, we need to creat it and its matrix
        {
            auto it1 = realStateOffsets.find(mstate1);
            auto it2 = realStateOffsets.find(mstate2);

            if(it1 != realStateOffsets.end() && it2 != realStateOffsets.end())// case where all of two ms are real DOF (non-mapped)
            {
                if (globalMatrix)
                {
                    r2.matrix = globalMatrix;
                    r2.offRow = it1->second;
                    r2.offCol = it2->second;

                    if(m_doPrintInfo)/////////////////////////////////////////////////////////
                    {
                        if (r2.matrix != nullptr)
                        {
                            msg_info() <<  "Giving Interaction Stiffness Matrix ["
                                    << r2.matrix->rowSize() << "." << r2.matrix->colSize()
                                    <<"] at offset ("<<r2.offRow <<","<< r2.offCol<<")  for interaction : "
                                    <<mstate1->getName()<<"["<<mstate1->getMatrixSize()<<"] --- "<<mstate2->getName()<<"["<<mstate2->getMatrixSize()<<"]" ;
                        }else{
                            msg_warning() << "Giving nullptr matrix  for interaction "
                                    << " for interaction : "
                                    << mstate1->getName()<<"["<<mstate1->getMatrixSize()<<"] --- "<<mstate2->getName()<<"["<<mstate2->getMatrixSize()<<"]" ;
                        }
                    }
                }
            }
            else //case where at least one ms is a mapped
            {
                defaulttype::BaseMatrix* m = createMatrixImpl(mstate1,mstate2,m_doPrintInfo);
                r2.matrix = m;
                r2.offRow = 0;
                r2.offCol = 0;
                //when creating an matrix, it dont have to be added before
                assert(interactionStiffnessBloc.find(pairMS) == interactionStiffnessBloc.end());

                if(m_doPrintInfo)/////////////////////////////////////////////////////////
                {
                    if (r2.matrix != nullptr)
                    {
                        msg_info() <<   "Giving Interaction Stiffness Matrix ["
                                << r2.matrix->rowSize() << "." << r2.matrix->colSize()
                                << "] at offset ("<<r2.offRow <<","<< r2.offCol<<")  for interaction : "
                                << mstate1->getName()<< "[" <<mstate1->getMatrixSize()<< "] --- " <<mstate2->getName()<<"[" <<mstate2->getMatrixSize()<<"]" ;
                    }
                    else
                    {
                        msg_info() << "Giving nullptr matrix  for interaction "
                                << " for interaction : "
                                << mstate1->getName()<<"["<<mstate1->getMatrixSize()<<"] --- "<<mstate2->getName()<<"["<<mstate2->getMatrixSize()<<"]" ;
                    }
                }
            }

            interactionStiffnessBloc[pairMS]=r2;
        }// end of the interaction is not added, we need to creat it and its matrix

    }//end of case where state1 # state2


    if(m_doPrintInfo && r2.matrix == nullptr)
    {
        msg_warning() << "nullptr matrix found for interaction " << mstate1->getName()<<" --- "<<mstate2->getName() ;
    }

    return r2;
}

void DefaultMultiMatrixAccessor::computeGlobalMatrix()
{
    const int lastMappingId = (int)mappingList.size() - 1;
    for(int id=lastMappingId; id>=0; --id)
    {
        sofa::core::BaseMapping* m_mapping = mappingList[id];
        const BaseMechanicalState* instate  = const_cast<const BaseMechanicalState*>(m_mapping->getMechFrom()[0]);
        const BaseMechanicalState* outstate  = const_cast<const BaseMechanicalState*>(m_mapping->getMechTo()[0]);

        const defaulttype::BaseMatrix* matrixJ = m_mapping->getJ();
        const auto nbR_J = matrixJ->rowSize();
        const auto nbC_J = matrixJ->colSize();

        if(m_doPrintInfo)/////////////////////////////////////////////////////////
        {
            msg_info() << "[CONTRIBUTION] of " << id << "-th mapping named : " << m_mapping->getName()
                    << " with fromModel : "<< instate->getName()
                    << " and toModel : "   << outstate->getName() <<" -- with matrix J["<<nbR_J <<"."<<nbC_J <<"]" ;
        }


        //for toModel -----------------------------------------------------------
        if(mappedMatrices.find(outstate) == mappedMatrices.end())
        {
            if(m_doPrintInfo)
            {
                msg_info() << "	[Propa.Stiff] WARNING toModel : "<< outstate->getName()<< " dont have stiffness matrix" ;
            }
        }

        if(diagonalStiffnessBloc.find(outstate) != diagonalStiffnessBloc.end())
        {
            //=================================
            //           K11 += Jt * K22 * J
            //=================================
            MatrixRef K1 = this->getMatrix(instate);
            MatrixRef K2 = this->getMatrix(outstate);

            const auto offset1  = K1.offset;
            const auto offset2  = K2.offset;
            const Index sizeK1 = Index(K1.matrix->rowSize() - offset1);
            const Index sizeK2 = Index(K2.matrix->rowSize() - offset2);

            if(m_doPrintInfo)/////////////////////////////////////////////////////////
            {
                msg_info() << "	[Propa.Stiff] propagating stiffness of : "<< outstate->getName()<< " to stifness "<<instate->getName()<< msgendl;
                msg_info() <<"	                    **multiplication of "
                        <<" K1["<<sizeK1<<"."<<sizeK1<< "]("<< offset1<<","<<offset1 <<  ")   =  "
                        <<" Jt["<<nbC_J<<"."<<nbR_J<< "] * "
                        <<" K2["<<sizeK2<<"."<<sizeK2<<"]("<< offset2<<","<<offset2 <<  ") * "
                        <<" J["<<nbR_J<<"."<<nbC_J<< "]"
                        << msgendl;
            }

            // Matrix multiplication  K11 += Jt * K22 * J
            for(Index i1 =0; i1 < sizeK1 ; ++i1)
            {
                for(Index j1 =0 ; j1 < sizeK1 ; ++j1)
                {
                    double Jt_K2_J_i1j1 = 0;

                    for(Index i2 =0 ; i2 < sizeK2 ; ++i2)
                    {
                        for(Index j2 =0 ; j2 < sizeK2 ; ++j2)
                        {
                            const double K2_i2j2 = (double) K2.matrix->element(offset2 + i2, offset2 + j2);
                            for(Index k2=0 ; k2 < sizeK2 ; ++k2)
                            {
                                const double Jt_i1k2 = (double) matrixJ->element( i1 , k2 ) ;
                                const double  J_k2j1 = (double) matrixJ->element( k2 , j1 ) ;

                                Jt_K2_J_i1j1 += Jt_i1k2 * K2_i2j2  * J_k2j1;
                            }
                        }
                    }
                    K1.matrix->add(offset1 + i1 , offset1 + j1 , Jt_K2_J_i1j1);
                }
            }
            // Matrix multiplication  K11 += Jt * K22 * J
        }

        std::vector<std::pair<const BaseMechanicalState*, const BaseMechanicalState*> > interactionList;
        for (std::map< std::pair<const BaseMechanicalState*, const BaseMechanicalState*>, InteractionMatrixRef >::iterator it = interactionStiffnessBloc.begin(), itend = interactionStiffnessBloc.end(); it != itend; ++it)
        {
            if(it->first.first == outstate || it->first.second == outstate )
            {
                interactionList.push_back(it->first);
            }
        }


        const size_t nbInteraction = interactionList.size();
        for(size_t i=0; i< nbInteraction; i++)
        {

            if(m_doPrintInfo)/////////////////////////////////////////////////////////
            {
                msg_info() << "	[Propa.Interac.Stiff] detected interaction between toModel : "<<interactionList[i].first->getName()
                        << " and : " <<interactionList[i].second->getName() ;
            }
            //                   |       |
            //                 MS1     MS2
            //                  /      /
            //                map    inter
            //                   \   /
            //                   MS3/
            //
            //           K_11 += Jt * K_33 * J
            //           I_12 += Jt * I_32
            //           I_21 +=      I_23 * J
            //

            if(interactionList[i].first == outstate)
            {
                InteractionMatrixRef I_32 = this->getMatrix( outstate , interactionList[i].second);
                InteractionMatrixRef I_12 = this->getMatrix( instate  , interactionList[i].second);
                //===========================
                //          I_12 += Jt * I_32
                //===========================
                const Index offR_I_12  = I_12.offRow;                       //      row offset of I12 matrix
                const Index offC_I_12  = I_12.offCol;                       //    colum offset of I12 matrix
                const Index nbR_I_12   = I_12.matrix->rowSize() - offR_I_12;//number of rows   of I12 matrix
                const Index nbC_I_12   = I_12.matrix->colSize() - offC_I_12;//number of colums of I12 matrix

                const Index offR_I_32  = I_32.offRow;                     //      row offset of I32 matrix
                const Index offC_I_32  = I_32.offCol;                     //    colum offset of I32 matrix
                const Index nbR_I_32 = I_32.matrix->rowSize() - offR_I_32;//number of rows   of I32 matrix
                const Index nbC_I_32 = I_32.matrix->colSize() - offC_I_32;//number of colums of I32 matrix


                if(m_doPrintInfo)/////////////////////////////////////////////////////////
                {
                    msg_info() << "	[Propa.Interac.Stiff] propagating interaction "
                            <<outstate->getName() << "--" <<interactionList[i].second->getName()
                            <<"  to : "
                            <<instate->getName() << "--" <<interactionList[i].second->getName() ;

                    msg_info() <<"	                    **multiplication of "
                            <<" I12["<<nbR_I_12<<"."<<nbC_I_12<< "]("<< offR_I_12<<","<<offC_I_12 <<  ")  =  "
                            <<" Jt["<<nbC_J<<"."<<nbR_J<< "]  *  "
                            <<" I32["<<nbR_I_32<<"."<<nbC_I_32<<"]("<< offR_I_32<<","<<offC_I_32 <<  ")" ;
                }

                // Matrix multiplication   I_12 += Jt * I_32
                for(Index _i = 0; _i < nbR_I_12 ; _i++)
                {
                    for(Index _j = 0; _j < nbC_I_12 ; _j++)
                    {
                        double Jt_I32_ij = 0;
                        for(Index _k = 0; _k < nbR_I_32 ; _k++)
                        {
                            const double Jt_ik    = (double) matrixJ->element( _k, _i ) ;
                            const double  I_32_kj = (double) I_32.matrix->element( offR_I_32 + _k, offC_I_32+_j) ;

                            Jt_I32_ij += Jt_ik  *  I_32_kj;
                        }
                        I_12.matrix->add(offR_I_12 + _i , offC_I_12 +  _j , Jt_I32_ij);
                    }
                }// Matrix multiplication   I_12 += Jt * I_32
            }

            if(interactionList[i].second == outstate)
            {
                InteractionMatrixRef I_21 = this->getMatrix(interactionList[i].first,instate);
                InteractionMatrixRef I_23 = this->getMatrix(interactionList[i].first,outstate);
                //=========================================
                //          I_21 +=      I_23 * J
                //=========================================
                const Index offR_I_21  = I_21.offRow;
                const Index offC_I_21  = I_21.offCol;
                const Index nbR_I_21 = I_21.matrix->rowSize() - offR_I_21;
                const Index nbC_I_21 = I_21.matrix->colSize() - offC_I_21;

                const Index offR_I_23  = I_23.offRow;
                const Index offC_I_23  = I_23.offCol;
                const Index nbR_I_23   = I_23.matrix->rowSize() - offR_I_23;
                const Index nbC_I_23   = I_23.matrix->colSize() - offC_I_23;

                if(m_doPrintInfo)/////////////////////////////////////////////////////////
                {
                    msg_info() << "	[Propa.Interac.Stiff] propagating interaction "
                            <<interactionList[i].first->getName()<< "--" <<outstate->getName()
                            <<" to : "<<interactionList[i].first->getName()<< "--" <<instate->getName() ;

                    msg_info() <<"	                    **multiplication of "
                            <<" I_21["<<nbR_I_21<<"."<<nbC_I_21<<"]("<< offR_I_21<<","<<offC_I_21 <<  ")  =  "
                            <<" I23["<<nbR_I_23<<"."<<nbC_I_23<< "]("<< offR_I_23<<","<<offC_I_23 <<  ")  *  "
                            <<" J["<<nbR_J<<"."<<nbC_J<<"]" ;
                }


                // Matrix multiplication  I_21 +=  I_23 * J
                for(Index _i = 0; _i < nbR_I_21 ; _i++)
                {
                    for(Index _j = 0; _j < nbC_I_21 ; _j++)
                    {
                        double I23_J_ij = 0;
                        for(Index _k = 0; _k < nbC_I_23 ; _k++)
                        {
                            const double I_23_ik = (double) I_23.matrix->element( offR_I_23 + _i, offC_I_23+_k) ;
                            const double J_kj    = (double) matrixJ->element( _k, _j ) ;

                            I23_J_ij += I_23_ik  * J_kj ;
                        }
                        I_21.matrix->add(offR_I_21 + _i , offC_I_21 + _j , I23_J_ij);
                    }
                }// Matrix multiplication  I_21 +=  I_23 * J

            }

            //after propagating the interaction, we remove the older interaction
            interactionStiffnessBloc.erase( interactionStiffnessBloc.find(interactionList[i]) );
            if(m_doPrintInfo)
            {
                msg_info() << "	--[Propa.Interac.Stiff] remove interaction of : "<<interactionList[i].first->getName()
                        << " and : " << interactionList[i].second->getName() ;
            }

        }//end of interaction loop

    }//end of mapping loop
}

defaulttype::BaseMatrix* DefaultMultiMatrixAccessor::createMatrix(const sofa::core::behavior::BaseMechanicalState* mstate1, const sofa::core::behavior::BaseMechanicalState* mstate2)
{
    return createMatrixImpl(mstate1, mstate2, false) ;
}

defaulttype::BaseMatrix* DefaultMultiMatrixAccessor::createMatrixImpl(const sofa::core::behavior::BaseMechanicalState* mstate1, const sofa::core::behavior::BaseMechanicalState* mstate2, bool doPrintInfo)
{
    // The auxiliar interaction matrix is added if and only if at least one of two state is not real state
    //assert(! (realStateOffsets.find(mstate1) != realStateOffsets.end() && realStateOffsets.find(mstate2) != realStateOffsets.end()) );
    FullMatrix<SReal>* m = new FullMatrix<SReal>;
    if(mstate1 == mstate2)
    {
        m->resize( mstate1->getMatrixSize(),mstate1->getMatrixSize());

        if(doPrintInfo)/////////////////////////////////////////////////////////
        {
            msg_info("DefaultMultiMatrixAccessor") << "			++ Creating matrix["<< m->rowSize() <<"."<< m->colSize() <<"]   for mapped state " << mstate1->getName() << "[" << mstate1->getMatrixSize()<<"]" ;
        }
    }
    else
    {
        m->resize( mstate1->getMatrixSize(),mstate2->getMatrixSize() );

        if(doPrintInfo)/////////////////////////////////////////////////////////
        {
            msg_info("DefaultMultiMatrixAccessor") << "			++ Creating interraction matrix["<< m->rowSize() <<"x"<< m->colSize()
                      << "] for interaction " << mstate1->getName() << "[" << mstate1->getMatrixSize()
                      << "] --- "             << mstate2->getName() << "[" << mstate2->getMatrixSize()<<"]" ;
        }
    }

    return m;
}

#ifdef SOFA_SUPPORT_CRS_MATRIX
//TODO separating in other file
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CRSMultiMatrixAccessor::addMechanicalMapping(sofa::core::BaseMapping* mapping)
{
    const sofa::defaulttype::BaseMatrix* jmatrix = mapping->getJ();

    if ((jmatrix != nullptr) && (mapping->isMechanical()) && (mapping->areMatricesMapped()))
    {
        const BaseMechanicalState* mappedState  = const_cast<const BaseMechanicalState*>(mapping->getMechTo()[0]);
        defaulttype::BaseMatrix* mappedstiffness;
        mappedstiffness = mapping->createMappedMatrix(mappedState,mappedState,&CRSMultiMatrixAccessor::createMatrix);
        mappedMatrices[mappedState]=mappedstiffness;

        mappingList.push_back(mapping);

        if(m_doPrintInfo)/////////////////////////////////////////////////////////
        {
            std::cout << "Mapping Visitor : adding validated MechanicalMapping " << mapping->getName()
                    << " with J["<< jmatrix->rowSize()<<"."<<jmatrix->colSize()<<"]" <<std::endl;
        }
    }
    else
    {
        std::cout << "	-- Warning DefaultMultiMatrixAccessor : mapping " << mapping->getName()<<" do not build matrices " << std::endl;
    }
}

defaulttype::BaseMatrix* CRSMultiMatrixAccessor::createMatrix(const sofa::core::behavior::BaseMechanicalState* mstate1, const sofa::core::behavior::BaseMechanicalState* mstate2)
{
    // The auxiliar interaction matrix is added if and only if at least one of two state is not real state
    //assert(! (realStateOffsets.find(mstate1) != realStateOffsets.end() && realStateOffsets.find(mstate2) != realStateOffsets.end()) );

    int nbDOFs1  = mstate1->getSize();
    int dofSize1 = mstate1->getDerivDimension();//getMatrixBlockSize();
    int nbDOFs2  = mstate2->getSize();
    int dofSize2 = mstate2->getDerivDimension();//getMatrixBlockSize();
    //int elementsize = globalMatrix->getElementSize();

    if (mstate1 == mstate2)
    {
        if(m_doPrintInfo)/////////////////////////////////////////////////////////
        {
            std::cout << "			++ Creating matrix Mapped Mechanical State  : "<< mstate1->getName()
                    <<" associated to K["<< mstate1->getMatrixSize() <<"x"<< mstate1->getMatrixSize() << "] in the format _"
                    << nbDOFs1 << "x"<< nbDOFs1 <<"_ of blocs _["
                    << dofSize1 << "x"<< dofSize1 <<"]_"
                    <<std::endl;
        }
        return createBlocSparseMatrix(dofSize1,dofSize1,sizeof(SReal) /*elementsize*/,nbDOFs1,nbDOFs1,MULTIMATRIX_VERBOSE);

    }
    else
    {

        if(m_doPrintInfo)/////////////////////////////////////////////////////////
        {
            std::cout << "			++ Creating matrix Interaction: "
                    << mstate1->getName() <<" -- "<< mstate2->getName()
                    <<" associated to K["<< mstate1->getMatrixSize() <<"x"<< mstate2->getMatrixSize() << "] in the format _"
                    << nbDOFs1 << "x"<< nbDOFs2 <<"_ of blocs _["
                    << dofSize1 << "x"<< dofSize2 <<"]_"
                    <<std::endl;
        }

        return createBlocSparseMatrix(dofSize1,dofSize2,sizeof(SReal) /*elementsize*/,nbDOFs1,nbDOFs2,MULTIMATRIX_VERBOSE);
    }
}

void CRSMultiMatrixAccessor::computeGlobalMatrix()
{
    if(m_doPrintInfo)/////////////////////////////////////////////////////////
    {
        std::cout << "==========================     VERIFICATION BLOC MATRIX FORMATS    ========================" <<std::endl << std::endl;

        for (std::map< const BaseMechanicalState*, MatrixRef >::iterator it = diagonalStiffnessBloc.begin(), itend = diagonalStiffnessBloc.end(); it != itend; ++it)
        {
            std::cout << " Mechanical State  : "<< it->first->getName()
                    <<" associated to K["<< it->second.matrix->rowSize() <<"x"<< it->second.matrix->colSize()<< "] in the format _"
                    << it->second.matrix->bRowSize() << "x"<< it->second.matrix->bColSize() <<"_ of blocs _["
                    << it->second.matrix->getBlockRows() << "x"<< it->second.matrix->getBlockCols() <<"]_"
                    <<std::endl;
        }

        std::map< std::pair<const BaseMechanicalState*, const BaseMechanicalState*>, InteractionMatrixRef >::iterator itBegin = interactionStiffnessBloc.begin();
        std::map< std::pair<const BaseMechanicalState*, const BaseMechanicalState*>, InteractionMatrixRef >::iterator itEnd = interactionStiffnessBloc.end();

        while(itBegin != itEnd)
        {
            std::cout << " Interaction: "
                    << itBegin->first.first->getName() <<" -- "<< itBegin->first.second->getName()
                    <<" associated to K["<< itBegin->second.matrix->rowSize() <<"x"<< itBegin->second.matrix->colSize()<< "] in the format _"
                    << itBegin->second.matrix->bRowSize() << "x"<< itBegin->second.matrix->bColSize() <<"_ of blocs _["
                    << itBegin->second.matrix->getBlockRows() << "x"<< itBegin->second.matrix->getBlockCols() <<"]_"
                    <<std::endl;

            ++itBegin;
        }

        const int lastMappingId = mappingList.size() - 1;
        for(int id=lastMappingId; id>=0; --id)
        {
            sofa::core::BaseMapping* m_mapping = mappingList[id];
            const BaseMechanicalState* instate  = const_cast<const BaseMechanicalState*>(m_mapping->getMechFrom()[0]);
            const BaseMechanicalState* outstate  = const_cast<const BaseMechanicalState*>(m_mapping->getMechTo()[0]);
            const defaulttype::BaseMatrix* matrixJ = m_mapping->getJ();

            std::cout << "  "<<id<< "-th Mapping : " <<m_mapping->getName()<< " associated to matrix J["
                    << matrixJ->rowSize() <<"x"<< matrixJ->colSize()<< "] in the format _"
                    << matrixJ->bRowSize() << "x"<< matrixJ->bColSize() <<"_ of blocs _["
                    << matrixJ->getBlockRows() << "x"<< matrixJ->getBlockCols() <<"]_"
                    <<std::endl;

            std::cout << "			inState  : "<< instate->getName()
                    <<" associated to K11["<< instate->getMatrixSize() <<"x"<< instate->getMatrixSize() << "] in the format _"
                    << instate->getSize() << "x"<< instate->getSize() <<"_ of blocs _["
                    << instate->getDerivDimension() << "x"<< instate->getDerivDimension() <<"]_"
                    <<std::endl;

            std::cout << "			outState  : "<< outstate->getName()
                    <<" associated to K11["<< outstate->getMatrixSize() <<"x"<< outstate->getMatrixSize() << "] in the format _"
                    << outstate->getSize() << "x"<< outstate->getSize() <<"_ of blocs _["
                    << outstate->getDerivDimension() << "x"<< outstate->getDerivDimension() <<"]_"
                    <<std::endl;
        }
        std::cout <<std::endl << "=======================     CONTRIBUTION CONTRIBUTION CONTRIBUTION     ======================" <<std::endl << std::endl;
    }


    const int lastMappingId = mappingList.size() - 1;
    for(int id=lastMappingId; id>=0; --id)
    {
        sofa::core::BaseMapping* m_mapping = mappingList[id];
        const BaseMechanicalState* instate  = const_cast<const BaseMechanicalState*>(m_mapping->getMechFrom()[0]);
        const BaseMechanicalState* outstate  = const_cast<const BaseMechanicalState*>(m_mapping->getMechTo()[0]);

        const defaulttype::BaseMatrix* matrixJ = m_mapping->getJ();
        const unsigned int nbR_J = matrixJ->rowSize();
        const unsigned int nbC_J = matrixJ->colSize();

        if(m_doPrintInfo)/////////////////////////////////////////////////////////
        {
            std::cout << "[CONTRIBUTION] of " << id << "-th mapping named : " << m_mapping->getName()
                    << " with fromModel : "<< instate->getName()
                    << " and toModel : "   << outstate->getName() <<" -- with matrix J["<<nbR_J <<"."<<nbC_J <<"]" <<std::endl;
        }


        //for toModel -----------------------------------------------------------
        if(mappedMatrices.find(outstate) == mappedMatrices.end())
        {
            if(m_doPrintInfo)
                std::cout << "	[Propa.Stiff] WARNING toModel : "<< outstate->getName()<< " dont have stiffness matrix"<<std::endl;
        }



        if(diagonalStiffnessBloc.find(outstate) != diagonalStiffnessBloc.end())
        {
            //=================================
            //           K11 += Jt * K22 * J
            //=================================
            MatrixRef K1 = this->getMatrix(instate);
            MatrixRef K2 = this->getMatrix(outstate);

            const unsigned int offset1  = K1.offset;
            const unsigned int offset2  = K2.offset;
            const unsigned int sizeK1 = K1.matrix->rowSize() - offset1;
            const unsigned int sizeK2 = K2.matrix->rowSize() - offset2;

            defaulttype::BaseMatrix* matrixJJ =const_cast<defaulttype::BaseMatrix*>(matrixJ);

            int JblocRsize     = matrixJJ->getBlockRows();
            int JblocCsize     = matrixJJ->getBlockCols();
            int JnbBlocCol     = matrixJJ->bColSize();

            int K2blocCsize    = K2.matrix->getBlockCols();
            int K2nbBlocCol   = K2.matrix->bColSize();

            int JelementSize   = matrixJJ->getElementSize();
            int MelementSize   = K2.matrix->getElementSize();

            // creating a tempo matrix  tempoMatrix
            defaulttype::BaseMatrix* tempoMatrix = createBlocSparseMatrix(JblocCsize, K2blocCsize, MelementSize, JnbBlocCol, K2nbBlocCol,MULTIMATRIX_VERBOSE);
            // Matrix multiplication  tempoMatrix += Jt * K22
            opAddMulJTM(tempoMatrix,   matrixJJ, K2.matrix,       0,0      , JblocRsize, JblocCsize , K2blocCsize, JelementSize,MelementSize,MULTIMATRIX_VERBOSE);
            // Matrix multiplication  K11         += tempoMatrix * J
            opAddMulMJ( K1.matrix  ,tempoMatrix,  matrixJJ, offset1,offset1, JblocCsize, K2blocCsize, JblocCsize , JelementSize,MelementSize,MULTIMATRIX_VERBOSE);

            delete tempoMatrix;



            if(m_doPrintInfo)/////////////////////////////////////////////////////////
            {
                std::cout << "	[Propa.Stiff] propagating stiffness of : "<< outstate->getName()<< " to stifness "<<instate->getName()<<std::endl;
                std::cout <<"	                    **multiplication of "
                        <<" K1["<<sizeK1<<"."<<sizeK1<< "]("<< offset1<<","<<offset1 <<  ")   =  "
                        <<" Jt["<<nbC_J<<"."<<nbR_J<< "] * "
                        <<" K2["<<sizeK2<<"."<<sizeK2<<"]("<< offset2<<","<<offset2 <<  ") * "
                        <<" J["<<nbR_J<<"."<<nbC_J<< "]"
                        <<std::endl;
            }
        }

        std::vector<std::pair<const BaseMechanicalState*, const BaseMechanicalState*> > interactionList;
        for (std::map< std::pair<const BaseMechanicalState*, const BaseMechanicalState*>, InteractionMatrixRef >::iterator it = interactionStiffnessBloc.begin(), itend = interactionStiffnessBloc.end(); it != itend; ++it)
        {
            if(it->first.first == outstate || it->first.second == outstate )
            {
                interactionList.push_back(it->first);
            }
            ++it;
        }


        const unsigned nbInteraction = interactionList.size();
        for(unsigned i=0; i< nbInteraction; i++)
        {

            if(m_doPrintInfo)/////////////////////////////////////////////////////////
            {
                std::cout << "	[Propa.Interac.Stiff] detected interaction between toModel : "<<interactionList[i].first->getName()
                        << " and : " <<interactionList[i].second->getName()<<std::endl;
            }
            //                   |       |
            //                 MS1     MS2
            //                  /      /
            //                map    inter
            //                   \   /
            //                   MS3/
            //
            //           K_11 += Jt * K_33 * J
            //           I_12 += Jt * I_32
            //           I_21 +=      I_23 * J
            //

            if(interactionList[i].first == outstate)
            {
                InteractionMatrixRef I_32 = this->getMatrix( outstate , interactionList[i].second);
                InteractionMatrixRef I_12 = this->getMatrix( instate  , interactionList[i].second);
                //===========================
                //          I_12 += Jt * I_32
                //===========================
                const unsigned int offR_I_12  = I_12.offRow;                       //      row offset of I12 matrix
                const unsigned int offC_I_12  = I_12.offCol;                       //    colum offset of I12 matrix
                const unsigned int nbR_I_12   = I_12.matrix->rowSize() - offR_I_12;//number of rows   of I12 matrix
                const unsigned int nbC_I_12   = I_12.matrix->colSize() - offC_I_12;//number of colums of I12 matrix

                const unsigned int offR_I_32  = I_32.offRow;                     //      row offset of I32 matrix
                const unsigned int offC_I_32  = I_32.offCol;                     //    colum offset of I32 matrix
                const unsigned int nbR_I_32 = I_32.matrix->rowSize() - offR_I_32;//number of rows   of I32 matrix
                const unsigned int nbC_I_32 = I_32.matrix->colSize() - offC_I_32;//number of colums of I32 matrix


                if(m_doPrintInfo)/////////////////////////////////////////////////////////
                {
                    std::cout << "	[Propa.Interac.Stiff] propagating interaction "
                            <<outstate->getName() << "--" <<interactionList[i].second->getName()
                            <<"  to : "
                            <<instate->getName() << "--" <<interactionList[i].second->getName()<<std::endl;

                    std::cout <<"	                    **multiplication of "
                            <<" I12["<<nbR_I_12<<"."<<nbC_I_12<< "]("<< offR_I_12<<","<<offC_I_12 <<  ")  =  "
                            <<" Jt["<<nbC_J<<"."<<nbR_J<< "]  *  "
                            <<" I32["<<nbR_I_32<<"."<<nbC_I_32<<"]("<< offR_I_32<<","<<offC_I_32 <<  ")" <<std::endl;
                }

                // Matrix multiplication   I_12 += Jt * I_32

                defaulttype::BaseMatrix* matrixJJ =const_cast<defaulttype::BaseMatrix*>(matrixJ);
                int JblocRsize     = matrixJJ->getBlockRows();
                int JblocCsize     = matrixJJ->getBlockCols();
                int I_32_blocCsize = I_32.matrix->getBlockCols();
                int JelementSize   = matrixJJ->getElementSize();
                int MelementSize   = I_32.matrix->getElementSize();

                opAddMulJTM(I_12.matrix,matrixJJ,I_32.matrix,offR_I_12,offC_I_12, JblocRsize, JblocCsize,I_32_blocCsize,JelementSize,MelementSize,MULTIMATRIX_VERBOSE);
            }

            if(interactionList[i].second == outstate)
            {
                InteractionMatrixRef I_21 = this->getMatrix(interactionList[i].first,instate);
                InteractionMatrixRef I_23 = this->getMatrix(interactionList[i].first,outstate);
                //=========================================
                //          I_21 +=      I_23 * J
                //=========================================
                const unsigned int offR_I_21  = I_21.offRow;
                const unsigned int offC_I_21  = I_21.offCol;
                const unsigned int nbR_I_21 = I_21.matrix->rowSize() - offR_I_21;
                const unsigned int nbC_I_21 = I_21.matrix->colSize() - offC_I_21;

                const unsigned int offR_I_23  = I_23.offRow;
                const unsigned int offC_I_23  = I_23.offCol;
                const unsigned int nbR_I_23   = I_23.matrix->rowSize() - offR_I_23;
                const unsigned int nbC_I_23   = I_23.matrix->colSize() - offC_I_23;

                if(m_doPrintInfo)/////////////////////////////////////////////////////////
                {
                    std::cout << "	[Propa.Interac.Stiff] propagating interaction "
                            <<interactionList[i].first->getName()<< "--" <<outstate->getName()
                            <<" to : "<<interactionList[i].first->getName()<< "--" <<instate->getName()<<std::endl;

                    std::cout <<"	                    **multiplication of "
                            <<" I_21["<<nbR_I_21<<"."<<nbC_I_21<<"]("<< offR_I_21<<","<<offC_I_21 <<  ")  =  "
                            <<" I23["<<nbR_I_23<<"."<<nbC_I_23<< "]("<< offR_I_23<<","<<offC_I_23 <<  ")  *  "
                            <<" J["<<nbR_J<<"."<<nbC_J<<"]" <<std::endl;
                }


                // Matrix multiplication  I_21 +=  I_23 * J
                defaulttype::BaseMatrix* matrixJJ =const_cast<defaulttype::BaseMatrix*>(matrixJ);
                int I_23_blocRsize = I_23.matrix->getBlockRows();
                int I_23_blocCsize = I_23.matrix->getBlockCols();
                int JblocCsize     = matrixJJ->getBlockCols();
                int JelementSize   = matrixJJ->getElementSize();
                int MelementSize   = I_23.matrix->getElementSize();

                opAddMulMJ(I_21.matrix,I_23.matrix,matrixJJ,offR_I_21,offC_I_21, I_23_blocRsize,I_23_blocCsize,JblocCsize,JelementSize,MelementSize,MULTIMATRIX_VERBOSE);

            }

            //after propagating the interaction, we remove the older interaction
            interactionStiffnessBloc.erase( interactionStiffnessBloc.find(interactionList[i]) );
            if(m_doPrintInfo)
            {
                std::cout << "	--[Propa.Interac.Stiff] remove interaction of : "<<interactionList[i].first->getName()
                        << " and : " << interactionList[i].second->getName()<<std::endl;
            }

        }//end of interaction loop

    }//end of mapping loop
}

#endif


} // namespace sofa::component::linearsolver
