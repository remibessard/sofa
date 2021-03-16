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
#include <SofaBaseLinearSolver/CGLinearSolver.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/helper/AdvancedTimer.h>


namespace sofa::component::linearsolver
{

/// Linear system solver using the conjugate gradient iterative algorithm
template<class TMatrix, class TVector>
CGLinearSolver<TMatrix,TVector>::CGLinearSolver()
    : f_maxIter( initData(&f_maxIter,(unsigned)25,"iterations","Maximum number of iterations of the Conjugate Gradient solution") )
    , f_tolerance( initData(&f_tolerance,(SReal)1e-5,"tolerance","Desired accuracy of the Conjugate Gradient solution (ratio of current residual norm over initial residual norm)") )
    , f_smallDenominatorThreshold( initData(&f_smallDenominatorThreshold,(SReal)1e-5,"threshold","Minimum value of the denominator in the conjugate Gradient solution") )
    , f_warmStart( initData(&f_warmStart,false,"warmStart","Use previous solution as initial solution") )
    , f_verbose( initData(&f_verbose,false,"verbose","Dump system state at each iteration") )
    , f_graph( initData(&f_graph,"graph","Graph of residuals at each iteration") )
{
    f_graph.setWidget("graph");
    f_maxIter.setRequired(true);
    f_tolerance.setRequired(true);
    f_smallDenominatorThreshold.setRequired(true);
    this->isIterativeSolver = true;
}

template<class TMatrix, class TVector>
void CGLinearSolver<TMatrix,TVector>::init()
{
    if(f_verbose.getValue())
    {
        this->f_printLog.setValue(true);
    }

    if(f_tolerance.getValue() < 0.0)
    {
        msg_warning() << "'tolerance' must be a positive value" << msgendl
                      << "default value used: 1e-5";
        f_tolerance.setValue(1e-5);
    }
    if(f_smallDenominatorThreshold.getValue() < 0.0)
    {
        msg_warning() << "'threshold' must be a positive value" << msgendl
                      << "default value used: 1e-5";
        f_smallDenominatorThreshold.setValue(1e-5);
    }

    timeStepCount = 0;
    equilibriumReached = false;
}

template<class TMatrix, class TVector>
void CGLinearSolver<TMatrix,TVector>::reinit()
{
    if(f_verbose.getValue())
    {
        this->f_printLog.setValue(true);
    }
}

template<class TMatrix, class TVector>
void CGLinearSolver<TMatrix,TVector>::resetSystem()
{
    f_graph.beginEdit()->clear();
    f_graph.endEdit();

    Inherit::resetSystem();
}

template<class TMatrix, class TVector>
void CGLinearSolver<TMatrix,TVector>::setSystemMBKMatrix(const sofa::core::MechanicalParams* mparams)
{
    sofa::helper::ScopedAdvancedTimer("CG-setSystemMBKMatrix");
    Inherit::setSystemMBKMatrix(mparams);
}

/// Solve Mx=b
template<class TMatrix, class TVector>
void CGLinearSolver<TMatrix,TVector>::solve(Matrix& M, Vector& x, Vector& b)
{
#ifdef SOFA_DUMP_VISITOR_INFO
    simulation::Visitor::printComment("ConjugateGradient");
#endif

#ifdef SOFA_DUMP_VISITOR_INFO
    simulation::Visitor::printNode("VectorAllocation");
#endif

    const core::ExecParams* params = core::ExecParams::defaultInstance();
    typename Inherit::TempVectorContainer vtmp(this, params, M, x, b);
    Vector& p = *vtmp.createTempVector();
    Vector& q = *vtmp.createTempVector();
    Vector& r = *vtmp.createTempVector();

    const bool verbose  = f_verbose.getValue();
    double rho, rho_1=0, alpha, beta;


    msg_info_when(verbose) << "b = " << b ;


    /// Compute the initial residual r
    if( f_warmStart.getValue() )
    {
        r = M * x;
        r.eq( b, r, -1.0 );   // initial residual r = b - Ax;
    }
    else
    {
        x.clear();
        r = b; // initial residual
    }

    /// Compute the norm of the right-hand-side vector b
    double normb = b.norm();


    std::map < std::string, sofa::helper::vector<SReal> >& graph = *f_graph.beginEdit();
    sofa::helper::vector<SReal>& graph_error = graph[std::string("Error")];
    graph_error.clear();
    sofa::helper::vector<SReal>& graph_den = graph[std::string("Denominator")];
    graph_den.clear();
    graph_error.push_back(1);
    unsigned nb_iter = 0;
    const char* endcond = "iterations";

    sofa::helper::AdvancedTimer::stepBegin("CG-Solve");

#ifdef SOFA_DUMP_VISITOR_INFO
    simulation::Visitor::printCloseNode("VectorAllocation");
#endif

    if(normb != 0.0)
    {
        for( nb_iter=1; nb_iter<=f_maxIter.getValue(); nb_iter++ )
        {
#ifdef SOFA_DUMP_VISITOR_INFO
            std::ostringstream comment;
            if (simulation::Visitor::isPrintActivated())
            {
                comment << "Iteration_" << nb_iter;
                simulation::Visitor::printNode(comment.str());
            }
#endif

            /// Compute p = r^2
            rho = r.dot(r);

            /// Compute the error from the norm of ρ and b
            double normr = sqrt(rho);
            double err = normr/normb;

            graph_error.push_back(err);


            /// Break condition = TOLERANCE criterion regarding the error err is reached
            if (err <= f_tolerance.getValue())
            {
                /// Tolerance met at first step, tolerance value might not be relevant
                if(nb_iter == 1 && timeStepCount == 0)
                {
                    msg_warning() << "tolerance reached at first iteration of CG" << msgendl
                                  << "Check the 'tolerance' data field, you might decrease it";
                }
                else
                {
                    if(nb_iter == 1 && !equilibriumReached)
                    {
                        msg_info() << "Equilibrium reached regarding tolerance";
                        equilibriumReached = true;
                    }
                    if(nb_iter > 1)
                    {
                        equilibriumReached = false;
                    }

                    endcond = "tolerance";
                    if( verbose )
                    {
                        msg_info() << "error = " << err <<", tolerance = " << f_tolerance.getValue();
                    }

#ifdef SOFA_DUMP_VISITOR_INFO
                    if (simulation::Visitor::isPrintActivated())
                        simulation::Visitor::printCloseNode(comment.str());
#endif
                    break;
                }
            }


            /// Compute the value of p, conjugate with x
            if( nb_iter==1 )    // FIRST step
                p = r;
            else                // ALL other steps
            {
                beta = rho / rho_1;

                /// Update p = p*beta + r;
                cgstep_beta(params, p,r,beta);
            }

            if( verbose )
            {
                msg_info() << "p : " << p;
            }

            /// Compute the matrix-vector product : M p
            q = M*p;

            if( verbose )
            {
                msg_info() << "q = M p : " << q;
            }

            /// Compute the denominator : p M p
            double den = p.dot(q);

            graph_den.push_back(den);

            if(den != 0.0)
            {
                /// Break condition = THRESHOLD criterion regarding the denominator is reached (but do at least one iteration)
                if (fabs(den) <= f_smallDenominatorThreshold.getValue())
                {
                    /// Threshold met at first step, threshold value might not be relevant
                    if(nb_iter == 1 && timeStepCount == 0)
                    {
                        msg_warning() << "denominator threshold reached at first iteration of CG" << msgendl
                                      << "Check the 'threshold' data field, you might decrease it";
                    }
                    else
                    {
                        if(nb_iter == 1 && !equilibriumReached)
                        {
                            msg_info() << "Equilibrium reached regarding threshold";
                            equilibriumReached = true;
                        }
                        if(nb_iter > 1)
                        {
                            equilibriumReached = false;
                        }

                        endcond = "threshold";
                        if( verbose )
                        {
                            msg_info() << "den = " << den <<", smallDenominatorThreshold = " << f_smallDenominatorThreshold.getValue() <<", err = " << err;
                        }

#ifdef SOFA_DUMP_VISITOR_INFO
                    if (simulation::Visitor::isPrintActivated())
                        simulation::Visitor::printCloseNode(comment.str());
#endif
                        break;
                    }
                }


                /// Compute the coefficient α for the conjugate direction
                alpha = rho/den;

                /// End of the CG step : update x and r
                cgstep_alpha(params, x,r,p,q,alpha);

                if( verbose )
                {
                    msg_info() << "den = " << den << ", alpha = " << alpha << ", x = " << x << ", r = " << r;
                }
            }
            else
            {
                msg_warning() << "den = 0.0, break the iterations";
                break;
            }

            rho_1 = rho;

#ifdef SOFA_DUMP_VISITOR_INFO
            if (simulation::Visitor::isPrintActivated())
                simulation::Visitor::printCloseNode(comment.str());
#endif
        }
    }
    // Case no forces applied, b=0
    else
    {
        endcond = "null norm of vector b";

        // If first step : check the value of threshold
        if( timeStepCount==0 )
        {
            p = r;
            q = M*p;
            double den = p.dot(q);

            if(den != 0.0)
            {
                if (fabs(den) <= f_smallDenominatorThreshold.getValue())
                {
                    msg_warning() << "denominator threshold reached at first iteration of CG" << msgendl
                                  << "Check the 'threshold' data field, you might decrease it";
                }
            }
            else
            {
                msg_info() << "no way to check the validity of : tolerance and threshold value";
            }
        }
    }

    sofa::helper::AdvancedTimer::stepEnd("CG-Solve");

    f_graph.endEdit();
    timeStepCount ++;

    sofa::helper::AdvancedTimer::valSet("CG iterations", nb_iter);

    dmsg_info() << "solve, nbiter = "<<nb_iter<<" stop because of "<<endcond;
    dmsg_info_when( verbose ) <<"solve, solution = "<< x ;

    vtmp.deleteTempVector(&p);
    vtmp.deleteTempVector(&q);
    vtmp.deleteTempVector(&r);
}

template<class TMatrix, class TVector>
inline void CGLinearSolver<TMatrix,TVector>::cgstep_beta(const core::ExecParams* /*params*/, Vector& p, Vector& r, SReal beta)
{
    // p = p*beta + r
    p *= beta;
    p += r;
}

template<class TMatrix, class TVector>
inline void CGLinearSolver<TMatrix,TVector>::cgstep_alpha(const core::ExecParams* /*params*/, Vector& x, Vector& r, Vector& p, Vector& q, SReal alpha)
{
    // x = x + alpha p
    x.peq(p,alpha);

    // r = r - alpha q
    r.peq(q,-alpha);
}

} // namespace sofa::component::linearsolver
