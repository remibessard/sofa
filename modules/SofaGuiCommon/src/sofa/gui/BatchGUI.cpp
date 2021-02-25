/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program. If not, see <http://www.gnu.org/licenses/>.              *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include "BatchGUI.h"



#include <sofa/simulation/Simulation.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/simulation/UpdateContextVisitor.h>
#include <sofa/helper/system/thread/CTime.h>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/program_options.hpp>
#include <sofa/core/objectmodel/IdleEvent.h>

#include <QCoreApplication>
#include <QApplication>

#include <sofa/simulation/events/SimulationStartEvent.h>
using sofa::simulation::SimulationStartEvent;

namespace sofa
{

namespace gui
{

const signed int BatchGUI::DEFAULT_NUMBER_OF_ITERATIONS = 1000;
signed int BatchGUI::nbIter = BatchGUI::DEFAULT_NUMBER_OF_ITERATIONS;
std::string BatchGUI::nbIterInp="";
BatchGUI::BatchGUI()
    : groot(nullptr)
	, m_animated(false)
{
}

BatchGUI::~BatchGUI()
{
}

class QSOFAApplication : public QCoreApplication
{
public:
	QSOFAApplication(int &argc, char ** argv)
		: QCoreApplication(argc, argv)
	{
		QCoreApplication::setOrganizationName("Sofa Consortium");
		QCoreApplication::setOrganizationDomain("sofa");
		QCoreApplication::setApplicationName("runSofa");
	}

protected:
	bool event(QEvent *event) override
	{
		switch (event->type())
		{
		default:
			return QCoreApplication::event(event);
		}
	}
};
QCoreApplication* application = nullptr;

int BatchGUI::mainLoop()
{
	std::cout << "in mainloop" << std::endl;

	if (groot)
    {   
        if (nbIter != -1)
        {   
            msg_info("BatchGUI") << "Computing " << nbIter << " iterations." << msgendl;
        }
        else
        {
            msg_info("BatchGUI") << "Computing infinite iterations." << msgendl;
        }

        sofa::helper::AdvancedTimer::begin("Animate");
        sofa::simulation::getSimulation()->animate(groot.get());

		currentSimulation()->getContext()->setAnimate(false);
        msg_info("BatchGUI") << "Processing." << sofa::helper::AdvancedTimer::end("Animate", groot.get()) << msgendl;
        sofa::simulation::Visitor::ctime_t rtfreq = sofa::helper::system::thread::CTime::getRefTicksPerSec();
        sofa::simulation::Visitor::ctime_t tfreq = sofa::helper::system::thread::CTime::getTicksPerSec();
        sofa::simulation::Visitor::ctime_t rt = sofa::helper::system::thread::CTime::getRefTime();
        sofa::simulation::Visitor::ctime_t t = sofa::helper::system::thread::CTime::getFastTime();
          
        signed int i = 1; //one simulatin step is animated above  
       

        while (i <= nbIter || nbIter == -1)
        {
			sofa::core::objectmodel::IdleEvent hb;
			if (groot)
			{
				groot->propagateEvent(core::ExecParams::defaultInstance(), &hb);
			}

			if (nbIter == -1 &&!m_animated)
				continue;

            if (i != nbIter)
            {
				application->processEvents();
                sofa::helper::ScopedAdvancedTimer("Animate");
                sofa::helper::AdvancedTimer::begin("Animate");
                sofa::simulation::getSimulation()->animate(groot.get());
                sofa::helper::AdvancedTimer::end("Animate");
            }

            if ( i == nbIter /*|| (nbIter == -1 && i%1000 == 0)*/ )
            {
                t = sofa::helper::system::thread::CTime::getFastTime()-t;
                rt = sofa::helper::system::thread::CTime::getRefTime()-rt;

                msg_info("BatchGUI") << i << " iterations done in " << ((double)t)/((double)tfreq) << " s ( " << (((double)tfreq)*i)/((double)t) << " FPS)." << msgendl;
                msg_info("BatchGUI") << i << " iterations done in " << ((double)rt)/((double)rtfreq) << " s ( " << (((double)rtfreq)*i)/((double)rt) << " FPS)." << msgendl;
                
                if (nbIter == -1) // Additional message for infinite iterations
                {
                     msg_info("BatchGUI") << "Press Ctrl + C (linux)/ Command + period (mac) to stop " << msgendl;
                }
            }

            i++;
        }
        
    }
    return 0;
}

void BatchGUI::redraw()
{
}

int BatchGUI::closeGUI()
{
    delete this;
    return 0;
}

void BatchGUI::playpauseGUI(bool value)
{
	if (groot)
		groot->getContext()->setAnimate(value);
	m_animated = value;
}

void BatchGUI::switchPlaypauseGUI()
{
	playpauseGUI(!m_animated);
}

void BatchGUI::setScene(sofa::simulation::Node::SPtr groot, const char* filename, bool )
{
    this->groot = groot;
    this->filename = (filename?filename:"");

    resetScene();
}


void BatchGUI::resetScene()
{
    sofa::simulation::Node* root = currentSimulation();

    if ( root )
    {
        root->setTime(0.);
        simulation::getSimulation()->reset ( root );

        sofa::simulation::UpdateSimulationContextVisitor(sofa::core::ExecParams::defaultInstance()).execute(root);
    }
}

void BatchGUI::startDumpVisitor()
{
#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Node* root = currentSimulation();
    if (root)
    {
        m_dumpVisitorStream.str("");
        sofa::simulation::Visitor::startDumpVisitor(&m_dumpVisitorStream, root->getTime());
    }
#endif
}

void BatchGUI::stopDumpVisitor()
{
#ifdef SOFA_DUMP_VISITOR_INFO
    sofa::simulation::Visitor::stopDumpVisitor();
    m_dumpVisitorStream.flush();
    m_dumpVisitorStream.str("");
#endif
}

sofa::simulation::Node* BatchGUI::currentSimulation()
{
    return groot.get();
}


BaseGUI* BatchGUI::CreateGUI(const char* name, sofa::simulation::Node::SPtr groot, const char* filename)
{
	QCoreApplication coreApplication();
	int  *argc = new int;
	char **argv = new char*[2];
	*argc = 1;
	argv[0] = strdup(BaseGUI::GetProgramName());
	argv[1] = nullptr;
	application = new QSOFAApplication(*argc, argv);

	//force locale to Standard C
	//(must be done immediatly after the QApplication has been created)
	QLocale locale(QLocale::C);
	QLocale::setDefault(locale);

    BatchGUI::mGuiName = name;
    BatchGUI* gui = new BatchGUI();
    gui->setScene(groot, filename);
    return gui;
}

int BatchGUI::RegisterGUIParameters(ArgumentParser* argumentParser)
{
    argumentParser->addArgument(
                boost::program_options::value<std::string>()
                ->notifier(setNumIterations),
                "nbIter,n",
                "(only batch) Number of iterations of the simulation"
                );
    //Parses the string and passes it to setNumIterations as argument
    return 0;
}

} // namespace gui

} // namespace sofaa
