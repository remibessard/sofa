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
#pragma once

#include <sofa/gui/BaseGUI.h>
#include <sofa/simulation/fwd.h>
#include <sofa/helper/ArgumentParser.h>
#include <string>

#include <QtCore/QObject>
#include <QMainWindow>

using sofa::helper::ArgumentParser;

namespace sofa::gui
{

class SOFA_SOFAGUICOMMON_API BatchGUI : public BaseGUI
{

public:

    /// @name methods each GUI must implement
    /// @{

    BatchGUI();

    void setScene(sofa::simulation::NodeSPtr groot, const char* filename="", bool temporaryFile=false) override;

    void resetScene();

    int mainLoop() override;
    void redraw() override;
    int closeGUI() override;

	void playpauseGUI(bool value);
	void switchPlaypauseGUI();

    static void setNumIterations(const std::string& nbIterInp) 
    {
        size_t inpLen= nbIterInp.length();
       
        if (nbIterInp == "infinite")
        {
            nbIter = -1;
        }
        else if (inpLen)
        {
            nbIter = std::stoi(nbIterInp);
        }
        else
        {
            nbIter = DEFAULT_NUMBER_OF_ITERATIONS;
        }
        
    }
    sofa::simulation::Node* currentSimulation() override;

    /// @}

    /// @name registration of each GUI
    /// @{

    static BaseGUI* CreateGUI(const char* name, sofa::simulation::NodeSPtr groot = nullptr, const char* filename = nullptr);
    static int RegisterGUIParameters(ArgumentParser* argumentParser);


    static const signed int DEFAULT_NUMBER_OF_ITERATIONS;
    /// @}

protected:
    /// The destructor should not be called directly. Use the closeGUI() method instead.
    ~BatchGUI() override;

    void startDumpVisitor();
    void stopDumpVisitor();

    std::ostringstream m_dumpVisitorStream;

    sofa::simulation::NodeSPtr groot;
    std::string filename;
    static signed int nbIter;
    static std::string nbIterInp;

	bool m_animated;
};

} // namespace sofa::gui
