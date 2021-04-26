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
#include <SofaLoader/config.h>

#include <sofa/core/objectmodel/BaseData.h>
#include <sofa/core/loader/MeshLoader.h>

#include <iostream>
#include <cstdio>
#include <sstream>

#include <tinyxml.h>

#include <SofaLoader/BaseVTKReader.h>
using sofa::component::loader::BaseVTKReader;

//namespace sofa::component::loader::basevtkreader
//{
//    class BaseVTKReader ;
//}

namespace sofa::component::loader
{
using basevtkreader::BaseVTKReader ;


class SOFA_SOFALOADER_API LegacyVTKReader : public BaseVTKReader
{
public:
    bool readFile(const char* filename) override;
};

class SOFA_SOFALOADER_API XMLVTKReader : public BaseVTKReader
{
public:
    bool readFile(const char* filename) override;
protected:
    bool loadUnstructuredGrid(TiXmlHandle datasetFormatHandle);
    bool loadPolydata(TiXmlHandle datasetFormatHandle);
    bool loadRectilinearGrid(TiXmlHandle datasetFormatHandle);
    bool loadStructuredGrid(TiXmlHandle datasetFormatHandle);
    bool loadStructuredPoints(TiXmlHandle datasetFormatHandle);
    bool loadImageData(TiXmlHandle datasetFormatHandle);
    BaseVTKDataIO* loadDataArray(TiXmlElement* dataArrayElement, int size, std::string type);
    BaseVTKDataIO* loadDataArray(TiXmlElement* dataArrayElement, int size);
    BaseVTKDataIO* loadDataArray(TiXmlElement* dataArrayElement);
};

/// Format doc: http://www.vtk.org/VTK/img/file-formats.pdf
/// http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html
class SOFA_SOFALOADER_API MeshVTKLoader : public sofa::core::loader::MeshLoader
{
    friend class ExtendedVTKLoader;
public:
    SOFA_CLASS(MeshVTKLoader,sofa::core::loader::MeshLoader);

public:
    core::objectmodel::BaseData* pointsData;
    core::objectmodel::BaseData* edgesData;
    core::objectmodel::BaseData* trianglesData;
    core::objectmodel::BaseData* quadsData;
    core::objectmodel::BaseData* tetrasData;
    core::objectmodel::BaseData* hexasData;

    bool doLoad() override;

protected:
    enum VTKFileType { NONE, LEGACY, XML };

    MeshVTKLoader();

    VTKFileType detectFileType(const char* filename);

    BaseVTKReader* reader;
    bool setInputsMesh();
    bool setInputsData();

    void doClearBuffers() override;
};

} /// namespace sofa::component::loader
