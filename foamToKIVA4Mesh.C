/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    foamToStarMesh

Description
    Reads an OpenFOAM mesh and writes a pro-STAR (v4) bnd/cel/vrt format.

Usage
    \b foamToStarMesh [OPTION]

    Options:
      - \par -noBnd
        Suppress writing the \c .bnd file

      - \par -scale \<factor\>
        Specify an alternative geometry scaling factor.
        The default is \b 1000 (scale \em [m] to \em [mm]).

Note
    The cellTable information available in the files
    \c constant/cellTable and \c constant/polyMesh/cellTableId
    will be used if available. Otherwise the cellZones are used when
    creating the cellTable information.

See also
    Foam::cellTable, Foam::meshWriter and Foam::meshWriters::STARCD

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "KIVA4MeshWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "read OpenFOAM mesh and write a pro-STAR (v4) bnd/cel/vrt format"
    );
    argList::noParallel();
    timeSelector::addOptions();

    argList::addOption
    (
        "scale",
        "factor",
        "geometry scaling factor - default is 1000 ([m] to [mm])"
    );
    argList::addBoolOption
    (
        "noBnd",
        "suppress writing a boundary (.bnd) file"
    );

    #include "setRootCase.H"  
    //https://github.com/OpenFOAM/OpenFOAM-4.x/blob/master/src/OpenFOAM/include/setRootCase.H
    #include "createTime.H"   
    //https://github.com/OpenFOAM/OpenFOAM-4.x/blob/master/src/OpenFOAM/include/createTime.H
    
    instantList timeDirs = timeSelector::select0(runTime, args);

    fileName exportName = "kiva4grid";//meshWriter::defaultMeshName;
    if (args.optionFound("case"))
    {
        exportName += '-' + args.globalCaseName();
    }

    // default: do not rescale// from [m] to [mm]
    scalar scaleFactor = 1;
    if (args.optionReadIfPresent("scale", scaleFactor))
    {
        if (scaleFactor <= 0)
        {
            scaleFactor = 1;
        }
    }
     
    //createMesh.H hace lo mismo pero para una fvmesh
    #include "createMesh.H"    //Crea la malla leyendo los archivos de OFoam
    //https://github.com/OpenFOAM/OpenFOAM-4.x/blob/master/src/OpenFOAM/include/createPolyMesh.H
   


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        #include "getTimeIndex.H"
        //Archivo incluido en la misma carpeta que este código fuente. No confundir con otro!!

        fvMesh::readUpdateState state = mesh.readUpdate();

        if (!timeI || state != fvMesh::UNCHANGED) //Ni idea de por qué mira que la malla no haya cambiado
        {
            meshWriters::STARCD writer(mesh, scaleFactor);

            if (args.optionFound("noBnd"))
            {
                //writer.noBoundary(); //Deshabilitado, no sé qué hace
            }
 
            fileName meshName(exportName);
            if (state != fvMesh::UNCHANGED)
            {
                //meshName += '_' + runTime.timeName(); //Deshabilitado, no sé qué hace
            }

            writer.write(meshName);
        }

        Info<< nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
