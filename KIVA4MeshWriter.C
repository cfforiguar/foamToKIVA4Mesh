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

\*---------------------------------------------------------------------------*/

#include "KIVA4MeshWriter.H"

#include "Time.H"
#include "SortableList.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::meshWriters::STARCD::defaultBoundaryName =
    "Default_Boundary_Region";

const Foam::label Foam::meshWriters::STARCD::foamToStarFaceAddr[4][6] =
{
    { 4, 5, 2, 3, 0, 1 },     // 11 = pro-STAR hex
    { 0, 1, 4, 5, 2, -1 },    // 12 = pro-STAR prism
    { 5, 4, 2, 0, -1, -1 },   // 13 = pro-STAR tetra
    { 0, 4, 3, 5, 2, -1 }     // 14 = pro-STAR pyramid
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshWriters::STARCD::findDefaultBoundary() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label id = -1;

    // find Default_Boundary_Region if it exists
    forAll(patches, patchi)
    {
        if (defaultBoundaryName == patches[patchi].name())
        {
            id = patchi;
            break;
        }
    }
    return id;
}


void Foam::meshWriters::STARCD::getCellTable()
{
    // read constant/polyMesh/propertyName
    IOList<label> ioList
    (
        IOobject
        (
            "cellTableId",
            mesh_.time().constant(),
            polyMesh::meshSubDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    bool useCellZones = false;
    cellTableId_.setSize(mesh_.nCells(), -1);

    // get information from constant/polyMesh/cellTableId if possible
    if (ioList.headerOk())
    {
        if (ioList.size() == mesh_.nCells())
        {
            cellTableId_.transfer(ioList);

            if (cellTable_.empty())
            {
                Info<< "no cellTable information available" << endl;
            }
        }
        else
        {
            WarningInFunction
                << ioList.objectPath() << " has incorrect number of cells "
                << " - use cellZone information"
                << endl;

            ioList.clear();
            useCellZones = true;
        }
    }
    else
    {
        useCellZones = true;
    }


    if (useCellZones)
    {
        if (cellTable_.empty())
        {
            Info<< "created cellTable from cellZones" << endl;
            cellTable_ = mesh_;
        }

        // track if there are unzoned cells
        label nUnzoned = mesh_.nCells();

        // get the cellZone <-> cellTable correspondence
        Info<< "matching cellZones to cellTable" << endl;

        forAll(mesh_.cellZones(), zoneI)
        {
            const cellZone& cZone = mesh_.cellZones()[zoneI];
            if (cZone.size())
            {
                nUnzoned -= cZone.size();

                label tableId = cellTable_.findIndex(cZone.name());
                if (tableId < 0)
                {
                    dictionary dict;

                    dict.add("Label", cZone.name());
                    dict.add("MaterialType", "fluid");
                    tableId = cellTable_.append(dict);
                }

                forAll(cZone, i)
                {
                    cellTableId_[cZone[i]] = tableId;
                }
            }
        }

        if (nUnzoned)
        {
            dictionary dict;

            dict.add("Label", "__unZonedCells__");
            dict.add("MaterialType", "fluid");
            label tableId = cellTable_.append(dict);

            forAll(cellTableId_, i)
            {
                if (cellTableId_[i] < 0)
                {
                    cellTableId_[i] = tableId;
                }
            }
        }
    }
}


void Foam::meshWriters::STARCD::writeHeader(Ostream& os, const char* filetype)
{
    os  << "PROSTAR_" << filetype << nl
        << 4000
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << " " << 0
        << endl;
}


void Foam::meshWriters::STARCD::writePoints(const fileName& prefix) const
{
    OFstream os(prefix + ".vrt");
    writeHeader(os, "VERTEX");

    // Set the precision of the points data to 10
    os.precision(10);

    // force decimal point for Fortran input
    os.setf(std::ios::showpoint);

    const pointField& points = mesh_.points();

    Info<< "Writing " << os.name() << " : "
        << points.size() << " points" << endl;

    forAll(points, ptI)
    {
        // convert [m] -> [mm]
        os  << scaleFactor_ * points[ptI].x() << " "
            << scaleFactor_ * points[ptI].y() << " "
            << scaleFactor_ * points[ptI].z() << nl;
    }
    os.flush();

}


void Foam::meshWriters::STARCD::writeCells(const fileName& prefix) const
{
    OFstream os(prefix + ".cel");
    writeHeader(os, "CELL");

    // this is what we seem to need
    // map foam cellModeller index -> star shape
    Map<label> shapeLookupIndex;
    shapeLookupIndex.insert(hexModel->index(), 11);
    shapeLookupIndex.insert(prismModel->index(), 12);
    shapeLookupIndex.insert(tetModel->index(), 13);
    shapeLookupIndex.insert(pyrModel->index(), 14);

    const cellShapeList& shapes = mesh_.cellShapes();
    const cellList& cells  = mesh_.cells();
    const faceList& faces  = mesh_.faces();
    const labelList& owner = mesh_.faceOwner();

    Info<< "Writing " << os.name() << " : "
        << cells.size() << " cells" << endl;

    forAll(cells, cellId)
    {
        label tableId = cellTableId_[cellId];
        label materialType  = 1;        // 1(fluid)
        if (cellTable_.found(tableId)) //UN: Esto puede ser interesante para lo del tipo de celda
        {
            const dictionary& dict = cellTable_[tableId];
            if (dict.found("MaterialType"))
            {
                word matType;
                dict.lookup("MaterialType") >> matType;
                if (matType == "solid")
                {
                    materialType = 2;
                }

            }
        }

        const cellShape& shape = shapes[cellId];
        label mapIndex = shape.model().index();

        // a registered primitive type
        if (shapeLookupIndex.found(mapIndex))
        {
            label shapeId = shapeLookupIndex[mapIndex];//UN: esto puede servir para lo de imprimir pirámides y demás
            const labelList& vrtList = shapes[cellId];

/*Impresión de diferentes características de la celda
            os  << cellId + 1               // # de celda
                << " " << shapeId
                << " " << vrtList.size()    // 
                << " " << tableId           // 
                << " " << materialType;     // Tipo de material
*/
            // primitives have <= 8 vertices, but prevent overrun anyhow
            // indent following lines for ease of reading

         //UN: acá imprime las celdas
         //Oden del hexaedro:
             os 
                << " " << vrtList[1] +1
                << " " << vrtList[2] +1
                << " " << vrtList[3] +1
                << " " << vrtList[0] +1
                << " " << vrtList[5] +1
                << " " << vrtList[6] +1
                << " " << vrtList[7] +1
                << " " << vrtList[4] +1
                << nl;
             //UN:Pirámides y tetraedros quedan pendientes mientras se confirma que el método del hexaedro funciona

        }
        else // treat as general polyhedral
        {
          //UN: general polyhedral is not supported
        }
    }
}


void Foam::meshWriters::STARCD::writeBoundary(const fileName& prefix) const
{
    OFstream os(prefix + ".bnd");
    writeHeader(os, "BOUNDARY");

    const cellShapeList& shapes = mesh_.cellShapes();
    const cellList& cells  = mesh_.cells();
    const faceList& faces  = mesh_.faces();
    const labelList& owner = mesh_.faceOwner();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // this is what we seem to need
    // these MUST correspond to foamToStarFaceAddr
    //
    Map<label> faceLookupIndex;
    faceLookupIndex.insert(hexModel->index(), 0);
    faceLookupIndex.insert(prismModel->index(), 1);
    faceLookupIndex.insert(tetModel->index(), 2);
    faceLookupIndex.insert(pyrModel->index(), 3);

    Info<< "Writing " << os.name() << " : "
        << (mesh_.nFaces() - patches[0].start()) << " boundaries" << endl;


    label defaultId = findDefaultBoundary();

    forAll(cells, celli)
    {
      const labelList& cFaces  = cells[celli];
      os
        << celli << " " << 4.0 << " ";   //tipo de celda
      forAll(cFaces, cFacei)
      {
        label cellFaceId = findIndex(cFaces, cFacei);
        os
          << cellFaceId << " ";    //tipo de cara left
/*
          << " " <<    //tipo de cara front
          << " " <<    //tipo de cara bottom
          << " " <<    //tipo de cara right
          << " " <<    //tipo de cara derrire
          << " " <<;    //tipo de cara top
*/
      }
        os
          << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshWriters::STARCD::STARCD
(
    const polyMesh& mesh,
    const scalar scaleFactor
)
:
    meshWriter(mesh, scaleFactor)
{
    boundaryRegion_.readDict(mesh_);
    cellTable_.readDict(mesh_);
    getCellTable();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshWriters::STARCD::~STARCD()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshWriters::STARCD::rmFiles(const fileName& baseName) const
{
    rm(baseName + ".vrt");
    rm(baseName + ".cel");
    rm(baseName + ".bnd");
    rm(baseName + ".inp");
}


bool Foam::meshWriters::STARCD::write(const fileName& meshName) const
{
    fileName baseName(meshName);

    if (baseName.empty())
    {
        baseName = meshWriter::defaultMeshName;

        if
        (
            mesh_.time().timeName() != "0"
         && mesh_.time().timeName() != mesh_.time().constant()
        )
        {
            baseName += "_" + mesh_.time().timeName();
        }
    }

    rmFiles(baseName);
    writePoints(baseName);
    writeCells(baseName);

    if (writeBoundary_)
    {
        writeBoundary(baseName);
    }

    return true;
}


// ************************************************************************* //
