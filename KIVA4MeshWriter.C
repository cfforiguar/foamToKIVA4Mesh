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
    { 0, 2, 4, 1, 3, 5 },     // 11 = KIVA4 hex
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
    Info<< "Writing " << os.name() << " : "
        << (mesh_.nFaces()) << " caras" << endl;
    const wordList lpatches=patches.names();
    os << lpatches[0]<< " 0 " << lpatches[1]<< " 1 " << lpatches[2]<< " 2 "
       << lpatches[3]<< " 3 " << lpatches[4]<< " 4 " << lpatches[5]<< " 5 " << endl;

 cellList cellsNew=cells;

    int cellCt = 0;
    forAll(cells, celli)
    {
      const labelList& cFaces  = cells[celli];
      const cellShape& shape = shapes[celli];
      label mapIndex = shape.model().index();
      int cFacesiCt = 0;
      if (faceLookupIndex.found(mapIndex))
      {
        cFacesiCt = 0;
        forAll(cFaces, cFacesi)
        {
          label cellFaceId = findIndex(cFaces, cFacesi);
          const faceList sFaces = shape.faces();
          forAll(sFaces, sFacei)
          {
             if (faces[cFaces[cFacesi]] == sFaces[sFacei])
             {
//                os << faces[cFaces[cFacesi]] << "° faces" << endl;
                os << sFaces[sFacei] << "° sFaces" << endl;
                os << shape << "° SHAPE" << endl;
                os << cFacesi << "° cFacesi#" << endl;
                os << sFacei << "° sFaces#" << endl;
//                os << mapIndex << "°mapIndex" << endl;
//                cellFaceId=sFacei;
//                mapIndex = faceLookupIndex[mapIndex];

//                os << cFacesi << " °face" << endl;
//                os << sFacei << " °mface" << endl;
                cellFaceId = foamToStarFaceAddr[0][sFacei];
                //mtx[labels de caras de la celda][labels de caras del parche]//
//                os << cellFaceId << " cellFaceId" << endl;
//                os << (cFacesiCt==cFacesi) << "° id ok?" << endl;
                os << cells[celli][cFacesi] << "° OLD" << endl;
                os << cells[celli][cellFaceId] << "° NEW" << endl;
                cellsNew[celli][cellFaceId]=cells[celli][cFacesi];
             }
            
          }
          cFacesiCt++;
        }
      }
      cellCt++;
    }







/*

 forAll(faces, facesi)
 {
   label cellId = owner[facesi];
   const labelList& cFaces  = cells[cellId];
   const cellShape& shape = shapes[cellId];
 
   label cellFaceId = findIndex(cFaces, facesi);

   label mapIndex = shape.model().index();
   if (faceLookupIndex.found(mapIndex))
   {
     const faceList sFaces = shape.faces();
     forAll(sFaces, sFacei)
     {
//
//       cells[celli][cFacei]=
//
       if (faces[facesi] == sFaces[sFacei])
       {
          cellFaceId = sFacei;
          mapIndex = faceLookupIndex[mapIndex];
          cellFaceId = foamToStarFaceAddr[mapIndex][cellFaceId];
          //mtx[labels de caras de la celda][labels de caras del parche]//
          cellsNew[cellId][cellFaceId]=0;//patches.whichPatch(facesi);
         break;
       }
      }
    }
 }
*/
    forAll(cells, celli)
    {
      const labelList& cFaces  = cellsNew[celli];
//      const faceList& faces = cellFaces_[celli];
      
      /////***************
      //const faceList& cFacesB = cells[celli].faces();
      /////***************
      os
//        << celli << " " //label de la celda
        << 10 << " "   //tipo de celda
        << cFaces << " "    //labels de las caras de la celda
        << "";
      //forAll(cFaces, cFacei)
      //{
         os
           << patches.whichPatch(cFaces[0]) << " "    //tipo de cara left
           << patches.whichPatch(cFaces[1]) << " "    //tipo de cara front
           << patches.whichPatch(cFaces[2]) << " "    //tipo de cara bottom
           << patches.whichPatch(cFaces[3]) << " "    //tipo de cara right
           << patches.whichPatch(cFaces[4]) << " "    //tipo de cara derriere
           << patches.whichPatch(cFaces[5]) << " "    //tipo de cara top
           << "||"
           << owner[cFaces[0]] << " "    //tipo de cara left
           << owner[cFaces[1]] << " "    //tipo de cara front
           << owner[cFaces[2]] << " "    //tipo de cara bottom
           << owner[cFaces[3]] << " "    //tipo de cara right
           << owner[cFaces[4]] << " "    //tipo de cara derriere
           << owner[cFaces[5]] << " "    //tipo de cara top
           << "||"
/*           << cFaces[0] << " "    //tipo de cara left
           << cFaces[1] << " "    //tipo de cara front
           << cFaces[2] << " "    //tipo de cara bottom
           << cFaces[3] << " "    //tipo de cara right
           << cFaces[4] << " "    //tipo de cara derriere
           << cFaces[5] << " "    //tipo de cara derriere
           << "||"
           << cellsNew  << " "
           << "||"
           << cells  << " "
           
           << "||"
           << owner << " "    //
           << "|| cFaces "
           << cFaces << " "    //
           << "|| cells "
           << cells << " "    //
           << "|| cells[celli][0] "
           << cells[celli][0] << " "    //
           << "|| faces[cells[celli][5]] "
           << faces[cells[celli][5]] << " "    //


          //# celda, # bicho correspondiente a la cara

           
/*
           << patches.whichPatch(cFaces[0]) << " "    //tipo de cara left
           << patches.whichPatch(cFaces[2]) << " "    //tipo de cara front
           << patches.whichPatch(cFaces[4]) << " "    //tipo de cara bottom
           << patches.whichPatch(cFaces[1]) << " "    //tipo de cara right
           << patches.whichPatch(cFaces[3]) << " "    //tipo de cara derriere
           << patches.whichPatch(cFaces[5]) << " "    //tipo de cara top
*/
           << ""; 
       //}
        os
          << endl;
          
//        const faceListList& ebrio=mesh_.cellFaces();
        
        const cellShape& shape = shapes[celli];
        //      Info<< "cell " << cellId + 1 << " face " << facei
        //          << " == " << faces[facei]
        //          << " is index " << cellFaceId << " from " << cFaces;

        // Unfortunately, the order of faces returned by
        //   primitiveMesh::cells() is not necessarily the same
        //   as defined by primitiveMesh::cellShapes()
        // Thus, for registered primitive types, do the lookup ourselves.
        // Finally, the cellModel face number is re-mapped to the
        // STAR-CD local face number

        label mapIndex = shape.model().index();
        
        label cellFaceId;// = findIndex(cFaces, facei); //No sé qué hace eso de finIndex, pero sirve para crear el tipo de label que se necesita

/*

        // a registered primitive type
        if (faceLookupIndex.found(mapIndex))
        {
            const faceList sFaces = shape.faces();
            forAll(sFaces, sFacei) //caras del polihedro std
            {
              forAll(cFaces, cFacei)//caras de la celda
              {
//              if (Cara forma celda sFaces(sFacei) == Cara celda face (???)
                if (sFacei == cFaces[cFacei])
                {
                    cellFaceId = sFacei;          
                    mapIndex = faceLookupIndex[mapIndex];
                    cellFaceId = foamToStarFaceAddr[mapIndex][cellFaceId];
                    os
                      << patches.whichPatch(cFaces[cellFaceId]) << " "    //tipo de cara left
                      << "";
                }
              }
            }
         }


*/




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
