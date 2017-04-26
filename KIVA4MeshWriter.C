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

bool  Foam::operator==(const labelList& a, const labelList& b)
{
   List<bool> fnd(a.size(), false);
 
     forAll(b, bI)
     {
        label curLabel = b[bI];

         bool found = false;
 
         forAll(a, aI)
         {
             if (a[aI] == curLabel)
             {
                 found = true;
                 fnd[aI] = true;
                 break;
             }
         }
         if (!found)
         {
             return false;
         }
     }
 
     // check if all faces on a were marked
     bool result = true;
 
     forAll(fnd, aI)
     {
         result = (result && fnd[aI]);
     }
     return result;
}


void Foam::meshWriters::STARCD::writeHeader(Ostream& os, const int nPoints, const int nCells)
{
    os  << "foamToKIVA4Mesh: By Universidad Nacional de Colombia"  << endl;
    
    os  << nCells << " " << nPoints << endl;
}


void Foam::meshWriters::STARCD::writePoints(Ostream& os) const
{
    // Set the precision of the points data to 10
    os.precision(10);

    // force decimal point for Fortran input
    os.setf(std::ios::showpoint);
    
    //Foam points    
    const pointField& points = mesh_.points();
    const cellList& cells  = mesh_.cells();

    //Shaped points
    const cellShapeList& sCells = mesh_.cellShapes();

    const cellShape& sCell = sCells[0];//[cellId];//UN: OJO, esto se deb adaptar para otras formas elementales
    const pointField& sPoints = sCell.points(points);
    

    Info<< "Writing " << os.name() << " : "
        << points.size() << " points" << endl
        << sCell.model() << "modelo" << endl
        << "";

/* //LABEL: TESTS 
    os  << "Mapeados a shape"<< nl;
    forAll(sPoints, ptI)
    {
        // convert [m] -> [mm]
        os
            << ptI << " "
            << scaleFactor_ * sPoints[ptI].x() << " "
            << scaleFactor_ * sPoints[ptI].y() << " "
            << scaleFactor_ * sPoints[ptI].z() << nl;
    }
*/ //LABEL: TESTS
    os  << "Sin mapear"<< nl;

    forAll(points, ptI)
    {
        // convert [m] -> [mm]
        os
            << ptI << " "
            << scaleFactor_ * points[ptI].x() << " "
            << scaleFactor_ * points[ptI].y() << " "
            << scaleFactor_ * points[ptI].z() << nl;
    }

    os.flush();

}


void Foam::meshWriters::STARCD::writeCells(Ostream& os) const
{

    // this is what we seem to need
    // map foam cellModeller index -> star shape
    Map<label> shapeLookupIndex;
    shapeLookupIndex.insert(hexModel->index(), 11);
    shapeLookupIndex.insert(prismModel->index(), 12);
    shapeLookupIndex.insert(tetModel->index(), 13);
    shapeLookupIndex.insert(pyrModel->index(), 14);

    const cellShapeList& sCells = mesh_.cellShapes();
    const cellList& cells  = mesh_.cells();
    const faceList& faces  = mesh_.faces();
    const labelList& owner = mesh_.faceOwner();

    Info<< "Writing " << os.name() << " : "
        << cells.size() << " cells" << endl;

    cellList posFaces = cells;//Me guarda la lista de caras que usaré para top, bottom, etc
    forAll(cells, cellId)
    {
        const cellShape& sCell = sCells[cellId];
        label mapIndex = sCell.model().index();
/* //LABEL: TESTS
        os 
            << sCell  << " sCell"
            << nl;


        const labelList& sFaces = sCell.meshFaces(faces, cells[cellId]);
        os 
            << sFaces  << " sFaces"
            << nl;
            
*/ //LABEL: TESTS

/* //LABEL: TESTS
        const edgeList& edges = mesh_.edges();
        const labelList sEgdes=sCell.meshEdges(edges,cells[cellId]); 
        os 
            << sEgdes  << " sEgdes"
            << nl;

          // a registered primitive type
//        if (shapeLookupIndex.found(mapIndex))
//        {
//            label shapeId = shapeLookupIndex[mapIndex];//UN: esto puede servir para lo de imprimir pirámides y demás

*/ //LABEL: TESTS
    //Foam points    
    const pointField& points = mesh_.points();

    //Shaped points

    const pointField& sPoints = sCell.points(points);
    
            const labelList& sCellPoints = sCell;
            /* //LABEL: TESTS
            os 
                << " " << sPoints[sCellPoints[0]] << points[sCellPoints[0]]
                << " " << sCellPoints[1]
                << " " << sCellPoints[2]
                << " " << sCellPoints[3]
                << " " << sCellPoints[4]
                << " " << sCellPoints[5]
                << " " << sCellPoints[6]
                << " " << sCellPoints[7]
                << nl;
            */ //LABEL: TESTS    
             const labelList& cellPoints = mesh_.cellPoints(cellId);//Puntos de la celda local
            /* //LABEL: TESTS 
             os 
                << " " << sPoints[cellPoints[0]] << points[cellPoints[0]]
                << " " << cellPoints[1]
                << " " << cellPoints[2]
                << " " << cellPoints[3]
                << " " << cellPoints[4]
                << " " << cellPoints[5]
                << " " << cellPoints[6]
                << " " << cellPoints[7]
                << nl;
             */ //LABEL: TESTS  
                
    const faceList& rawModelFaces = sCell.model().modelFaces();
//        os << rawModelFaces  << " rawModelFaces" << nl; //Caras fictis que no se usan en el programa
        
    const faceList 	slistOfFaces = sCell.model().faces(sCellPoints); 
//        os << slistOfFaces  << " slistOfFaces" << nl; //Caras según los ptos sin mapear
//        os << sCellPoints  << " sCellPoints" << nl;
        
    const faceList 	listOfFaces = sCell.model().faces(cellPoints);
//        os << listOfFaces  << " listOfFaces" << nl; //Caras fictis que no se usan en el programa
//        os << cellPoints  << " cellPoints" << nl;
        
        
        
    const labelList 	LabelOfFaces = sCell.meshFaces(faces,cells[0]);
/* //LABEL: TESTS    
        os << LabelOfFaces  << " LabelOfFaces" << nl;
        os << cells  << " cells" << nl;
        os << LabelOfFaces[0]  << " LabelOfFaces[0]" << nl;
        os << faces[LabelOfFaces[0]]  << " faces[LabelOfFaces[0]] " << nl;
        
    os << sCell  << " sCell" << nl;
    os << cells  << " cells" << nl;
    os << sCell.faces()  << " sCell.faces()" << nl; //Caras según los ptos sin mapear
    os << faces  << " faces= cells_.faces()" << nl; //Caras en el orden de OFOAM
    os << sCell.meshFaces(faces,cells[0])  << " sCell.meshFaces(faces,cells[0])" << nl;
//    os << " " << "sPoints" << sPoints << nl;
//    os << " " << "points" << points << nl;

    os << " " << "sCell.faces()[0][0]" << sCell.faces()[0][0] << nl;
    os << " " << "sCell.faces()[0] " << sCell.faces()[0] << nl;
*/ //LABEL: TESTS
    const vectorField& facesCentres = mesh_.faceCentres(); //Return cell centres as volVectorField. More...   
    const vectorField& cellsCentres = mesh_.cellCentres(); //Return cell centres as volVectorField. More...       
    vectorField distFaceCentres = facesCentres;
    
    //UN: Para las demás caras, usar label sCellnFaces=sCell.nFaces() 
    //UN:    Esto funciona para cualquier cara pero lo hace leeeento. Es mejor separar por tipos de caras y luego hacer un loop para cada tipo
    float test[6] {};
//    forAll(cells, cellId)
    const labelList& cFaces  = cells[cellId];
    os << " " << "posFaces " << posFaces  << nl;
    forAll (cFaces,cFacesi)//Caras de la respectiva celda
    {
      distFaceCentres[cFacesi]=facesCentres[cFacesi]-cellsCentres[cellId];
      os << " " << "distFaceCentres [" << cFacesi << "] " << distFaceCentres[cFacesi] << nl;
      //OFOAM: Left, Right, Front, Back, Bottom, Top
      //KIVA: ???? Left, Right, Front, Back, Bottom, Top
      if (test[0]>distFaceCentres[cFacesi].x()) posFaces[cellId][0]=cFaces[cFacesi];//x min = posFaces[cellId][0]
      if (test[1]<distFaceCentres[cFacesi].x()) posFaces[cellId][1]=cFaces[cFacesi];//x max =       .
      if (test[2]>distFaceCentres[cFacesi].y()) posFaces[cellId][2]=cFaces[cFacesi];//y min =       .
      if (test[3]<distFaceCentres[cFacesi].y()) posFaces[cellId][3]=cFaces[cFacesi];//y max =       .
      if (test[4]>distFaceCentres[cFacesi].z()) posFaces[cellId][4]=cFaces[cFacesi];//z min =       .
      if (test[5]<distFaceCentres[cFacesi].z()) posFaces[cellId][5]=cFaces[cFacesi];//z max =       .
      //UN: ESTE ALGORITMO FALLARÍA SI ALGÚN ELEMENTO ESTÁ A EXACTAMENTE 45° PUES AL MENOS 2 CARAS TENDRÍAN VALORES IDÉNTICOS EN CADA COORDENADA
      os << " " << "cFaces[" << cFacesi << "] " << cFaces[cFacesi] << nl;
      label cellFaceId = findIndex(cFaces, cFacesi);//UN: ¿¿¿¿pasa el índice local de cara cFacesi a índice global de caras mesh.faces()????
      os << "cellFaceId " << cellFaceId << nl;
    }
    
    //Tomar los ptos de cada celda. ¿cell[0....7]?
 
    //Usar dichos ptos para acceder a los fPoints respectivos
    //Hallar la posición de dichos pts en la celda
    
      
    int lCaras [8][3] = {//véase 5.1.3: Cell Shapes:
                         //https://cfd.direct/openfoam/user-guide/mesh-description/
                        {4, 0, 2},//Caras que intersecan el pto 0 según manual oFoam 
                        {4, 1, 2},//Caras que intersecan el pto 1 según manual oFoam 
                        {4, 1, 3},//Caras que intersecan el pto 2 según manual oFoam 
                        {4, 0, 3},//Caras que intersecan el pto 3 según manual oFoam 
                        {5, 0, 2},//Caras que intersecan el pto 4 según manual oFoam 
                        {5, 1, 2},//Caras que intersecan el pto 5 según manual oFoam 
                        {5, 1, 3},//Caras que intersecan el pto 6 según manual oFoam 
                        {5, 0, 3},//Caras que intersecan el pto 7 según manual oFoam 
                        };
    
    const labelListList& fPoints=mesh_.pointFaces();//Quiero que sean los fPoints de cada celda, en vez de todos
      os << "fPoints " << fPoints << nl;
    
    const  labelListList& cPoints = mesh_.cellPoints();//Labels de los puntos de las celdas
    labelListList scPoints = cPoints;//Puntos de la celda modelo (shapeada)
      os << "cPoints " << cPoints << nl;
      
    forAll(cPoints, celli)      //Mirar c/u de las celdas
    {
      //Ubica las caras de cada cela en el orden dado en las intersecciones de lCaras
      //De este modo que se pueden ubicar los puntos en posFaces
      List <labelList> pointFacesShape(0);
      for (int i=0;i<=7;i++)//Caras de la respectiva celda
      {
        List <label> tempLabelList(0);//No idea on how to clean a list
        for (int j=0;j<=2;j++)//Caras de la respectiva celda
        {
          label tempLabel= posFaces[celli][lCaras[i][j]];
          tempLabelList.append(tempLabel);
        }
        pointFacesShape.append(tempLabelList);
      }
      
      //Encuentra la ubicación de los puntos al comparar la tabla de caras intersectadas
      // con las caras reales del volumen construido por openFOam
      const labelList& cIPoints = cPoints[celli];
      os << "cIPoints " << cIPoints << nl;
      forAll (cIPoints,cIPointsi)  //Puntos de la celda actual
      {
//        os << "cIPoints[" << cIPointsi << "]= " << cIPoints[cIPointsi] << nl;
        label curPoint = cIPoints[cIPointsi]; //Current point's label
//        os << "curPoint " << lPoint << nl;
        for (int j=0;j<=7;j++)     //Comparar con puntos de la celda modelo
        {
            if (fPoints[curPoint]==pointFacesShape[j]) scPoints[celli][j]=curPoint;
        }
      }
    }    
    os << "scPoints " << scPoints << nl;

    //UN:Pirámides y tetraedros quedan pendientes mientras se confirma que el método del hexaedro funciona

//        }
/*        else // treat as general polyhedral
        {
          //UN: general polyhedral is not supported
        }
        */
    }
}


void Foam::meshWriters::STARCD::writeBoundary(Ostream& os) const
{

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

/*//Recupera y escribe los nombres de los parches 
    const wordList lpatches=patches.names();
    os << lpatches[0]<< " 0 " << lpatches[1]<< " 1 " << lpatches[2]<< " 2 "
       << lpatches[3]<< " 3 " << lpatches[4]<< " 4 " << lpatches[5]<< " 5 " << endl;
*/

 cellList cellsNew=cells;


    //      Info<< "cell " << cellId + 1 << " face " << facei
    //          << " == " << faces[facei]
    //          << " is index " << cellFaceId << " from " << cFaces;

    // Unfortunately, the order of faces returned by
    //   primitiveMesh::cells() is not necessarily the same
    //   as defined by primitiveMesh::cellShapes()
    // Thus, for registered primitive types, do the lookup ourselves.
    // Finally, the cellModel face number is re-mapped to the
    // STAR-CD local face number

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
/*
                os << sFaces[sFacei] << "° sFaces" << endl;
                os << shape << "° SHAPE" << endl;
                os << cFacesi << "° cFacesi#" << endl;
                os << sFacei << "° sFaces#" << endl;
*/
//                os << mapIndex << "°mapIndex" << endl;
//                cellFaceId=sFacei;
//                mapIndex = faceLookupIndex[mapIndex];

//                os << cFacesi << " °face" << endl;
//                os << sFacei << " °mface" << endl;
                cellFaceId = foamToStarFaceAddr[0][sFacei];
                //mtx[labels de caras de la celda][labels de caras del parche]//
//                os << cellFaceId << " cellFaceId" << endl;
//                os << (cFacesiCt==cFacesi) << "° id ok?" << endl;
//                os << cells[celli][cFacesi] << "° OLD" << endl;
//                os << cells[celli][cellFaceId] << "° NEW" << endl;
                cellsNew[celli][cellFaceId]=cells[celli][cFacesi];
             }
            
          }
          cFacesiCt++;
        }
      }
      cellCt++;
    }


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
//        << cFaces << " "    //labels de las caras de la celda
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
/*
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
/*       
        os
          << endl;
*/
          
//        const faceListList& ebrio=mesh_.cellFaces();
        
        const cellShape& shape = shapes[celli];

        label mapIndex = shape.model().index();
        
        label cellFaceId;// = findIndex(cFaces, facei); //No sé qué hace eso de finIndex, pero sirve para crear el tipo de label que se necesita






         os
           << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshWriters::STARCD::STARCD
(
    const fvMesh& mesh,
    const scalar scaleFactor
)
:
    meshWriter(mesh, scaleFactor)
{
    boundaryRegion_.readDict(mesh_);
    cellTable_.readDict(mesh_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshWriters::STARCD::~STARCD()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshWriters::STARCD::rmFiles(const fileName& baseName) const
{
    rm(baseName);
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
    OFstream os(baseName);
    const pointField& points = mesh_.points();
    const cellList& cells  = mesh_.cells();
    writeHeader(os, points.size(), cells.size());
    writePoints(os);
    writeCells(os);
    
    
    if (writeBoundary_)
    {
        writeBoundary(os);
    }
    //UNPending "structured mesh" part
    os 
      << "0" << endl
      << "0" << endl;
    return true;
}


// ************************************************************************* //
