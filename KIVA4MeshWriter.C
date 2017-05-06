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


Foam::cellList Foam::meshWriters::STARCD::rePosFaces(Ostream& os) const //Me guarda la lista de caras que usaré para top, bottom, etc
{
    const cellList& cells  = mesh_.cells();
    cellList posFaces = cells;//Me guarda la lista de caras que usaré para top, bottom, etc
    const vectorField& facesCentres = mesh_.faceCentres(); //Return cell centres as volVectorField. More...   
    const vectorField& cellsCentres = mesh_.cellCentres(); //Return cell centres as volVectorField. More...       
    vectorField distFaceCentres = facesCentres;
    //<DEBUG>
    const faceList& faces  = mesh_.faces();
    const vectorField& 	vFaces = mesh_.faceAreas();
    const labelList& owner = mesh_.faceOwner();
    os << "vFaces "<< vFaces << nl;
    os << "Owner "<< owner << nl;
    //<\DEBUG>

    forAll (cells, cellId)
    {
  //    os << "facesCentres"<< facesCentres << nl;
  //        os << "cellsCentres"<< cellsCentres << nl;            
      //UN: Para las demás caras, usar label sCellnFaces=sCell.nFaces() 
      //UN:    Esto funciona para cualquier cara pero lo hace leeeento. Es mejor separar por tipos de caras y luego hacer un loop para cada tipo
      double  test[6] {-1.0e-6,-1.0e-6,-1.0e-6,1.0e-6,1.0e-6,1.0e-6};
             
  //    forAll(cells, cellId)
      const labelList& cFaces  = cells[cellId];
  //    os << " " << "posFaces " << posFaces  << nl;
    
      forAll (cFaces,cFacesi)//Caras de la respectiva celda
      {
        label curFace = cFaces[cFacesi]; //Current faces's labels     
        distFaceCentres[curFace]=facesCentres[curFace]-cellsCentres[cellId];
        vector vFace=vFaces[curFace];
        vFace /= mag(vFace);
        os << "cellId "<< cellId << nl;
        os << "  vFace "<< vFace << nl;
        if (owner[curFace]!=cellId){
          os << "     Owned "<< cellId << nl;
          vFace=vFace*-1.0;
          os << "     vFace New"<< vFace << nl;
        }
        //<DEBUG>
//        os << "curFacePoints "<< faces[curFace] << nl;
        //<\DEBUG>
//        os << "distFaceCentres[curFace]"<< distFaceCentres[curFace] << nl;
  /* //LABEL     //OFOAM: Left, Right, Front, Back, Bottom, Top
        if (test[0]>vFace.x()) posFaces[cellId][0]=curFace;//x min = posFaces[cellId][0]
        if (test[1]<vFace.x()) posFaces[cellId][1]=curFace;//x max =       .
        if (test[2]>vFace.y()) posFaces[cellId][2]=curFace;//y min =       .
        if (test[3]<vFace.y()) posFaces[cellId][3]=curFace;//y max =       .
        if (test[4]>vFace.z()) posFaces[cellId][4]=curFace;//z min =       .
        if (test[5]<vFace.z()) posFaces[cellId][5]=curFace;//z max =       .
  */ //LABEL
        
   //LABEL     //KIVA: left 0, front 1, bottom 2, right 3, derrire 4, top 5
        if (test[0]>vFace.x()){
          test[0]=vFace.x();
          posFaces[cellId][0]=curFace;//x min = posFaces[cellId][0]
          os << "-x= " << vFace.x() << nl;
        }
        if (test[1]>vFace.y()){
          test[1]=vFace.y();
          posFaces[cellId][1]=curFace;//y min =       .
          os << "-y= " <<  vFace.y() << nl;
        }
        if (test[2]>vFace.z()){
          test[2]=vFace.z();
          posFaces[cellId][2]=curFace;//z min =       .
          os << "-z= " <<  vFace.z() << nl;
        }
        if (test[3]<vFace.x()){
          test[3]=vFace.x();
          posFaces[cellId][3]=curFace;//x max =       .
          os << "x= " <<  vFace.x() << nl;
        }
        if (test[4]<vFace.y()){
          test[4]=vFace.y();
          posFaces[cellId][4]=curFace;//y max =       .
          os << "y= " << vFace.y() << nl;
        }
        if (test[5]<vFace.z()){
          test[5]=vFace.z();
          posFaces[cellId][5]=curFace;//z max =       .
          os << "z= " << vFace.z() << nl;
        }
   //LABEL
        
        //UN: ESTE ALGORITMO FALLARÍA SI ALGÚN ELEMENTO ESTÁ A EXACTAMENTE 45° PUES AL MENOS 2 CARAS TENDRÍAN VALORES IDÉNTICOS EN CADA COORDENADA
        label cellFaceId = findIndex(cFaces, cFacesi);//UN: ¿¿¿¿pasa el índice local de cara cFacesi a índice global de caras mesh.faces()????
      }
    }
    os << posFaces << endl;
    return posFaces;
}

bool  Foam::meshWriters::STARCD::isContained (const labelList& a, const labelList& b) const
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
     /*
         if (!found)
         {
             return false;
         }
     */
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
    os  << "Sin mapear"<< nl;
*/ //LABEL: TESTS

    forAll(points, ptI)
    {
        // convert [m] -> [mm]
        os
 //           << ptI << " "
            << scaleFactor_ * points[ptI].x() << " "
            << scaleFactor_ * points[ptI].y() << " "
            << scaleFactor_ * points[ptI].z() << nl;
    }

    os.flush();

}


void Foam::meshWriters::STARCD::writeCells(Ostream& os,const cellList& posFaces) const
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
        
    forAll (cells, cellId)
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

    //Tomar los ptos de cada celda. ¿cell[0....7]?
 
    //Usar dichos ptos para acceder a los fPoints respectivos
    //Hallar la posición de dichos pts en la celda
    }
    int lCaras [8][3] = 
     //-----Según KIVA4
                        {//véase 5.1.3: Cell Shapes:
                         //https://cfd.direct/openfoam/user-guide/mesh-description/
                        {2, 1, 3},//Caras que intersecan el pto 0 según manual 
                        {2, 3, 4},//Caras que intersecan el pto 1 según manual 
                        {2, 4, 0},//Caras que intersecan el pto 2 según manual
                        {2, 0, 1},//Caras que intersecan el pto 3 según manual
                        {5, 1, 3},//Caras que intersecan el pto 4 según manual
                        {5, 3, 4},//Caras que intersecan el pto 5 según manual
                        {5, 4, 0},//Caras que intersecan el pto 6 según manual
                        {5, 0, 1},//Caras que intersecan el pto 7 según manual
                        };
     //-----Según KIVA4

    /* //-----Según oFoam
                        {//véase 5.1.3: Cell Shapes:
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
    */ //-----Según oFoam

    const labelListList& fPoints=mesh_.pointFaces();//Quiero que sean los fPoints de cada celda, en vez de todos
    
    const  labelListList& cPoints = mesh_.cellPoints();//Labels de los puntos de las celdas
    labelListList scPoints = cPoints;//Puntos de la celda modelo (shapeada)
    
    forAll(scPoints, celli)      //Mirar c/u de las celdas
    {
      labelList& scPointsCell=scPoints[celli];
      forAll(scPointsCell,scPointsCelli) 
        scPointsCell[scPointsCelli]=-1;
    }
    
    
    
    forAll(cPoints, celli)      //Mirar c/u de las celdas
    {
      //Ubica las caras de cada celda en el orden dado en las intersecciones de lCaras
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
//      os << "pointFacesShape"<< pointFacesShape << nl;
      //Encuentra la ubicación de los puntos al comparar la tabla de caras intersectadas
      // con las caras reales del volumen construido por openFOam
      const labelList& cIPoints = cPoints[celli];
      forAll (cIPoints,cIPointsi)  //Puntos de la celda actual
      {
//        os << "cIPoints[" << cIPointsi << "]= " << cIPoints[cIPointsi] << nl;
        label curPoint = cIPoints[cIPointsi]; //Current point's label
//        os << "curPoint " << lPoint << nl;
        for (int j=0;j<=7;j++)     //Comparar con puntos de la celda modelo
        {
//            os << "fPoints[" <<curPoint << "]=" << fPoints[curPoint]<< nl;
            if (isContained(pointFacesShape[j],fPoints[curPoint]))
                scPoints[celli][j]=curPoint+1;
        }
      }
    }
    
    forAll(cells, celli)
    {
         os
           << scPoints[celli][0] << " "
           << scPoints[celli][1] << " "
           << scPoints[celli][2] << " "
           << scPoints[celli][3] << " "
           << scPoints[celli][4] << " "
           << scPoints[celli][5] << " "
           << scPoints[celli][6] << " "
           << scPoints[celli][7] << " "
           << endl;
    }

    //UN:Pirámides y tetraedros quedan pendientes mientras se confirma que el método del hexaedro funciona

//        }
/*        else // treat as general polyhedral
        {
          //UN: general polyhedral is not supported
        }
        */
    
}


void Foam::meshWriters::STARCD::writeBoundary(Ostream& os, const cellList& posFaces) const
{
    const cellList& cells  = mesh_.cells();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const wordList patchNames=patches.names();
    labelList swPatches(0);
    wordList patchChange;
    labelList  patchCodes (19,0);
    cellList newPatches=posFaces;   
    
    Info<< "Writing " << os.name() << " : "
        << (mesh_.nFaces() - patches[0].start()) << " boundaries" << endl;
    Info<< "Writing " << os.name() << " : "
        << (mesh_.nFaces()) << " caras" << endl;

    patchCodes[ 0]=10  ; patchChange.append(word("moving"));
    patchCodes[ 1]=11  ; patchChange.append(word("v1bottom"));
    patchCodes[ 2]=12  ; patchChange.append(word("v1top"));
    patchCodes[ 3]=13  ; patchChange.append(word("v2bottom"));
    patchCodes[ 4]=14  ; patchChange.append(word("v2top"));
    patchCodes[ 5]=15  ; patchChange.append(word("v3bottom"));
    patchCodes[ 6]=16  ; patchChange.append(word("v3top"));
    patchCodes[ 7]=17  ; patchChange.append(word("v4bottom"));
    patchCodes[ 8]=18  ; patchChange.append(word("v4top"));
    patchCodes[ 9]=20  ; patchChange.append(word("solid"));
    patchCodes[10]=21  ; patchChange.append(word("solidh"));
    patchCodes[11]=30  ; patchChange.append(word("axis"));
    patchCodes[12]=40  ; patchChange.append(word("fluid"));
    patchCodes[13]=50  ; patchChange.append(word("periodf"));
    patchCodes[14]=60  ; patchChange.append(word("periodd"));
    patchCodes[15]=70  ; patchChange.append(word("inflow"));
    patchCodes[16]=80  ; patchChange.append(word("outflow"));
    patchCodes[17]=90  ; patchChange.append(word("presin"));
    patchCodes[18]=100 ; patchChange.append(word("presout"));


    forAll(patchNames, patchNamesi)//Seleciono parches válidos para ahorrarme unos ciclos
    {
      forAll(patchChange, patchChangei)
      {
        if (patchNames[patchNamesi]==patchChange[patchChangei]){
            swPatches.append(patchCodes[patchChangei]);
/* //Debug
             os
               << "swPatches[" << patchNamesi << "]= " << swPatches[patchNamesi] << endl
               << "patchCodes[" << patchChangei << "]= " << patchCodes[patchChangei] << endl
               << "";
*/ //Debug
        }
      }
    }



    forAll(cells, celli)//Reasigno los números de parche
    {
      const cell curCell=posFaces[celli];
      
      forAll(curCell, curCelli)
      {
        const label currPatch = patches.whichPatch(curCell[curCelli]);
        
        if (currPatch==-1){
          newPatches[celli][curCelli]=patchCodes[12];/////caso donde el parche = -1, entonces es fluido
          continue;
        }
        if (currPatch>-1)
              newPatches[celli][curCelli]=swPatches[currPatch]; //Acá hago la asignación de los demás parches
      }
    }
    
   
    forAll(cells, celli)
    {

/*       
       os
         << 10 << "patches "   //tipo de celda
         << patches.whichPatch(posFaces[celli][0]) << " "
         << patches.whichPatch(posFaces[celli][1]) << " "
         << patches.whichPatch(posFaces[celli][2]) << " "
         << patches.whichPatch(posFaces[celli][3]) << " "
         << patches.whichPatch(posFaces[celli][4]) << " "
         << patches.whichPatch(posFaces[celli][5]) << " "
         << endl;
*/         
       os
         << 10 << " "   //tipo de celda
         << newPatches[celli][0] << " "
         << newPatches[celli][1] << " "
         << newPatches[celli][2] << " "
         << newPatches[celli][3] << " "
         << newPatches[celli][4] << " "
         << newPatches[celli][5] << " "
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
    cellList posFaces=rePosFaces(os);
    writeCells(os,posFaces);
    
    
    if (writeBoundary_)
    {
        writeBoundary(os,posFaces);
    }
    //UNPending "structured mesh" part
    os 
      << "0" << endl
      << "0" << endl;
    return true;
}


// ************************************************************************* //
