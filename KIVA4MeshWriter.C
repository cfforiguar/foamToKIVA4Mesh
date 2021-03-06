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


Foam::cellList Foam::meshWriters::STARCD::rePosFaces(Ostream& os) const //Me guarda la lista de caras que usaré para top, bottom, etc
{
    const cellList& cells  = mesh_.cells();
    cellList posFaces = cells;//Me guarda la lista de caras que usaré para top, bottom, etc
    const vectorField& facesCentres = mesh_.faceCentres(); //Return cell centres as volVectorField. More...   
    const vectorField& cellsCentres = mesh_.cellCentres(); //Return cell centres as volVectorField. More...       
    //<DEBUG>
    const faceList& faces  = mesh_.faces();
    const vectorField& 	vFaces = mesh_.faceAreas();
    const labelList& owner = mesh_.faceOwner();
//    os << "vFaces "<< vFaces << nl;
//    os << "Owner "<< owner << nl;
    //<\DEBUG>

    forAll (cells, cellId)
    {
      const labelList& cFaces  = cells[cellId];
      forAll (cFaces,cFacesi)//Caras de la respectiva celda
      {        
        posFaces[cellId][cFacesi]=-1;
      }
    }
//        os << "posFaces: Inicial "<< posFaces << nl;

   //UN: Creo variables que clonaré más adelante para guardar los datos de las caras de la celda
    List <labelList> dummyLabelListList(0);
    labelList dummyLabelList(2);dummyLabelList[0]=-1;dummyLabelList[1]=-1;
    List <vector> dummyVectorList(0);
    vector dummyVector=facesCentres[0]; dummyVector[0]=-9;dummyVector[1]=-9;dummyVector[2]=-9;
    for (int i=0;i<6;i++)//Caras de la respectiva celda
    {
      dummyLabelListList.append(dummyLabelList);
      dummyVectorList.append(dummyVector);
    }

    forAll (cells, cellId)
    {
  //    os << "facesCentres"<< facesCentres << nl;
  //        os << "cellsCentres"<< cellsCentres << nl;            
      //UN: Para las demás caras, usar label sCellnFaces=sCell.nFaces() 
      //UN:    Esto funciona para cualquier cara pero lo hace leeeento. Es mejor separar por tipos de caras y luego hacer un loop para cada tipo
      double  test[6] {0.0,0.0,0.0,0.0,0.0,0.0};
      const labelList& cFaces  = cells[cellId];
  //    os << " " << "posFaces " << posFaces  << nl;
      List <vector> cuFaceVector=dummyVectorList;
      List <labelList> cOppoFaces=dummyLabelListList;
      for (label iFaces=0;iFaces<=5;iFaces++){//Caras de la respectiva celda
        label curFace = cFaces[iFaces]; //Current faces's labels     
        label	curOppFace=cells[cellId].opposingFaceLabel(curFace, faces);//Opposite of current faces's labels  
        cOppoFaces[iFaces][0]=curFace;
        cOppoFaces[iFaces][1]=curOppFace;
        vector uOppFaces =facesCentres[curFace]-facesCentres[curOppFace]; uOppFaces /= mag(uOppFaces);
        cuFaceVector[iFaces]=uOppFaces;
      }      

      double max=0.0;      
      int posMax[6][3]=  //[Place][Axis]
          {{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}};
      //Busque las componentes mayores de dirección de caras
      for (int Axis=0;Axis<=2;Axis++)
      {   //Encuentre los valores máximos/eje
      
        bool NotRanked[]={1,1,1,1,1,1}; //Lleva la cuenta de los vectores ya ranqueados
        for (int Place=0;Place<6;Place++)
        {  //Saque los 3 primeros ejes por magnitud
          int maxLocation=0;
          bool flag=0;
          for (int i=0;i<6;i++)
          {
            if (NotRanked[i])
            {  //Si aún no ha sido clasificado el valor, úselo como máximo
              max=cuFaceVector[i][Axis];
              maxLocation=i;
              break;
            }
          }
          for (int currVectori=0;currVectori<=5;currVectori++)
          {//Hace la comparación
            if (NotRanked[currVectori])
            {  //Si aún no ha sido clasificado el valor, úselo
              double currValue=cuFaceVector[currVectori][Axis];
              if (max <= currValue)//UN: Si no funciona, pedir que sea un valor > 0
              {  //Si supera el valor máximo hasta la fecha, guárdelo
                maxLocation=currVectori;
                max=currValue;
                flag=1;
              }
            }
          }
          if (flag){
            posMax[Place][Axis]=maxLocation;
            NotRanked[maxLocation]=0;
            flag=0;
          }
        }
/*        
        os << "NotRanked " << nl; 
        for(int q=0;q<6;q++){
          os << " "<< NotRanked[q]; 

        }
*/

      }


/*      
      os << "cuFaceVector "<< cuFaceVector << nl; 
      os << " " << "posMax "  << nl;
      for (int i =0;i<6;i++){
        for (int j =0;j<3;j++){
          os << " " << posMax[i][j];
        }
        os << nl;
      }
*/      
      
//      os << " cellId " << cellId << nl;
      
      //Mirar si hay ejes repetidos
      //Y vetar los ejes no repetdos
      bool Repetida[3]={1,1,1};
      int escogidos[3]={-1,-1,-1};
      bool cFallida=0;
      for (int Axis=0;Axis<3;Axis++)
      { 
        bool flag=0;       
        if(posMax[0][Axis]==-1) {cFallida=1;os << "ERROR: Vectores unitarios no determinados "  << nl;};
        for (int Axis2=0;Axis2<3;Axis2++)//Caras de la respectiva celda
        {        
          if(isContained(cOppoFaces[posMax[0][Axis]],cOppoFaces[posMax[0][Axis2]])
              && Axis!=Axis2) 
          {
            flag=1;
          }
        }
        Repetida[Axis]=flag;
        if (!flag){
          escogidos[Axis]=posMax[0][Axis];
        }
      }
/*      
      os << " Repetida";
      for (int j =0;j<3;j++){
        os << " " << Repetida[j] ;
      } os << nl;
      os << " escogidos";
      for (int j =0;j<3;j++){
        os << " " << escogidos[j] ;
      } os << nl;      
*/      
      
      //Ranquear los ejes
      //Asignarlos si no hay repetidos
      //Si hay repetidos
      //-Tomar el eje donde el repetido es mayor
      //  -Vetar el repetido mayor
      //-Tomar otro eje y asignarle el siguiente en el ranking que no este ya escogido
      
      
      //Ranquear los ejes repetidos
      bool  NotRanked[]={1,1,1,1,1,1}; //Lleva la cuenta de los vectores ya ranqueados
      for (int i=0;i<3;i++)  {NotRanked[i]=Repetida[i];}
      
      int axisOrder[3]={-1,-1,-1};
      for (int Place=0;Place<3;Place++)
      {  //Saque los 3 primeros ejes por magnitud
        int maxLocation=0;
        bool flag=0;
        for (int Axis=0;Axis<3;Axis++)
        {
          if (NotRanked[Axis])
          {  //Escoja otro valor no ranqueado para clasificarlo
            max=cuFaceVector[posMax[0][Axis]][Axis];
            maxLocation=Axis;
            break;
          }
        }
        for (int Axis=0;Axis<=2;Axis++)
        {//Hace la comparación
          if (NotRanked[Axis])
          {  //Si aún no ha sido clasificado el valor, úselo
            double currValue=cuFaceVector[posMax[0][Axis]][Axis];
            if (max<=currValue)
            {  //Si supera el valor máximo hasta la fecha, guárdelo
              maxLocation=Axis;
              max=currValue;
              flag=1;
            }
          }
        }
        if (flag){
          axisOrder[Place]=maxLocation;
          NotRanked[maxLocation]=0;
          flag=0;
        }
      }

/*      
      os << " axisOrder";
      for (int j =0;j<3;j++){
        os << " " << axisOrder[j] ;
      } os << nl;
*/
      
      for (int Place=0;Place<3;Place++)
      {
        if (axisOrder[Place]>-1){
          for (int cont=0;cont<6;cont++)
          {
            int Candidate=posMax[cont][axisOrder[Place]];
            bool axisFlag[3]={1,1,1};
            for(int i=0;i<3;i++){
              if (escogidos[i]>-1){
                axisFlag[i]=!isContained(cOppoFaces[Candidate],cOppoFaces[escogidos[i]]);
              }
            }
            if (axisFlag[0] && axisFlag[1] && axisFlag[2])
            {
              escogidos[axisOrder[Place]]=Candidate;
              break;
            }
          }
        }        
      }

/*
      os << " escogidos";
      for (int j =0;j<3;j++){
        os << " " << escogidos[j] ;
      } os << nl;

      os << cOppoFaces << nl;
*/      
      //LABEL     //KIVA: left 0, front 1, bottom 2, right 3, derrire 4, top 5
      posFaces[cellId][0]=cOppoFaces[escogidos[0]] [1];//x min = posFaces[cellId][0]
      posFaces[cellId][3]=cOppoFaces[escogidos[0]] [0];//x max =       .
      
      posFaces[cellId][1]=cOppoFaces[escogidos[1]] [1];//y min =       .
      posFaces[cellId][4]=cOppoFaces[escogidos[1]] [0];//y max =       .
      
      posFaces[cellId][2]=cOppoFaces[escogidos[2]] [1];//z min =       .
      posFaces[cellId][5]=cOppoFaces[escogidos[2]] [0];//z max =       .

      
    }
/*
    os << "posFaces: Final " << posFaces << endl;
*/


    forAll (cells, cellId)
    {
      const labelList& cFaces  = cells[cellId];
      bool cFallida=0;
      forAll (cFaces,cFacesi)//Caras de la respectiva celda
      {        
        if(posFaces[cellId][cFacesi]==-1){cFallida=1;os << "ERROR: CELDA INDETERMINADA: "  << nl;};
        forAll (cFaces,cFacesii)//Caras de la respectiva celda
        {        
          if(posFaces[cellId][cFacesi]==posFaces[cellId][cFacesii]
              && cFacesi!=cFacesii)
               {cFallida=1;os << "ERROR: CARAS REPETIDAS: "  << nl;};
        }
      }
        if(cFallida==1){ 
          os << "cellId " << cellId << nl;
          os << "posFaces " << posFaces[cellId] << nl;
          os << "nodeLabels " << cells[cellId] << nl;
                 /*          forAll (cFaces,cFacesi)//Caras de la respectiva celda
          {        
            if(posFaces[cellId][cFacesi]==-1) bool cFallida=1;
          }      */
        }
    }
//    os << "posFaces " << posFaces << nl;
    return posFaces;
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
    
    
    forAll(patchNames, patchNamesi)//Seleciono parches presentes para ahorrarme unos ciclos
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

    //Reasigno los números de parche
    forAll(cells, celli)
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
    
    labelList swZones(0);
    const cellZoneMesh&  cellZones=mesh_.cellZones();
    const wordList zoneNames=cellZones.names();
    labelList  cellZCodes (8,0);
    wordList cellZChange;
    labelList NewCZones(cells.size(),0);
    
//    os << "cellZones " << cellZones << endl;
//    os << "cellZones[0][0]" << cellZones[0][0] << endl;//cellLabels      List<label> 
//    os << "test " << cellZones.whichZone(patchCodes[ 0]) <<endl;    
//    os << "zoneNames " << zoneNames << endl;
    
    cellZCodes[ 0]=0   ; cellZChange.append(word("inactive"));
    cellZCodes[ 1]=10  ; cellZChange.append(word("squish"));
    cellZCodes[ 2]=11  ; cellZChange.append(word("bowl"));
    cellZCodes[ 3]=14  ; cellZChange.append(word("dome"));
    cellZCodes[ 4]=20  ; cellZChange.append(word("port1"));
    cellZCodes[ 5]=30  ; cellZChange.append(word("port2"));
    cellZCodes[ 6]=40  ; cellZChange.append(word("port3"));
    cellZCodes[ 7]=50  ; cellZChange.append(word("port4"));



    forAll(zoneNames, zoneNamesi)//Seleciono sólo zonas presentes para ahorrarme unos ciclos
    {
      forAll(cellZChange, cellZChangei)
      {
        if (zoneNames[zoneNamesi]==cellZChange[cellZChangei]){
            swZones.append(cellZCodes[cellZChangei]);
/* //Debug
             os
               << "swZones[" << zoneNamesi << "]= " << swZones[zoneNamesi] << endl
               << "cellZCodes[" << cellZChangei << "]= " << cellZCodes[cellZChangei] << endl
               << "";
*/  //Debug
        }
      }
    }
    

    forAll(cells, celli)//Reasigno los números zona de OpenFoam a Kiva4
    {
      NewCZones[celli]=swZones[cellZones.whichZone(celli)];
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
         << NewCZones [celli]    << " "   //tipo de celda
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
