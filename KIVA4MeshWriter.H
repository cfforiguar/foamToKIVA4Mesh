/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::meshWriters::STARCD

Description
    Writes polyMesh in pro-STAR (v4) bnd/cel/vrt format

    The cellTableId and cellTable information are used (if available).
    Otherwise the cellZones are used (if available).

SourceFiles
    STARCDMeshWriter.C

\*---------------------------------------------------------------------------*/

#ifndef KIVA4MeshWriter.H
#define KIVA4MeshWriter.H

#include "meshWriter.H"
#include "IOstream.H"
#include "fvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace meshWriters
{

/*---------------------------------------------------------------------------*\
                    Class meshWriters::STARCD Declaration
\*---------------------------------------------------------------------------*/

class STARCD
:
    public meshWriter
{
    // Private Data

        static const char* defaultBoundaryName;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        STARCD(const STARCD&);

        //- Disallow default bitwise assignment
        void operator=(const STARCD&);

        //- Pro-STAR 4+ header format
        static void writeHeader(Ostream&, const int nPoints, const int ncells);

        //- Write points
        void writePoints(Ostream&) const;

        //- Write cells
        void writeCells(Ostream&,const cellList&) const;

        //- Write boundary
        void writeBoundary(Ostream&,const cellList&) const;
        
        bool isContained (const labelList&, const labelList&) const;
        
        cellList rePosFaces(Ostream&) const;

public:

    // Static data members

        //- Face addressing from OpenFOAM faces -> pro-STAR faces
        static const label foamToStarFaceAddr[4][6];


    // Constructors

        //- Open a file for writing
        STARCD
        (
            const fvMesh&,
            const scalar scaleFactor = 1.0
        );


    //- Destructor
    virtual ~STARCD();


    // Member Functions

        // Edit

            //- Remove STAR-CD files for the baseName
            void rmFiles(const fileName& baseName) const;


        // Write

            //- Write volume mesh
            virtual bool write
            (
                const fileName& meshName = fileName::null
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meshWriters
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
