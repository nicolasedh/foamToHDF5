/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::exportToHDF5

Description

SourceFiles
    exportToHDF5.C


\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
//#include <mpi.h>
#include "hdf5.h"
#include <vector>
#include "runTimeSelectionTables.H"

#ifndef exportToHDF5_H
#define exportToHDF5_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class exportToHDF5 Declaration
\*---------------------------------------------------------------------------*/

class exportToHDF5
{
protected:
    // Private data

        //- Name of this set of writeRegisteredObject
        word name_;

        //- Refererence to Db
        const objectRegistry& obr_;

private:
    // Private data


        //- Description of data_
        word filename_;
        hid_t file_id_;
        bool updateMesh_;
        dictionary dict_;

    // Private Member Functions
        //- Disallow default bitwise copy construct
        //exportToHDF5(const exportToHDF5&){}

        //- Disallow default bitwise assignment
        void operator=(const exportToHDF5&);

        //Open a file for apending
        void openFile();
        void closeFile();

        //write an array on all procs to the group
        // "prepathXXpostpath" where XX is the processor number
        //if ncols>1 then array is 2D ie array[0][0] second index is col index
        // currently int or double arrays are supported.
        template <class T> void writeArrayOnAllProcs
        (
                std::string prepath,
                std::string postpath,
                int nelems,
                int ncols,
                const T * arraydata
        );
        //4 helper functions for writeArrayOnAllProcs
        //(parts of the function that need to be specific to the type)
        hid_t customH5Dcreate(char *name, hid_t spaceid, hid_t plinkcreate,
                              hid_t plistDCreate, double unused=0.0);
        hid_t customH5Dcreate(char *name, hid_t spaceid, hid_t plinkcreate,
                              hid_t plistDCreate, int unused=0);
        void customH5Dwrite(hid_t datasetid, hid_t plist_id, const double * arraydata);
        void customH5Dwrite(hid_t datasetid, hid_t plist_id, const int * arraydata);

        //! Write an attribute to all procs.
        //! Only one int is supported at the moment.
        template<class T> void writeAttributeOnGroupOnAllProcs
        (
                word prepath,
                word postpath,
                word name,
                const T data
        );
        //!helper function that writes on one proc only.
        //!remember all procs have to write all attributes
        //!it's not this function responsability
        template<class T> void writeAttributeOnGroup
        (
                word location,
                word attrname,
                const T data
        );
        //4 helper functions for writeAttributeOnAllProcs
        //(parts of the function that need to be specific to the type)

        //return a reference to the mesh
        const fvMesh& mymesh();

        //write faces, points,owners and neighbour
        void writeMesh();

        //! writes patch information, name, startFace, nFaces
        void writeBoundaries(word location);

        void writeFields();
        template <class cmpType> void writeVolField(GeometricField<cmpType, fvPatchField, volMesh> field, int ncomp=1);



public:

    // Static data members

        //- Static data staticData

    TypeName("exportToHDF5")
    // Constructors

        //- Construct null
        //exportToHDF5();

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        exportToHDF5
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool loadFromFiles = false
        );
        //- Construct from Istream
        exportToHDF5(Istream&);



    // Selectors

        //- Select null constructed
        //static autoPtr<exportToHDF5> New();


    //- Destructor
    virtual ~exportToHDF5();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    //- Return name of the writeRegisteredObject
    virtual const word& name() const
    {
        return name_;
    }

    //- Read the writeRegisteredObject data
    virtual void read(const dictionary&);

    //- Execute, currently does nothing
    virtual void execute()
    {}

    //- Execute at the final time-loop, currently does nothing
    virtual void end();

    //- Called when time was set at the end of the Time::operator++
    virtual void timeSet();

    //- Write the writeRegisteredObject
    virtual void write();

    //- Update for changes of mesh
    virtual void updateMesh(const mapPolyMesh&)
    {}

    //- Update for changes of mesh
    virtual void movePoints(const polyMesh&)
    {}

    // Friend Functions

    // Friend Operators
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
