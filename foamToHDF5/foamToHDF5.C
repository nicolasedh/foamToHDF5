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

Application
    foamToCGNS

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "argList.H"
#include "IOobjectList.H"
#include "ListOps.H"

//#include "hdf5.h"
//#include "hdf5_hl.h"
#include "exportToHDF5.H"

#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:
int main(int argc, char *argv[])
{
    Info << "starting app" << endl;

    timeSelector::addOptions();
    argList::removeOption("noFunctionObjects");
    argList::removeOption("newTimes");
    #include "addRegionOption.H"
//    argList::noParallel();
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of fields to be reconstructed. Eg, '(U T p)'."
        "Regular expressions not currently supported"
    );
    argList::addOption("file","word","Filname to save to. Default result.foamcgns");
    argList::addNote("Convert mesh and solution files to custom CGNS file. "
                     "Only volScalarFields and volVectorFields are supported.");

    #include "setRootCase.H"
    #include "createTime.H"


    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"


    //Create List of selected fields
    // we need to know the type in order to read so devide
    // the list of fields into separate fields for each type.
    // then read the fields separatly
    IOobjectList objects(mesh, runTime.timeName());
    IOobjectList scalarFields = objects.lookupClass("volScalarField");
    IOobjectList vectorFields = objects.lookupClass("volVectorField");
    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }
    wordList fields;
    Info << "Selecting fields ";
    forAllConstIter(IOobjectList, scalarFields, fieldIter)
    {
        if
        (
            selectedFields.empty()
         || selectedFields.found(fieldIter()->name())
        )
        {
            fields.append(fieldIter()->name());
            Info << fieldIter()->name() << " ";
        }
    }
    forAllConstIter(IOobjectList, vectorFields, fieldIter)
    {
        if
        (
            selectedFields.empty()
         || selectedFields.found(fieldIter()->name())
        )
        {
            fields.append(fieldIter()->name());
            Info << fieldIter()->name() << " ";
        }
    }
    Info << endl << endl;

    word filename = args.optionLookupOrDefault<word>("file","results.h5foam");


    dictionary dict;
    dict.add("filename",filename);
    dict.add("fields",fields);
    dict.add("region",regionName);


    //label in my compilation is WM_LABEL_SIZE (32)
/////////////////////////////////////////////////////////////
    exportToHDF5 exporter("hejhej",runTime,dict);

    forAll(timeDirs,timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        exporter.timeSet();
        Info << "Time = " << runTime.timeName() << endl;
        //Read all selected fields into db
        forAllConstIter(IOobjectList, scalarFields, fieldIter)
        {
            //don't read fields user didn't ask for
            if(findIndex(fields,fieldIter()->name()) < 0 )
                    continue;
            static autoPtr<volScalarField> scalar;
            scalar.reset
            (
                 new volScalarField
                 (
                    IOobject
                    (
                        fieldIter()->name(),
                        runTime.timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE,
                        true
                    ),
                    mesh
                )
            );
            scalar().checkIn();
        }


        forAllConstIter(IOobjectList, vectorFields, fieldIter)
        {
            //don't read fields user didn't ask for
            if(findIndex(fields,fieldIter()->name()) < 0 )
                    continue;
            static autoPtr<volVectorField> vec;
            vec.reset
            (
                 new volVectorField
                 (
                    IOobject
                    (
                        fieldIter()->name(),
                        runTime.timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE,
                        true
                    ),
                    mesh
                 )
           );
           vec().checkIn();
        }
        exporter.write();

    }

    exporter.end();
////////////////////////////////////////////////////////////////
    Info<< endl << "End\n" << endl;

    return 0;

}
// ************************************************************************* //
