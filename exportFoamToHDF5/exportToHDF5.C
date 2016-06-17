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

\*---------------------------------------------------------------------------*/


#include "exportToHDF5.H"
#include "IOobjectList.H"
#include "genericFvPatch.H"


namespace Foam
{
defineTypeNameAndDebug(exportToHDF5, 0);
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::exportToHDF5::openFile()
{
    //2do: more convinent to clear file when openening.
    //     chane H5F_ACC_TRUNC to H5F_ACC_RDWR | H5F_ACC_CREAT
    //     both open Read&write and create if the file doesn exist
    Info << "Opening hdf5:" << filename_ << endl;

    if(Pstream::parRun())
    {
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
        file_id_ = H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);
    }
    else
        file_id_ = H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

}

void Foam::exportToHDF5::closeFile()
{

//    if(Pstream::parRun())
//        MPI_Barrier(MPI_COMM_WORLD);

    H5Fclose(file_id_);
    file_id_ =-1;
}

const fvMesh& Foam::exportToHDF5::mymesh()
{
    word region = dict_.lookupOrDefault<word>("region",fvMesh::defaultRegion);
    return obr_.time().lookupObject<fvMesh>(region);
}

template <class T> void Foam::exportToHDF5::writeArrayOnAllProcs
(
        std::__cxx11::string prepath,
        std::__cxx11::string postpath,
        int nelems,
        int ncols,
        const T *arraydata
)
{

    //distribute list of number of elements
    List<label> nItemsPerProc(Pstream::nProcs());
    nItemsPerProc[Pstream::myProcNo()] = nelems;
    Pstream::gatherList(nItemsPerProc);
    Pstream::scatterList(nItemsPerProc);

    hsize_t dims[2];
    //dims[0] = nelems; //will be set later
    dims[1] = ncols;
    int RANK = 1;
    if(ncols>1)
        RANK=2;

    //While writing parallel all procs must create all datasets
    char name[1000];
    for(int pi=0;pi<Pstream::nProcs();pi++)
    {
        dims[0] = nItemsPerProc[pi];
        //create space
        hid_t spaceid = H5Screate_simple(RANK,dims,NULL);


        hid_t plinkcreate = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(plinkcreate, 1);

        std::sprintf(name,"%s%d%s",prepath.c_str(),pi,postpath.c_str());
        //create data set
        hid_t plistDCreate = H5Pcreate(H5P_DATASET_CREATE);
        hid_t setid = customH5Dcreate(name, spaceid,
                                 plinkcreate, plistDCreate,T());
        H5Sclose(spaceid);
        H5Dclose(setid);
        H5Pclose(plistDCreate);
        H5Pclose(plinkcreate);
    }

    std::sprintf(name,"%s%d%s",prepath.c_str(),Pstream::myProcNo(),postpath.c_str());

    hid_t datasetid= H5Dopen2(file_id_,name,H5P_DEFAULT);

    if(Pstream::parRun())
    {

        hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
//        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        if (nelems > 0)
            customH5Dwrite(datasetid, plist_id,arraydata);
        H5Pclose(plist_id);
    }
    else
        customH5Dwrite(datasetid, H5P_DEFAULT,arraydata);
//        H5Dwrite(datasetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
//                        arraydata);

    H5Dclose(datasetid);
}

hid_t Foam::exportToHDF5::customH5Dcreate
(
        char *name,
        hid_t spaceid,
        hid_t plinkcreate,
        hid_t plistDCreate,
        double unused
)
{
    return H5Dcreate2(file_id_,name, H5T_NATIVE_DOUBLE, spaceid,
                             plinkcreate, plistDCreate, H5P_DEFAULT);
}
hid_t Foam::exportToHDF5::customH5Dcreate
(
        char * name,
        hid_t spaceid,
        hid_t plinkcreate,
        hid_t plistDCreate,
        int unused
)
{
    return H5Dcreate2(file_id_,name, H5T_NATIVE_INT, spaceid,
                             plinkcreate, plistDCreate, H5P_DEFAULT);
}
void Foam::exportToHDF5::customH5Dwrite(hid_t datasetid,hid_t plist_id, const double *arraydata)
{
    H5Dwrite(datasetid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plist_id,
                        arraydata);
}

void Foam::exportToHDF5::customH5Dwrite(hid_t datasetid, hid_t plist_id, const int *arraydata)
{
    H5Dwrite(datasetid, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_id,
                        arraydata);
}

template <class T> void Foam::exportToHDF5::writeAttributeOnGroupOnAllProcs
(
        word prepath,
        word postpath,
        word name,
        const T data
)
{
    //2do: currently only takes 1 element of data.
    //     could be extended by list<T> instead.
    //     see how we write patches.
    //distribute list of number of elements
    List<T> procData(Pstream::nProcs());
    procData[Pstream::myProcNo()] = data;
    Pstream::gatherList(procData);
    Pstream::scatterList(procData);

    char location[100];
    for(int pi=0;pi<Pstream::nProcs();pi++)
    {

        std::sprintf(location,"%s%d%s",prepath.data(),pi,postpath.data());
        writeAttributeOnGroup(location,name,data);
/*
        //H5Gopen can fail in which case we create the group. Suppress error messages
        //Save old error handler
//        herr_t (*old_func)(void*);
        H5E_auto2_t  old_func;
        void * old_client_data;
        H5Eget_auto2(H5E_DEFAULT,&old_func, &old_client_data);
        // Turn off error handling
        H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
        hid_t gid  = H5Gopen(file_id_,location,H5P_DEFAULT);
        //Restore previous error handler
        H5Eset_auto2(H5E_DEFAULT,old_func, old_client_data);
        //create group if it doesn't exist
        if (gid <0 )
        {
            hid_t plinkcreate = H5Pcreate(H5P_LINK_CREATE);
            H5Pset_create_intermediate_group(plinkcreate, 1);
            gid = H5Gcreate2(file_id_,location,plinkcreate,H5P_DEFAULT,H5P_DEFAULT);
            H5Pclose(plinkcreate);
        }
        hid_t spid = H5Screate_simple(RANK,dims,NULL);
        hid_t atid = H5Acreate2 (gid, name.data(), H5T_NATIVE_INT, spid,
                                 H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(atid,H5T_NATIVE_INT,&procData[pi]);
        H5Aclose(atid);
        H5Sclose(spid);
        H5Gclose(gid);
        */
    }
}
template <class T> void Foam::exportToHDF5::writeAttributeOnGroup
(
        word location,
        word attrname,
        const T data
)
{
    int RANK=1;
    hsize_t dims[1]={1};

    //H5Gopen can fail in which case there is no groupw.
    //Suppress error messages and create the group

    //Save old error handler
    H5E_auto2_t  old_func;
    void * old_client_data;
    H5Eget_auto2(H5E_DEFAULT,&old_func, &old_client_data);
    // Turn off error handling
    H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
    hid_t gid  = H5Gopen(file_id_,location.c_str(),H5P_DEFAULT);
    //Restore previous error handler
    H5Eset_auto2(H5E_DEFAULT,old_func, old_client_data);
    //create group if it doesn't exist
    if (gid <0 )
    {
        hid_t plinkcreate = H5Pcreate(H5P_LINK_CREATE);
        H5Pset_create_intermediate_group(plinkcreate, 1);
        gid = H5Gcreate2(file_id_,location.c_str(),plinkcreate,H5P_DEFAULT,H5P_DEFAULT);
        H5Pclose(plinkcreate);
    }
    hid_t spid = H5Screate_simple(RANK,dims,NULL);
    hid_t atid = H5Acreate2 (gid, attrname.data(), H5T_NATIVE_INT, spid,
                             H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(atid,H5T_NATIVE_INT,&data);
    H5Aclose(atid);
    H5Sclose(spid);
    H5Gclose(gid);
}

void Foam::exportToHDF5::writeMesh()
{
    if(! updateMesh_)
        return;
    //write points
    std::string location = "/constant/region0/polyMesh/processor";
    const pointField pointsF = mymesh().points();
    double points[pointsF.size()][3];
    forAll(pointsF,pI)
    {
        points[pI][0]=pointsF[pI].x();
        points[pI][1]=pointsF[pI].y();
        points[pI][2]=pointsF[pI].z();
    }
    writeArrayOnAllProcs(location,"/points",pointsF.size(),3,&points[0][0]);

    //Write faces (collect them to a 1D array)
    // F0nP F0p0 F0p1 ... F1nP F1p0 F1p1 ...
    std::vector<label> pointsOfFaces;
    pointsOfFaces.reserve(4*mymesh().faces().size());
    forAll(mymesh().faces(),fi)
    {
        face f = mymesh().faces()[fi];
        pointsOfFaces.push_back(f.size());
        //loop over points in face
        forAll(f,fpi)
        {
            pointsOfFaces.push_back(f[fpi]);
        }
    }
    //2do write attribute here nFaces

    writeArrayOnAllProcs
            (
                location,
                "/faces",
                pointsOfFaces.size(),
                1,
                pointsOfFaces.data()
            );

    //Write owner
//    Info << mesh().faceOwner() << endl << mesh().faceOwner().size();

    writeArrayOnAllProcs
            (
                location,
                "/owner",
                mymesh().faceOwner().size(),
                1,
                mymesh().faceOwner().cdata()
            );
    writeArrayOnAllProcs
            (
                location,
                "/neighbour",
                mymesh().faceNeighbour().size(),
                1,
                mymesh().faceNeighbour().cdata()
            );
    writeAttributeOnGroupOnAllProcs(location,"","nFaces",mymesh().nFaces());//redundant = owners.size
    writeAttributeOnGroupOnAllProcs(location,"","nCells",mymesh().nCells());

    //Write boundaries
    writeBoundaries(location);

    updateMesh_ = false;
}

void Foam::exportToHDF5::writeBoundaries(word location)
{
    //HDF5 demands that all procs write all attributes.
    //all procs have all ordinary patches but not procboudnaries.
    //so create a list of wordList (size nProcs). Where the wordList
    //containsthe bc of each proc.
    List<wordList> allPatches(Pstream::nProcs());
    allPatches[Pstream::myProcNo()] = mymesh().boundaryMesh().names();
    Pstream::gatherList(allPatches);
    Pstream::scatterList(allPatches);

    //now we need to do the same for all attributes we wish to write
    // i.e. startFace and nFaces
    List<labelList> allStartFaces(Pstream::nProcs());
    List<labelList> allNFaces(Pstream::nProcs());
    //create local labelLists of startFace and end Face
    labelList localStartFace(mymesh().boundaryMesh().size(),0);
    labelList localNFace(mymesh().boundaryMesh().size(),0);
    forAll(mymesh().boundaryMesh(),pI)
    {
        const polyPatch * patch = mymesh().boundaryMesh()(pI);
        localStartFace[pI] = patch->start();
        localNFace[pI] = patch->faceCentres().size();
    }

    allStartFaces[Pstream::myProcNo()]=localStartFace;
    allNFaces[Pstream::myProcNo()]=localNFace;
    Pstream::gatherList(allStartFaces);
    Pstream::scatterList(allStartFaces);
    Pstream::gatherList(allNFaces);
    Pstream::scatterList(allNFaces);

    char plocation[100];
    forAll(allPatches,lstI)
    {
        wordList procPatches = allPatches[lstI];
        labelList procStartFaces = allStartFaces[lstI];
        labelList procNFaces = allNFaces[lstI];
        forAll(procPatches,pI)
        {
            sprintf
            (
                plocation,
                "%s%d/boundary/%s",
                location.c_str(),
                lstI, //same as procnr
                procPatches[pI].c_str()
            );//plocation could be /constant/polyMesh/processorXX/boundary/patchName
            writeAttributeOnGroup(plocation,"startFace",procStartFaces[pI]);
            writeAttributeOnGroup(plocation,"nFaces",procNFaces[pI]);
        }
    }

}

void Foam::exportToHDF5::writeFields()
{

    //2do: what if list is empty, output all or fatal error and dump toc?
    wordList fields = dict_.lookupOrDefault<wordList>("fields",wordList());
    forAll(fields, fI)
    {
        word fieldName=fields[fI];
        if(! mymesh().HashTable::found(fieldName))
        {
            Warning << "Did not find field "
                          << fieldName << endl;
            continue;
        }

//        char prelocation[100];
//        std::sprintf(prelocation,"/%s/region0/processor",mymesh().time().timeName().c_str());
//        char postlocation[100];
//        std::sprintf(postlocation,"/%s/internalField",fieldName.c_str());
        word fieldType=mymesh().find(fieldName)()->type();

        if(fieldType=="volScalarField")
        {
            volScalarField scalarf = mymesh().lookupObject<volScalarField>(fieldName);
            writeVolField(scalarf,1);
        }
        else if (fieldType == "volVectorField")
        {
            volVectorField vec = mymesh().lookupObject<volVectorField>(fieldName);
            writeVolField(vec,3);
        }
    }

}
template <class cmpType> void Foam::exportToHDF5::writeVolField
(
        GeometricField<cmpType, fvPatchField, volMesh> field,
        int ncomp
)
{//2do pass const &Type instead
    char prelocation[500];
    std::sprintf(prelocation,"/%s/region0/processor",mymesh().time().timeName().c_str());
    char postlocation[500];
    std::sprintf(postlocation,"/%s/internalField",field.name().c_str());
    char location[1010];
    //but all compontents in a "2D" array
    double arr2D[field.size()*ncomp];
    for(int i=0;i<ncomp;i++)
    {
        tmp<Field<double> > cmpt = field.internalField().component(i);
//internalField().component(i)
        forAll(field,cI)
                arr2D[cI*ncomp+i] = cmpt()[cI];

    }
    writeArrayOnAllProcs(prelocation,postlocation,field.size(),ncomp,arr2D);

    //---------------------write patches
    //HDF5 demands that all procs write all attributes.
    //all procs have all ordinary patches but not procboundaries.
    //so create a list of wordList (size nProcs). Where the wordList
    //containsthe bc of each proc.
    List<wordList> allPatches(Pstream::nProcs());
    allPatches[Pstream::myProcNo()] = mymesh().boundaryMesh().names();
    Pstream::gatherList(allPatches);
    Pstream::scatterList(allPatches);

    //now we tell all procs of the number of faces (and get the same)
    List<labelList> allNFaces(Pstream::nProcs());
    //create local labelLists of startFace and end Face
    labelList localNFace(mymesh().boundaryMesh().size(),0);
    forAll(mymesh().boundaryMesh(),pI)
    {
        const polyPatch * patch = mymesh().boundaryMesh()(pI);
        localNFace[pI] = patch->faceCentres().size();
    }

    allNFaces[Pstream::myProcNo()]=localNFace;

    Pstream::gatherList(allNFaces);
    Pstream::scatterList(allNFaces);


    //create data sets on all procs on global
    //patchlist, looping over procs
    hsize_t dims[2];
    dims[1] = ncomp;
    int RANK = 1;
    if(ncomp>1)
        RANK=2;
    forAll(allPatches,lstI)
    {
        wordList procPatches = allPatches[lstI];
        labelList procNFaces = allNFaces[lstI];
        forAll(procPatches,pI)
        {
            word name = procPatches[pI];
            sprintf
            (
                postlocation,
                "/%s/%s/",
                field.name().c_str(),
                name.c_str()
            );

            dims[0] = procNFaces[pI];
            hid_t spaceid = H5Screate_simple(RANK,dims,NULL);
            hid_t plinkcreate = H5Pcreate(H5P_LINK_CREATE);
            H5Pset_create_intermediate_group(plinkcreate, 1);
            std::sprintf(location,"%s%d%s",prelocation,lstI,postlocation);
            hid_t plistDCreate = H5Pcreate(H5P_DATASET_CREATE);
            hid_t setid = customH5Dcreate
                    (
                        location,
                        spaceid,
                        plinkcreate,
                        plistDCreate,//2do speed up by not copying
                        field.internalField().component(0)()[0]//just to know int or double
                    );
            H5Sclose(spaceid);
            H5Dclose(setid);
            H5Pclose(plistDCreate);
            H5Pclose(plinkcreate);
        }
    }

    //Write the patch data here. Can be done locally
    forAll(field.boundaryField(),bfI)
    {
        word name = mymesh().boundaryMesh()[bfI].name();
        word type = mymesh().boundaryMesh()[bfI].type();

        label size = field.boundaryField()(bfI)->size();
        dims[0] = size;

        //skip empty patches since they misbehave. Have no data?
//        if(type == "empty" || type == "processor" ||dims[0] <1)
        if(size<1)
            continue;

        sprintf(postlocation,"/%s/%s/",field.name().c_str(),name.c_str());
        std::sprintf
        (
            location,
            "%s%d%s",
            prelocation,
            Pstream::myProcNo(),
            postlocation
        );

        //Open the data set
        hid_t datasetid= H5Dopen2(file_id_,location,H5P_DEFAULT);
        //convert data to 2D vector
        double vec2Dpatch[size][ncomp];
        for(int i=0;i<ncomp;i++)
        {
            tmp<Field<double> > cmpt = field.boundaryField()[bfI].component(i);
            forAll(field.boundaryField()[bfI],faceI)
                vec2Dpatch[faceI][i] = cmpt->operator [](faceI);
        }

        if(Pstream::parRun())
        {
            hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    //        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
//            if (size > 0)
            customH5Dwrite(datasetid, plist_id,&vec2Dpatch[0][0]);
            H5Pclose(plist_id);
        }
        else
            customH5Dwrite(datasetid, H5P_DEFAULT,&vec2Dpatch[0][0]);

        H5Dclose(datasetid);
    }

}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::exportToHDF5::exportToHDF5
(
        const word &name,
        const objectRegistry & obr,
        const dictionary & dict,
        const bool loadFromFiles
)
:
          name_(name),
          obr_(obr),
          updateMesh_(true)

{
    //name_ is the dictname and filename is the hdf5 file
    filename_ = dict.lookupOrDefault<word>("filename",name+".h5foam");
    dict_=dict;
    openFile();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::exportToHDF5::~exportToHDF5()
{
    closeFile();
    Pout << "I'm melting" << endl;

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::exportToHDF5::write()
{
    //2do: open one file for each time and close
    writeMesh();
    writeFields();

    // flush file in case we crash later on
    H5Fflush(file_id_, H5F_SCOPE_GLOBAL);

}

void Foam::exportToHDF5::timeSet()
{
//    Info << "time set was called" << obr_.time().timeName()<< endl;
}

void Foam::exportToHDF5::read(const dictionary& dict)
{//2do if filename changed close and open file
  dict_ = dict;
}

void Foam::exportToHDF5::end()
{
    Pout << "ending now" << endl;
//    closeFile();
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::exportToHDF5::operator=(const exportToHDF5& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::exportToHDF5::operator=(const Foam::exportToHDF5&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
/**
 * To test if we are missing a close:
 *     ssize_t norphans;
    norphans = H5Fget_obj_count(fp, H5F_OBJ_ALL);
    if (norphans > 1) { // expect 1 for the file we have not closed
        int i;
        H5O_info_t info;
        char name[64];
        hid_t * objects = calloc(norphans, sizeof(hid_t));
        H5Fget_obj_ids(fp, H5F_OBJ_ALL, -1, objects);
        for (i=0; i<norphans; i++) {
            H5Oget_info(objects[i], &info);
            H5Iget_name(objects[i], name, 64);
            printf("%d of %zd things still open: %d with name %s of type %d",
                  i, norphans, objects[i], name, info.type);
        }
        free(objects);
    }
 *
  */
