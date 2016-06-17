/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkH5FoamReader.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include <string.h>
#include <sstream>
#include <locale>

//#include "vtkByteSwap.h"
//#include "vtkCellArray.h"
//#include "vtkCellData.h"
#include "vtkErrorCode.h"
//#include "vtkFloatArray.h"
//#include "vtkIncrementalPointLocator.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
//#include "vtkMergePoints.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkUnstructuredGrid.h"
#include "vtkIdList.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkGenericCell.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkCellData.h"
#include "vtkArray.h"
#include "vtkDoubleArray.h"
#include "vtkMath.h"
#include "vtkCollection.h"


#include "vtkH5FoamReader.h"
#include "H5FDmpio.h"
#include "mpi.h"

#define CREATE_P(CLASS,NAME)  \
    vtkSmartPointer<CLASS> NAME = vtkSmartPointer<CLASS>::New();

vtkStandardNewMacro(vtkH5FoamReader);

//2do:
//   * läsa patches
//      ^ strukturera block. fundera lite. Måste nog ha procs som parent.
//      ^ men man kan så klart sköta valet av fields/patches själv
//   * testa ett litet större fall med ~30 procs mot of:s läsare
//   * parallellt
//vtkCxxSetObjectMacro(vtkH5FoamReader, Locator, vtkIncrementalPointLocator);

//------------------------------------------------------------------------------
// Construct object 
vtkH5FoamReader::vtkH5FoamReader()
{
  this->FileName = NULL;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  updateMesh = true;
  this->times = vtkSmartPointer<vtkDoubleArray>::New();
  this->timesNames = vtkSmartPointer<vtkVariantArray>::New();
  this->file_id = -1;
}

//------------------------------------------------------------------------------
vtkH5FoamReader::~vtkH5FoamReader()
{
    closeFile();
  //this->SetFileName(0);
}

int vtkH5FoamReader::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector)
{
    if (!this->FileName || strlen(this->FileName) == 0)
      {
      vtkErrorMacro("FileName has to be specified!");
      this->SetErrorCode(vtkErrorCode::NoFileNameError);
      return 0;
      }
    openFile();

    findTimes();
    vtkIdType nsteps = times->GetNumberOfTuples();
    double * trange = new double[2];
    trange[0]=times->GetValue(0);
    trange[1]=times->GetValue(nsteps-1);
    vtkInformation * outInfo = outputVector->GetInformationObject(0);
    outInfo->Set(
                vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
                times->GetPointer(0), nsteps
                );
    outInfo->Set(
                vtkStreamingDemandDrivenPipeline::TIME_RANGE(),
                trange,2);

    //2do create a list of fields (in current time only?)
    //2do create a list of blocks i.e. patches,zones regions etc
    closeFile();

    return 1;
}

//------------------------------------------------------------------------------
void vtkH5FoamReader::openFile()
{
    if(false)
    {
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
        file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, plist_id);
        H5Pclose(plist_id);
    }
    else
        file_id = H5Fopen(this->FileName,H5F_ACC_RDONLY,H5P_DEFAULT);
    //2do error if HDF5 don't like the file
    //2do write attribute so we know it's our file type.

}

void vtkH5FoamReader::closeFile()
{
    if(file_id > 0)
        H5Fclose(file_id);
    file_id=-1;

}

int vtkH5FoamReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  output = vtkMultiBlockDataSet::SafeDownCast(
              outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // All of the data in the first piece.
  if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
    {
        return 0;
    }

  if (!this->FileName || *this->FileName == 0)
    {
        vtkErrorMacro(<<"A FileName must be specified.");
        this->SetErrorCode(vtkErrorCode::NoFileNameError);
        return 0;
    }

    openFile();
    readMesh();
    double requestedTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
    vtkIdType timeIndex = findClosestSolutionIndex(requestedTime);

    readFields(timesNames->GetValue(timeIndex));
    closeFile();

  return 1;
}

void vtkH5FoamReader::readMesh()
{
    //2do get region name
    vtkStdString meshLoc="/constant/region0/polyMesh";
    hid_t gid_pMesh = H5Gopen2(file_id,meshLoc.c_str(),H5P_DEFAULT);
    //get num procs
    H5G_info_t ginfo;
    H5Gget_info(gid_pMesh,&ginfo);
    vtkSmartPointer<vtkVariantArray> procs = vtkSmartPointer<vtkVariantArray>::New();
    findProcs(meshLoc.c_str(),procs);
    int nprocs=procs->GetNumberOfTuples();
    char procLoc[1000];
    for(int pI=0;pI<procs->GetNumberOfTuples();pI++)
    {
        sprintf(procLoc,"%s/%s",meshLoc.c_str(),procs->GetValue(pI).ToString().c_str());
        int procN = vtkVariant(procs->GetValue(pI).ToString().substr(9)).ToInt();
        readProcMesh(procLoc,procN);
    }

    H5Gclose(gid_pMesh);


}
void vtkH5FoamReader::readProcMesh(char * location,int procN)
{
    hid_t gid_proc = H5Gopen2(file_id,location,H5P_DEFAULT);
    //2do: decompose poly cells!
    //      vtkPolyhedron::DecomposeAPolyhedronCell()
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnstructuredGrid> ugrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    readPoints(gid_proc,"points",points);
    vtkDoubleArray * tst = vtkDoubleArray::New();
    int ncols;
    readDoubleArray(gid_proc,"points",tst,&ncols);

    ugrid->SetPoints(points);

    vtkSmartPointer<vtkIntArray> faces = vtkSmartPointer<vtkIntArray>::New();
    readArray(gid_proc,"faces",faces.GetPointer());
    vtkIdType nfaces = faces->GetNumberOfTuples();

    vtkSmartPointer<vtkIntArray> owner = vtkSmartPointer<vtkIntArray>::New();
    readArray(gid_proc,"owner",owner.GetPointer());
    vtkIdType nowner = owner->GetNumberOfTuples();

    vtkSmartPointer<vtkIntArray> neighbour = vtkSmartPointer<vtkIntArray>::New();
    readArray(gid_proc,"neighbour",neighbour.GetPointer());
    vtkIdType nneighbour = neighbour->GetNumberOfTuples();

    H5Gclose(gid_proc);
    //Count the number of cells
    //which is just max of max(owner) max(neighbour), last cell could own no faces
    vtkIdType nCells = std::max(owner->GetRange()[1],neighbour->GetRange()[1]) +1;

    //count number of faces per cell
    vtkSmartPointer<vtkIntArray> nFacesPerCell = vtkSmartPointer<vtkIntArray>::New();
    nFacesPerCell->SetNumberOfTuples(nCells);
    for(vtkIdType i=0;i<nCells;i++)
        nFacesPerCell->SetValue(i,0);

    nFacesPerCell->SetNumberOfTuples(nCells);
    for(vtkIdType i=0;i<nneighbour;i++)
    {
        vtkIdType cI = neighbour->GetValue(i);
        nFacesPerCell->SetValue(cI,nFacesPerCell->GetValue(cI)+1);
    }
    for(vtkIdType i=0;i<nowner;i++)
    {
        vtkIdType cI = owner->GetValue(i);
        nFacesPerCell->SetValue(cI,nFacesPerCell->GetValue(cI)+1);
    }


    //Create face connectivity
    vtkSmartPointer<vtkPolyData> faceCon = vtkSmartPointer<vtkPolyData>::New();
    faceCon->SetPoints(points);
    faceCon->Allocate(nowner);//allocate actuall number of faces
    // from faces: (NpointsInface1 f1p1, f1p2 .. ... , NpointsFaceN,fNp1, fNp2 ..)
    vtkIdList * facepnts = vtkIdList::New();
    for(vtkIdType i=0;i<nfaces;)
    {
        int fnp = faces->GetValue(i); //nr points if current face
        i++;

        facepnts->SetNumberOfIds(fnp);
        for(vtkIdType j=0;j<fnp;j++,i++)
            facepnts->SetId(j,faces->GetValue(i));
        faceCon->InsertNextCell(VTK_POLYGON,facepnts);

    }
    facepnts->Delete();

    //This is maybe a bit memory intensive
    //but I havn't figured out a better way
    //store the cell connectivies, i.e.
    //list for each cell. The list defines a polyhedra
    //(nFaces nPointsFace0 Face0Point0 Face0Point1 ... nPointsFace1 ...)
    //Since owner might not be sorted we have to create these lists
    std::vector<vtkIdList*> cellCons;
    cellCons.reserve(nCells);
    //initialize all lists and set nFaces
    for(vtkIdType cI=0;cI<nCells;cI++)
    {
        //initialisera och lägg in nFacesPerCell
        vtkIdList * ccons  = vtkIdList::New();
        ccons->InsertNextId(nFacesPerCell->GetValue(cI));
        cellCons[cI] = ccons;
    }

    // add the current face to the owner cell
    for (vtkIdType fI=0;fI<nowner;fI++)
    {
        vtkIdType cI = owner->GetValue(fI);
        vtkIdList * ccons  = cellCons[cI];
        vtkCell * face = faceCon->GetCell(fI);
        ccons->InsertNextId(face->GetNumberOfPoints());
        for(vtkIdType i=0;i<face->GetNumberOfPoints();i++)
            ccons->InsertNextId(face->GetPointId(i));

    }
    //add the current face to the neigbour cell
    for (vtkIdType fI=0;fI<nneighbour;fI++)
    {
        vtkIdType cI = neighbour->GetValue(fI);
        vtkIdList * ccons  = cellCons[cI];
        vtkCell * face = faceCon->GetCell(fI);
        ccons->InsertNextId(face->GetNumberOfPoints());
        for(vtkIdType i=0;i<face->GetNumberOfPoints();i++)
            ccons->InsertNextId(face->GetPointId(i));

    }
    //Now all cells have a list which is used
    //to define a cell. Then delete the list for clean up
    for(vtkIdType cI=0; cI < nCells;cI++)
    {
        vtkIdList * ccon = cellCons[cI];
        ugrid->InsertNextCell(VTK_POLYHEDRON,ccon);
        ccon->Delete();
    }

    output->SetBlock(procN,ugrid);
    char name[100];
    output->GetMetaData(procN)->Set(vtkCompositeDataSet::NAME(),
    name);


}

void vtkH5FoamReader::readArray(const hid_t gid_proc, const char *name, vtkIntArray *ids)
{
    //2do gör specifik funktion för att läsa INT/DOUBLE
    hid_t did = H5Dopen2(gid_proc,name,H5P_DEFAULT);
    hid_t data_type =  H5Dget_type(did);
    H5T_class_t dclas = H5Tget_class(data_type);

    hid_t spid = H5Dget_space(did);

    hid_t rank = H5Sget_simple_extent_ndims(spid);

    hsize_t dims[rank];
    H5Sget_simple_extent_dims(spid, dims, NULL);
    hsize_t nIds = dims[0];
    ids->SetNumberOfTuples(nIds);


    hid_t prop = H5Dget_create_plist(did);
    //2do if chunked ?
    hid_t memspace = H5Screate_simple(rank,dims,NULL);
    H5Dread(did, H5T_NATIVE_INT, memspace, spid,
                 H5S_ALL, ids->GetPointer(0));//not does not work with vtkIdList

    H5Sclose(memspace);
    H5Pclose(prop);
    H5Sclose(spid);
    H5Dclose(did);

}

void vtkH5FoamReader::readDoubleArray
(
        const hid_t gid_proc,
        const char *name,
        vtkDoubleArray *data,
        int *ncols
)
{
    hid_t did = H5Dopen2(gid_proc,name,H5P_DEFAULT);
    hid_t data_type =  H5Dget_type(did);
    H5T_class_t dclass = H5Tget_class(data_type);
    if(dclass != H5T_FLOAT)
        cout << "Tried to read " << name << "as a double array " << endl;

    hid_t spid = H5Dget_space(did);

    hid_t rank = H5Sget_simple_extent_ndims(spid);

    hsize_t dims[rank];
    H5Sget_simple_extent_dims(spid, dims, NULL);
    hsize_t nrows = dims[0];


    int ncc;
    if (rank==1)
        ncc=1;
    else if(rank ==2)
        ncc=dims[1];
    *ncols=ncc;
    //handle rank >2? Naaa
    data->SetNumberOfComponents(ncc);
    data->SetNumberOfTuples(nrows);

    double tmpdata[nrows*ncc];
    hid_t prop = H5Dget_create_plist(did);
    //2do if chunked ?
    hid_t memspace = H5Screate_simple(rank,dims,NULL);
    H5Dread(did, H5T_NATIVE_DOUBLE, memspace, spid,
                 H5S_ALL, tmpdata);

    for(int i =0;i<nrows;i++)
    {
        for(int j=0;j<ncc;j++)
            data->SetComponent(i,j,tmpdata[i*ncc+j]);
    }
    H5Sclose(memspace);
    H5Pclose(prop);
    H5Sclose(spid);
    H5Dclose(did);

}

void vtkH5FoamReader::readPoints(const hid_t gid_proc,const char *name,vtkPoints * points)
{
    //2do remove and replace with read double or at least use that function.

    hid_t did = H5Dopen2(gid_proc,name,H5P_DEFAULT);

    hid_t spid = H5Dget_space(did);

    hid_t rank = H5Sget_simple_extent_ndims(spid);

    hsize_t dims[rank];
    H5Sget_simple_extent_dims(spid, dims, NULL);
    hsize_t npoints = dims[0];
    points->SetNumberOfPoints(dims[0]);

    //2do error if rank != 2
    double tmp[dims[1]*dims[0]];

    hid_t prop = H5Dget_create_plist(did);
    //2do if chunked ?
    hid_t memspace = H5Screate_simple(rank,dims,NULL);

    H5Dread(did, H5T_NATIVE_DOUBLE, memspace, spid,
                 H5S_ALL, tmp);

    H5Sclose(memspace);
    H5Pclose(prop);
    H5Sclose(spid);
    H5Dclose(did);

    double pnt[3];
    //store all points in vtkPoints
    for(vtkIdType i=0;i<npoints;i++)
    {
        for(int j=0;j<3;j++)
            pnt[j] = tmp[i*3+j];
        points->SetPoint(i,pnt);
    }
}

void vtkH5FoamReader::findProcs(const char *location,vtkVariantArray * procNames)
{
//2do sortera numeriskt?
    hid_t gid_pMesh = H5Gopen2(file_id,location,H5P_DEFAULT);
    //get num procs
    H5G_info_t ginfo;
    H5Gget_info(gid_pMesh,&ginfo);
    int nprocs=0;char name[1000];
    for (int i = 0;i< ginfo.nlinks;i++)
    {
        H5Gget_objname_by_idx(gid_pMesh,i,name,1000);
        vtkStdString strname(name);
        int pos = strname.find("processor");
        if(pos==strname.length())
                continue;
        int procN=vtkVariant(strname.substr(pos+9)).ToInt() ;
        nprocs++;

        procNames->InsertNextValue(vtkVariant(name));
    }
    H5Gclose(gid_pMesh);
}

void vtkH5FoamReader::readFields(vtkVariant timeName)
{
    char location[1000];
    sprintf(location,"/%s/region0",timeName.ToString().c_str());
    //read all procs
    //2do get region name
    vtkSmartPointer<vtkVariantArray> procs = vtkSmartPointer<vtkVariantArray>::New();
    findProcs(location,procs);
    for(vtkIdType pI=0;pI<procs->GetNumberOfTuples();pI++)
    {
        char ploc[1000];
        sprintf(ploc,"%s/%s",location,procs->GetVariantValue(pI).ToString().c_str());
        //This is == pI for serial runs but not parallel
        int procN = vtkVariant(procs->GetValue(pI).ToString().substr(9)).ToInt();
        readProcFields(ploc,procN);

    }


}
void vtkH5FoamReader::readProcFields(char *location, int procId)
{
    hid_t gid = H5Gopen2(file_id,location,H5P_DEFAULT);
    H5G_info_t ginfo;
    H5Gget_info(gid,&ginfo);
    char fieldName[1000];
    char internalFieldLoc[1000];
    for (int i = 0;i< ginfo.nlinks;i++)
    {
        H5Gget_objname_by_idx(gid,i,fieldName,1000);
        sprintf(internalFieldLoc,"%s/%s/internalField",location,fieldName);
        //2do fixa vid patches,
        vtkUnstructuredGrid * ugrid = vtkUnstructuredGrid::SafeDownCast(this->GetOutput()->GetBlock(procId));
        vtkIdType ncells = ugrid->GetNumberOfCells();
        vtkSmartPointer<vtkDoubleArray> field = vtkSmartPointer<vtkDoubleArray>::New();
        field->SetName(fieldName);
        int ncols=0;
        readDoubleArray(gid,internalFieldLoc,field,&ncols);
        ugrid->GetCellData()->AddArray(field);

    }
    H5Gclose(gid);
}

int vtkH5FoamReader::findTimes()
{
    //LIST TIMES
        char name[1000];
        hid_t gid = H5Gopen2(file_id,"/",H5P_DEFAULT);
        H5G_info_t info;
        H5Gget_info(gid,&info);

        for (int i = 0;i< info.nlinks;i++)
        {
            H5Gget_objname_by_idx(gid,i,name,1000);
            vtkVariant vname(name);
            bool isDouble;
            double t = vname.ToDouble(&isDouble);
            if(isDouble)
            {
                times->InsertNextValue(t);
                timesNames->InsertNextValue(vname);
            }
        }
        H5Gclose(gid);


  return times->GetNumberOfTuples();
}

vtkIdType vtkH5FoamReader::findClosestSolutionIndex(double Time)
{
    double dt=1e24;
    vtkIdType ti=0;//time index
    for(vtkIdType i=0;i<times->GetNumberOfTuples();i++)
    {
        if(dt > std::abs(Time - times->GetValue(i) ) )
        {
            ti=i;
            dt=std::abs(Time-times->GetValue(i));
        }
    }

  return ti;
}


//------------------------------------------------------------------------------
void vtkH5FoamReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: "
     <<(this->FileName ? this->FileName : "(none)") << "\n";

}
