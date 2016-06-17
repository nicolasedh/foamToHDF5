/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkH5FoamReader.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkH5FoamReader - read ASCII or binary stereo lithography files
// .SECTION Description
// vtkH5FoamReader is a source object that reads ASCII or binary stereo
// lithography files (.stl files). The FileName must be specified to
// vtkH5FoamReader. The object automatically detects whether the file is
// ASCII or binary.
//
// .stl files are quite inefficient since they duplicate vertex
// definitions. By setting the Merging boolean you can control whether the
// point data is merged after reading. Merging is performed by default,
// however, merging requires a large amount of temporary storage since a
// 3D hash table must be constructed.

// .SECTION Caveats
// Binary files written on one system may not be readable on other systems.
// vtkSTLWriter uses VAX or PC byte ordering and swaps bytes on other systems.

#ifndef vtkH5FoamReader_h
#define vtkH5FoamReader_h

#include "vtkIOGeometryModule.h" // For export macro
//#include "vtkAbstractPolyDataReader.h"
#include "vtkMultiBlockDataSetAlgorithm.h"
#include "vtkSmartPointer.h"
#include "vtkVariant.h"
#include "vtkVariantArray.h"
#include <map>
#include <vector>
#include "hdf5.h"

class vtkCellArray;
class vtkFloatArray;
class vtkIncrementalPointLocator;
class vtkPoints;
class vtkDoubleArray;
class vtkIntArray;
class vtkUnstructuredGrid;

typedef std::vector<std::string> stringVec;
typedef std::map<std::string,stringVec> componentMap;


class VTKIOGEOMETRY_EXPORT vtkH5FoamReader : public vtkMultiBlockDataSetAlgorithm
{
public:
  vtkTypeMacro(vtkH5FoamReader,vtkMultiBlockDataSetAlgorithm);
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct object 
  static vtkH5FoamReader *New();



protected:
  vtkH5FoamReader();
  ~vtkH5FoamReader();
  char * FileName;

  int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **vtkNotUsed(inputVector), vtkInformationVector *outputVector);
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
  //Data
  hid_t file_id;
  bool updateMesh;
//  hid_t ncells,npoints,nfaces;
  vtkSmartPointer<vtkDoubleArray> times;
  vtkSmartPointer<vtkVariantArray> timesNames;
  vtkSmartPointer<vtkUnstructuredGrid> ugrid;
  vtkMultiBlockDataSet * output;


  //Member functions
  void openFile();
  void closeFile();
  void readMesh();
  void readProcMesh(char * location, int procN);
  void readPoints(const hid_t gid_proc, const char *name, vtkPoints *points);
  //can currently only read vtkIntArray or
  void readArray(const hid_t gid_proc, const char *name, vtkIntArray * ids);
  void readDoubleArray(const hid_t gid_proc, const char *name, vtkDoubleArray * data,int *ncols);
  void readFields(vtkVariant timeName);
  void readProcFields(char * location,int procId);
  void findProcs(const char * location, vtkVariantArray *procNames);
  vtkIdType findClosestSolutionIndex(double Time);
  int findTimes();
  vtkH5FoamReader(const vtkH5FoamReader&);  // Not implemented.
  void operator=(const vtkH5FoamReader&);  // Not implemented.
};

#endif
