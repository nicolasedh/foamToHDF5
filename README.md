# foamToHDF5
functionObject that saves OpenFOAM results to a HDF5 file. It also includes a reader for ParaView.


###########################Install notes ################################

In order to compile one needs first to compile hdf5. I found it
quite a pain to compile parallel hdf5 together with paraview
so instead I compiled the reader with a custom build of hdf5 and
ignored the one that comes with vtk. 

############## Compile hdf5 ########################
#first compile hdf5, this is how I did it:
#download version 1.8.16 or git. Then

cd hdf5_git_src
#this is for git version only, as it doesn't include
#the configure script
./autogen.sh

CC=mpicc CXX=mpic++ ./configure \
   --disable-fortran \
   --disable-hl \
   --prefix=$WM_PROJECT_INST_DIR/installs/hdf5 \
   --enable-parallel

#incompatibel med parallel
#   --enable-cxx \

make -j 4
make install

############## compile the writer and reader #####
Just run the ./Allwmake script and it might work.
I've assumed that hdf5 was installed in the same
location as above.

Finaly to load the paraview plugin reader by
starting paraview and gooing to tools->manage plugins
Locate the libH5FoamReader.so in the lib dir of paraview
Load it and tick auto load on start up

