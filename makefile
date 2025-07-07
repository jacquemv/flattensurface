# the links more likely to break
LIB_CGAL = /usr/lib64/libSFCGAL.so.2
LIB_LAPACK = /usr/lib64/liblapack.so.3
INC_EIGEN3 = /usr/include/eigen3

CC = gcc -O3
CPP = g++ -O3
OPENMESH_TARGET = openmesh/dgpc
OPENMESH_LIBS = -lOpenMeshCore -lOpenMeshTools
CGAL_TARGET = cgal/parameterize
CGAL_OPTS = -DCGAL_EIGEN3_ENABLED -I$(INC_EIGEN3) $(LIB_CGAL) -frounding-math -lgmp
SPHERMAP_TARGET = sphermap/spharm_local_smoothing
SPHERMAP_LIBS = -lm $(LIB_LAPACK)
TARGETS = $(OPENMESH_TARGET) $(CGAL_TARGET) $(SPHERMAP_TARGET)

compile: $(TARGETS)

dependencies: install_cgal install_openmesh install_sphermap

clean:
	rm -f $(TARGETS)

# compile executables

$(OPENMESH_TARGET): openmesh/dgpc.cpp
	$(CPP) -o $(OPENMESH_TARGET) openmesh/dgpc.cpp $(OPENMESH_LIBS)

$(CGAL_TARGET): cgal/parameterize.cpp
	$(CPP) $(CGAL_OPTS) cgal/parameterize.cpp -o $(CGAL_TARGET)

$(SPHERMAP_TARGET): sphermap/spharm_local_smoothing.c
	$(CC) $(SPHERMAP_LIBS) sphermap/spharm_local_smoothing.c -o $(SPHERMAP_TARGET) 

# install dependencies on Fedora linux

install_cgal:
	sudo dnf install SFCGAL CGAL-devel eigen3-devel

install_openmesh:
	sudo dnf install OpenMesh OpenMesh-devel

install_sphermap:
	sudo dnf install octave lapack
