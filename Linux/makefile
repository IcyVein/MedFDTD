objects = MedFDTD.o computeFraction.o computeMassAveragedSAR.o \
		 computeType1Voxel.o computeType2Voxel.o               \
		 computeType3Voxel.o

libs = -lmpi

MedFDTD : $(objects)
		g++ -o MedFDTD.exe $(objects) -L./ $(libs)

MedFDTD.o : 1D_FDTD.H 3D_FDTD_DEFINES.h BUILDOBJECTS.H         \
			COMPUTE.H EXTENSIONS.H INITIALIZE_DATA_AND_FILES.H \
			POWERSOURCE.H SETUP.H WRITEFIELD.H                 \
			computeMassAveragedSAR.h mpi.h mpio.h

.PHONY : clean
clean : 
	-rm $(objects)

