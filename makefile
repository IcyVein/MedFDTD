objects = MedFDTD.o MASS_AVERAGE_SAR.o WHOLE_BODY_AVERAGE_SAR.o

libs = -lmpi

MedFDTD : $(objects)
		g++ -o MedFDTD.exe $(objects) -L./ $(libs)

MedFDTD.o : 3D_FDTD_DEFINES.H BUILDOBJECTS.H	        	   \
			COMPUTE.H EXTENSIONS.H INITIALIZE_DATA_AND_FILES.H \
			POWERSOURCE.H  SETUP.H WRITEFIELD.H	dispersion.h   \
			MASS_AVERAGE_SAR.H mpi.h mpio.h

.PHONY : clean
clean : 
	-del $(objects)

