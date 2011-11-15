//Version : 1.01
//@Author: wangmeng, yudong
//		   Britton Chance Center for Biomedical Photonics
//		   Wuhan National Laboratory for Optoelectronics
//		   Huazhong University of Science and Technology (HUST)
//Last modified: 2011.9.28.

/*
 * Header files (Libraries to be included)
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>

using namespace std;

#define MPICH_SKIP_MPICXX
#include "mpi.h"
#pragma comment (lib,"mpi.lib")

#include "3D_FDTD_DEFINES.h"
#include "1D_FDTD.h"
#include "INITIALIZE_DATA_AND_FILES.h"
#include "SETUP.h"
#include "BUILDOBJECTS.h"
#include "POWERSOURCE.h"
#include "COMPUTE.h"
#include "WRITEFIELD.h"

#include "computeMassAveragedSAR.h"

#include "EXTENSIONS.h"

int main(int argc, char* argv[])
{
	char project_path_name[MAX_SIZE_OF_PATH];
	strcpy(project_path_name, path_proj);
	strcat(strcat(project_path_name,proj_name),".txt");
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	if (myrank == 0)
	{
		if (newProject(project_path_name))
		{
			printf("Creat project : %s fail.\n", project_path_name);
			fflush(stdout);
			return (FILE_ERROR);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (openProject(project_path_name))
	{
		printf("Open project : %s fail.\n", project_path_name);
		fflush(stdout);
		return(FILE_ERROR);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	initializeFile();

	MPI_Barrier(MPI_COMM_WORLD);
	initializePart1();

	MPI_Barrier(MPI_COMM_WORLD);
	setUp();

	MPI_Barrier(MPI_COMM_WORLD);
	setUpCPML();

	MPI_Barrier(MPI_COMM_WORLD);
	initializePart2();

	MPI_Barrier(MPI_COMM_WORLD);
	if (nprocs > 1)/* parallel */
	{
		if (myrank == 0)
		{
			printf("In Parallel.\n");
			fflush(stdout);
		}
		if (abcNo == 2)
		{
			printf("The absorbing boundary No. %d is invalid, use PML now.\n", abcNo);
		}
		compute();
	}
	else if (abcNo == 1)
		computeOneCPU();
	else if (abcNo == 2)
		computeOneCPU_Mur2();
	else
	{
		printf("The absorbing boundary No. %d is not found, use PML now.\n", abcNo);
		computeOneCPU();
	}

	printf("CPU%d/%d is quit successfully.\n", myrank, nprocs);
	fflush(stdout);
	MPI_Finalize();

	return (SUCCESS);
}

/*****************************************************************************************/
/*
 * Function name: openProject
 * Description: open a project file
 * Parameters: 
 *			 proj_path
 * Return: 
 *			 1				Successful
 *			 0				fail
 */
int openProject(char* proj_path)
{
	int i;
	FILE* fp_proj;
	if (myrank == 0)
	{
		printf("Loading Project : %s\n", proj_path);
		fflush(stdout);
	}
	if (fp_proj = fopen(proj_path, "r"))
	{
		fscanf(fp_proj, "nMax=%d\n", &nMax);

		fscanf(fp_proj, "_spaceX=%d,_spaceY=%d,_spaceZ=%d\n", &_spaceX, &_spaceY, &_spaceZ);
		fscanf(fp_proj, "abcNo=%d\n", &abcNo);
		fscanf(fp_proj, "thicknessOfPml=%d\n", &thicknessOfPml);
		fscanf(fp_proj, "padding=%d\n", &padding);
		fscanf(fp_proj, "mediaNum=%d\n", &mediaNum);

		fscanf(fp_proj, "_isource=%d,_jsource=%d,_ksource=%d\n", &_isource, &_jsource, &_ksource);
		fscanf(fp_proj, "port=%c\n", &port);
		fscanf(fp_proj, "sourceType=%d\n", &sourceType);
		fscanf(fp_proj, "waveForm=%d\n", &waveForm);
		fscanf(fp_proj, "amp=%lf\n", &amp);
		fscanf(fp_proj, "freq=%lf\n", &freq);
		fscanf(fp_proj, "t0=%d\n", &t0);
		fscanf(fp_proj, "pulse_width=%d\n", &pulse_width);

		fscanf(fp_proj, "#BEGIN# Field Save\n");
		fscanf(fp_proj, "save_plane_amount=%d\n", &save_plane_amount);
		fp_save_field_file = (field_file*) malloc (save_plane_amount * sizeof(field_file));
		for (i = 0; i < save_plane_amount; i++)
		{
			fscanf(fp_proj, "saveStart=%d,", &fp_save_field_file[i].sp.start);
			fscanf(fp_proj, "saveEnd=%d,", &fp_save_field_file[i].sp.end);
			fscanf(fp_proj, "saveStep=%d,", &fp_save_field_file[i].sp.step);
			fscanf(fp_proj, "savePlaneNo=%d,", &fp_save_field_file[i].sp.plane_no);
			fscanf(fp_proj, "slice=%d\n", &fp_save_field_file[i].sp.slice);
		}
		fscanf(fp_proj, "#END# Field Save\n");
#ifdef _SAR
		fscanf(fp_proj, "#BEGIN# SAR Save\n");

		fscanf(fp_proj, "#BEGIN# Mass Averaged SAR Save\n");
		fscanf(fp_proj, "Amount of Mass Averaged SAR %d\n", &nXgSAR);
		if (myrank)
		{
			XgSAR = (double*) calloc (nXgSAR, sizeof(double));
		}
		for (i = 0; i<nXgSAR; ++i)
			fscanf(fp_proj, "Mass Averaged SAR %lf\n", &XgSAR[i]);
		fscanf(fp_proj, "#END# Mass Averaged SAR Save\n");

		fscanf(fp_proj, "#BEGIN# LocalSAR Save\n");
		fscanf(fp_proj, "save_localSAR_amount=%d\n", &save_localSAR_amount);

		pSAR = (localSAR*) calloc (save_localSAR_amount+1, sizeof(localSAR));

		for (i = 0; i < save_localSAR_amount; i++)
		{
			fscanf(fp_proj, "saveLocalSARStart=%d,", &pSAR[i].start);
			fscanf(fp_proj, "saveLocalSAREnd=%d,", &pSAR[i].end);
			fscanf(fp_proj, "saveLocalSARPlaneNo=%d,", &pSAR[i].plane_no);
			fscanf(fp_proj, "LocalSARslice=%d\n", &pSAR[i].slice);
		}
		pSAR[i].start = pSAR[i-1].start;
		pSAR[i].end = pSAR[i-1].end;
		pSAR[i].plane_no = pSAR[i-1].plane_no;
		pSAR[i].slice = pSAR[i-1].slice+1;
		fscanf(fp_proj, "#END# LocalSAR Save\n");

		fscanf(fp_proj, "#END# SAR Save\n");
#endif
		fclose(fp_proj);

		if (abcNo == 0)
		{
			thicknessOfPml = 0;
		}
		else if (abcNo == 2 && nprocs == 1)
		{
			thicknessOfPml = 0;
		}
		paddingX_1 = padding;
		paddingX_2 = padding;
		paddingY_1 = padding;
		paddingY_2 = padding;
		paddingZ_1 = padding;
		paddingZ_2 = padding;

		spaceX = _spaceX + paddingX_1 + paddingX_2 + 1;
		spaceY = _spaceY + paddingY_1 + paddingY_2 + 1;
		spaceZ = _spaceZ + paddingZ_1 + paddingZ_2 + 1;

		Imax = spaceX + 2 * thicknessOfPml + 2;
		Jmax = spaceY + 2 * thicknessOfPml + 2;
		Kmax = spaceZ + 2 * thicknessOfPml + 2;

		nx_procs = (int*) malloc (nprocs * sizeof(int));
		for (i = 0; i<nprocs; ++i)
			nx_procs[i] = Imax/nprocs;
		for (i = nprocs - 1; i>=nprocs - Imax%nprocs; --i)
			nx_procs[i]++;

		for (i = 0; i<myrank; ++i)
			is += nx_procs[i];
		ie = is + nx_procs[myrank]-1;
		_global_is = is;
		_global_ie = ie;
		if (nprocs > 1)
		{
			Imax = nx_procs[myrank] + 1;

			cout<<"CPU"<<myrank<<": "<<is<<" --> "<<ie<<endl;
		}

		isource = _isource + paddingX_1 + thicknessOfPml;
		jsource = _jsource + paddingY_1 + thicknessOfPml;
		ksource = _ksource + paddingZ_1 + thicknessOfPml;

		for (i = 0; i < save_plane_amount; ++i)
			switch (fp_save_field_file[i].sp.plane_no)
			{
				case 1 :
					fp_save_field_file[i].sp.slice = fp_save_field_file[i].sp.slice + paddingX_1 + thicknessOfPml;
					break;
				case 2 :
					fp_save_field_file[i].sp.slice = fp_save_field_file[i].sp.slice + paddingY_1 + thicknessOfPml;
					break;
				case 3 :
					fp_save_field_file[i].sp.slice = fp_save_field_file[i].sp.slice + paddingZ_1 + thicknessOfPml;
					break;
			}
#ifdef _SAR
		for (i = 0; i<save_localSAR_amount+1; ++i)
		{
			pSAR[i].slice = pSAR[i].slice + paddingX_1 + thicknessOfPml;
		}
#endif
		return SUCCESS;
	}
	else
		return FILE_ERROR;
}

/*****************************************************************************************/
/*
 * Function name: buildObject
 * Description: Creat geometry(cube¡¢antenna)
 * Parameters:
 * Return:
 */
void buildObject()/* %% Can be modified : build model here */
{
	if (flag <13)
	{
		_DIPOLE_ANTENNA _dipole_antenna;/* %% Can be modified : define a antenna */
		_dipole_antenna.direction = 1;/* %% Can be modified */
		_dipole_antenna.feed_x = _isource;/* %% Can be modified : feed point coordination in x-axis */
		_dipole_antenna.feed_y = _jsource;/* %% Can be modified : feed point coordination in y-axis */
		_dipole_antenna.feed_z = _ksource;/* %% Can be modified : feed point coordination in z-axis */
		_dipole_antenna.impedance = 0;/* %% Can be modified */
		_dipole_antenna.length_high = 11;/* %% Can be modified */
		_dipole_antenna.length_low  = 11;/* %% Can be modified */
//		creatDipoleAntenna(_dipole_antenna);/* %% Can be modified : create the antenna "_dipole_antenna" */
	}
	/* flag >=13 Will be executed before compute SAR */
//	buildPlane(5, 10, 15, 5, 10+10, 15+10, 0);/* %% Can be modified : a 2D plane, coordinate range from 1 to _spaceX,Y,Z , MKS units*/
	buildCuboid(24, 5, 20, 7, 21, 21, createNewMedia(maxMedia, 0.1, 81, 1e6));/* %% Can be modified : a 3D cuboid, coordinate range from 1 to _spaceX,Y,Z , MKS units*/
//	buildCuboid(31, 5, 20, 10, 21, 21, createNewMedia(maxMedia, 2, 200, 2e3));
}

/*****************************************************************************************/
/*
 * Function name: newProject
 * Description: creat a project file
 * Parameters:
 *			 proj_path
 * Return:
 *			 1				Successful
 *			 0				fail
 */
int newProject(char* proj_path)
{
	int i;
	FILE* fp_proj;
	if (isNewProj == 0)
		return -1;
	printf("Building new project...\n");
	fflush(stdout);

/* Start */
	nMax=2000;/* %% Can be modified : Max time step */
	/* Size of grids(x, y, z) */
	_spaceX=100;/* %% Can be modified */
	_spaceY=50;/* %% Can be modified */
	_spaceZ=50;/* %% Can be modified */
	abcNo = 1;/* %% Can be modified : type of absorbing boundary : 0--PEC; 1--PML; 2--Mur2(only used in serial) */
	thicknessOfPml=7;/* %% Can be modified */
	padding=10;/* %% Can be modified : space from PML to computational grid */
	mediaNum = 0;/* %% Can be modified : number of media type in media file */

	_isource=0;/* %% Can be modified : Location of source in x-axis: range 1 to _spaceX */
	_jsource=1;/* %% Can be modified : Location of source in y-axis: range 1 to _spaceY*/
	_ksource=0;/* %% Can be modified : Location of source in z-axis: range 1 to _spaceZ*/
	port='x';/* %% Can be modified : voltage direction, 'x'--Ex, 'y'--Ey, 'z'--Ez */
	sourceType=1;/* %% Can be modified : set source type, 0--Point source, 1--Plane source, 2--Dipole source */
	/* Wave form: -1--User, 0--Sine : freq; 1--Gauss : t0, pulse_width; 2--raised cosine : pulse_width;
	 *		       3--differential Gauss : t0, pulse_width; 4-- 3-truncated cosine : pulse_width; 
	 *		       5--modulation Gauss : freq, t0, pulse_width.
	 */
	waveForm=0;/* %% Can be modified */
	amp=1000;/* %% Can be modified */
	freq=1.8e9;/* %% Can be modified */
	t0 = 150;/* %% Can be modified */
	pulse_width = 120;/* %% Can be modified */

	save_plane_amount = 3;/* %% Can be modified : Amount of plane need to save */

	fp_save_field_file = (field_file*) malloc (save_plane_amount * sizeof(field_file));/* Defaule setting: save plane in the plane where have feed point */
	int temp[3];
	if (sourceType == 1)
	{
		temp[0] = _spaceZ/2;
		temp[1] = _spaceX/2;
		temp[2] = _spaceY/2;
	}
	else
	{
		temp[0] = _ksource;
		temp[1] = _isource;
		temp[2] = _jsource;
	}
	int temp_plane = 1;
	for (i = 0; i < save_plane_amount; i++)
	{
		fp_save_field_file[i].sp.start=1;/* %% Can be modified */
		fp_save_field_file[i].sp.end=nMax;/* %% Can be modified */
		fp_save_field_file[i].sp.step=1;/* %% Can be modified */
		fp_save_field_file[i].sp.plane_no=temp_plane;/* %% Can be modified : set stored plane, 1--xy, 2--yz, 3--xz*/

		fp_save_field_file[i].sp.slice=temp[temp_plane-1];/* %% Can be modified */
		if (temp_plane<3) temp_plane++;
	}

#ifdef _SAR
	nXgSAR = 2;/* %% Can be modified */
	XgSAR = (double*) calloc (nXgSAR, sizeof(double));
	/* %% Can be modified */
	XgSAR[0] = 1.0;
	XgSAR[1] = 10.0;

	if (nXgSAR)
		save_localSAR_amount = _spaceZ;
	else
		save_localSAR_amount = 0;/* %% Can be modified : if you need stored local SAR */

	pSAR = (localSAR*) calloc (save_localSAR_amount, sizeof(localSAR));

	dt = 1.0 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));

	for (i = 0; i<save_localSAR_amount; ++i)
	{
		pSAR[i].start = nMax-(int)(1/freq/dt)+1;/* %% Can be modified : a cycle */
		pSAR[i].end = nMax;/* %% Can be modified */
		pSAR[i].plane_no = 1;/* %% Can be modified */
		pSAR[i].slice = _ksource+i;/* %% Can be modified */
	}
	if (nXgSAR)
		for (i = 0; i<save_localSAR_amount; ++i)
			pSAR[i].slice = 1+i;
#endif

/*End*/

	if (fp_proj = fopen(proj_path, "w+"))/* Write a new project file */
	{
		fprintf(fp_proj, "nMax=%d\n", nMax);

		fprintf(fp_proj, "_spaceX=%d,_spaceY=%d,_spaceZ=%d\n", _spaceX, _spaceY, _spaceZ);
		fprintf(fp_proj, "abcNo=%d\n", abcNo);
		fprintf(fp_proj, "thicknessOfPml=%d\n", thicknessOfPml);
		fprintf(fp_proj, "padding=%d\n", padding);
		fprintf(fp_proj, "mediaNum=%d\n", mediaNum);

		fprintf(fp_proj, "_isource=%d,_jsource=%d,_ksource=%d\n", _isource, _jsource, _ksource);
		fprintf(fp_proj, "port=%c\n", port);
		fprintf(fp_proj, "sourceType=%d\n", sourceType);
		fprintf(fp_proj, "waveForm=%d\n", waveForm);
		fprintf(fp_proj, "amp=%lf\n", amp);
		fprintf(fp_proj, "freq=%lf\n", freq);
		fprintf(fp_proj, "t0=%d\n", t0);
		fprintf(fp_proj, "pulse_width=%d\n", pulse_width);

		fprintf(fp_proj, "#BEGIN# Field Save\n");
		fprintf(fp_proj, "save_plane_amount=%d\n", save_plane_amount);
		for (i = 0; i < save_plane_amount; i++)
		{
			fprintf(fp_proj, "saveStart=%d,", fp_save_field_file[i].sp.start);
			fprintf(fp_proj, "saveEnd=%d,", fp_save_field_file[i].sp.end);
			fprintf(fp_proj, "saveStep=%d,", fp_save_field_file[i].sp.step);
			fprintf(fp_proj, "savePlaneNo=%d,", fp_save_field_file[i].sp.plane_no);
			fprintf(fp_proj, "slice=%d\n", fp_save_field_file[i].sp.slice);
		}
		fprintf(fp_proj, "#END# Field Save\n");
#ifdef _SAR
		fprintf(fp_proj, "#BEGIN# SAR Save\n");

		fprintf(fp_proj, "#BEGIN# Mass Averaged SAR Save\n");
		fprintf(fp_proj, "Amount of Mass Averaged SAR %d\n", nXgSAR);
		for (i = 0; i<nXgSAR; ++i)
			fprintf(fp_proj, "Mass Averaged SAR %lf\n", XgSAR[i]);
		fprintf(fp_proj, "#END# Mass Averaged SAR Save\n");

		fprintf(fp_proj, "#BEGIN# LocalSAR Save\n");
		fprintf(fp_proj, "save_localSAR_amount=%d\n", save_localSAR_amount);
		for (i = 0; i < save_localSAR_amount; ++i)
		{
			fprintf(fp_proj, "saveLocalSARStart=%d,", pSAR[i].start);
			fprintf(fp_proj, "saveLocalSAREnd=%d,", pSAR[i].end);
			fprintf(fp_proj, "saveLocalSARPlaneNo=%d,", pSAR[i].plane_no);
			fprintf(fp_proj, "LocalSARslice=%d\n", pSAR[i].slice);
		}
		fprintf(fp_proj, "#END# LocalSAR Save\n");

		fprintf(fp_proj, "#END# SAR Save\n");
#endif

		fclose(fp_proj);
		return SUCCESS;
	}
	else
		return FILE_ERROR;
}

