//Version : 1.1
//@Author: wangmeng, yudong
//		   Britton Chance Center for Biomedical Photonics
//		   Wuhan National Laboratory for Optoelectronics
//		   Huazhong University of Science and Technology (HUST)
//Last modified: 2011.12.20.

/*
 * Header files (Libraries to be included)
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>

using namespace std;

#define MPICH_SKIP_MPICXX
#include "mpi.h"
#pragma comment (lib,"mpi.lib")

#include "3D_FDTD_DEFINES.H"
#include "INITIALIZE_DATA_AND_FILES.H"
#include "SETUP.H"
#include "BUILDOBJECTS.H"
#include "POWERSOURCE.H"
#include "COMPUTE.H"
#include "dispersion.h"
#include "WRITEFIELD.H"

#include "MASS_AVERAGE_SAR.H"

#include "EXTENSIONS.H"

#include "process.h"

int main(int argc, char** argv)
{
	int return_data = 0;
	char project_path_name[MAX_SIZE_OF_PATH];
	char config_path[MAX_SIZE_OF_PATH];
    time_t date;
    time(&date);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    cout<<"#"<<_getpid()<<endl;
    fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	if (argc == 1)
		strcpy(config_path, ".\\config.txt");
	else
		strcpy(config_path, argv[1]);
    FILE* fp_proj = fopen(config_path, "r");
    if (argc>2)
        strcpy(proj_name, argv[2]);
    else
        strcpy(proj_name, "output");
    if (myrank == 0)
    {
        printf("Loading setting file: %s\n", config_path);
    }
    fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	if (fp_proj == NULL)
	{
        if (myrank == 0)
        {
		    printf("fail, exit.\n");
        }
		fflush(stdout);
        fp_log = fopen(".\\MedFDTD_Error.log", "a");
		if (fp_log == NULL)
			fp_log = fopen(".\\MedFDTD_Error.log", "w");
        fprintf(fp_log, "MedFDTD can't find the configure file in %s, please check you input\n", config_path);
        fclose(fp_log);
		return(FILE_ERROR);
	}
    else
    {
        openProject(fp_proj);
		if (myrank == 0)
			printf("\n");
		fflush(stdout);
	    MPI_Barrier(MPI_COMM_WORLD);
    }
    if (myrank == 0)
    {
        fp_log = fopen(path_log, "a");
		if (fp_log == NULL)
			fp_log = fopen(path_log, "w");
        if (fp_log == NULL)
            printf("Create log file %s fail.\n", path_log);
        else
            printf("Log file : %s\n", path_log);
        printf("--------------------------\n");
        printf("Start MedFDTD at %s", ctime(&date));
        printf("--------------------------\n");
        fprintf(fp_log, "----------------------------------------------------\n");
        fprintf(fp_log, "Start MedFDTD at %s", ctime(&date));
        fprintf(fp_log, "----------------------------------------------------\n");
        fflush(stdout);
    }

	MPI_Barrier(MPI_COMM_WORLD);
	initializeFile();

	MPI_Barrier(MPI_COMM_WORLD);
	initializePart1();

	MPI_Barrier(MPI_COMM_WORLD);
	return_data = setUp();
	if(return_data == MEM_ERROR)
	{
		MPI_Finalize();
		exit(MEM_ERROR);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	setUpCPML();

	MPI_Barrier(MPI_COMM_WORLD);
	initializePart2();
    if (media_dispersion == 1)
        initializeDispersion();

	MPI_Barrier(MPI_COMM_WORLD);

	if (nprocs > 1)/* parallel */
	{
		if (myrank == 0)
		{
            printf("Total Memory requirements = %ld KB\n", mem_count / 1024);
			printf("In Parallel.\n");
			fprintf(fp_log, "In Parallel.\n");
            fprintf(fp_log, "Total Memory requirements = %ld KB\n", mem_count / 1024);
			fflush(stdout);
		}
		if (abcNo == 2 && myrank == 0)
		{
			printf("The absorbing boundary No. %d is invalid, use PML now.\n", abcNo);
			fprintf(fp_log, "The absorbing boundary No. %d is invalid, use PML now.\n", abcNo);
			fflush(stdout);
		}
        if (media_dispersion == 1)
            computeDispersion();
        else
		    compute();
	}
	else if (abcNo == 1)
		computeOneCPU();
	else if (abcNo == 2)
		computeOneCPU_Mur2();
	else
	{
        if (myrank == 0)
        {
		    printf("The absorbing boundary No. %d is not found, use PML now.\n", abcNo);
		    fprintf(fp_log, "The absorbing boundary No. %d is not found, use PML now.\n", abcNo);
			fflush(stdout);
        }
		computeOneCPU();
	}

    if (need_compute_temperature_rise)
    {
    	MPI_Barrier(MPI_COMM_WORLD);
        computeTemperatureRise(path_sar);
	    MPI_Barrier(MPI_COMM_WORLD);
    }

    if (myrank == 0)
    {
        time(&date);
        printf("--------------------------\n");
        printf("Complete at %s", ctime(&date));
        printf("--------------------------\n");
        fprintf(fp_log, "----------------------------------------------------\n");
        fprintf(fp_log, "Complete at %s", ctime(&date));
        fprintf(fp_log, "----------------------------------------------------\n");
        fflush(stdout);
        fclose(fp_log);
    }

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
int openProject(FILE* fp)
{
	int i;
	if (fp)
	{
        fscanf(fp, "<Path>\n");
        fscanf(fp, "path_save=%s\n", path_save);
        //fscanf(fp, "path_log=%s\n", path_log);
        sprintf(path_log, "%s", path_save);
        fscanf(fp, "<Time>\n");
        fscanf(fp, "nMax=%d\n", &nMax);
        fscanf(fp, "Convergence=%f\n", &convergenceTarget);
        fscanf(fp, "dt=%lf\n", &dt);

        fscanf(fp, "<Mesh>\n");
		fscanf(fp, "_spaceX=%d,_spaceY=%d,_spaceZ=%d\n", &_spaceX, &_spaceY, &_spaceZ);
        fscanf(fp, "dx=%lf,dy=%lf,dz=%lf\n", &dx, &dy, &dz);
        fscanf(fp, "padding=%d\n", &padding);
        if (dt == 0)
        {
			dt = 0.9999 / (C * sqrt(1.0/(dx * dx) + 1.0/(dy * dy)+ 1.0/(dz * dz)));
        }

        fscanf(fp, "<Absorbing boundary>\n");
		fscanf(fp, "abcNo=%d\n", &abcNo);
		fscanf(fp, "thicknessOfPml=%d\n", &thicknessOfPml);

        fscanf(fp, "<Power source>\n");
        fscanf(fp, "sourceType=%d\n", &sourceType);
        fscanf(fp, "sourceNum=%d\n", &sourceNum);
        for (i = 0; i<sourceNum; ++i)
        {
            int iSrc, jSrc, kSrc;
            fscanf(fp, "%d,%d,%d\n", &iSrc, &jSrc, &kSrc);
            _isource[i] = iSrc;
            _jsource[i] = jSrc;
            _ksource[i] = kSrc;
        }
		fscanf(fp, "port=%c\n", &port);

		fscanf(fp, "waveForm=%d\n", &waveForm);
		fscanf(fp, "amp=%lf\n", &amp);
        fscanf(fp, "freq=%lf\n", &freq);
	    fscanf(fp, "t0=%d\n", &t0);
	    fscanf(fp, "pulse_width=%d\n", &pulse_width);
        fscanf(fp, "pathSRC=%s\n", pathSRC);

        fscanf(fp, "<Model>\n");
        fscanf(fp, "-<Import model>\n");
        fscanf(fp, "path_data=%s\n", path_data);
        if (strcmp(path_data, "NULL\n"))
        {
            memset(path_data, 0, sizeof(path_data));
            strcpy(path_data, path_save);
        }
        fscanf(fp, "model_name=%s\n", model_name);
        fscanf(fp, "media_name=%s\n", media_name);
        fscanf(fp, "media_dispersion=%d\n", &media_dispersion);
		fscanf(fp, "mediaNum=%d\n", &mediaNum);
        fscanf(fp, "-<Build object>\n");
        fscanf(fp, "object_num=%d\n", &object_num);
        object_creat = (int**) calloc (object_num, sizeof(int*));
        for (i = 0; i<object_num; ++i)
        {
            object_creat[i] = (int*) calloc (8, sizeof(int));
        }
        double temp_media[3] = {0.0};
        for (i = 0; i<object_num; ++i)
        {
            fscanf(fp, "%d,%d,%d,%d,%d,%d,%d,%lf,%lf,%lf\n", &object_creat[i][0], 
                    &object_creat[i][1], &object_creat[i][2], &object_creat[i][3], 
                    &object_creat[i][4], &object_creat[i][5], &object_creat[i][6], 
                    &temp_media[0], &temp_media[1], &temp_media[2]);
            object_creat[i][7] = createNewMedia(maxMedia, temp_media[0], temp_media[1], temp_media[2]);
        }
        fscanf(fp, "-<Antenna>\n");
        fscanf(fp, "antenna_amount=%d\n", &antenna_amount);
        _dipole_antenna = (_DIPOLE_ANTENNA*) calloc (antenna_amount, sizeof(_DIPOLE_ANTENNA));
        for (i = 0; i<antenna_amount; ++i)
        {
            fscanf(fp, "antenna_direction=%d\n", &_dipole_antenna[i].direction);
            fscanf(fp, "antenna_feed_x=%d,antenna_feed_y=%d,antenna_feed_z=%d\n", &_dipole_antenna[i].feed_x, &_dipole_antenna[i].feed_y, &_dipole_antenna[i].feed_z);
            fscanf(fp, "antenna_impedance=%lf\n", &_dipole_antenna[i].impedance);
            fscanf(fp, "antenna_length_high=%d,antenna_length_low=%d\n", &_dipole_antenna[i].length_high, &_dipole_antenna[i].length_low);
        }
		fscanf(fp, "<Field Save>\n");
		fscanf(fp, "save_plane_amount=%d\n", &save_plane_amount);
		fp_save_field_file = (field_file*) calloc (save_plane_amount, sizeof(field_file));
		for (i = 0; i < save_plane_amount; i++)
		{
			fscanf(fp, "%d,", &fp_save_field_file[i].sp.start);
			fscanf(fp, "%d,", &fp_save_field_file[i].sp.end);
			fscanf(fp, "%d,", &fp_save_field_file[i].sp.step);
			fscanf(fp, "%d,", &fp_save_field_file[i].sp.plane_no);
			fscanf(fp, "%d\n", &fp_save_field_file[i].sp.slice);
		}
#ifdef _SAR
		fscanf(fp, "<SAR>\n");

		fscanf(fp, "-<Whole body SAR>\n");
		fscanf(fp, "Whole body SAR %d\n", &whole_body_sar);

		fscanf(fp, "-<Mass Averaged SAR>\n");
		fscanf(fp, "%d\n", &nXgSAR);
        sourceT = (int)(1.0/freq/dt)+1;
        if (sourceT<=0 || waveForm!=0)
            sourceT = nMax/100;
        if (convergenceTarget == 0)
            sourceT = nMax;
		XgSAR = (double*) calloc (nXgSAR, sizeof(double));
		for (i = 0; i<nXgSAR; ++i)
        {
			fscanf(fp, "%lf\n", &XgSAR[i]);
            XgSAR[i] = XgSAR[i]/1000; // g to kg
        }

        if (nXgSAR || whole_body_sar)
        {
            save_localSAR_amount = _spaceZ;
            pSAR = (localSAR*) calloc (save_localSAR_amount+1, sizeof(localSAR));

		    for (i = 0; i < save_localSAR_amount+1; i++)
    		{
				pSAR[i].length = sourceT;
				pSAR[i].plane_no = 1;
				pSAR[i].slice = i+1;
    		}
			fscanf(fp, "-<LocalSAR>\n");
			fscanf(fp, "save_localSAR_amount=0\n");
        }
        else
        {
    		fscanf(fp, "-<LocalSAR>\n");
	    	fscanf(fp, "save_localSAR_amount=%d\n", &save_localSAR_amount);

            if (convergenceTarget == 0)
		        pSAR = (localSAR*) calloc (save_localSAR_amount+1, sizeof(localSAR));
            else
		        pSAR = (localSAR*) calloc (_spaceZ+1, sizeof(localSAR));

		    for (i = 0; i < save_localSAR_amount; i++)
    		{
    			fscanf(fp, "saveLocalSARLength=%d,", &pSAR[i].length);
    			fscanf(fp, "saveLocalSARPlaneNo=%d,", &pSAR[i].plane_no);
    			fscanf(fp, "LocalSARslice=%d\n", &pSAR[i].slice);
    		}
    		pSAR[i].length = pSAR[i-1].length;
    		pSAR[i].plane_no = pSAR[i-1].plane_no;
    		pSAR[i].slice = pSAR[i-1].slice+1;
        }
#endif
        fscanf(fp, "<Temperature Rise>\n");
        fscanf(fp, "need compute temperature rise %d\n", &need_compute_temperature_rise);
        if (need_compute_temperature_rise)
        {
            fscanf(fp, "sar_file=%s\n", path_sar);
            fscanf(fp, "ds=%lf\n", &ds);
            fscanf(fp, "dt_T=%lf\n", &dt_T);
            fscanf(fp, "max time=%d\n", &nMaxT);
        }
		fclose(fp);

        sprintf(path_log, "%s%s.log", path_log, proj_name);

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

        //PML Layers(With PEC)
	    nxPML_1 = thicknessOfPml + 1;
	    nxPML_2 = thicknessOfPml + 1;
	    nyPML_1 = thicknessOfPml + 1;
	    nyPML_2 = thicknessOfPml + 1;
	    nzPML_1 = thicknessOfPml + 1;
    	nzPML_2 = thicknessOfPml + 1;

		nx_procs = (int*) malloc (nprocs * sizeof(int));
		for (i = 0; i<nprocs; ++i)
			nx_procs[i] = Imax/nprocs;
		for (i = nprocs - 1; i>=nprocs - Imax%nprocs; --i)
			nx_procs[i]++;

        int loadBalancingCell = 4;// 调整包含PML计算区域的负载
        int loadBalancingCellAll = 2*loadBalancingCell;
        int loadBalancingThreshold = 4;
        int computationalDomainEdge = thicknessOfPml+1 + padding +1;
        planeWaveIndex = 2;

		if (nprocs>=loadBalancingThreshold)
		{
			nx_procs[0] -= loadBalancingCell;
			nx_procs[nprocs-1] -= loadBalancingCell;
			for (i = 1; i<nprocs-1; ++i)
				nx_procs[i] += loadBalancingCellAll/(nprocs-2);
			for (i = 1; i<=loadBalancingCellAll%(nprocs-2); ++i)
				nx_procs[i]++;
		}

        loadBalancingCellAll = 0;
        if (nx_procs[0] < computationalDomainEdge+planeWaveIndex)//是CPU0和CPU_End都至少含有1个总场计算区域
        {
            loadBalancingCellAll += computationalDomainEdge+planeWaveIndex-nx_procs[0];
            nx_procs[0] = computationalDomainEdge+planeWaveIndex;
        }
        if (nx_procs[nprocs-1] < computationalDomainEdge)
        {
            nx_procs[nprocs-1] = computationalDomainEdge;
            loadBalancingCellAll += computationalDomainEdge-nx_procs[nprocs-1];
        }
        if (loadBalancingCellAll)
        {
            for (i = 1; i<nprocs-1; ++i)
		    	nx_procs[i] -= loadBalancingCellAll/(nprocs-2);
		    for (i = 1; i<=loadBalancingCellAll%(nprocs-2); ++i)
		    	nx_procs[i]--;
        }

		for (i = 0; i<myrank; ++i)
			is += nx_procs[i];
		ie = is + nx_procs[myrank]-1;
		_global_is = is;
		_global_ie = ie;

		if (nprocs > 1)
		{
			Imax = nx_procs[myrank] + 1;
		}

        setOutputRange();

        for (i = 0; i<sourceNum; ++i)
        {
		    isource[i] = _isource[i] + paddingX_1 + thicknessOfPml;
		    jsource[i] = _jsource[i] + paddingY_1 + thicknessOfPml;
		    ksource[i] = _ksource[i] + paddingZ_1 + thicknessOfPml;
        }

        int sourceNumTotal = sourceNum;
        sourceNum = 0;

        for (int nSrc = 0; nSrc<sourceNumTotal; ++nSrc) 
        {
	        if (is <= isource[nSrc] && isource[nSrc] <= ie)
	        {
	        	if (myrank != 0)
                {
	        		isource[sourceNum++] = isource[nSrc] - is + 1;
                    jsource[sourceNum-1] = jsource[nSrc];
                    ksource[sourceNum-1] = ksource[nSrc];
                }
                else 
                    sourceNum++;
	        }
	        else if (isource[nSrc] == ie+1 && (port == 'y' || port == 'z'))
	        {
	        	if (myrank != 0)
                {
	        		isource[sourceNum++] = isource[nSrc] - is + 1;
                    jsource[sourceNum-1] = jsource[nSrc];
                    ksource[sourceNum-1] = ksource[nSrc];
                }
                else 
                    sourceNum++;
	        }
	        else;
        }

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
 * Description: Creat geometry(cube、antenna)
 * Parameters:
 * Return:
 */
void buildObject()
{
    int i = 0;
	if (flag <13)
	{
        for (i = 0; i<antenna_amount; ++i)
        {
	        creatDipoleAntenna(_dipole_antenna[i]);
        }
	}
	/* flag >=13 Will be executed before compute SAR */
    for (i = 0; i<object_num; ++i)
    {
        switch(object_creat[i][0])
        {
            case 1:
                buildLine(object_creat[i][1], object_creat[i][2], object_creat[i][3],
                          object_creat[i][4], object_creat[i][5], object_creat[i][6],
                          object_creat[i][7]);
                break;
            case 2:
                buildPlane(object_creat[i][1], object_creat[i][2], object_creat[i][3],
                           object_creat[i][4], object_creat[i][5], object_creat[i][6],
                           object_creat[i][7]);
                break;
            case 3:
                buildCuboid(object_creat[i][1], object_creat[i][2], object_creat[i][3],
                            object_creat[i][4], object_creat[i][5], object_creat[i][6],
                            object_creat[i][7]);
                break;
        }
    }
}

void setOutputRange() {// 存储范围是半开半闭区间：[beg, end)
    if (myrank == 0) {
        outputRangeE[0] = nxPML_1;
        outputRangeE[1] = Imax - 1;
    }
    else if (myrank == nprocs - 1) {
        outputRangeE[0] = 1;
        outputRangeE[1] = Imax - nxPML_2;
    }
    else {
        outputRangeE[0] = 1;
        outputRangeE[1] = Imax;
    }
    outputRangeE[2] = nyPML_1;
    outputRangeE[3] = Jmax - nyPML_2;
    outputRangeE[4] = nzPML_1;
    outputRangeE[5] = Kmax - nzPML_2;
}

