/*
 * File Description£ºExtensions function
 * Data:
 *				Ex(i,j,k) : Ex, i,j,k; Range:[nxPML_1,Imax-nxPML_2),[nyPML_1,Jmax-nyPML_2),[nzPML_1,Kmax-nzPML_2);
 *				Ey(i,j,k) : Ey, Range: Ditto;
 *				Ez(i,j,k) : Ez, Range: Ditto.
 *				epsilon,sigma : epsilon = loadData3(path_eps, size, 1.0);
 *											sigma = loadData3(path_sigma, size, 0.0);
 */
void addFunctions()
{
	/* codes */
}

#ifdef _SAR
localSAR* initializeLocalSAR(localSAR* pSAR, int nSAR)
{
	/* Initialize localSAR */
	int i, n;

	int save_length = 0;

	if (myrank == 0)
		save_length = nx_procs[myrank]-nxPML_1-paddingX_1;
	else if (myrank == nprocs - 1)
		save_length = nx_procs[myrank]-nxPML_2-paddingX_2;
	else
		save_length = nx_procs[myrank];

	for (n = 0; n<nSAR; ++n)
	{
		pSAR[n].RMS_ex = (double**) calloc (save_length, sizeof(double*));
		pSAR[n].RMS_ey = (double**) calloc (save_length+1, sizeof(double*));
		pSAR[n].RMS_ez = (double**) calloc (save_length+1, sizeof(double*));
		for (i = 0; i<save_length+1; ++i)
		{
			pSAR[n].RMS_ey[i] = (double*) calloc (_spaceY+1, sizeof(double));
			pSAR[n].RMS_ez[i] = (double*) calloc (_spaceY+1, sizeof(double));
		}
		for (i = 0; i<save_length; ++i)
		{
			pSAR[n].RMS_ex[i] = (double*) calloc (_spaceY+1, sizeof(double));
		}
	}

	return pSAR;
}

int freeLocalSARRMS(localSAR* pSAR, int nSAR)
{
	int i, n;

	int save_length = 0;

	if (myrank == 0)
		save_length = nx_procs[myrank]-nxPML_1-paddingX_1;
	else if (myrank == nprocs - 1)
		save_length = nx_procs[myrank]-nxPML_2-paddingX_2;
	else
		save_length = nx_procs[myrank];
	for (n = 0; n<nSAR; ++n)
	{
		for (i = 0; i<save_length+1; ++i)
		{
			free(pSAR[n].RMS_ey[i]);
			free(pSAR[n].RMS_ez[i]);
		}
		for (i = 0; i<save_length; ++i)
		{
			free(pSAR[n].RMS_ex[i]);
		}
	}
	return 0;
}
int freeLocalSARData(localSAR* pSAR, int nSAR)
{
	int i, n;

	int save_length = 0;

	if (myrank == 0)
		save_length = nx_procs[myrank]-nxPML_1-paddingX_1;
	else if (myrank == nprocs - 1)
		save_length = nx_procs[myrank]-nxPML_2-paddingX_2;
	else
		save_length = nx_procs[myrank];

	for (n = 0; n<nSAR; ++n)
	{
		for (i = 0; i<save_length; ++i)
		{
			free(pSAR[n].localSARData[i]);
		}
	}
	free(pSAR);
	return 0;
}

void computeRMS(localSAR* pSAR, int timestep)
{
	int i, j, k=0;

	if( (*pSAR).start > timestep || (*pSAR).end < timestep)	return;

	int save_length;
	if (myrank == 0)
		save_length = nx_procs[myrank]-nxPML_1-paddingX_1;
	else if (myrank == nprocs - 1)
		save_length = nx_procs[myrank]-nxPML_2-paddingX_2;
	else
		save_length = nx_procs[myrank];

	int change_by_pml = 0;
	if (myrank == 0)
		change_by_pml = nxPML_1+paddingX_1;
	else
		change_by_pml = is;

	switch ((*pSAR).plane_no)
	{
	case 1:
		for (i = 0; i<save_length+1; ++i)
		{
			for (j = 0; j<_spaceY+1; ++j)
			{
				(*pSAR).RMS_ey[i][j] += Ey(i+change_by_pml, j+nyPML_1+paddingY_1, (*pSAR).slice-1)*
									    Ey(i+change_by_pml, j+nyPML_1+paddingY_1, (*pSAR).slice-1);
				(*pSAR).RMS_ez[i][j] += Ez(i+change_by_pml, j+nyPML_1+paddingY_1, (*pSAR).slice)*
									    Ez(i+change_by_pml, j+nyPML_1+paddingY_1, (*pSAR).slice);
			}
		}

		for (i = 0; i<save_length; ++i)
		{
			for (j = 0; j<_spaceY+1; ++j)
			{
				(*pSAR).RMS_ex[i][j] += Ex(i+change_by_pml, j+nyPML_1+paddingY_1, (*pSAR).slice-1)*
									    Ex(i+change_by_pml, j+nyPML_1+paddingY_1, (*pSAR).slice-1);
			}
		}
		break;
	}

	if (timestep == (*pSAR).end)
	{
		for (i = 0; i<save_length+1; ++i)
		{
			for (j = 0; j<_spaceY+1; ++j)
			{
				(*pSAR).RMS_ey[i][j] /= ((*pSAR).end-(*pSAR).start+1);
				(*pSAR).RMS_ez[i][j] /= ((*pSAR).end-(*pSAR).start+1);
			}
		}
		for (i = 0; i<save_length; ++i)
		{
			for (j = 0; j<_spaceY+1; ++j)
			{
				(*pSAR).RMS_ex[i][j] /= ((*pSAR).end-(*pSAR).start+1);
			}
		}
	}
}

void computeLocalSAR(localSAR sSAR1, localSAR sSAR2, double*** localSARData)
{
	int i, j;

	int save_length;
	if (myrank == 0)
		save_length = nx_procs[myrank]-nxPML_1-paddingX_1;
	else if (myrank == nprocs - 1)
		save_length = nx_procs[myrank]-nxPML_2-paddingX_2;
	else
		save_length = nx_procs[myrank];

	int pml_change_X = 1, pml_change_Y = thicknessOfPml + paddingY_1 + 1;
	if (myrank == 0)
		pml_change_X = thicknessOfPml + paddingX_1 + 1;

	for (i = 0; i<save_length; ++i)
	{
		for (j = 0; j<_spaceY; ++j)
		{
			if (fabs(media[modelDataX[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1]].rho)       > 1e-6 &&
				fabs(media[modelDataX[i+pml_change_X][j+pml_change_Y +1][sSAR1.slice-1]].rho)    > 1e-6 &&
				fabs(media[modelDataX[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1 +1]].rho)    > 1e-6 &&
				fabs(media[modelDataX[i+pml_change_X][j+pml_change_Y +1][sSAR1.slice-1 +1]].rho) > 1e-6 &&
				fabs(media[modelDataY[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1]].rho)       > 1e-6 &&
				fabs(media[modelDataY[i+pml_change_X +1][j+pml_change_Y][sSAR1.slice-1]].rho)    > 1e-6 &&
				fabs(media[modelDataY[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1 +1]].rho)    > 1e-6 &&
				fabs(media[modelDataY[i+pml_change_X +1][j+pml_change_Y][sSAR1.slice-1 +1]].rho) > 1e-6 &&
				fabs(media[modelDataZ[i+pml_change_X][j+pml_change_Y][sSAR1.slice]].rho)         > 1e-6 &&
				fabs(media[modelDataZ[i+pml_change_X+1][j+pml_change_Y][sSAR1.slice]].rho)       > 1e-6 &&
				fabs(media[modelDataZ[i+pml_change_X][j+pml_change_Y+1][sSAR1.slice]].rho)       > 1e-6 &&
				fabs(media[modelDataZ[i+pml_change_X+1][j+pml_change_Y+1][sSAR1.slice]].rho)     > 1e-6   )
			{
				(*localSARData)[i][j] = (sSAR1.RMS_ex[i][j]    * media[modelDataX[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1]].sigma      / media[modelDataX[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1]].rho      + 
										 sSAR1.RMS_ex[i][j+1]  * media[modelDataX[i+pml_change_X][j+pml_change_Y +1][sSAR1.slice-1]].sigma   / media[modelDataX[i+pml_change_X][j+pml_change_Y +1][sSAR1.slice-1]].rho   + 
										 sSAR2.RMS_ex[i][j]    * media[modelDataX[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1 +1]].sigma   / media[modelDataX[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1 +1]].rho   + 
										 sSAR2.RMS_ex[i][j+1]  * media[modelDataX[i+pml_change_X][j+pml_change_Y +1][sSAR1.slice-1 +1]].sigma/ media[modelDataX[i+pml_change_X][j+pml_change_Y +1][sSAR1.slice-1 +1]].rho )/4 +
										(sSAR1.RMS_ey[i][j]    * media[modelDataY[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1]].sigma      / media[modelDataY[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1]].rho      + 
										 sSAR1.RMS_ey[i+1][j]  * media[modelDataY[i+pml_change_X +1][j+pml_change_Y][sSAR1.slice-1]].sigma   / media[modelDataY[i+pml_change_X +1][j+pml_change_Y][sSAR1.slice-1]].rho   + 
										 sSAR2.RMS_ey[i][j]    * media[modelDataY[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1 +1]].sigma   / media[modelDataY[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1 +1]].rho   + 
										 sSAR2.RMS_ey[i+1][j]  * media[modelDataY[i+pml_change_X +1][j+pml_change_Y][sSAR1.slice-1 +1]].sigma/ media[modelDataY[i+pml_change_X +1][j+pml_change_Y][sSAR1.slice-1 +1]].rho )/4 +
										(sSAR1.RMS_ez[i][j]    * media[modelDataZ[i+pml_change_X][j+pml_change_Y][sSAR1.slice]].sigma        / media[modelDataZ[i+pml_change_X][j+pml_change_Y][sSAR1.slice]].rho        + 
										 sSAR1.RMS_ez[i+1][j]  * media[modelDataZ[i+pml_change_X+1][j+pml_change_Y][sSAR1.slice]].sigma      / media[modelDataZ[i+pml_change_X+1][j+pml_change_Y][sSAR1.slice]].rho      + 
										 sSAR1.RMS_ez[i][j+1]  * media[modelDataZ[i+pml_change_X][j+pml_change_Y+1][sSAR1.slice]].sigma      / media[modelDataZ[i+pml_change_X][j+pml_change_Y+1][sSAR1.slice]].rho      + 
										 sSAR1.RMS_ez[i+1][j+1]* media[modelDataZ[i+pml_change_X+1][j+pml_change_Y+1][sSAR1.slice]].sigma    / media[modelDataZ[i+pml_change_X+1][j+pml_change_Y+1][sSAR1.slice]].rho     )/4;
			}
			else
			{
				(*localSARData)[i][j] = 0.0;
			}
		}
	}
}

void writeLocalSAR(FILE* fp, double** localSARData)
{
	int i, j;

	int save_length;
	if (myrank == 0)
		save_length = nx_procs[myrank]-nxPML_1-paddingX_1;
	else if (myrank == nprocs - 1)
		save_length = nx_procs[myrank]-nxPML_2-paddingX_2;
	else
		save_length = nx_procs[myrank];

	for (i = 0; i<save_length; ++i)
	{
		for (j = 0; j<_spaceY-1; ++j)
		{
			fprintf(fp, "%f ", localSARData[i][j]);
		}
		fprintf(fp, "%f\n", localSARData[i][j]);
	}
}

int computeXgSAR(int nXgSAR)
{
	int i, j, k;
	double*** localSARData = NULL;

	int spaceSize[3] = {_spaceX, _spaceY, _spaceZ};
	int resolution[3] = {1, 1, 1};

	if (myrank)
		spaceSize[0] = nx_procs[myrank];
	if (myrank == nprocs - 1)
		spaceSize[0] -= thicknessOfPml + paddingX_2 + 2; /* Bug: may cause <0 */

	char path_log_file[MAX_SIZE_OF_PATH];
	if (myrank == 0)
	{
		strcpy(path_log_file, path_log);
		strcat(path_log_file, proj_name);
		strcat(path_log_file, ".log");
		if (!(fp_log = fopen(path_log_file, "w+")))
		{
			printf("Load file %s fail.\n", path_log);
			fflush(stdout);
		}
	}

	freeLocalSARRMS(pSAR, _spaceZ);
	localSARData = loadLocalSAR(spaceSize);
	freeLocalSARData(pSAR, _spaceZ);

	double* tempData1D;
	double*** tempData3D;
	if (myrank)
	{
		tempData1D = trans3DTo1D(localSARData, spaceSize);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Request request_send;

	MPI_Request request_recv; 
	MPI_Status status_recv;

	if (myrank)
	{
		MPI_Isend(tempData1D, spaceSize[0]*spaceSize[1]*spaceSize[2], MPI_DOUBLE, 0, 1000+myrank, MPI_COMM_WORLD, &request_send);
	}
	else
	{
		printf("Sending Data to root CPU.\n");
		fflush(stdout);
		int is_change = -(thicknessOfPml + paddingX_1 + 1);
		for (int rank = 1; rank<nprocs; ++rank)
		{
			if (rank)
				spaceSize[0] = nx_procs[rank];
			if (rank == nprocs - 1)
				spaceSize[0] -= thicknessOfPml + paddingX_2 + 2;
			tempData1D = initArray(spaceSize[0]*spaceSize[1]*spaceSize[2]);

			MPI_Irecv(tempData1D, spaceSize[0]*spaceSize[1]*spaceSize[2], MPI_DOUBLE, rank, 1000+rank, MPI_COMM_WORLD, &request_recv);
			MPI_Wait(&request_recv, &status_recv);

			tempData3D = trans1DTo3D(tempData1D, spaceSize);
			is_change += nx_procs[rank-1];

			for (i = 0; i<spaceSize[0]; ++i)
				for (j = 0; j<spaceSize[1]; ++j)
					for (k = 0; k<spaceSize[2]; ++k)
					{
						localSARData[i+is_change][j][k] = tempData3D[i][j][k];
					}
			double*** tempP = tempData3D;
			freeArray3(tempP, spaceSize[0], spaceSize[1], spaceSize[2]);
		}
		char path_sar[128];
		sprintf(path_sar, "%slocalsar.txt ", path_save);
		FILE* fp = fopen(path_sar, "w+");
		for (k = 0; k<spaceSize[2]; ++k)
				for (j = 0; j<spaceSize[1]; ++j)
				{
					for (i = 0; i<_spaceX; ++i)
					{
						fprintf(fp, "%.12e ", localSARData[i][j][k]);
					}
					fprintf(fp, "\n");
				}
		fclose(fp);
		printf("Sending Data End.\n");
		fflush(stdout);
	}
	spaceSize[0] = _spaceX;
	if (myrank == 0)
	{
		rho = initArray3(_spaceX, _spaceY, _spaceZ, 0.0);
		char*** modelData2 = loadData3(path_model, size);
		for(i = paddingX_1; i<spaceX-paddingX_2-1; ++i)
			for(j = paddingY_1; j<spaceY-paddingY_2-1; ++j)
				for(k = paddingZ_1; k<spaceZ-paddingZ_2-1; ++k)
					rho[i-paddingX_1][j-paddingY_1][k-paddingZ_1] = media[modelData2[i][j][k]].rho;
		freeArray3Char(modelData2, spaceX, spaceY, spaceZ);
		buildObject();
		for (i = 0; i<nXgSAR; ++i)
		{
			printf("\nCompute Mass Averaged SAR %f, %d SAR calculation leave.\n", XgSAR[i], nXgSAR-i-1);
			fflush(stdout);
			computeMassAveragedSAR(localSARData, rho, spaceSize, resolution, &XgSAR[i], path_save, fp_log);
		}
	}

	return 0;
}

double*** loadLocalSAR(int* spaceSize)
{
	int i, j, k;
	double*** localSARData = NULL;

	localSARData = (double***) calloc (spaceSize[0], sizeof(double**));
	for (i = 0; i<spaceSize[0]; ++i)
	{
		localSARData[i] = (double**) calloc (spaceSize[1], sizeof(double*));
		for (j = 0; j<spaceSize[1]; ++j)
			localSARData[i][j] = (double*) calloc (spaceSize[2], sizeof(double));
	}

	if (myrank == 0)
		spaceSize[0] = nx_procs[0] - thicknessOfPml - paddingX_1 - 1;
	for (i = 0; i<spaceSize[0]; ++i)
		for (j = 0; j<spaceSize[1]; ++j)
			for (k = 0; k<spaceSize[2]; ++k)
				localSARData[i][j][k] = pSAR[k].localSARData[i][j];
	if (myrank == 0)
		spaceSize[0] = _spaceX;
	return localSARData;
}

double* trans3DTo1D(double*** data, int* size)
{
	double* newData = (double*) calloc (size[0]*size[1]*size[2], sizeof(double));
	for (int i = 0; i<size[0]; ++i)
		for (int j = 0; j<size[1]; ++j)
			for (int k = 0; k<size[2]; ++k)
				newData[i*size[1]*size[2] + j*size[2] + k] = data[i][j][k];
	
	return newData;
}
double*** trans1DTo3D(double* data, int* size)
{
	double*** newData = initArray3(size[0], size[1], size[2], 0.0);
	for (int i = 0; i<size[0]; ++i)
		for (int j = 0; j<size[1]; ++j)
			for (int k = 0; k<size[2]; ++k)
				newData[i][j][k] = data[i*size[1]*size[2] + j*size[2] + k];

	return newData;
}

#endif

