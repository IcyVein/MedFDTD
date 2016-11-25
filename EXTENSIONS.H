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

#define T(i,j,k) T[(i)*size[1]*size[2]+(j)*size[2]+k]
int computeTemperatureRise(char* path_sar)
{
	char path_temp[MAX_SIZE_OF_PATH];
    sprintf(path_temp, "%stemperature_rise_%d.txt", path_save, myrank);

	int i, j, k;

	if (myrank == 0)
		cout<<"Load model data."<<endl;
	fflush(stdout);
	unsigned char*** model_w = loadData3T(path_model, size);
	if (myrank == 0)
		cout<<"Load sar data."<<endl;
	fflush(stdout);
	double*** sar_w = loadSar(path_sar, size);

	int* nx_procs = (int*) calloc (nprocs, sizeof(int));
	int is=0, ie=0;
	for (i = 0; i<nprocs; ++i)
		nx_procs[i] = size[0]/nprocs;
	for (i = nprocs - 1; i>=nprocs - size[0]%nprocs; --i)
		nx_procs[i]++;
	for (i = 0; i<myrank; ++i)
		is += nx_procs[i];
	ie = is + nx_procs[myrank]-1;
	size[0] = nx_procs[myrank];

	unsigned char*** model;
	double*** sar;
	if (myrank == 0)
	{
		model = initArray3Char(size[0]+1, size[1], size[2]);
		for (i = 0; i<size[0]+1; ++i)
			for (j = 0; j<size[1]; ++j)
				for (k = 0; k<size[2]; ++k)
					model[i][j][k] = model_w[i+is][j][k];
		sar = (double***) calloc (size[0]+1, sizeof(double**));
		for (i = 0; i<size[0]+1; ++i)
		{
			sar[i] = (double**) calloc (size[1], sizeof(double*));
			for (j = 0; j<size[1]; ++j)
			{
				sar[i][j] = (double*) calloc (size[2], sizeof(double));
				for (k = 0; k<size[2]; ++k)
					sar[i][j][k] = sar_w[i+is][j][k];
			}
		}
	}
	else if (myrank == nprocs-1)
	{
		model = initArray3Char(size[0]+1, size[1], size[2]);
		for (i = 0; i<size[0]+1; ++i)
			for (j = 0; j<size[1]; ++j)
				for (k = 0; k<size[2]; ++k)
					model[i][j][k] = model_w[i+is-1][j][k];
		sar = (double***) calloc (size[0]+1, sizeof(double**));
		for (i = 0; i<size[0]+1; ++i)
		{
			sar[i] = (double**) calloc (size[1], sizeof(double*));
			for (j = 0; j<size[1]; ++j)
			{
				sar[i][j] = (double*) calloc (size[2], sizeof(double));
				for (k = 0; k<size[2]; ++k)
					sar[i][j][k] = sar_w[i+is-1][j][k];
			}
		}
	}
	else
	{
		model = initArray3Char(size[0]+2, size[1], size[2]);
		for (i = 0; i<size[0]+2; ++i)
			for (j = 0; j<size[1]; ++j)
				for (k = 0; k<size[2]; ++k)
					model[i][j][k] = model_w[i+is-1][j][k];
		sar = (double***) calloc (size[0]+2, sizeof(double**));
		for (i = 0; i<size[0]+2; ++i)
		{
			sar[i] = (double**) calloc (size[1], sizeof(double*));
			for (j = 0; j<size[1]; ++j)
			{
				sar[i][j] = (double*) calloc (size[2], sizeof(double));
				for (k = 0; k<size[2]; ++k)
					sar[i][j][k] = sar_w[i+is-1][j][k];
			}
		}
	}

    MPI_Barrier(MPI_COMM_WORLD);
	if (myrank == 0)
		cout<<"Computing"<<endl;
	fflush(stdout);
	temperatureRise(path_temp, sar, model, size);
	if (myrank == 0)
		cout<<"Done."<<endl;
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);

	return 0;

}
int temperatureRise(char* path, double*** sar, unsigned char*** model, int* size)
{
	int i, j, k, n;

	double dx = ds;
	double dt = dt_T;
	int nmin = 1;
	int nmax = (int)(nMaxT/dt_T);
	double Tmax = nmax*dt;
	double ta = 0.0, tb = 0.0;
	double h = 10.5;
	double rho, c, K, B;
	int backup = 0;

    float* T = (float*) calloc ((size[0]+2)*size[1]*size[2], sizeof(float));//T[0]&T[max] for data exchage.
	for (i = 0; i<size[0]+2; ++i)
		for (j = 0; j<size[1]; ++j)
			for (k = 0; k<size[2]; ++k)
			{
				T(i, j, k) = ta;
			}
	if (myrank == 0)
	{
		for (i = 1; i<size[0]+2; ++i)
			for (j = 0; j<size[1]; ++j)
				for (k = 0; k<size[2]; ++k)
				{
					if (model[i-1][j][k])
					{
						T(i, j, k) = tb;
					}
				}
	}
	else if (myrank == nprocs - 1)
	{
		for (i = 0; i<size[0]+1; ++i)
			for (j = 0; j<size[1]; ++j)
				for (k = 0; k<size[2]; ++k)
				{
					if (model[i][j][k])
					{
						T(i, j, k) = tb;
					}
				}
	}
	else 
	{
		for (i = 0; i<size[0]+2; ++i)
			for (j = 0; j<size[1]; ++j)
				for (k = 0; k<size[2]; ++k)
				{
					if (model[i][j][k])
						T(i, j, k) = tb;
				}
	}
	/*if (backup)
	{
		if (myrank == 0)
		{
			printf("Loading backup %d\n", nmin/10000);
			fflush(stdout);
		}
		int return_val = loadBackUp(T, size, nmin/10000);
		if (myrank == 0 && return_val == 0)
		{
			printf("Load backup done.\n");
			fflush(stdout);
		}
	}*/

	time_t t = time(NULL);
	time_t tt;
	FILE* fp = fopen(path, "w+");
	char path_backup[MAX_SIZE_OF_PATH];
	FILE* fp_backup;

	MPI_Request request_send_r, request_send_l, request_recv_r, request_recv_l;
	MPI_Status status_send_r, status_send_l, status_recv_r, status_recv_l;

    MPI_Barrier(MPI_COMM_WORLD);
	for (n = nmin; n<=nmax; ++n)
	{
        MPI_Barrier(MPI_COMM_WORLD);
		for (i = 1; i<=size[0]; ++i)
		{
			for (j = 1; j<size[1]-1; ++j)
			{
				for (k = 1; k<size[2]-1; ++k)
				{
					if (model[i][j][k] != 0)
					{
						rho = media[model[i][j][k]].rho;
						c = media[model[i][j][k]].spec_heat;
						K = media[model[i][j][k]].K;
						B = media[model[i][j][k]].B;
						T(i, j, k) = T(i, j, k) + dt/rho/c*(
								K*((T(i+1, j, k)-2*T(i, j, k)+T(i-1, j, k))/dx/dx+(T(i, j+1, k)-2*T(i, j, k)+T(i, j-1, k))/dx/dx+(T(i, j, k+1)-2*T(i, j, k)+T(i, j, k-1))/dx/dx)+rho*sar[i][j][k]-B*(T(i, j, k)-tb) );
					}
					else 
					{
						rho = media[model[i][j][k]].rho;
						c = media[model[i][j][k]].spec_heat;
						K = media[model[i][j][k]].K;
						B = media[model[i][j][k]].B;
						h = 50.0;
						T(i, j, k) = (K*T(i, j, k) + ta*h*dx) / (K+h*dx);
					}
				}
			}
		}
        MPI_Barrier(MPI_COMM_WORLD);
        if (myrank == 0)
        {
            MPI_Isend(&T(size[0], 0, 0), size[1]*size[2], MPI_FLOAT, myrank+1, 100+myrank+1, MPI_COMM_WORLD, &request_send_r);
			MPI_Wait(&request_send_r, &status_send_r);
        }
        else if (myrank == nprocs-1)
        {
            MPI_Irecv(&T(0, 0, 0), size[1]*size[2], MPI_FLOAT, myrank-1, 100+myrank, MPI_COMM_WORLD, &request_recv_l);
			MPI_Wait(&request_recv_l, &status_recv_l);
        }
        else
        {
            MPI_Isend(&T(size[0], 0, 0), size[1]*size[2], MPI_FLOAT, myrank+1, 100+myrank+1, MPI_COMM_WORLD, &request_send_r);
			MPI_Wait(&request_send_r, &status_send_r);
            MPI_Irecv(&T(0, 0, 0), size[1]*size[2], MPI_FLOAT, myrank-1, 100+myrank, MPI_COMM_WORLD, &request_recv_l);
			MPI_Wait(&request_recv_l, &status_recv_l);
        }
        MPI_Barrier(MPI_COMM_WORLD);
		if (myrank == 0)
        {
            MPI_Irecv(&T(size[0]+1, 0, 0), size[1]*size[2], MPI_FLOAT, myrank+1, 200+myrank, MPI_COMM_WORLD, &request_recv_r);
			MPI_Wait(&request_recv_r, &status_recv_r);
        }
        else if (myrank == nprocs-1)
        {
            MPI_Isend(&T(1, 0, 0), size[1]*size[2], MPI_FLOAT, myrank-1, 200+myrank-1, MPI_COMM_WORLD, &request_send_l);
			MPI_Wait(&request_send_l, &status_send_l);
        }
        else
        {
            MPI_Isend(&T(1, 0, 0), size[1]*size[2], MPI_FLOAT, myrank-1, 200+myrank-1, MPI_COMM_WORLD, &request_send_l);
			MPI_Wait(&request_send_l, &status_send_l);
            MPI_Irecv(&T(size[0]+1, 0, 0), size[1]*size[2], MPI_FLOAT, myrank+1, 200+myrank, MPI_COMM_WORLD, &request_recv_r);
			MPI_Wait(&request_recv_r, &status_recv_r);
        }
		MPI_Barrier(MPI_COMM_WORLD);
		if (n%500 == 0 || n==nmax)
		{
			for (j = 0; j<size[1]; ++j)
			{
				for (i = 1; i<=size[0]-1; ++i)
					fprintf(fp, "%e ", T(i, j, 118));
				fprintf(fp, "%e\n", T(i, j, 118));
			}
			tt = time(NULL)-t;
			if (myrank == 0)
				printf("Step : %d/%d (%f), need %d/%f\n", n, nmax, (double)tt/(double)(n-nmin+1), tt, (double)tt/(double)(n-nmin+1)*(nmax-nmin+1));
			fflush(stdout);
		}
		/*if (n%10000==0 || n==nmax)
		{
			sprintf(path_backup, ".\\backup_%d_%d.txt", myrank, n/10000);
			fp_backup = fopen(path_backup, "w+");
			for (k = 0; k<size[2]; ++k)
				for (j = 0; j<size[1]; ++j)
				{
					for (i = 0; i<=size[0]; ++i)
						fprintf(fp_backup, "%e ", T(i, j, k));
					fprintf(fp_backup, "%e\n", T(i, j, k));
				}
			fclose(fp_backup);
		}*/
	}

	fclose(fp);

	return 0;
}

unsigned char*** loadData3T(char* path, int* size)
{
    int i;
	int ldi, ldj, ldk;
	int schedule = 0, schedule_past = 0; 
	unsigned char*** A;
	int temp_data;
	unsigned char temp_c;
	FILE *fp;
	time_t load_time;

	A = initArray3Char(size[0], size[1], size[2]);

	if (path == NULL)
	{
		printf("None model import\n");
		fflush(stdout);
		return (A);
	}

	if (!(fp = fopen(path,"r")))
	{
		printf("Open data file ");
		printf(path);
		printf(" fail.\n");
		return (A);
	}

	load_time = clock();

	for (ldk = 0; ldk < size[2]; ldk++)
	{
		for (ldj = 0; ldj < size[1]; ldj++)
		{
			for (ldi = 0; ldi < size[0] - 1; ldi++)
			{
				fscanf(fp, "%d ", &temp_data);
				A[ldi][ldj][ldk] = temp_data;
			}
			fscanf(fp, "%d", &temp_data);
			A[ldi][ldj][ldk] = temp_data;
			fscanf(fp, "%c", &temp_c);
		}
		schedule_past = schedule;
		schedule = (int)((ldk)/(double)(size[2]-1)*20);
		if(myrank == 0 && (clock() - load_time > 2000 || schedule-schedule_past))
		{
			printf("\rLoading model: |");
            for (i = 0; i<schedule; ++i)
			{
                printf("=");
			}
            for (i = schedule; i<20; ++i)
			{
                printf(" ");
			}
            printf("| %4.1f%%", (ldk)/(double)(size[2]-1)*100);
			fflush(stdout);
			load_time = clock();
		}
	}
	if (myrank == 0)
		printf("\n");
	fclose(fp);
	return (A);
}

double*** loadSar(char* path, int* size)
{
	int i, j, k;
	double temp;

	double*** sar = (double***) calloc (size[0], sizeof(double**));
	for (i = 0; i<size[0]; ++i)
	{
		sar[i] = (double**) calloc (size[1], sizeof(double*));
		for (j = 0; j<size[1]; ++j)
		{
			sar[i][j] = (double*) calloc (size[2], sizeof(double));
			for (k = 0; k<size[2]; ++k)
				sar[i][j][k] = 0.0;
		}
	}

	FILE* fp = fopen(path, "r");
	if (fp==NULL)
	{
		printf("%s is no found.\n", path);
		exit(1);
	}
	for (k = 0; k<size[2]; ++k)
	{
		for(j = 0; j<size[1]; ++j)
		{
			for (i = 0; i<size[0]-1; ++i)
			{
				fscanf(fp, "%lf ", &temp);
				sar[i][j][k] = temp;
			}
			fscanf(fp, "%lf\n", &temp);
			sar[i][j][k] = temp;
		}
	}
	fclose(fp);

	return sar;
}

int loadBackUp(float* T, int* size, int num)
{
	char filename[256];
	sprintf(filename, "backup_%d_%d.txt", myrank, num);
	FILE* fp = fopen(filename, "r");

	int i, j, k;
	float temp;

	for (k = 0; k<size[2]; ++k)
		for (j = 0; j<size[1]; ++j)
		{
			for (i = 0; i<size[0]; ++i)
			{
				fscanf(fp, "%e ", &temp);
				T(i, j, k) = temp;
			}
			fscanf(fp, "%e\n", &temp);
			T(i, j, k) = temp;
		}
	MPI_Status status_send_l, status_recv_r;
	MPI_Request request_send_l, request_recv_r;
	MPI_Barrier(MPI_COMM_WORLD);
	if (myrank == 0)
    {
        MPI_Irecv(&T(size[0]+1, 0, 0), size[1]*size[2], MPI_FLOAT, myrank+1, 200+myrank, MPI_COMM_WORLD, &request_recv_r);
		MPI_Wait(&request_recv_r, &status_recv_r);
    }
    else if (myrank == nprocs-1)
    {
        MPI_Isend(&T(1, 0, 0), size[1]*size[2], MPI_FLOAT, myrank-1, 200+myrank-1, MPI_COMM_WORLD, &request_send_l);
		MPI_Wait(&request_send_l, &status_send_l);
    }
    else
    {
        MPI_Isend(&T(1, 0, 0), size[1]*size[2], MPI_FLOAT, myrank-1, 200+myrank-1, MPI_COMM_WORLD, &request_send_l);
		MPI_Wait(&request_send_l, &status_send_l);
        MPI_Irecv(&T(size[0]+1, 0, 0), size[1]*size[2], MPI_FLOAT, myrank+1, 200+myrank, MPI_COMM_WORLD, &request_recv_r);
		MPI_Wait(&request_recv_r, &status_recv_r);
    }
	MPI_Barrier(MPI_COMM_WORLD);

	return 0;
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
		pSAR[n].RMS_ex = (float**) calloc (save_length, sizeof(float*));
		pSAR[n].RMS_ey = (float**) calloc (save_length+1, sizeof(float*));
		pSAR[n].RMS_ez = (float**) calloc (save_length+1, sizeof(float*));
        mem_count += (save_length*3+2) * sizeof(float*);
		for (i = 0; i<save_length+1; ++i)
		{
			pSAR[n].RMS_ey[i] = (float*) calloc (_spaceY+1, sizeof(float));
			pSAR[n].RMS_ez[i] = (float*) calloc (_spaceY+1, sizeof(float));
            mem_count += (_spaceY+1)*2*sizeof(float);
		}
		for (i = 0; i<save_length; ++i)
		{
			pSAR[n].RMS_ex[i] = (float*) calloc (_spaceY+1, sizeof(float));
            mem_count += (_spaceY+1)*sizeof(float);
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
            mem_count -= (_spaceY+1)*2*sizeof(float);
		}
		for (i = 0; i<save_length; ++i)
		{
			free(pSAR[n].RMS_ex[i]);
            mem_count -= (_spaceY+1)*sizeof(float);
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

//	if( (*pSAR).start > timestep || (*pSAR).end < timestep)	return;

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
/*
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
	}*/
}

float checkRMS()
{
    int i, j, k;
    float meanRMS = 0.0;

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
    	for (j = 0; j<_spaceY+1; ++j)
    	{
            for(k = 0; k<_spaceZ+1; ++k)
            {
                meanRMS += pSAR[k].RMS_ez[i][j];
            }
        }
    }

    if (save_length > 0)
        meanRMS = meanRMS/save_length/(_spaceY+1)/(_spaceZ+1)/sourceT;

    return meanRMS;
}

void resetRMS()
{
    int i, j, k;

	int save_length;
	if (myrank == 0)
		save_length = nx_procs[myrank]-nxPML_1-paddingX_1;
	else if (myrank == nprocs - 1)
		save_length = nx_procs[myrank]-nxPML_2-paddingX_2;
	else
		save_length = nx_procs[myrank];

    for (k = 0; k<_spaceZ+1; ++k)
	{
		for (i = 0; i<save_length+1; ++i)
		{
			for (j = 0; j<_spaceY+1; ++j)
			{
				pSAR[k].RMS_ey[i][j] = 0;
				pSAR[k].RMS_ez[i][j] = 0;
			}
		}

		for (i = 0; i<save_length; ++i)
		{
			for (j = 0; j<_spaceY+1; ++j)
			{
				pSAR[k].RMS_ex[i][j] = 0;
			}
		}
	}
}

void computeLocalSAR(localSAR sSAR1, localSAR sSAR2, float*** localSARData)
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
				(*localSARData)[i][j] = (
                                        (sSAR1.RMS_ex[i][j]    * media[modelDataX[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1]].sigma      / media[modelDataX[i+pml_change_X][j+pml_change_Y][sSAR1.slice-1]].rho      + 
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
										 sSAR1.RMS_ez[i+1][j+1]* media[modelDataZ[i+pml_change_X+1][j+pml_change_Y+1][sSAR1.slice]].sigma    / media[modelDataZ[i+pml_change_X+1][j+pml_change_Y+1][sSAR1.slice]].rho     )/4
                                        )/sSAR1.length;
			}
			else
			{
				(*localSARData)[i][j] = 0.0;
			}
		}
	}
}

void writeLocalSAR(FILE* fp, float** localSARData)
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
	float*** localSARData = NULL;

	int spaceSize[3] = {_spaceX, _spaceY, _spaceZ};

	if (myrank)
		spaceSize[0] = nx_procs[myrank];
	if (myrank == nprocs - 1)
		spaceSize[0] -= thicknessOfPml + paddingX_2 + 2; /* Bug: may cause <0 */
	if (spaceSize[0]<0)
		spaceSize[0] = 0;

	freeLocalSARRMS(pSAR, _spaceZ);
	localSARData = loadLocalSAR(spaceSize);
	freeLocalSARData(pSAR, _spaceZ);

	float* tempData1D;
	float*** tempData3D;
	if (myrank)
	{
		tempData1D = trans3DTo1D(localSARData, spaceSize);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Request request_send;

	MPI_Request request_recv; 
	MPI_Status status_recv;

	if (myrank && nprocs>1)
	{
		MPI_Isend(tempData1D, spaceSize[0]*spaceSize[1]*spaceSize[2], MPI_FLOAT, 0, 1000+myrank, MPI_COMM_WORLD, &request_send);
	}
	else if (nprocs>1)
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
			tempData1D = initArrayFloat(spaceSize[0]*spaceSize[1]*spaceSize[2]);

			MPI_Irecv(tempData1D, spaceSize[0]*spaceSize[1]*spaceSize[2], MPI_FLOAT, rank, 1000+rank, MPI_COMM_WORLD, &request_recv);
			MPI_Wait(&request_recv, &status_recv);

			tempData3D = trans1DTo3D(tempData1D, spaceSize);
			is_change += nx_procs[rank-1];

			for (i = 0; i<spaceSize[0]; ++i)
				for (j = 0; j<spaceSize[1]; ++j)
					for (k = 0; k<spaceSize[2]; ++k)
					{
						localSARData[i+is_change][j][k] = tempData3D[i][j][k];
					}
			float* tempP1 = tempData1D;
			free(tempP1);
			float*** tempP3 = tempData3D;
			freeArray3Float(tempP3, spaceSize[0], spaceSize[1], spaceSize[2]);
		}
		printf("Sending Data End.\n");
		fflush(stdout);
	}
	spaceSize[0] = _spaceX;
	if (myrank == 0 && nprocs>1)
	{
		rho = initArray3Float(_spaceX, _spaceY, _spaceZ, 0.0);
		unsigned char*** modelData2 = loadData3(path_model, size);
        unsigned char*** modelData3 = initArray3Char(size[0], size[1], size[2]);
		char path_sar[MAX_SIZE_OF_PATH];
		printf("Writting local sar.\n");
        fflush(stdout);
		sprintf(path_sar, "%slocalsar.txt ", path_save);
		FILE* fp_localsar = fopen(path_sar, "w+");
		for (k = 0; k<spaceSize[2]; ++k)
			for (j = 0; j<spaceSize[1]; ++j)
			{
				for (i = 0; i<_spaceX; ++i)
				{
					if (fabs(media[modelData2[i+paddingX_1][j+paddingY_1][k+paddingZ_1]].rho) < 1e-6)
						localSARData[i][j][k] = 0.0;
					fprintf(fp_localsar, "%e ", localSARData[i][j][k]);
				}
				fprintf(fp_localsar, "\n");
			}
		fclose(fp_localsar);
		for(i = paddingX_1; i<spaceX-paddingX_2-1; ++i)
			for(j = paddingY_1; j<spaceY-paddingY_2-1; ++j)
				for(k = paddingZ_1; k<spaceZ-paddingZ_2-1; ++k)
                {
                    modelData3[i-paddingX_1][j-paddingY_1][k-paddingZ_1] = modelData2[i][j][k];
					rho[i-paddingX_1][j-paddingY_1][k-paddingZ_1] = media[modelData2[i][j][k]].rho*dx*dy*dz;  /* Mass of cell */
                }
		freeArray3Char(modelData2, spaceX, spaceY, spaceZ);
        double rho2[MAX_NUM_OF_MEDIA];
        for (i = 0; i<maxMedia; ++i)
        {
            rho2[i] = media[i].rho*dx*dy*dz;
        }
		buildObject();
		if (whole_body_sar)
		{
			//compute_whole_body_average_sar(localSARData, rho, spaceSize, fp_log);
			compute_whole_body_average_sar(localSARData, modelData3, rho2, spaceSize, fp_log);
		}
		for (i = 0; i<nXgSAR; ++i)
		{
			printf("\nCompute mass averaged sar %fkg, %d sar calculation leave.\n", XgSAR[i], nXgSAR-i-1);
			fprintf(fp_log, "\nCompute mass averaged sar %fkg, %d sar calculation leave.\n", XgSAR[i], nXgSAR-i-1);
			fflush(stdout);
			compute_mass_average_sar(localSARData, rho, spaceSize, XgSAR[i], 0.2, path_save, fp_log);
		}
	}

	return 0;
}

float*** loadLocalSAR(int* spaceSize)
{
	int i, j, k;
	float*** localSARData = NULL;

	localSARData = (float***) calloc (spaceSize[0], sizeof(float**));
    mem_count += spaceSize[0]*sizeof(float**);
	for (i = 0; i<spaceSize[0]; ++i)
	{
		localSARData[i] = (float**) calloc (spaceSize[1], sizeof(float*));
        mem_count += spaceSize[1]*sizeof(float*);
		for (j = 0; j<spaceSize[1]; ++j)
        {
			localSARData[i][j] = (float*) calloc (spaceSize[2], sizeof(float));
            mem_count += spaceSize[2]*sizeof(float);
        }
	}

	if (myrank == 0 && nprocs>1)
		spaceSize[0] = nx_procs[0] - thicknessOfPml - paddingX_1 - 1;
	for (i = 0; i<spaceSize[0]; ++i)
		for (j = 0; j<spaceSize[1]; ++j)
			for (k = 0; k<spaceSize[2]; ++k)
				localSARData[i][j][k] = pSAR[k].localSARData[i][j];
	if (myrank == 0)
		spaceSize[0] = _spaceX;
	return localSARData;
}

float* trans3DTo1D(float*** data, int* size)
{
	float* newData = (float*) calloc (size[0]*size[1]*size[2], sizeof(float));
	for (int i = 0; i<size[0]; ++i)
		for (int j = 0; j<size[1]; ++j)
			for (int k = 0; k<size[2]; ++k)
				newData[i*size[1]*size[2] + j*size[2] + k] = data[i][j][k];
	
	return newData;
}
float*** trans1DTo3D(float* data, int* size)
{
	float*** newData = initArray3Float(size[0], size[1], size[2], 0.0);
	for (int i = 0; i<size[0]; ++i)
		for (int j = 0; j<size[1]; ++j)
			for (int k = 0; k<size[2]; ++k)
				newData[i][j][k] = data[i*size[1]*size[2] + j*size[2] + k];

	return newData;
}

#endif

