/*****************************************************************************************/
/*
 * Function name: initializeFile
 * Description: Initialize File name, path, and Creat File
 * Parameters: 
 * Return: 
 */
void initializeFile()
{
	int i;
	char temp[MAX_SIZE_OF_PATH];
	char save_path_ex[MAX_SIZE_OF_PATH], save_path_ey[MAX_SIZE_OF_PATH], save_path_ez[MAX_SIZE_OF_PATH];
	char save_path_hx[MAX_SIZE_OF_PATH], save_path_hy[MAX_SIZE_OF_PATH], save_path_hz[MAX_SIZE_OF_PATH];

	save_field_file_name = (char**) malloc (save_plane_amount * sizeof(char*));

	for (i = 0; i < save_plane_amount; i++)
	{
		//FILE to Save
		strcpy(save_path_ex, path_save);
		strcpy(save_path_ey, path_save);
		strcpy(save_path_ez, path_save);
		strcpy(save_path_hx, path_save);
		strcpy(save_path_hy, path_save);
		strcpy(save_path_hz, path_save);
		save_field_file_name[i] = (char*) malloc (MAX_SIZE_OF_PATH * sizeof(char));
		strcpy(save_field_file_name[i], proj_name);
		switch (fp_save_field_file[i].sp.plane_no)
		{
		case 1 :
			strcat(save_field_file_name[i], ".xy");
			int2str(fp_save_field_file[i].sp.slice - thicknessOfPml, temp);
			strcat(save_field_file_name[i], temp);
			break;
		case 2 :
			strcat(save_field_file_name[i], ".yz");
			int2str(fp_save_field_file[i].sp.slice - thicknessOfPml, temp);
			strcat(save_field_file_name[i], temp);
			break;
		case 3 :
			strcat(save_field_file_name[i], ".xz");
			int2str(fp_save_field_file[i].sp.slice - thicknessOfPml, temp);
			strcat(save_field_file_name[i], temp);
			break;
		default :
			break;
		}
		int2str(myrank, temp);

		fp_save_field_file[i].ex = fopen(strcat(strcat(strcat(save_path_ex, save_field_file_name[i]),"_ex.txt"),temp), "w+");
		fp_save_field_file[i].ey = fopen(strcat(strcat(strcat(save_path_ey, save_field_file_name[i]),"_ey.txt"),temp), "w+");
		fp_save_field_file[i].ez = fopen(strcat(strcat(strcat(save_path_ez, save_field_file_name[i]),"_ez.txt"),temp), "w+");
		fp_save_field_file[i].hx = fopen(strcat(strcat(strcat(save_path_hx, save_field_file_name[i]),"_hx.txt"),temp), "w+");
		fp_save_field_file[i].hy = fopen(strcat(strcat(strcat(save_path_hy, save_field_file_name[i]),"_hy.txt"),temp), "w+");
		fp_save_field_file[i].hz = fopen(strcat(strcat(strcat(save_path_hz, save_field_file_name[i]),"_hz.txt"),temp), "w+");
	}

#ifdef _SAR
	int n;
	pSAR = initializeLocalSAR(pSAR, save_localSAR_amount+1);
	for (n = 0; n<save_localSAR_amount; ++n)
	{
		pSAR[n].localSARData = (double**) calloc (spaceX, sizeof(double*));
		for (i = 0; i<spaceX; ++i)
			pSAR[n].localSARData[i] = (double*) calloc (spaceY, sizeof(double));

		char temp_SAR[MAX_SIZE_OF_PATH];
		char save_path_localSAR[MAX_SIZE_OF_PATH];
		strcpy(save_path_localSAR, path_save);
		strcat(save_path_localSAR, proj_name);
		strcat(save_path_localSAR, ".xy");
		int2str(pSAR[n].slice - thicknessOfPml, temp_SAR);
		strcat(save_path_localSAR, temp_SAR);
		strcat(save_path_localSAR, ".sar");
		int2str(myrank, temp_SAR);
		pSAR[n].fp = fopen(strcat(save_path_localSAR, temp_SAR), "w+");
	}
#endif

}

/*****************************************************************************************/
/*
 * Function name: initializePart1
 * Description: Initialize parameters, init Arrays
 * Parameters: 
 * Return: 
 */
void initializePart1()
{
	if (myrank == 0)
	{
		printf("Initializing Array...\n");
		fflush(stdout);
	}
	muO = 4.0 * pi * 1.0E-7;
	epsO = 1.0 / (C * C * muO);

	//PML Layers(With PEC)
	nxPML_1 = thicknessOfPml + 1;
	nxPML_2 = thicknessOfPml + 1;
	nyPML_1 = thicknessOfPml + 1;
	nyPML_2 = thicknessOfPml + 1;
	nzPML_1 = thicknessOfPml + 1;
	nzPML_2 = thicknessOfPml + 1;

	//Dynamic memory allocation
	if (myrank != nprocs-1 && nprocs > 1)
		Imax = Imax + 1;
	else if (nprocs == 1)
		Imax = Imax - 1;

	CAx = initArray3(Imax, Jmax-1, Kmax-1, 1.0);
	CBx = initArray3(Imax, Jmax-1, Kmax-1, 0.0);

	CAy = initArray3(Imax, Jmax-1, Kmax-1, 1.0);
	CBy = initArray3(Imax, Jmax-1, Kmax-1, 0.0);

	CAz = initArray3(Imax, Jmax-1, Kmax-1, 1.0);
	CBz = initArray3(Imax, Jmax-1, Kmax-1, 0.0);

	ID1_EpsR = initArray3(Imax+1, Jmax, Kmax, 1.0);/*多LOAD一格以读取真实的ie+1格系数*/
	ID1_SigR = initArray3(Imax+1, Jmax, Kmax, 0.0);

	ID2_EpsR = initArray3(Imax+1, Jmax, Kmax, 1.0);
	ID2_SigR = initArray3(Imax+1, Jmax, Kmax, 0.0);

	ID3_EpsR = initArray3(Imax+1, Jmax, Kmax, 1.0);
	ID3_SigR = initArray3(Imax+1, Jmax, Kmax, 0.0);

#ifdef _SAR
	ID1_RhoR = initArray3(Imax+1, Jmax, Kmax, 0.0);
	ID2_RhoR = initArray3(Imax+1, Jmax, Kmax, 0.0);
	ID3_RhoR = initArray3(Imax+1, Jmax, Kmax, 0.0);
#endif

	psi_Ezx_1 = initArray3(nxPML_1, Jmax, Kmax, 0.0);

	psi_Ezx_2 = initArray3(nxPML_2, Jmax, Kmax, 0.0);

	psi_Hyx_1 = initArray3(nxPML_1-1, Jmax, Kmax, 0.0);

	psi_Hyx_2 = initArray3(nxPML_2-1, Jmax, Kmax, 0.0);

	psi_Ezy_1 = initArray3(Imax, nyPML_1, Kmax, 0.0);

	psi_Ezy_2 = initArray3(Imax, nyPML_2, Kmax, 0.0);

	psi_Hxy_1 = initArray3(Imax, nyPML_1-1, Kmax, 0.0);

	psi_Hxy_2 = initArray3(Imax, nyPML_2-1, Kmax, 0.0);

	psi_Hxz_1 = initArray3(Imax, Jmax-1, nzPML_1-1, 0.0);

	psi_Hxz_2 = initArray3(Imax, Jmax-1, nzPML_2-1, 0.0);

	psi_Hyz_1 = initArray3(Imax, Jmax, nzPML_1-1, 0.0);

	psi_Hyz_2 = initArray3(Imax, Jmax, nzPML_2-1, 0.0);

	psi_Exz_1 = initArray3(Imax, Jmax, nzPML_1, 0.0);

	psi_Exz_2 = initArray3(Imax, Jmax, nzPML_2, 0.0);

	psi_Eyz_1 = initArray3(Imax, Jmax-1, nzPML_1, 0.0);

	psi_Eyz_2 = initArray3(Imax, Jmax-1, nzPML_2, 0.0);

	psi_Hzx_1 = initArray3(nxPML_1-1, Jmax-1, Kmax-1, 0.0);

	psi_Hzx_2 = initArray3(nxPML_2-1, Jmax-1, Kmax-1, 0.0);

	psi_Eyx_1 = initArray3(nxPML_1, Jmax-1, Kmax-1, 0.0);

	psi_Eyx_2 = initArray3(nxPML_2, Jmax-1, Kmax-1, 0.0);

	psi_Hzy_1 = initArray3(Imax, nyPML_1-1, Kmax-1, 0.0);

	psi_Hzy_2 = initArray3(Imax, nyPML_2-1, Kmax-1, 0.0);

	psi_Exy_1 = initArray3(Imax, nyPML_1, Kmax-1, 0.0);

	psi_Exy_2 = initArray3(Imax, nyPML_2, Kmax-1, 0.0);

	den_ex = initArray(Imax);

	den_hx = initArray(Imax);

	den_ey = initArray(Jmax-1);

	den_hy = initArray(Jmax-1);

	den_ez = initArray(Kmax-1);

	den_hz = initArray(Kmax-1);

	be_x_1 = initArray(nxPML_1);

	ce_x_1 = initArray(nxPML_1);

	alphae_x_PML_1 = initArray(nxPML_1);

	sige_x_PML_1 = initArray(nxPML_1);

	kappae_x_PML_1 = initArray(nxPML_1);

	bh_x_1 = initArray(nxPML_1-1);

	ch_x_1 = initArray(nxPML_1-1);

	alphah_x_PML_1 = initArray(nxPML_1-1);

	sigh_x_PML_1 = initArray(nxPML_1-1);

	kappah_x_PML_1 = initArray(nxPML_1-1);

	be_x_2 = initArray(nxPML_2);

	ce_x_2 = initArray(nxPML_2);

	alphae_x_PML_2 = initArray(nxPML_2);


	sige_x_PML_2 = initArray(nxPML_2);


	kappae_x_PML_2 = initArray(nxPML_2);


	bh_x_2 = initArray(nxPML_2-1);


	ch_x_2 = initArray(nxPML_2-1);

	alphah_x_PML_2 = initArray(nxPML_2-1);

	sigh_x_PML_2 = initArray(nxPML_2-1);

	kappah_x_PML_2 = initArray(nxPML_2-1);

	be_y_1 = initArray(nyPML_1);

	ce_y_1 = initArray(nyPML_1);

	alphae_y_PML_1 = initArray(nyPML_1);

	sige_y_PML_1 = initArray(nyPML_1);

	kappae_y_PML_1 = initArray(nyPML_1);

	bh_y_1 = initArray(nyPML_1-1);

	ch_y_1 = initArray(nyPML_1-1);

	alphah_y_PML_1 = initArray(nyPML_1-1);

	sigh_y_PML_1 = initArray(nyPML_1-1);

	kappah_y_PML_1 = initArray(nyPML_1-1);

	be_y_2 = initArray(nyPML_2);

	ce_y_2 = initArray(nyPML_2);

	alphae_y_PML_2 = initArray(nyPML_2);

	sige_y_PML_2 = initArray(nyPML_2);

	kappae_y_PML_2 = initArray(nyPML_2);

	bh_y_2 = initArray(nyPML_2-1);

	ch_y_2 = initArray(nyPML_2-1);

	alphah_y_PML_2 = initArray(nyPML_2-1);

	sigh_y_PML_2 = initArray(nyPML_2-1);

	kappah_y_PML_2 = initArray(nyPML_2-1);

	be_z_1 = initArray(nzPML_1);

	ce_z_1 = initArray(nzPML_1);

	alphae_z_PML_1 = initArray(nzPML_1);

	sige_z_PML_1 = initArray(nzPML_1);

	kappae_z_PML_1 = initArray(nzPML_1);

	bh_z_1 = initArray(nzPML_1-1);

	ch_z_1 = initArray(nzPML_1-1);

	alphah_z_PML_1 = initArray(nzPML_1-1);

	sigh_z_PML_1 = initArray(nzPML_1-1);

	kappah_z_PML_1 = initArray(nzPML_1-1);

	be_z_2 = initArray(nzPML_2);

	ce_z_2 = initArray(nzPML_2);

	alphae_z_PML_2 = initArray(nzPML_2);

	sige_z_PML_2 = initArray(nzPML_2);

	kappae_z_PML_2 = initArray(nzPML_2);

	bh_z_2 = initArray(nzPML_2-1);

	ch_z_2 = initArray(nzPML_2-1);

	alphah_z_PML_2 = initArray(nzPML_2-1);

	sigh_z_PML_2 = initArray(nzPML_2-1);

	kappah_z_PML_2 = initArray(nzPML_2-1);

	if (myrank != nprocs-1 && nprocs > 1)
		Imax = Imax - 1;
	else if (nprocs == 1)
		Imax = Imax + 1;

	if(myrank == 0)
	{
		printf("Done.\n");
		fflush(stdout);
	}
}

/*****************************************************************************************/
/*
 * Function name: initializePart2
 * Description: Init Field Arrays
 * Parameters: 
 * Return: 
 */
void initializePart2()
{
	if(myrank == 0)
	{
		printf("Initializing Field Arrays...");
		fflush(stdout);
	}
	if (myrank == 0 && nprocs > 1)
	{
		Ez = initArray((Imax+1) * Jmax * Kmax);

		Ey = initArray((Imax+1) * (Jmax-1) * (Kmax-1));

		Ex = initArray(Imax * Jmax * (Kmax-1));

		Hx = initArray((Imax+1) * (Jmax-1) * Kmax);

		Hy = initArray((Imax+1) * Jmax * Kmax);

		Hz = initArray((Imax+1) * (Jmax-1) * (Kmax-1));

#ifdef _POYNTING
		preEz = initArray((Imax+1) * Jmax * Kmax);

		preEy = initArray((Imax+1) * (Jmax-1) * (Kmax-1));

		preEx = initArray(Imax * Jmax * (Kmax-1));
#endif
	}
	else if (myrank == nprocs - 1 && nprocs > 1)
	{
		Ez = initArray(Imax * Jmax * Kmax);

		Ey = initArray(Imax * (Jmax-1) * (Kmax-1));

		Ex = initArray(Imax * Jmax * (Kmax-1));

		Hx = initArray(Imax * (Jmax-1) * Kmax);

		Hy = initArray((Imax+1) * Jmax * Kmax);

		Hz = initArray((Imax+1) * (Jmax-1) * (Kmax-1));

#ifdef _POYNTING
		preEz = initArray(Imax * Jmax * Kmax);

		preEy = initArray(Imax * (Jmax-1) * (Kmax-1));

		preEx = initArray(Imax * Jmax * (Kmax-1));
#endif
	}
	else if (nprocs > 1)
	{
		Ez = initArray((Imax+1) * Jmax * Kmax);

		Ey = initArray((Imax+1) * (Jmax-1) * (Kmax-1));

		Ex = initArray(Imax * Jmax * (Kmax-1));

		Hx = initArray((Imax+1) * (Jmax-1) * Kmax);

		Hy = initArray((Imax+1) * Jmax * Kmax);

		Hz = initArray((Imax+1) * (Jmax-1) * (Kmax-1));

#ifdef _POYNTING
		preEz = initArray((Imax+1) * Jmax * Kmax);

		preEy = initArray((Imax+1) * (Jmax-1) * (Kmax-1));

		preEx = initArray(Imax * Jmax * (Kmax-1));
#endif
	}
	else/* One CPU */
	{
		Ez = initArray(Imax * Jmax * Kmax);

		Ey = initArray(Imax * (Jmax-1) * (Kmax-1));

		Ex = initArray((Imax-1) * Jmax * (Kmax-1));

		Hx = initArray(Imax * (Jmax-1) * Kmax);

		Hy = initArray((Imax-1) * Jmax * Kmax);

		Hz = initArray((Imax-1) * (Jmax-1) * (Kmax-1));

		/* Mur2 */
		Ez_n_1 = initArray(Imax * Jmax * Kmax);

		Ey_n_1 = initArray(Imax * (Jmax-1) * (Kmax-1));

		Ex_n_1 = initArray((Imax-1) * Jmax * (Kmax-1));

		Ez_n = initArray(Imax * Jmax * Kmax);

		Ey_n = initArray(Imax * (Jmax-1) * (Kmax-1));

		Ex_n = initArray((Imax-1) * Jmax * (Kmax-1));

#ifdef _POYNTING
		preEz = initArray(Imax * Jmax * Kmax);

		preEy = initArray(Imax * (Jmax-1) * (Kmax-1));

		preEx = initArray((Imax-1) * Jmax * (Kmax-1));
#endif
	}
	if(myrank == 0)
	{
		printf("Done.\n");
		fflush(stdout);
	}
}


void int2str(int i, char *s)
{
	sprintf(s,"%d",i);
}

/*****************************************************************************************/
/*
 * Function name: initArray
 * Description: Init a 1D Array
 * Parameters: 
 *			 sizeOfArray
 * Return: 
 *			 Point to Array
 */
double* initArray(int sizeOfArray)
{
	double* newArray = (double*) calloc (sizeOfArray, sizeof(double));
	if (newArray == NULL)
		return (0);
	else
	{
		mem_count += sizeOfArray*sizeof(double);
		return (newArray);
	}
}

/*****************************************************************************************/
/*
 * Function name: initArray3
 * Description: Init a 3D Array
 * Parameters: 
 *			 imax、jmax、kmax
 *			 defaultValue
 * Return: 
 *			 Point to Array
 */
double*** initArray3(int imax, int jmax, int kmax, double defaultValue)
{
	int i, j, k;
	double*** data;
	try
	{
		data = (double ***) calloc (imax, sizeof(double **));
		if (data == NULL)
			throw MEM_ERROR;
		for(i = 0; i < imax; i++)
		{
			data[i] = (double **) calloc (jmax, sizeof(double *));
			if (data[i] == NULL)
				throw MEM_ERROR;

			for(j = 0; j < jmax; j++)
			{
				data[i][j] = (double *) calloc (kmax, sizeof(double));
				if (data[i][j] == NULL)
					throw MEM_ERROR;

				for(k = 0; k < kmax; k++)
				{
					data[i][j][k] = defaultValue;
				}
			}
		}
	}
	catch(int &er)
	{
		if(er == MEM_ERROR)
			cout<<"ERROR : memory allocate fail."<<endl;
	}
	mem_count += imax*jmax*kmax*sizeof(double)+imax*jmax*sizeof(double*)+imax*sizeof(double**);
	return data;
}

void freeArray3(double*** data, int imax, int jmax, int kmax)
{
	int i, j;
	for(i = 0; i < imax; i++)
	{
		for(j = 0; j < jmax; j++)
		{
			free(data[i][j]);
		}
		free(data[i]);
	}
	free(data);
	mem_count -= imax*jmax*kmax*sizeof(double)+imax*jmax*sizeof(double*)+imax*sizeof(double**);
}

