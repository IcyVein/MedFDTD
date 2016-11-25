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
    mem_count += save_plane_amount*sizeof(char*);

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
        mem_count += MAX_SIZE_OF_PATH*sizeof(char);
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
	int save_length;
    if (convergenceTarget == 0)
	    pSAR = initializeLocalSAR(pSAR, save_localSAR_amount+1);
    else
	    pSAR = initializeLocalSAR(pSAR, _spaceZ+1);

	if (myrank == 0)
		save_length = nx_procs[myrank]-nxPML_1-paddingX_1;
	else if (myrank == nprocs - 1)
		save_length = nx_procs[myrank]-nxPML_2-paddingX_2;
	else
		save_length = nx_procs[myrank];

	for (n = 0; n<save_localSAR_amount; ++n)
	{
		pSAR[n].localSARData = (float**) calloc (save_length, sizeof(float*));
        mem_count += save_length*sizeof(float*);
		for (i = 0; i<save_length; ++i)
        {
			pSAR[n].localSARData[i] = (float*) calloc (spaceY, sizeof(float));
            mem_count += spaceY * sizeof(float);
        }

		char temp_SAR[MAX_SIZE_OF_PATH];
		char save_path_localSAR[MAX_SIZE_OF_PATH];
		strcpy(save_path_localSAR, path_save);
		strcat(save_path_localSAR, proj_name);
		strcat(save_path_localSAR, ".xy");
		int2str(pSAR[n].slice - thicknessOfPml, temp_SAR);
		strcat(save_path_localSAR, temp_SAR);
		strcat(save_path_localSAR, ".sar");
		int2str(myrank, temp_SAR);
		//pSAR[n].fp = fopen(strcat(save_path_localSAR, temp_SAR), "w+"); // .sar FILE
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

	//Dynamic memory allocation
	if (myrank != nprocs-1 && nprocs > 1)
		Imax = Imax + 1;
	else if (nprocs == 1)
		Imax = Imax - 1;

    fprintf(fp_mem, "Part1_modelData_0:\t%d\n", mem_count);
	modelDataX = initArray3Char(Imax+1, Jmax, Kmax);/*多LOAD一格以读取真实的ie+1格系数*/

	modelDataY = initArray3Char(Imax+1, Jmax, Kmax);

	modelDataZ = initArray3Char(Imax+1, Jmax, Kmax);
    fprintf(fp_mem, "Part1_modelData_1:\t%d\n", mem_count);

	psi_Ezx_1 = initArray3Float(nxPML_1, Jmax, Kmax, 0.0);
    fprintf(fp_mem, "Part1_Psi_1:\t%d\n", mem_count);

	psi_Ezx_2 = initArray3Float(nxPML_2, Jmax, Kmax, 0.0);

	psi_Hyx_1 = initArray3Float(nxPML_1-1, Jmax, Kmax, 0.0);

	psi_Hyx_2 = initArray3Float(nxPML_2-1, Jmax, Kmax, 0.0);

	psi_Ezy_1 = initArray3Float(Imax, nyPML_1, Kmax, 0.0);

	psi_Ezy_2 = initArray3Float(Imax, nyPML_2, Kmax, 0.0);

	psi_Hxy_1 = initArray3Float(Imax, nyPML_1-1, Kmax, 0.0);

	psi_Hxy_2 = initArray3Float(Imax, nyPML_2-1, Kmax, 0.0);

	psi_Hxz_1 = initArray3Float(Imax, Jmax-1, nzPML_1-1, 0.0);

	psi_Hxz_2 = initArray3Float(Imax, Jmax-1, nzPML_2-1, 0.0);

	psi_Hyz_1 = initArray3Float(Imax, Jmax, nzPML_1-1, 0.0);

	psi_Hyz_2 = initArray3Float(Imax, Jmax, nzPML_2-1, 0.0);

	psi_Exz_1 = initArray3Float(Imax, Jmax, nzPML_1, 0.0);

	psi_Exz_2 = initArray3Float(Imax, Jmax, nzPML_2, 0.0);

	psi_Eyz_1 = initArray3Float(Imax, Jmax-1, nzPML_1, 0.0);

	psi_Eyz_2 = initArray3Float(Imax, Jmax-1, nzPML_2, 0.0);

	psi_Hzx_1 = initArray3Float(nxPML_1-1, Jmax-1, Kmax-1, 0.0);

	psi_Hzx_2 = initArray3Float(nxPML_2-1, Jmax-1, Kmax-1, 0.0);

	psi_Eyx_1 = initArray3Float(nxPML_1, Jmax-1, Kmax-1, 0.0);

	psi_Eyx_2 = initArray3Float(nxPML_2, Jmax-1, Kmax-1, 0.0);

	psi_Hzy_1 = initArray3Float(Imax, nyPML_1-1, Kmax-1, 0.0);

	psi_Hzy_2 = initArray3Float(Imax, nyPML_2-1, Kmax-1, 0.0);

	psi_Exy_1 = initArray3Float(Imax, nyPML_1, Kmax-1, 0.0);

	psi_Exy_2 = initArray3Float(Imax, nyPML_2, Kmax-1, 0.0);
    fprintf(fp_mem, "Part1_Psi_2:\t%d\n", mem_count);

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
		printf("Initializing Field Arrays...\n");
		fflush(stdout);
	}
	if (myrank == 0 && nprocs > 1)
	{
		maxE = new vectorD3(Imax+1, vectorD2(Jmax+1, vectorD1(Kmax+1, 0)));

		Ez = initArrayFloat((Imax+1) * Jmax * Kmax);

		Ey = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));

		Ex = initArrayFloat(Imax * Jmax * (Kmax-1));

		Hx = initArrayFloat((Imax+1) * (Jmax-1) * Kmax);

		Hy = initArrayFloat((Imax+1) * Jmax * Kmax);

		Hz = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));

#ifdef _POYNTING
		preEz = initArray((Imax+1) * Jmax * Kmax);

		preEy = initArray((Imax+1) * (Jmax-1) * (Kmax-1));

		preEx = initArray(Imax * Jmax * (Kmax-1));
#endif
	}
	else if (myrank == nprocs - 1 && nprocs > 1)
	{
		maxE = new vectorD3(Imax+1, vectorD2(Jmax+1, vectorD1(Kmax+1, 0)));

		Ez = initArrayFloat(Imax * Jmax * Kmax);

		Ey = initArrayFloat(Imax * (Jmax-1) * (Kmax-1));

		Ex = initArrayFloat(Imax * Jmax * (Kmax-1));

		Hx = initArrayFloat(Imax * (Jmax-1) * Kmax);

		Hy = initArrayFloat((Imax+1) * Jmax * Kmax);

		Hz = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));

#ifdef _POYNTING
		preEz = initArray(Imax * Jmax * Kmax);

		preEy = initArray(Imax * (Jmax-1) * (Kmax-1));

		preEx = initArray(Imax * Jmax * (Kmax-1));
#endif
	}
	else if (nprocs > 1)
	{
		maxE = new vectorD3(Imax+1, vectorD2(Jmax+1, vectorD1(Kmax+1, 0)));

		Ez = initArrayFloat((Imax+1) * Jmax * Kmax);

		Ey = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));

		Ex = initArrayFloat(Imax * Jmax * (Kmax-1));

		Hx = initArrayFloat((Imax+1) * (Jmax-1) * Kmax);

		Hy = initArrayFloat((Imax+1) * Jmax * Kmax);

		Hz = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));

#ifdef _POYNTING
		preEz = initArray((Imax+1) * Jmax * Kmax);

		preEy = initArray((Imax+1) * (Jmax-1) * (Kmax-1));

		preEx = initArray(Imax * Jmax * (Kmax-1));
#endif
	}
	else/* One CPU */
	{
		Ez = initArrayFloat(Imax * Jmax * Kmax);

		Ey = initArrayFloat(Imax * (Jmax-1) * (Kmax-1));

		Ex = initArrayFloat((Imax-1) * Jmax * (Kmax-1));

		Hx = initArrayFloat(Imax * (Jmax-1) * Kmax);

		Hy = initArrayFloat((Imax-1) * Jmax * Kmax);

		Hz = initArrayFloat((Imax-1) * (Jmax-1) * (Kmax-1));

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
}

void freeFDTDData()
{
	free(Ex);
	free(Ey);
	free(Ez);
	free(Hx);
	free(Hy);
	free(Hz);

	if (myrank != nprocs-1 && nprocs > 1)
		Imax = Imax + 1;
	else if (nprocs == 1)
		Imax = Imax - 1;

	freeArray3Float(psi_Ezx_1, nxPML_1, Jmax, Kmax);

	freeArray3Float(psi_Ezx_2, nxPML_2, Jmax, Kmax);

	freeArray3Float(psi_Hyx_1, nxPML_1-1, Jmax, Kmax);

	freeArray3Float(psi_Hyx_2, nxPML_2-1, Jmax, Kmax);

	freeArray3Float(psi_Ezy_1, Imax, nyPML_1, Kmax);

	freeArray3Float(psi_Ezy_2, Imax, nyPML_2, Kmax);

	freeArray3Float(psi_Hxy_1, Imax, nyPML_1-1, Kmax);

	freeArray3Float(psi_Hxy_2, Imax, nyPML_2-1, Kmax);

	freeArray3Float(psi_Hxz_1, Imax, Jmax-1, nzPML_1-1);

	freeArray3Float(psi_Hxz_2, Imax, Jmax-1, nzPML_2-1);

	freeArray3Float(psi_Hyz_1, Imax, Jmax, nzPML_1-1);

	freeArray3Float(psi_Hyz_2, Imax, Jmax, nzPML_2-1);

	freeArray3Float(psi_Exz_1, Imax, Jmax, nzPML_1);

	freeArray3Float(psi_Exz_2, Imax, Jmax, nzPML_2);

	freeArray3Float(psi_Eyz_1, Imax, Jmax-1, nzPML_1);

	freeArray3Float(psi_Eyz_2, Imax, Jmax-1, nzPML_2);

	freeArray3Float(psi_Hzx_1, nxPML_1-1, Jmax-1, Kmax-1);

	freeArray3Float(psi_Hzx_2, nxPML_2-1, Jmax-1, Kmax-1);

	freeArray3Float(psi_Eyx_1, nxPML_1, Jmax-1, Kmax-1);

	freeArray3Float(psi_Eyx_2, nxPML_2, Jmax-1, Kmax-1);

	freeArray3Float(psi_Hzy_1, Imax, nyPML_1-1, Kmax-1);

	freeArray3Float(psi_Hzy_2, Imax, nyPML_2-1, Kmax-1);

	freeArray3Float(psi_Exy_1, Imax, nyPML_1, Kmax-1);

	freeArray3Float(psi_Exy_2, Imax, nyPML_2, Kmax-1);

	free(den_ex);

	free(den_hx);

	free(den_ey);

	free(den_hy);

	free(den_ez);

	free(den_hz);

	free(be_x_1);

	free(ce_x_1);

	free(bh_x_1);

	free(ch_x_1);

	free(be_x_2);

	free(ce_x_2);

	free(bh_x_2);

	free(ch_x_2);

	free(be_y_1);

	free(ce_y_1);

	free(bh_y_1);

	free(ch_y_1);

	free(be_y_2);

	free(ce_y_2);

	free(bh_y_2);

	free(ch_y_2);

	free(be_z_1);

	free(ce_z_1);

	free(bh_z_1);

	free(ch_z_1);

	free(be_z_2);

	free(ce_z_2);

	free(bh_z_2);

	free(ch_z_2);

	if (myrank != nprocs-1 && nprocs > 1)
		Imax = Imax - 1;
	else if (nprocs == 1)
		Imax = Imax + 1;
}


void int2str(int i, char *s)
{
	sprintf(s,"%d",i);
}

/*****************************************************************************************/
/*
 * Function name: initArrayFloat
 * Description: Init a 1D Array
 * Parameters: 
 *			 sizeOfArray
 * Return: 
 *			 Point to Array
 */
float* initArrayFloat(int sizeOfArray)
{
	float* newArray;
	try
	{
		newArray = (float*) calloc (sizeOfArray, sizeof(float));
		if (newArray == NULL)
			throw MEM_ERROR;
	}
	catch (int &er)
	{
		if (er == MEM_ERROR)
		{
			printf("ERROR: memory allocate fail in initArray float.\n");
			exit(MEM_ERROR);
		}
		fflush(stdout);
	}
	mem_count += sizeOfArray*sizeof(float);
	return (newArray);
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
	double* newArray;
	try
	{
		newArray = (double*) calloc (sizeOfArray, sizeof(double));
		if (newArray == NULL)
			throw MEM_ERROR;
	}
	catch (int &er)
	{
		if (er == MEM_ERROR)
		{
			printf("ERROR: memory allocate fail in initArray.\n");
			exit(MEM_ERROR);
		}
		fflush(stdout);
	}
	mem_count += sizeOfArray*sizeof(double);
	return (newArray);
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
	catch (int &er)
	{
		if (er == MEM_ERROR)
		{
			printf("ERROR: memory allocate fail in initArray3.\n");
			exit(MEM_ERROR);
		}
		fflush(stdout);
	}
	mem_count += imax*jmax*kmax*sizeof(double)+imax*jmax*sizeof(double*)+imax*sizeof(double**);
	return data;
}

/*****************************************************************************************/
/*
 * Function name: initArray3Float
 * Description: Init a 3D Array
 * Parameters: 
 *			 imax、jmax、kmax
 *			 defaultValue
 * Return: 
 *			 Point to Array
 */
float*** initArray3Float(int imax, int jmax, int kmax, float defaultValue)
{
	int i, j, k;
	float*** data;
	try
	{
		data = (float ***) calloc (imax, sizeof(float **));
		if (data == NULL)
			throw MEM_ERROR;
		for(i = 0; i < imax; i++)
		{
			data[i] = (float **) calloc (jmax, sizeof(float *));
			if (data[i] == NULL)
				throw MEM_ERROR;

			for(j = 0; j < jmax; j++)
			{
				data[i][j] = (float *) calloc (kmax, sizeof(float));
				if (data[i][j] == NULL)
					throw MEM_ERROR;

				for(k = 0; k < kmax; k++)
				{
					data[i][j][k] = defaultValue;
				}
			}
		}
	}
	catch (int &er)
	{
		if (er == MEM_ERROR)
		{
			printf("ERROR: memory allocate fail in initArray3 float.\n");
			exit(MEM_ERROR);
		}
		fflush(stdout);
	}
	mem_count += imax*jmax*kmax*sizeof(float)+imax*jmax*sizeof(float*)+imax*sizeof(float**);
	return data;
}

/*
 * Function name: initArray3Char
 * Description: Init a 3D char Array
 * Parameters: 
 *			 imax、jmax、kmax
 * Return: 
 *			 Point to Array
 */
unsigned char*** initArray3Char(int imax, int jmax, int kmax)
{
	int i, j, k;
	unsigned char*** data;
	try
	{
		data = (unsigned char***) calloc (imax, sizeof(unsigned char**));
		if (data == NULL)
			throw MEM_ERROR;
		for(i = 0; i < imax; i++)
		{
			data[i] = (unsigned char**) calloc (jmax, sizeof(unsigned char*));
			if (data[i] == NULL)
				throw MEM_ERROR;

			for(j = 0; j < jmax; j++)
			{
				data[i][j] = (unsigned char*) calloc (kmax, sizeof(unsigned char));
				if (data[i][j] == NULL)
					throw MEM_ERROR;

				for(k = 0; k < kmax; k++)
				{
					data[i][j][k] = 0; /* vaccum */
				}
			}
		}
	}
	catch (int &er)
	{
		if (er == MEM_ERROR)
		{
			printf("ERROR: memory allocate fail in initArray3 char.\n");
			exit(MEM_ERROR);
		}
		fflush(stdout);
	}
	mem_count += imax*jmax*kmax*sizeof(unsigned char)+imax*jmax*sizeof(unsigned char*)+imax*sizeof(unsigned char**);
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

void freeArray3Float(float*** data, int imax, int jmax, int kmax)
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
	mem_count -= imax*jmax*kmax*sizeof(float)+imax*jmax*sizeof(float*)+imax*sizeof(float**);
}

void freeArray3Char(unsigned char*** data, int imax, int jmax, int kmax)
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
	mem_count -= imax*jmax*kmax*sizeof(unsigned char)+imax*jmax*sizeof(unsigned char*)+imax*sizeof(unsigned char**);
}
/*
void mycolloc(int size, int sizeoftype)
{
	MEMORY_NUM += size * sizeoftype;
	return calloc(size, sizeoftype);
}
*/

