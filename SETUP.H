int setUp()
{
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  COMPUTING FIELD UPDATE EQUATION COEFFICIENTS
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	DA = 1.0;
	DB = dt / muO;
	media[0].epsilon = 1.0; /* Vaccum */
    media[0].epsilon_inf = 1.0;
    media[0].delay= 0.0;
	media[0].sigma = 0.0;
	media[0].rho = 0.0;
	media[1].epsilon = 0.0; /* PEC */
    media[1].epsilon_inf = 0.0;
    media[1].delay= 0.0;
	media[1].sigma = 0.0;
	media[1].rho = 0.0;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  LOAD MODEL AND COMPUTE CA,CB
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	size[0] = _spaceX;
	size[1] = _spaceY;
	size[2] = _spaceZ;

	//FILE to Load
    /*
	strcpy(path_model, path_data);
	strcat(path_model, model_name);
	strcpy(path_media, path_data);
	strcat(path_media, media_name);*/
    sprintf(path_model, "%s%s", path_data, model_name);
    sprintf(path_media, "%s%s", path_data, media_name);
	
    if (media_dispersion == 1)
    {
        if (loadMediaData_Dispersion(path_media, mediaNum) == MEM_ERROR)
            return MEM_ERROR;
    }
    else
    {
	    if (loadMediaData(path_media, mediaNum) == MEM_ERROR)
            return MEM_ERROR;
    }
	try
	{
		modelData = loadData3(path_model, size);
	}
	catch(int &er)
	{
		exit(MEM_ERROR);
	}

	CA[0] = 1.0;	/* vaccum */
	CB[0] = dt/epsO;
    CB_PML[0] = dt/epsO;
	CA[1] = 0.0;	/* PEC */
	CB[1] = 0.0;
	CB_PML[1] = 0.0;

	if (nprocs == 1)
		loadModelOneCPU();
	else
		loadModel();
	freeArray3Char(modelData, spaceX, spaceY, spaceZ);

	buildObject();

	if (myrank != nprocs-1 && nprocs>1)
		Imax = Imax + 1;
	else if (nprocs == 1)
		Imax = Imax - 1;
/*	freeArray3Char(modelDataX, Imax+1, Jmax, Kmax);
	freeArray3Char(modelDataY, Imax+1, Jmax, Kmax);
	freeArray3Char(modelDataZ, Imax+1, Jmax, Kmax);*/
	if (myrank != nprocs-1 && nprocs>1)
		Imax = Imax - 1;
	else if (nprocs == 1)
		Imax = Imax + 1;

	if (waveForm == -1)
		power = loadSRC(pathSRC);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  PML parameters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	sig_x_max = 0.75 * (0.8 * (m + 1) / (dx * sqrt(muO / (epsO * epsR))));
	sig_y_max = 0.75 * (0.8 * (m + 1) / (dy * sqrt(muO / (epsO * epsR))));
	sig_z_max = 0.75 * (0.8 * (m + 1) / (dz * sqrt(muO / (epsO * epsR))));
	alpha_x_max = 0.24;
	alpha_y_max = alpha_x_max;
	alpha_z_max = alpha_x_max;
	kappa_x_max = 15;
	kappa_y_max = kappa_x_max;
	kappa_z_max = kappa_x_max;
	if (myrank == 0)
	{
		printf("Time step = %e\n", dt);
		printf("Number of steps = %d\n", nMax);
		printf("Size of Simulation region = %d X %d X %d\n", spaceX, spaceY, spaceZ);
		printf("Total Simulation time = %e Seconds\n", nMax * dt);
		fprintf(fp_log, "Time step = %e\n", dt);
		fprintf(fp_log, "Number of steps = %d\n", nMax);
		fprintf(fp_log, "Size of Simulation region = %d X %d X %d\n", spaceX, spaceY, spaceZ);
		fprintf(fp_log, "Total Simulation time = %e Seconds\n", nMax * dt);
		fflush(stdout);
	}
	return SUCCESS;
 }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//SET CPML PARAMETERS IN EACH DIRECTION
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void setUpCPML()
{
	int i, j, k, ii, jj, kk;

	if (abcNo == 1)
	{

	for(i = 0; i < nxPML_1; ++i)
	{
		sige_x_PML_1[i] = sig_x_max * pow(((nxPML_1 - 1 - i)
						/ (nxPML_1 - 1.0)), m);
		alphae_x_PML_1[i] = alpha_x_max * pow(((double)i
						  / (nxPML_1 - 1.0)), ma);
		kappae_x_PML_1[i] = 1.0 + (kappa_x_max - 1.0)
						  * pow((nxPML_1 - 1- i) / (nxPML_1 - 1.0), m);
		be_x_1[i] = exp(-(sige_x_PML_1[i] / kappae_x_PML_1[i]
				  + alphae_x_PML_1[i]) * dt / epsO);

		if ((sige_x_PML_1[i] == 0.0) &&
			(alphae_x_PML_1[i] == 0.0) && (i == nxPML_1 - 1))
		{
			ce_x_1[i] = 0.0;

		}
		else
		{
			ce_x_1[i] = sige_x_PML_1[i] * (be_x_1[i] - 1.0)
					  / (sige_x_PML_1[i] + kappae_x_PML_1[i] * alphae_x_PML_1[i])
					  / kappae_x_PML_1[i];
		}
	}

	for(i = 0; i < nxPML_1 - 1; ++i)
	{
		sigh_x_PML_1[i] = sig_x_max * pow(((nxPML_1 - 1 - i - 0.5)
						/ (nxPML_1-1.0)), m);
		alphah_x_PML_1[i] = alpha_x_max * pow(((i + 1 -0.5)
						  / (nxPML_1-1.0)), ma);
		kappah_x_PML_1[i] = 1.0 + (kappa_x_max - 1.0)
						  * pow(((nxPML_1 - 1 - i - 0.5) / (nxPML_1 - 1.0)), m);
		bh_x_1[i] = exp(-(sigh_x_PML_1[i] / kappah_x_PML_1[i]
				  + alphah_x_PML_1[i]) * dt / epsO);
		ch_x_1[i] = sigh_x_PML_1[i] * (bh_x_1[i] - 1.0)
				  / (sigh_x_PML_1[i] + kappah_x_PML_1[i] * alphah_x_PML_1[i])
				  / kappah_x_PML_1[i];
	}

	for(i = 0; i < nxPML_2; ++i)
	{
		sige_x_PML_2[i] = sig_x_max * pow(((nxPML_2 - 1 - i)
						/ (nxPML_2 - 1.0)), m);
		alphae_x_PML_2[i] = alpha_x_max * pow(((double)i
						  / (nxPML_2 - 1.0)), ma);
		kappae_x_PML_2[i] = 1.0 + (kappa_x_max - 1.0)
						  * pow((nxPML_2 - 1- i) / (nxPML_2 - 1.0), m);
		be_x_2[i] = exp(-(sige_x_PML_2[i] / kappae_x_PML_2[i]
				  + alphae_x_PML_2[i]) * dt / epsO);

		if ((sige_x_PML_2[i] == 0.0) &&
			(alphae_x_PML_2[i] == 0.0) && (i == nxPML_2 - 1))
		{
			ce_x_2[i] = 0.0;

		}
		else
		{
			ce_x_2[i] = sige_x_PML_2[i] * (be_x_2[i] - 1.0)
					  / (sige_x_PML_2[i] + kappae_x_PML_2[i] * alphae_x_PML_2[i])
					  / kappae_x_PML_2[i];
		}

	}

	for(i = 0; i < nxPML_2 - 1; ++i)
	{
		sigh_x_PML_2[i] = sig_x_max * pow(((nxPML_2 - 1 - i - 0.5)
						/ (nxPML_2-1.0)), m);
		alphah_x_PML_2[i] = alpha_x_max * pow(((i + 1 -0.5)
						  / (nxPML_2-1.0)), ma);
		kappah_x_PML_2[i] = 1.0 + (kappa_x_max - 1.0)
						  * pow(((nxPML_2 - 1 - i - 0.5) / (nxPML_2 - 1.0)), m);
		bh_x_2[i] = exp(-(sigh_x_PML_2[i] / kappah_x_PML_2[i]
				  + alphah_x_PML_2[i]) * dt / epsO);
		ch_x_2[i] = sigh_x_PML_2[i] * (bh_x_2[i] - 1.0)
				  / (sigh_x_PML_2[i] + kappah_x_PML_2[i] * alphah_x_PML_2[i])
				  / kappah_x_PML_2[i];
	}

	for(j = 0; j < nyPML_1; ++j)
	{
		sige_y_PML_1[j] = sig_y_max * pow(((nyPML_1 - 1 - j)
						/ (nyPML_1 - 1.0)), m);
		alphae_y_PML_1[j] = alpha_y_max * pow(((double)j
						  / (nyPML_1 - 1.0)), ma);
		kappae_y_PML_1[j] = 1.0 + (kappa_y_max - 1.0)
						  * pow((nyPML_1 - 1- j) / (nyPML_1 - 1.0), m);
		be_y_1[j] = exp(-(sige_y_PML_1[j] / kappae_y_PML_1[j]
				  + alphae_y_PML_1[j]) * dt / epsO);

		if ((sige_y_PML_1[j] == 0.0) &&
			 (alphae_y_PML_1[j] == 0.0) && (j == nyPML_1 - 1))
		{
			ce_y_1[j] = 0.0;

		}
		else
		{
			ce_y_1[j] = sige_y_PML_1[j] * (be_y_1[j] - 1.0)
					  / (sige_y_PML_1[j] + kappae_y_PML_1[j] * alphae_y_PML_1[j])
					  / kappae_y_PML_1[j];
		}
	}

	for(j = 0; j < nyPML_1 - 1; ++j)
	{
		sigh_y_PML_1[j] = sig_y_max * pow(((nyPML_1 - 1 - j - 0.5)
						/ (nyPML_1-1.0)), m);
		alphah_y_PML_1[j] = alpha_y_max * pow(((j + 1 -0.5)
						  / (nyPML_1-1.0)), ma);
		kappah_y_PML_1[j] = 1.0 + (kappa_y_max - 1.0)
						  * pow(((nyPML_1 - 1 - j - 0.5) / (nyPML_1 - 1.0)), m);
		bh_y_1[j] = exp(-(sigh_y_PML_1[j] / kappah_y_PML_1[j]
				  + alphah_y_PML_1[j]) * dt / epsO);
		ch_y_1[j] = sigh_y_PML_1[j] * (bh_y_1[j] - 1.0)
				  / (sigh_y_PML_1[j] + kappah_y_PML_1[j] * alphah_y_PML_1[j])
				  / kappah_y_PML_1[j];
	}

	for(j = 0; j < nyPML_2; ++j)
	{
		sige_y_PML_2[j] = sig_y_max * pow(((nyPML_2 - 1 - j)
						/ (nyPML_2 - 1.0)), m);
		alphae_y_PML_2[j] = alpha_y_max * pow(((double)j
						  / (nyPML_2 - 1.0)), ma);
		kappae_y_PML_2[j] = 1.0 + (kappa_y_max - 1.0)
						  * pow((nyPML_2 - 1- j) / (nyPML_2 - 1.0), m);
		be_y_2[j] = exp(-(sige_y_PML_2[j] / kappae_y_PML_2[j]
				  + alphae_y_PML_2[j]) * dt / epsO);

		if ((sige_y_PML_2[j] == 0.0) &&
			(alphae_y_PML_2[j] == 0.0) && (j == nyPML_2 - 1))
		{
			ce_y_2[j] = 0.0;
		}
		else
		{
			ce_y_2[j] = sige_y_PML_2[j] * (be_y_2[j] - 1.0)
					  / (sige_y_PML_2[j] + kappae_y_PML_2[j] * alphae_y_PML_2[j])
					  / kappae_y_PML_2[j];
		}
	}

	for(j = 0; j < nyPML_2 - 1; ++j)
	{
		sigh_y_PML_2[j] = sig_y_max * pow(((nyPML_2 - 1 - j - 0.5)
						/ (nyPML_2-1.0)), m);
		alphah_y_PML_2[j] = alpha_y_max * pow(((j + 1 -0.5)
						  / (nyPML_2-1.0)), ma);
		kappah_y_PML_2[j] = 1.0 + (kappa_y_max - 1.0)
						  * pow(((nyPML_2 - 1 - j - 0.5) / (nyPML_2 - 1.0)), m);
		bh_y_2[j] = exp(-(sigh_y_PML_2[j] / kappah_y_PML_2[j]
				  + alphah_y_PML_2[j]) * dt / epsO);
		ch_y_2[j] = sigh_y_PML_2[j] * (bh_y_2[j] - 1.0)
				  / (sigh_y_PML_2[j] + kappah_y_PML_2[j] * alphah_y_PML_2[j])
				  / kappah_y_PML_2[j];
	}

	for(k = 0; k < nzPML_1; ++k)
	{
		sige_z_PML_1[k] = sig_z_max * pow(((nzPML_1 - 1 - k)
						/ (nzPML_1 - 1.0)), m);
		alphae_z_PML_1[k] = alpha_z_max * pow(((double)k
						  / (nzPML_1 - 1.0)), ma);
		kappae_z_PML_1[k] = 1.0 + (kappa_z_max - 1.0)
						  * pow((nzPML_1 - 1- k) / (nzPML_1 - 1.0), m);
		be_z_1[k] = exp(-(sige_z_PML_1[k] / kappae_z_PML_1[k]
				  + alphae_z_PML_1[k]) * dt / epsO);

		if ((sige_z_PML_1[k] == 0.0) &&
			(alphae_z_PML_1[k] == 0.0) && (k == nzPML_1 - 1))
		{
			ce_z_1[k] = 0.0;
		}
		else
		{
			ce_z_1[k] = sige_z_PML_1[k] * (be_z_1[k] - 1.0)
					  / (sige_z_PML_1[k] + kappae_z_PML_1[k] * alphae_z_PML_1[k])
					  / kappae_z_PML_1[k];
		}
	}

	for(k = 0; k < nzPML_1 - 1; ++k)
	{
		sigh_z_PML_1[k] = sig_z_max * pow(((nzPML_1 - 1 - k - 0.5)
						/ (nzPML_1-1.0)), m);
		alphah_z_PML_1[k] = alpha_z_max * pow(((k + 1 -0.5)
						  / (nzPML_1-1.0)), ma);
		kappah_z_PML_1[k] = 1.0 + (kappa_z_max - 1.0)
						  * pow(((nzPML_1 - 1 - k - 0.5) / (nzPML_1 - 1.0)), m);
		bh_z_1[k] = exp(-(sigh_z_PML_1[k] / kappah_z_PML_1[k]
				  + alphah_z_PML_1[k]) * dt / epsO);
		ch_z_1[k] = sigh_z_PML_1[k] * (bh_z_1[k] - 1.0)
				  / (sigh_z_PML_1[k] + kappah_z_PML_1[k] * alphah_z_PML_1[k])
				  / kappah_z_PML_1[k];
	}

	for(k = 0; k < nzPML_2; ++k)
	{
		sige_z_PML_2[k] = sig_z_max * pow(((nzPML_2 - 1 - k)
						/ (nzPML_2 - 1.0)), m);
		alphae_z_PML_2[k] = alpha_z_max * pow(((double)k
						  / (nzPML_2 - 1.0)), ma);
		kappae_z_PML_2[k] = 1.0 + (kappa_z_max - 1.0)
						  * pow((nzPML_2 - 1- k) / (nzPML_2 - 1.0), m);
		be_z_2[k] = exp(-(sige_z_PML_2[k] / kappae_z_PML_2[k]
				  + alphae_z_PML_2[k]) * dt / epsO);

		if ((sige_z_PML_2[k] == 0.0) &&
			(alphae_z_PML_2[k] == 0.0) && (k == nzPML_2 - 1))
		{
			ce_z_2[k] = 0.0;
		}
		else
		{
			ce_z_2[k] = sige_z_PML_2[k] * (be_z_2[k] - 1.0)
					  / (sige_z_PML_2[k] + kappae_z_PML_2[k] * alphae_z_PML_2[k])
					  / kappae_z_PML_2[k];
		}
	}

	for(k = 0; k < nzPML_2 - 1; ++k)
	{
		sigh_z_PML_2[k] = sig_z_max * pow(((nzPML_2 - 1 - k - 0.5)
						/ (nzPML_2-1.0)), m);
		alphah_z_PML_2[k] = alpha_z_max * pow(((k + 1 -0.5)
						  / (nzPML_2-1.0)), ma);
		kappah_z_PML_2[k] = 1.0 + (kappa_z_max - 1.0)
						  * pow(((nzPML_2 - 1 - k - 0.5) / (nzPML_2 - 1.0)), m);
		bh_z_2[k] = exp(-(sigh_z_PML_2[k] / kappah_z_PML_2[k]
				  + alphah_z_PML_2[k]) * dt / epsO);
		ch_z_2[k] = sigh_z_PML_2[k] * (bh_z_2[k] - 1.0)
				  / (sigh_z_PML_2[k] + kappah_z_PML_2[k] * alphah_z_PML_2[k])
				  / kappah_z_PML_2[k];
	}

	}/* End if (abcNo == 1) */

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  DENOMINATORS FOR FIELD UPDATES
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	if (myrank != nprocs-1 && nprocs > 1)
		Imax = Imax+1;
	else if (nprocs == 1)
		Imax = Imax-1;

	ii = nxPML_2 - 2;

	for(i = 0; i < Imax; ++i)
	{
		if (i < nxPML_1-1 && myrank == 0)
		{
			den_hx[i] = 1.0 / (kappah_x_PML_1[i] * dx);
		}
		else if (i >= Imax - nxPML_2 + (!(nprocs-1)) && myrank == nprocs - 1)
		{
			den_hx[i] = 1.0 / (kappah_x_PML_2[ii] * dx);
			ii = ii - 1;
		}
		else
		{
			den_hx[i] = 1.0 / dx;
		}
	}

	jj = nyPML_2 - 2;

	for(j = 0; j < Jmax-1; ++j)
	{
		if (j < nyPML_1-1)
		{
			den_hy[j] = 1.0 / (kappah_y_PML_1[j] * dy);
		}
		else if (j >= Jmax - nyPML_2)
		{
			den_hy[j] = 1.0 / (kappah_y_PML_2[jj] * dy);
			jj = jj - 1;
		}
		else
		{
			den_hy[j] = 1.0 / dy;
		}
	}

	kk = nzPML_2 - 2;

	for (k = 1; k < Kmax - 1; ++k)
	{
		if (k < nzPML_1)
		{
			den_hz[k] = 1.0 / (kappah_z_PML_1[k-1] * dz);
		}
		else if (k >= Kmax - nzPML_2)
		{
			den_hz[k] = 1.0 / (kappah_z_PML_2[kk] * dz);
			kk = kk - 1;
		}
		else
		{
			den_hz[k] = 1.0 / dz;
		}
	}

	ii = nxPML_2 - 1;

	for(i = 0; i < Imax; ++i)
	{
		if (i < nxPML_1 && myrank == 0)
		{
		 den_ex[i] = 1.0 / (kappae_x_PML_1[i] * dx);
		}
		else if (i >= Imax - nxPML_2 + (!(nprocs-1)) && myrank == nprocs - 1)
		{
			den_ex[i] = 1.0 / (kappae_x_PML_2[ii] * dx);
			ii = ii - 1;
		}
		else
		{
			den_ex[i] = 1.0 / dx;
		}
	}

	jj = nyPML_2 - 1;

	for(j = 0; j < Jmax - 1; ++j)
	{
		if (j < nyPML_1)
		{
			den_ey[j] = 1.0 / (kappae_y_PML_1[j] * dy);
		}
		else if (j >= Jmax - nyPML_2)
		{
			den_ey[j] = 1.0 / (kappae_y_PML_2[jj] * dy);
			jj = jj - 1;
		}
		else
		{
			den_ey[j] = 1.0 / dy;
		}
	}

	kk = nzPML_2 - 1;

	for(k = 0; k < Kmax - 1; ++k)
	{
		if (k < nzPML_1 && abcNo == 1)
		{
			den_ez[k] = 1.0 / (kappae_z_PML_1[k] * dz);
		}
		else if (k >= Kmax - 1 - nzPML_2 && abcNo == 1)
		{
			den_ez[k] = 1.0 / (kappae_z_PML_2[kk] * dz);
			kk = kk - 1;
		}
		else
		{
			den_ez[k] = 1.0 / dz;
		}
	}
	if (myrank != nprocs-1 && nprocs>1)
		Imax = Imax-1;
	else if (nprocs == 1)
		Imax = Imax+1;

	free(alphae_x_PML_1);
	free(alphae_x_PML_2);
	free(alphae_y_PML_1);
	free(alphae_y_PML_2);
	free(alphae_z_PML_1);
	free(alphae_z_PML_2);
	free(alphah_x_PML_1);
	free(alphah_x_PML_2);
	free(alphah_y_PML_1);
	free(alphah_y_PML_2);
	free(alphah_z_PML_1);
	free(alphah_z_PML_2);

	free(kappae_x_PML_1);
	free(kappae_x_PML_2);
	free(kappae_y_PML_1);
	free(kappae_y_PML_2);
	free(kappae_z_PML_1);
	free(kappae_z_PML_2);
	free(kappah_x_PML_1);
	free(kappah_x_PML_2);
	free(kappah_y_PML_1);
	free(kappah_y_PML_2);
	free(kappah_z_PML_1);
	free(kappah_z_PML_2);

	free(sige_x_PML_1);
	free(sige_x_PML_2);
	free(sige_y_PML_1);
	free(sige_y_PML_2);
	free(sige_z_PML_1);
	free(sige_z_PML_2);
	free(sigh_x_PML_1);
	free(sigh_x_PML_2);
	free(sigh_y_PML_1);
	free(sigh_y_PML_2);
	free(sigh_z_PML_1);
	free(sigh_z_PML_2);
}

void yeeCubeCap (int I, int J, int K, unsigned char mediaNo)
{
    if (fabs(media[mediaNo].epsilon-1)<1e-6)
    	modelDataX[I][J][K-1] = mediaNo;
    else if (fabs(media[mediaNo].epsilon-2)<1e-6)
	    modelDataY[I][J][K-1] = mediaNo;
    else if (fabs(media[mediaNo].epsilon-3)<1e-6)
	    modelDataZ[I][J][K] = mediaNo;
}

void yeeCube (int I, int J, int K, unsigned char mediaNo)
{
	K=K-1;
	//set face 1 (for EX)
	modelDataX[I][J][K] = mediaNo;
	modelDataX[I][J][K + 1] = mediaNo;
	modelDataX[I][J + 1][K + 1] = mediaNo;
	modelDataX[I][J + 1][K] = mediaNo;

	//set face 2 (for EY)
	modelDataY[I][J][K] = mediaNo;
	modelDataY[I + 1][J][K] = mediaNo;
	modelDataY[I + 1][J][K + 1] = mediaNo;
	modelDataY[I][J][K + 1] = mediaNo;

	K=K+1;
	//set face 3 (for EZ)
	modelDataZ[I][J][K] = mediaNo;
	modelDataZ[I + 1][J][K] = mediaNo;
	modelDataZ[I + 1][J + 1][K] = mediaNo;
	modelDataZ[I][J + 1][K] = mediaNo;
}

/*****************************************************************************************/
/*
 * Function name: loadSRC
 * Description: Load *.src Files(xFDTD)
 * Parameters: 
 *			 path
 * Return: 
 *			 Point to Array, Save Row 2--end in *.src
 */
double* loadSRC(char* path)
{
	int ldi;
	double temp, *A;
	FILE *fp;
	if (!(fp = fopen(path,"r"))/* && myrank == 0*/)
	{
		printf("Open data file ");
		printf(path);
		printf(" fail.\n");
		if (myrank == 0)
        {
            fprintf(fp_log, "Open data file ");
		    fprintf(fp_log, path);
		    fprintf(fp_log, " fail.\n");
        }
        fflush(stdout);
		return (NULL);
	}
	fscanf(fp,"%d\n", &sizeSRC);
	A = (double *)malloc(sizeSRC * sizeof(double));
    mem_count += sizeSRC*sizeof(double);
	for (ldi = 0; ldi < sizeSRC; ldi++)
	{
		fscanf(fp, "%lf\t%lf\n", &A[ldi], &temp);
	}

	fclose(fp);
	return (A);
}

/*****************************************************************************************/
/*
 * Function name: loadData3
 * Description: Load a 3D Array
 * Parameters: 
 *			 path
 *			 size
 *			 defauleValue
 * Return: 
 *			 Point to 3D Array
 */
unsigned char*** loadData3(char* path, int* size)
{
    int i;
	int ldi, ldj, ldk;
	int schedule = 0, schedule_past = 0; 
	unsigned char*** A;
	int temp_data;
	unsigned char temp_c;
	FILE *fp;
	time_t load_time;

	A = initArray3Char(size[0] + paddingX_1 + paddingX_2 + 1,
						size[1] + paddingY_1 + paddingY_2 + 1,
						size[2] + paddingZ_1 + paddingZ_2 + 1);

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
		printf(" fail in CPU%d.\n", myrank);
		if (myrank == 0)
        {
            fprintf(fp_log, "Can't open model file in %s\n", path);
        }
        fflush(stdout);
		return (A);
	}

	load_time = clock();

	for (ldk = paddingZ_1; ldk < size[2] + paddingZ_1; ldk++)
	{
		for (ldj = paddingY_1; ldj < size[1] + paddingY_1; ldj++)
		{
			for (ldi = paddingX_1; ldi < size[0] + paddingX_1 - 1; ldi++)
			{
				fscanf(fp, "%d ", &temp_data);
				A[ldi][ldj][ldk] = temp_data;
			}
			fscanf(fp, "%d", &temp_data);
			A[ldi][ldj][ldk] = temp_data;
			fscanf(fp, "%c", &temp_c);
			if (temp_c != '\n')
			{/*
				printf("Model file error in (%d, %d, %d), end with %d char\n", ldi-paddingX_1, ldj-paddingY_1, ldk-paddingZ_1, (int)temp_c);
				throw MEM_ERROR;
			*/}
		}
		schedule_past = schedule;
		schedule = (int)((ldk-paddingZ_1)/(double)(size[2]-1)*20);
		if(clock() - load_time > 2000 || schedule-schedule_past)
		{
			if (myrank == 0)
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
                printf("| %4.1f%%", (ldk-paddingZ_1)/(double)(size[2]-1)*100);
				fflush(stdout);
			}
			load_time = clock();
		}
	}

	if (myrank == 0)
	{
		printf("\n");
		fprintf(fp_log, "Loading model from %s\n", path);
		fflush(stdout);
	}
	fclose(fp);
	return (A);
}

int loadMediaData(char* path, int size)
{
	int i;
	int id;
	double eps, sig, rho, spec_heat, K, B;
	FILE* fp = fopen(path, "r");
	if (fp == NULL)	
	{
		printf("Open data file %s fail in CPU%d\n", path, myrank);
		if (myrank == 0)
	        fprintf(fp_log, "Can't open media file in %s\n", path);
		fflush(stdout);
		return FILE_ERROR;
	}
	for (i = 0; i<size; ++i)
	{
		fscanf(fp, "%d %lf %lf %lf", &id, &sig, &eps, &rho);
		if (id>=maxMedia) maxMedia = id+1;// BUILDOBJECTS\createNewMedia函数需要的数据
        if (need_compute_temperature_rise)
        {
            fscanf(fp, " %lf %lf %lf", &spec_heat, &K, &B);
        }
        fscanf(fp, "\n");
        if (id == 0 || id == 1)
        {
            printf("Media id cannot be 0(vaccum) or 1(pec), quit.\n");
            fflush(stdout);
            fprintf(fp_log, "Media id cannot be 0(vaccum) or 1(pec), quit.\n");
            return MEM_ERROR;
        }
        if (id == 201) // 电容器，数据格式为："id 电容值 极化方向 0\n"
        {
            media[i].sigma = sig;// 电容值C
		    media[id].epsilon = eps;// 极化方向1——3 ——> x——z
		    media[id].rho= rho;
            CA[id] = 1.0;
            double tempPolar = 0.0;
            if (fabs(eps-1)<1e-6)
                tempPolar = dx/dy/dz;
            else if (fabs(eps-2)<1e-6)
                tempPolar = dy/dx/dz;
            else if (fabs(eps-3)<1e-6)
                tempPolar = dz/dx/dy;
            else
            {
                printf("电容器极性设置错误，请重新设置（1-->x, 2-->y, 3-->z），当前设置为%lf\n", eps);
                fflush(stdout);
                return MEM_ERROR;
            }
            CA[id] = 1.0;
            CB[id] = dt/epsO/(1.0 + sig*tempPolar/epsO);
        }
		media[id].sigma= sig;
		media[id].epsilon = eps;
		media[id].rho= rho;
        if (need_compute_temperature_rise)
        {
			media[id].spec_heat = spec_heat;
            media[id].K = K;
            media[id].B = B;
        }
		CA[id] = ((1.0 - sig * dt / (2.0 * eps * epsO))
				/ (1.0 + sig * dt / (2.0 * eps * epsO)));
		CB[id] = (dt / (eps * epsO))
				/ (1.0 + sig * dt / (2.0 * eps * epsO));
	}
	if (maxMedia > MAX_NUM_OF_MEDIA && myrank == 0)
	{
		printf("Too many media!\n");
		fprintf(fp_log, "Too many media!\n");
		fflush(stdout);
		return MEM_ERROR;
	}
    if (need_compute_temperature_rise)
	{
		media[0].rho = 1.4128;
		media[0].spec_heat = 1012;
		media[0].K = 0.023;
		media[0].B = 0.0;
	}

	return SUCCESS;
}

/*****************************************************************************************/
/*
 * Function name: loadModel
 * Description: Load Model data to Yee Grids.
 * Parameters: 
 * Return: 
 */
void loadModel()
{
	int i, j, k;
	int change_by_is = 0;
	if (myrank == 0)
		change_by_is = -is;
	else
		change_by_is = -is + 1;
	if (myrank != nprocs - 1)
		ie = ie+1;
	if (myrank != 0)
		is = is-1;
	for(i = nxPML_1; i < nxPML_1+paddingX_1+_spaceX+paddingX_2; i++)
	{
		if(is<=i && i<=ie)
		{
			for(j = nyPML_1; j < Jmax - nyPML_2; j++)
			{
				for(k = nzPML_1; k < Kmax - nzPML_2; k++)
				{
					if (modelData[i - nxPML_1][j - nyPML_1][k - nzPML_1] == 0) /* vaccum */
						continue;
                    if (modelData[i - nxPML_1][j - nyPML_1][k - nzPML_1] == 201) // 电容器
                    {
                        yeeCubeCap(i+change_by_is, j, k, modelData[i - nxPML_1][j - nyPML_1][k - nzPML_1]);
                    }
                    else
					    yeeCube(i+change_by_is, j, k, modelData[i - nxPML_1][j - nyPML_1][k - nzPML_1]);
				}
			}
		}
	}
	if (myrank != nprocs - 1)
		ie = ie-1;
	if (myrank != 0)
		is = is+1;
}
void loadModelOneCPU()
{
	int i, j, k;
	for(i = nxPML_1; i < Imax - nxPML_2; i++)
	{

		for(j = nyPML_1; j < Jmax - nyPML_2; j++)
		{

			for(k = nzPML_1; k < Kmax - nzPML_2; k++)
			{
				if (modelData[i - nxPML_1][j - nyPML_1][k - nzPML_1] == 0)
					continue;
				yeeCube(i, j, k, modelData[i - nxPML_1][j - nyPML_1][k - nzPML_1]); 
			}
		}
	}
}

