void setUp()
{
	//Time step
	dt = 1.0 / (C * sqrt(1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz)));

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  COMPUTING FIELD UPDATE EQUATION COEFFICIENTS
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	DA = 1.0;
	DB = dt / muO;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  LOAD MODEL AND COMPUTE CA,CB
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	size[0] = _spaceX;
	size[1] = _spaceY;
	size[2] = _spaceZ;

	//FILE to Load
	strcpy(path_eps, path_data);
	strcat(path_eps, data_name_eps);
	strcpy(path_sigma, path_data);
	strcat(path_sigma, data_name_sigma);
	strcpy(path_rho, path_data);
	strcat(path_rho, data_name_rho);
	
	epsilon = loadData3(path_eps, size, 1.0);
	sigma = loadData3(path_sigma, size, 0.0);
#ifdef _SAR
	rho = loadData3(path_rho, size, 0.0);
#endif

	if (nprocs == 1)
		loadModelOneCPU();
	else
		loadModel();
	freeArray3(epsilon, spaceX, spaceY, spaceZ);
	freeArray3(sigma, spaceX, spaceY, spaceZ);

	buildObject();
	loadCACB();

	if (myrank != nprocs-1 && nprocs>1)
		Imax = Imax + 1;
	else if (nprocs == 1)
		Imax = Imax - 1;
	freeArray3(ID1_EpsR, Imax+1, Jmax, Kmax);
	freeArray3(ID2_EpsR, Imax+1, Jmax, Kmax);
	freeArray3(ID3_EpsR, Imax+1, Jmax, Kmax);
#ifndef _SAR
	freeArray3(ID1_SigR, Imax+1, Jmax, Kmax);
	freeArray3(ID2_SigR, Imax+1, Jmax, Kmax);
	freeArray3(ID3_SigR, Imax+1, Jmax, Kmax);
#endif
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
		printf("Total Memory requirements = %ld KB\n", mem_count / 1024);
		fflush(stdout);
	}
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

void loadCACB()
{
	int i, j, k;
	double epsr, sigr;
	if(myrank != nprocs-1 && nprocs >1)
		Imax = Imax+1;
	else if (nprocs == 1)
		Imax = Imax-1;//i range from 0(1) to Imax-1
	for(i = 0; i < Imax; ++i)
	{
		for(j = 1; j < Jmax-1; ++j)
		{
			for(k = 0; k < Kmax-1; ++k)
			{
					epsr = ID1_EpsR[i][j][k];
					sigr = ID1_SigR[i][j][k];

					if (epsr < -0.1)
					{
						CAx[i][j][k] = 0;
						CBx[i][j][k] = 0;
						continue;
					}
				CAx[i][j][k] = (1.0 - sigr * dt / (2.0 * epsr * epsO)) /
								(1.0 + sigr * dt / (2.0 * epsr * epsO));
				CBx[i][j][k] = (dt / (epsr * epsO)) /
								(1.0 + sigr * dt / (2.0 * epsr * epsO));
			}
		}
	}

	for(k = 0; k < Kmax-1; ++k)
	{
		for(i = 1; i < Imax; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
					epsr = ID2_EpsR[i][j][k];
					sigr = ID2_SigR[i][j][k];

					if (epsr < -0.1)
					{
						CAy[i][j][k] = 0;
						CBy[i][j][k] = 0;
						continue;
					}
				CAy[i][j][k] = (1.0 - sigr * dt / (2.0 * epsr * epsO)) /
								(1.0 + sigr * dt / (2.0 * epsr * epsO));
				CBy[i][j][k] = (dt / (epsr * epsO)) /
								(1.0 + sigr * dt / (2.0 * epsr * epsO));
			}
		}
	}
	for(k = 1; k < Kmax-1; ++k)
	{
		for(i = 1; i < Imax; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
					epsr = ID3_EpsR[i][j][k];
					sigr = ID3_SigR[i][j][k];

					if (epsr < -0.1)
					{
						CAz[i][j][k] = 0;
						CBz[i][j][k] = 0;
						continue;
					}
				CAz[i][j][k] = (1.0 - sigr * dt / (2.0 * epsr * epsO)) /
								(1.0 + sigr * dt / (2.0 * epsr * epsO));
				CBz[i][j][k] = (dt / (epsr * epsO)) /
								(1.0 + sigr * dt / (2.0 * epsr * epsO));
			}
		}
	}
	if(myrank != nprocs-1 && nprocs>1)
		Imax = Imax-1;
	else if (nprocs == 1)
		Imax = Imax+1;
}

void yeeCube (int I, int J, int K, double EpsR, double SigR, double RhoR)
{
	K=K-1;
	//set face 1 (for EX)
	ID1_EpsR[I][J][K] = EpsR;
	ID1_EpsR[I][J][K + 1] = EpsR;
	ID1_EpsR[I][J + 1][K + 1] = EpsR;
	ID1_EpsR[I][J + 1][K] = EpsR;
	ID1_SigR[I][J][K] = SigR;
	ID1_SigR[I][J][K + 1] = SigR;
	ID1_SigR[I][J + 1][K + 1] = SigR;
	ID1_SigR[I][J + 1][K] = SigR;

	//set face 2 (for EY)
	ID2_EpsR[I][J][K] = EpsR;
	ID2_EpsR[I + 1][J][K] = EpsR;
	ID2_EpsR[I + 1][J][K + 1] = EpsR;
	ID2_EpsR[I][J][K + 1] = EpsR;
	ID2_SigR[I][J][K] = SigR;
	ID2_SigR[I + 1][J][K] = SigR;
	ID2_SigR[I + 1][J][K + 1] = SigR;
	ID2_SigR[I][J][K + 1] = SigR;

	K=K+1;
	//set face 3 (for EZ)
	ID3_EpsR[I][J][K] = EpsR;
	ID3_EpsR[I + 1][J][K] = EpsR;
	ID3_EpsR[I + 1][J + 1][K] = EpsR;
	ID3_EpsR[I][J + 1][K] = EpsR;
	ID3_SigR[I][J][K] = SigR;
	ID3_SigR[I + 1][J][K] = SigR;
	ID3_SigR[I + 1][J + 1][K] = SigR;
	ID3_SigR[I][J + 1][K] = SigR;

#ifdef _SAR
	K = K-1;
	ID1_RhoR[I][J][K] = RhoR;
	ID1_RhoR[I][J][K + 1] = RhoR;
	ID1_RhoR[I][J + 1][K + 1] = RhoR;
	ID1_RhoR[I][J + 1][K] = RhoR;

	ID2_RhoR[I][J][K] = RhoR;
	ID2_RhoR[I + 1][J][K] = RhoR;
	ID2_RhoR[I + 1][J][K + 1] = RhoR;
	ID2_RhoR[I][J][K + 1] = RhoR;

	K = K+1;
	ID3_RhoR[I][J][K] = RhoR;
	ID3_RhoR[I + 1][J][K] = RhoR;
	ID3_RhoR[I + 1][J + 1][K] = RhoR;
	ID3_RhoR[I][J + 1][K] = RhoR;
#endif
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
		fflush(stdout);
		return (NULL);
	}
	fscanf(fp,"%d\n", &sizeSRC);
	A = (double *)malloc(sizeSRC * sizeof(double));
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
double*** loadData3(char* path, int* size, double defauleValue)
{
	int ldi, ldj, ldk;
	double*** A;
	FILE *fp;
	time_t load_time;

	A = initArray3(size[0] + paddingX_1 + paddingX_2 + 1,
				   size[1] + paddingY_1 + paddingY_2 + 1,
				   size[2] + paddingZ_1 + paddingZ_2 + 1,
				   defauleValue);

	if (path == NULL)
		return (A);

	if (!(fp = fopen(path,"r")))
	{
		printf("Open data file ");
		printf(path);
		printf(" fail in CPU%d.\n", myrank);
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
				fscanf(fp, "%lf ", &A[ldi][ldj][ldk]);
			}
			fscanf(fp, "%lf\n", &A[ldi][ldj][ldk]);
		}
		if(clock() - load_time > 2000)
		{
			if (myrank == 0)
			{
				printf("\rLoading %s : %4.2f%%", path, (ldk-paddingZ_1)/(double)(size[2]-1)*100);
				fflush(stdout);
			}
			load_time = clock();
		}
	}

	if (myrank == 0)
	{
		printf("\rLoading %s : 100.00%%\n", path);
		fflush(stdout);
	}
	fclose(fp);
	return (A);
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
					if (fabs(epsilon[i - nxPML_1][j - nyPML_1][k - nzPML_1] - 1) < 1e-6 &&
						fabs(sigma[i - nxPML_1][j - nyPML_1][k - nzPML_1]) < 1e-6)
						continue;
						yeeCube(i+change_by_is, j, k, epsilon[i - nxPML_1][j - nyPML_1][k - nzPML_1], 
													  sigma[i - nxPML_1][j - nyPML_1][k - nzPML_1],
													  rho[i - nxPML_1][j - nyPML_1][k - nzPML_1]);
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
				if (fabs(epsilon[i - nxPML_1][j - nyPML_1][k - nzPML_1] - 1) < 1e-6 && 
					fabs(sigma[i - nxPML_1][j - nyPML_1][k - nzPML_1]) < 1e-6)
					continue;
				yeeCube(i, j, k, epsilon[i - nxPML_1][j - nyPML_1][k - nzPML_1], 
								 sigma[i - nxPML_1][j - nyPML_1][k - nzPML_1],
								 rho[i - nxPML_1][j - nyPML_1][k - nzPML_1]);
			}
		}
	}
}

