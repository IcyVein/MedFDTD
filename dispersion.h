/*****************************************************************************************/
/*
 * Function name: compute
 * Description: The main compute function
 * Parameters: 
 * Return
 */
void computeDispersion()
{
	int i, n;
	double tStart, tCurrent;
	double tStepStart, tStepEnd;
	double tStep=0.0;
	double tNear20[20] = {0.0};
	int dn = 10;
	int is_have_powersource = 0;
    char path_check[MAX_SIZE_OF_PATH];
    sprintf(path_check, "%scheckRMS_%d.txt", path_save, myrank);
    checkData = 0, checkDataPast = 0;

	/* Find the feed point */
    for (int nSrc = 0; nSrc<sourceNum; ++nSrc) {
	if (is <= isource[nSrc] && isource[nSrc] <= ie)
	{
		is_have_powersource = 1;
		int powersource_rank = myrank;
		if (myrank != 0)
			isource[nSrc] = isource[nSrc] - is + 1;
	}
	else if (isource[nSrc] == ie+1 && (port == 'y' || port == 'z'))
	{
		is_have_powersource = 1;
		int powersource_rank = myrank;
		if (myrank != 0)
			isource[nSrc] = isource[nSrc] - is + 1;
	}
	else;
    }

	if (myrank == 0)
	{
		is = 0;
		ie = Imax - 2;
	}
	else
	{
		is = 1;
		ie = Imax-1;
	}

#ifdef _POYNTING
	double poyingt = 0;
	double reducePoyingt = 0;
	double aveRadiationPower = 0;
	FILE* fp_POYNTING;

	if (myrank == 0)
	{
		char path_p[MAX_SIZE_OF_PATH];
		sprintf(path_p, "%spoynting.txt ", path_save);
		fp_POYNTING = fopen(path_p,"w+");
	}

	rad_region.xStart = 2 + paddingX_1 + thicknessOfPml;
	rad_region.xEnd = -2 +_spaceX + paddingX_1 + thicknessOfPml;
	rad_region.yStart = 2 + paddingY_1 + thicknessOfPml;
	rad_region.yEnd = -2 +_spaceY + paddingY_1 + thicknessOfPml;
	rad_region.zStart = 2 + paddingZ_1 + thicknessOfPml;
	rad_region.zEnd = -2 +_spaceZ + paddingZ_1 + thicknessOfPml;
	rad_region.computeX_1 = 1;
	rad_region.computeX_2 = 1;

	if(_global_is <= rad_region.xStart && rad_region.xEnd <= _global_ie )
	{
		if (myrank != 0)
		{
			rad_region.xStart -= _global_is-1;
			rad_region.xEnd -= _global_is-1;
		}
	}
	else if(_global_is <= rad_region.xStart && rad_region.xStart <= _global_ie && _global_ie <= rad_region.xEnd)
	{
		if (myrank != 0)
		{
			rad_region.xStart -= _global_is-1;
			
		}
		rad_region.xEnd = ie;
		rad_region.computeX_2 = 0;
	}
	else if(rad_region.xStart <= _global_is && _global_is <= rad_region.xEnd && _global_ie <= rad_region.xEnd)
	{
		if (myrank != 0)
		{
			rad_region.xStart = is;
			rad_region.xEnd = ie;
			rad_region.computeX_1 = 0;
			rad_region.computeX_2 = 0;
		}
	}
	else if(rad_region.xStart <= _global_is && _global_is <= rad_region.xEnd && rad_region.xEnd <= _global_ie)
	{
		if (myrank != 0)
		{
			rad_region.xStart = is;
			rad_region.xEnd -= _global_is-1;
			rad_region.computeX_1 = 0;
		}
	}
	else
	{
		rad_region.xStart = -1;
		rad_region.xEnd = -1;
		rad_region.computeX_1 = 0;
		rad_region.computeX_2 = 0;
	}
#endif

    float plane_wave_power = 0.0;

	MPI_Request request_sendy_r, request_sendy_l;
	MPI_Request request_sendz_r, request_sendz_l;
	MPI_Request request_recvy_r, request_recvy_l;
	MPI_Request request_recvz_r, request_recvz_l;

	MPI_Status status_recvy_r, status_recvy_l;
	MPI_Status status_recvz_r, status_recvz_l;
	MPI_Status status_sendy_r, status_sendy_l;
	MPI_Status status_sendz_r, status_sendz_l;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  BEGIN TIME STEP
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if (myrank == 0)
	{
		printf("Begin time-stepping...\n");
		fflush(stdout);
	}
	tStart = MPI_Wtime();
	time_t t[8] = {0}, t0;
	for(n = 1; n <= nMax; ++n)
	{
		tStepStart = MPI_Wtime();
		if(n % dn == 0 && myrank == 0)
		{
			if(n <= 20)
            {
				printf("Step: %d, %fsec/step, %d/%d sec left/total.\n",
						n, tStep * 20/n, (int)((nMax-n+1) * tStep) * 20/n, (int)((nMax-n+1) * tStep - tStart + tCurrent) * 20/n );
            }
			else
            {
				printf("Step: %d, %fsec/step, %d/%d sec left/total.\n",
						n, tStep, (int)((nMax-n+1) * tStep), (int)((nMax-n+1) * tStep - tStart + tCurrent));
            }
			fflush(stdout);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		t[6] += time(NULL)-t0;

	//-----------------------------------------------------------

		/* H */
		t0 = time(NULL);
		if (0<myrank && myrank<nprocs-1)
			computeFieldH_Dispersion();
		else if (myrank == 0)
			computeFieldH_0_Dispersion();
		else /* myrank == nprocs - 1 */
			computeFieldH_nprocsSub1_Dispersion();
		t[0] += time(NULL)-t0;

		/* H-PML */
		t0 = time(NULL);
		if (0<myrank && myrank<nprocs-1)
			computePMLH_Dispersion();
		else if (myrank == 0)
			computePMLH_0_Dispersion();
		else /* myrank == nprocs - 1 */
			computePMLH_nprocsSub1_Dispersion();
		t[1] += time(NULL)-t0;

        if( sourceType == 1/*If Plane Wave*/)
		{
			powerSourcePlaneWaveH(n);/* -1/2 - dx/2/C/dt */
		}

        /* Recv */
        if (0<myrank && myrank<nprocs-1)
		{
			MPI_Irecv(&Hy(ie+1,0,0), Jmax * Kmax, MPI_FLOAT, myrank+1, 300+myrank, MPI_COMM_WORLD, &request_recvy_r);
			MPI_Irecv(&Hy(is-1,0,0), Jmax * Kmax, MPI_FLOAT, myrank-1, 100+myrank, MPI_COMM_WORLD, &request_recvy_l);
			MPI_Irecv(&Hz(is-1,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank-1, 200+myrank, MPI_COMM_WORLD, &request_recvz_l);	
			MPI_Irecv(&Hz(ie+1,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank+1, 400+myrank, MPI_COMM_WORLD, &request_recvz_r);
		}
		else if (myrank == 0)
		{
			MPI_Irecv(&Hy(ie+1,0,0), Jmax * Kmax, MPI_FLOAT, myrank+1, 300+myrank, MPI_COMM_WORLD, &request_recvy_r);
			MPI_Irecv(&Hz(ie+1,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank+1, 400+myrank, MPI_COMM_WORLD, &request_recvz_r);
		}
		else /* myrank == nprocs - 1 */
		{
			MPI_Irecv(&Hy(is-1,0,0), Jmax * Kmax, MPI_FLOAT, myrank-1, 100+myrank, MPI_COMM_WORLD, &request_recvy_l);
			MPI_Irecv(&Hz(is-1,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank-1, 200+myrank, MPI_COMM_WORLD, &request_recvz_l);
		}
		t[2] += time(NULL)-t0;
		MPI_Barrier(MPI_COMM_WORLD);

		/* Send */
		t0 = time(NULL);
		if (0<myrank && myrank<nprocs-1)
		{
			MPI_Isend(&Hy(is,0,0), Jmax * Kmax, MPI_FLOAT, myrank-1, 300+myrank-1, MPI_COMM_WORLD, &request_sendy_l);
			MPI_Isend(&Hy(ie,0,0), Jmax * Kmax, MPI_FLOAT, myrank+1, 100+myrank+1, MPI_COMM_WORLD, &request_sendy_r);
			MPI_Isend(&Hz(is,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank-1, 400+myrank-1, MPI_COMM_WORLD, &request_sendz_l);
			MPI_Isend(&Hz(ie,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank+1, 200+myrank+1, MPI_COMM_WORLD, &request_sendz_r);
		}
		else if (myrank == 0)
		{
			MPI_Isend(&Hy(ie,0,0), Jmax * Kmax, MPI_FLOAT, myrank+1, 100+myrank+1, MPI_COMM_WORLD, &request_sendy_r);
			MPI_Isend(&Hz(ie,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank+1, 200+myrank+1, MPI_COMM_WORLD, &request_sendz_r);
		}
		else /* myrank == nprocs - 1 */
		{
			MPI_Isend(&Hy(is,0,0), Jmax * Kmax, MPI_FLOAT, myrank-1, 300+myrank-1, MPI_COMM_WORLD, &request_sendy_l);
			MPI_Isend(&Hz(is,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank-1, 400+myrank-1, MPI_COMM_WORLD, &request_sendz_l);
		}

		//MPI_Barrier(MPI_COMM_WORLD);

		/* E */
		t0 = time(NULL);
		if (0<myrank && myrank<nprocs-1)
		{
			MPI_Wait(&request_recvy_l, &status_recvy_l);
			MPI_Wait(&request_sendy_l, &status_sendy_l);

			MPI_Wait(&request_recvz_l, &status_recvz_l);
			MPI_Wait(&request_sendz_l, &status_sendz_l);

			computeFieldE_Dispersion();

			MPI_Wait(&request_recvy_r, &status_recvy_r);
			MPI_Wait(&request_sendy_r, &status_sendy_r);
			MPI_Wait(&request_recvz_r, &status_recvz_r);
			MPI_Wait(&request_sendz_r, &status_sendz_r);
			computeFieldE_right_Dispersion();
		}
		else if (myrank == 0)
		{
			MPI_Wait(&request_recvy_r, &status_recvy_r);
			MPI_Wait(&request_sendy_r, &status_sendy_r);
			MPI_Wait(&request_recvz_r, &status_recvz_r);
			MPI_Wait(&request_sendz_r, &status_sendz_r);
			computeFieldE_0_Dispersion();
		}
		else /* myrank == nprocs - 1 */
		{
			MPI_Wait(&request_recvy_l, &status_recvy_l);
			MPI_Wait(&request_sendy_l, &status_sendy_l);
			MPI_Wait(&request_recvz_l, &status_recvz_l);
			MPI_Wait(&request_sendz_l, &status_sendz_l);

			computeFieldE_nprocsSub1_Dispersion();
		}
		t[4] += time(NULL)-t0;

		/* E-PML */
		t0 = time(NULL);
		if (0<myrank && myrank<nprocs-1)
		{
			computePMLE_Dispersion();
		}
		else if (myrank == 0)
		{
			computePMLE_0_Dispersion();
		}
		else /* myrank == nprocs - 1 */
			computePMLE_nprocsSub1_Dispersion();
		t[5] += time(NULL)-t0;

		t0 = time(NULL);
		//MPI_Barrier(MPI_COMM_WORLD);
		t[6] += time(NULL)-t0;

	//-----------------------------------------------------------
	//   Apply a point source (Hard)
	//-----------------------------------------------------------

		t0 = time(NULL);
		if(is_have_powersource || sourceType == 1/*If Plane Wave*/)
		{
			powerSource(n);
		}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  WRITE TO OUTPUT FILES
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		saveData(n);
		MPI_Barrier(MPI_COMM_WORLD);
#ifdef _SAR
		for (i = 0; i<_spaceZ+1; ++i)
		{
			computeRMS(&pSAR[i], n);
		}
        if (n%(sourceT) == 0 && n+sourceT<=nMax)
        {
            checkDataPast = checkData;
            checkData = checkRMS();
            if (checkData != 0 && checkData != checkDataPast)
                convergence = 10.0*log10(fabs(checkData-checkDataPast)/checkData);
            else if (checkData == checkDataPast)
                convergence = convergenceTarget;
            else
                convergence = 0.0;
            resetRMS();
            fflush(stdout);
            MPI_Allreduce(&convergence, &convergenceCurr, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
            if (convergenceCurr <= convergenceTarget)
            {
                nMax = nMax < (n+sourceT) ? nMax:(n+sourceT); /* nMax = min(nMax, n+sourceT) */
                if (myrank == 0)
                {
                    printf("Convergence: %2.1f <= %2.1fdB, set max time step: %d.\n", convergenceCurr, convergenceTarget, nMax);
                    fprintf(fp_log, "Convergence: %2.1f <= %2.1fdB, set max time step: %d.\n", convergenceCurr, convergenceTarget, nMax);
                    fflush(stdout);
                }
            }
            else if (myrank == 0)
            {
                printf("Convergence: %2.1f / %2.1fdB.\n", convergenceCurr, convergenceTarget);
                fprintf(fp_log, "Convergence: %2.1f / %2.1fdB.\n", convergenceCurr, convergenceTarget);
                fflush(stdout);
            }
        }
        if (n == (nMax-sourceT))
            resetRMS();
#endif

#ifdef _POYNTING
		poyingt = radiationPower(rad_region);

		MPI_Reduce(&poyingt, &reducePoyingt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if (myrank == 0)
			fprintf(fp_POYNTING, "%lf\n", reducePoyingt);
		
		if(n>nMax-(int)(1/(freq*dt)) && n<=nMax && myrank == 0)
		{
			aveRadiationPower += reducePoyingt;
		}
#endif

		tCurrent = MPI_Wtime();
		tStepEnd = MPI_Wtime();
		tNear20[n%20] = tStepEnd - tStepStart;
		tStep = 0;
		for (i = 0; i<20; ++i)
			tStep += tNear20[i];
		tStep /= 20;
		t[7] += time(NULL)-t0;
	}//  END TIME STEP
	freeFDTDData();

	for (i = 0; i < save_plane_amount; i++)
	{
		fclose(fp_save_field_file[i].ex);
		fclose(fp_save_field_file[i].ey);
		fclose(fp_save_field_file[i].ez);
		fclose(fp_save_field_file[i].hx);
		fclose(fp_save_field_file[i].hy);
		fclose(fp_save_field_file[i].hz);
	}
	flag = 13; /* End FDTD Iteration */
#ifdef _SAR
	for (i = 0; i<save_localSAR_amount; ++i)
	{
		computeLocalSAR(pSAR[i], pSAR[i+1], &pSAR[i].localSARData);
	}

	for (i = 0; i<save_localSAR_amount; ++i)
	{
		//writeLocalSAR(pSAR[i].fp, pSAR[i].localSARData);
		//fclose(pSAR[i].fp);
	}
#endif

	freeArray3Char(modelDataX, Imax+1, Jmax, Kmax);
	freeArray3Char(modelDataY, Imax+1, Jmax, Kmax);
	freeArray3Char(modelDataZ, Imax+1, Jmax, Kmax);

#ifdef _SAR
	if (nXgSAR || whole_body_sar)
	{
		computeXgSAR(nXgSAR);
	}
#endif

#ifdef _POYNTING
	if (myrank == 0)
	{
		aveRadiationPower /= (int)(1/(freq*dt));
		printf("Averaged radiation power = %eW\n", aveRadiationPower);
		fclose(fp_POYNTING);
	}
#endif
   /* 
    if (myrank == 0)
    {
        int j, k;
        char path_EMax[MAX_SIZE_OF_PATH];
        sprintf(path_EMax, "%s\\EMax.txt", path_save);
        FILE* fp_EMax = fopen(path_EMax, "w+");
        if (fp_EMax == NULL)
            cout<<"Open EMax.txt fail."<<endl;
        else
        {
	        for(k = nzPML_1; k < Kmax-nzPML_2; ++k)
	        {
	        	for(j = nyPML_1; j < Jmax-nyPML_2-1; ++j)
	        	{
                    fprintf(fp_EMax, "%lf\t", maxEz(0, j, k));
	        	}
                fprintf(fp_EMax, "\n");
	        }
            fclose(fp_EMax);
        }
        char path_EMin[MAX_SIZE_OF_PATH];
        sprintf(path_EMin, "%s\\EMin.txt", path_save);
        FILE* fp_EMin = fopen(path_EMin, "w+");

        if (fp_EMin == NULL)
            cout<<"Open EMin.txt fail."<<endl;
        else
        {
	        for(k = nzPML_1; k < Kmax-nzPML_2; ++k)
	        {
	        	for(j = nyPML_1; j < Jmax-nyPML_2-1; ++j)
	        	{
                    fprintf(fp_EMin, "%lf\t", minEz(0, j, k));
	        	}
                fprintf(fp_EMin, "\n");
	        }
            fclose(fp_EMin);
        }
    }*/
    if (sourceType == 1 && myrank == 0)
    {
        plane_wave_power = amp*amp/2/sqrt(muO/epsO);
        printf("Power of plane wave : %eW/m2\n", plane_wave_power);
        fprintf(fp_log, "Power of plane wave : %eW/m2\n", plane_wave_power);
        fflush(stdout);
    }
	if (myrank == 0)
	{
		tCurrent = MPI_Wtime();
		printf("All computational tasks are completed.\n");
		printf("Time consumption : %d second, %d time step.\n", (int)(tCurrent - tStart), nMax);
		fprintf(fp_log, "Time consumption : %d second, %d time step.\n", (int)(tCurrent - tStart), nMax);
		fflush(stdout);
	}
}

/*****************************************************************************************/
/*
 * Function name：computeFieldX_Y
 * Description: compute, X = E or H
 *					     Y = 0 or NULL or nprocsSub1 --> rank = 0，middle，nprocs-1
 * Parameters: 
 * Return
 */
void computeFieldH_0_Dispersion()
{
	int i, j, k;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Hx
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i <= ie+1; ++i)
	{
		for (j = 0; j < Jmax-1; ++j)
		{
			for (k = 1; k < Kmax-1; ++k)
			{
				Hx(i,j,k) = DA * Hx(i,j,k) + DB *
					((Ez(i,j,k) - Ez(i,j+1,k)) * den_hy[j]  +
					(Ey(i,j,k) - Ey(i,j,k-1)) * den_hz[k] );
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Hy
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i <= ie; ++i)
	{
		for(j = 0; j < Jmax-1; ++j)
		{
			for(k = 1; k < Kmax-1; ++k)
			{
				Hy(i,j,k) = DA * Hy(i,j,k) + DB *
					((Ez(i+1,j,k) - Ez(i,j,k)) * den_hx[i] +
					(Ex(i,j,k-1) - Ex(i,j,k)) * den_hz[k] );
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Hz
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i <= ie; ++i)
	{
		for(j = 0; j< Jmax-1; ++j)
		{
			for(k = 0; k < Kmax-1; ++k)
			{
				Hz(i,j,k) = DA * Hz(i,j,k) + DB
							* ((Ey(i,j,k) - Ey(i+1,j,k)) * den_hx[i] +
							(Ex(i,j+1,k) - Ex(i,j,k)) * den_hy[j]);
			}
		}
	}
}

void computeFieldH_Dispersion()
{
	int i, j, k;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Hx
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i <= ie+1; ++i)
	{
		for (j = 0; j < Jmax-1; ++j)
		{
			for (k = 1; k < Kmax-1; ++k)
			{
				Hx(i,j,k) = DA * Hx(i,j,k) + DB *
					((Ez(i,j,k) - Ez(i,j+1,k)) * den_hy[j]  +
					(Ey(i,j,k) - Ey(i,j,k-1)) * den_hz[k] );
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Hy
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i <= ie; ++i)
	{
		for(j = 0; j < Jmax-1; ++j)
		{
			for(k = 1; k < Kmax-1; ++k)
			{
				Hy(i,j,k) = DA * Hy(i,j,k) + DB *
					((Ez(i+1,j,k) - Ez(i,j,k)) * den_hx[i] +
					(Ex(i,j,k-1) - Ex(i,j,k)) * den_hz[k] );
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Hz
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i <= ie; ++i)
	{
		for(j = 0; j< Jmax-1; ++j)
		{
			for(k = 0; k < Kmax-1; ++k)
			{
				Hz(i,j,k) = DA * Hz(i,j,k) + DB
							* ((Ey(i,j,k) - Ey(i+1,j,k)) * den_hx[i] +
							(Ex(i,j+1,k) - Ex(i,j,k)) * den_hy[j]);
			}
		}
	}
}

void computeFieldH_nprocsSub1_Dispersion()
{
	int i, j, k;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Hx
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i < ie; ++i)//ie(max)= Imax -1
	{
		for (j = 0; j < Jmax-1; ++j)
		{
			for (k = 1; k < Kmax-1; ++k)
			{
				Hx(i,j,k) = DA * Hx(i,j,k) + DB *
					((Ez(i,j,k) - Ez(i,j+1,k)) * den_hy[j]  +
					(Ey(i,j,k) - Ey(i,j,k-1)) * den_hz[k] );
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Hy
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i < ie; ++i)
	{
		for(j = 0; j < Jmax-1; ++j)
		{
			for(k = 1; k < Kmax-1; ++k)
			{
				Hy(i,j,k) = DA * Hy(i,j,k) + DB *
					((Ez(i+1,j,k) - Ez(i,j,k)) * den_hx[i] +
					(Ex(i,j,k-1) - Ex(i,j,k)) * den_hz[k] );
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Hz
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i < ie; ++i)
	{
		for(j = 0; j< Jmax-1; ++j)
		{
			for(k = 0; k < Kmax-1; ++k)
			{
				Hz(i,j,k) = DA * Hz(i,j,k) + DB
							* ((Ey(i,j,k) - Ey(i+1,j,k)) * den_hx[i] +
							(Ex(i,j+1,k) - Ex(i,j,k)) * den_hy[j]);
			}
		}
	}
}

void computePMLH_0_Dispersion()
{
	int i, j, k;
	int ii, i2, jj, j2, kk, k2;
	for(i = is; i <= ie+1; ++i)
	{
		jj = nyPML_2 - 2;
		j2 = Jmax - nyPML_2;
		for(j = 0; j < nyPML_1-1; ++j)
		{
			for (k = 1; k < Kmax-1; ++k)
			{
				psi_Hxy_1[i][j][k] = bh_y_1[j] * psi_Hxy_1[i][j][k]
								   + ch_y_1[j] * (Ez(i,j,k) - Ez(i,j+1,k)) / dy;
				Hx(i,j,k) = Hx(i,j,k) + DB * psi_Hxy_1[i][j][k];

				psi_Hxy_2[i][jj][k] = bh_y_2[jj] * psi_Hxy_2[i][jj][k]
									+ ch_y_2[jj] * (Ez(i,j2,k) - Ez(i,j2+1,k)) / dy;
				Hx(i,j2,k) = Hx(i,j2,k) + DB * psi_Hxy_2[i][jj][k];
			}
			--jj;
			++j2;
		}
	}

	for(i = is; i <= ie+1; ++i)
	{
		for(j = 0; j < Jmax-1; ++j)
		{
			kk = nzPML_2 - 2;
			k2 = Kmax - nzPML_2;
			for(k = 1; k < nzPML_1; ++k)
			{
				psi_Hxz_1[i][j][k-1] = bh_z_1[k-1] * psi_Hxz_1[i][j][k-1]
									 + ch_z_1[k-1] * (Ey(i,j,k) - Ey(i,j,k-1)) / dz;
				Hx(i,j,k) = Hx(i,j,k) + DB * psi_Hxz_1[i][j][k-1];

				psi_Hxz_2[i][j][kk] = bh_z_2[kk] * psi_Hxz_2[i][j][kk]
									+ ch_z_2[kk] * (Ey(i,j,k2) - Ey(i,j,k2-1)) / dz;
				Hx(i,j,k2) = Hx(i,j,k2) + DB * psi_Hxz_2[i][j][kk];
				--kk;
				++k2;
			}
		}
	}

	ii = nxPML_2 - 2;
	i2 = Imax - nxPML_2;
	for(i = 0; i < nxPML_1-1; ++i)
	{
		for(j = 0; j < Jmax-1; ++j)
		{
			for (k = 1; k < Kmax-1; ++k)
			{
				psi_Hyx_1[i][j][k] = bh_x_1[i] * psi_Hyx_1[i][j][k]
						+ ch_x_1[i] * (Ez(i+1,j,k) - Ez(i,j,k)) / dx;
				Hy(i,j,k) = Hy(i,j,k) + DB * psi_Hyx_1[i][j][k];
			}
		}
		--ii;
		++i2;
	}

	for(i = is; i <= ie; ++i)
	{
		for(j = 0; j < Jmax-1; ++j)
		{
			kk = nzPML_2 - 2;
			k2 = Kmax - nzPML_2;
			for(k = 1; k < nzPML_1; ++k)
			{
				psi_Hyz_1[i][j][k-1] = bh_z_1[k-1] * psi_Hyz_1[i][j][k-1]
									 + ch_z_1[k-1] * (Ex(i,j,k-1) - Ex(i,j,k)) / dz;
				Hy(i,j,k) = Hy(i,j,k) + DB * psi_Hyz_1[i][j][k-1];

				psi_Hyz_2[i][j][kk] = bh_z_2[kk] * psi_Hyz_2[i][j][kk]
									+ ch_z_2[kk] * (Ex(i,j,k2-1) - Ex(i,j,k2)) / dz;
				Hy(i,j,k2) = Hy(i,j,k2) + DB * psi_Hyz_2[i][j][kk];
				--kk;
				++k2;
			}
		}
	}

	ii = nxPML_2 - 2;
	i2 = Imax - nxPML_2;
	for(i = 0; i < nxPML_1-1; ++i)
	{
		for(j = 0; j < Jmax-1; ++j)
		{
			for(k = 0; k < Kmax-1; ++k)
			{
				psi_Hzx_1[i][j][k] = bh_x_1[i] * psi_Hzx_1[i][j][k]
								   + ch_x_1[i] * (Ey(i,j,k) - Ey(i+1,j,k)) / dx;
				Hz(i,j,k) = Hz(i,j,k) + DB * psi_Hzx_1[i][j][k];
			}
		}
		--ii;
		++i2;
	}

	for(i = is; i <= ie; ++i)
	{
		jj = nyPML_2 - 2;
		j2 = Jmax - nyPML_2;
		for(j = 0; j < nyPML_1-1; ++j)
		{
			for(k = 0; k < Kmax-1; ++k)
			{
				psi_Hzy_1[i][j][k] = bh_y_1[j] * psi_Hzy_1[i][j][k]
								   + ch_y_1[j] * (Ex(i,j+1,k) - Ex(i,j,k)) / dy;
				Hz(i,j,k) = Hz(i,j,k) + DB* psi_Hzy_1[i][j][k];
				psi_Hzy_2[i][jj][k] = bh_y_2[jj] * psi_Hzy_2[i][jj][k]
									+ ch_y_2[jj] * (Ex(i,j2+1,k) - Ex(i,j2,k)) / dy;
				Hz(i,j2,k) = Hz(i,j2,k) + DB * psi_Hzy_2[i][jj][k];
			}
			--jj;
			++j2;
		}
	}
}

void computePMLH_Dispersion()
{
	int i, j, k;
	int jj, j2, kk, k2;
	for(i = is; i <= ie+1; ++i)
		{
			jj = nyPML_2 - 2;
			j2 = Jmax - nyPML_2;
			for(j = 0; j < nyPML_1-1; ++j)
			{
				for (k = 1; k < Kmax-1; ++k)
				{
					psi_Hxy_1[i][j][k] = bh_y_1[j] * psi_Hxy_1[i][j][k]
									   + ch_y_1[j] * (Ez(i,j,k) - Ez(i,j+1,k)) / dy;
					Hx(i,j,k) = Hx(i,j,k) + DB * psi_Hxy_1[i][j][k];

					psi_Hxy_2[i][jj][k] = bh_y_2[jj] * psi_Hxy_2[i][jj][k]
										+ ch_y_2[jj] * (Ez(i,j2,k) - Ez(i,j2+1,k)) / dy;
					Hx(i,j2,k) = Hx(i,j2,k) + DB * psi_Hxy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}
		for(i = is; i <= ie+1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				kk = nzPML_2 - 2;
				k2 = Kmax - nzPML_2;
				for(k = 1; k < nzPML_1; ++k)
				{
					psi_Hxz_1[i][j][k-1] = bh_z_1[k-1] * psi_Hxz_1[i][j][k-1]
										 + ch_z_1[k-1] * (Ey(i,j,k) - Ey(i,j,k-1)) / dz;
					Hx(i,j,k) = Hx(i,j,k) + DB * psi_Hxz_1[i][j][k-1];

					psi_Hxz_2[i][j][kk] = bh_z_2[kk] * psi_Hxz_2[i][j][kk]
										+ ch_z_2[kk] * (Ey(i,j,k2) - Ey(i,j,k2-1)) / dz;
					Hx(i,j,k2) = Hx(i,j,k2) + DB * psi_Hxz_2[i][j][kk];
					--kk;
					++k2;
				}
			}
		}

		for(i = is; i <= ie; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				kk = nzPML_2 - 2;
				k2 = Kmax - nzPML_2;
				for(k = 1; k < nzPML_1; ++k)
				{
					psi_Hyz_1[i][j][k-1] = bh_z_1[k-1] * psi_Hyz_1[i][j][k-1]
										 + ch_z_1[k-1] * (Ex(i,j,k-1) - Ex(i,j,k)) / dz;
					Hy(i,j,k) = Hy(i,j,k) + DB * psi_Hyz_1[i][j][k-1];

					psi_Hyz_2[i][j][kk] = bh_z_2[kk] * psi_Hyz_2[i][j][kk]
										+ ch_z_2[kk] * (Ex(i,j,k2-1) - Ex(i,j,k2)) / dz;
					Hy(i,j,k2) = Hy(i,j,k2) + DB * psi_Hyz_2[i][j][kk];
					--kk;
					++k2;
				}
			}
		}

		for(i = is; i <= ie; ++i)
		{
			jj = nyPML_2 - 2;
			j2 = Jmax - nyPML_2;
			for(j = 0; j < nyPML_1-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					psi_Hzy_1[i][j][k] = bh_y_1[j] * psi_Hzy_1[i][j][k]
									   + ch_y_1[j] * (Ex(i,j+1,k) - Ex(i,j,k)) / dy;
					Hz(i,j,k) = Hz(i,j,k) + DB*  psi_Hzy_1[i][j][k];
					psi_Hzy_2[i][jj][k] = bh_y_2[jj] * psi_Hzy_2[i][jj][k]
										+ ch_y_2[jj] * (Ex(i,j2+1,k) - Ex(i,j2,k)) / dy;
					Hz(i,j2,k) = Hz(i,j2,k) + DB * psi_Hzy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}
}

void computePMLH_nprocsSub1_Dispersion()
{
	int i, j, k;
	int ii, i2, jj, j2, kk, k2;
	for(i = is; i < Imax-1; ++i)
		{
			jj = nyPML_2 - 2;
			j2 = Jmax - nyPML_2;
			for(j = 0; j < nyPML_1-1; ++j)
			{
				for (k = 1; k < Kmax-1; ++k)
				{
					psi_Hxy_1[i][j][k] = bh_y_1[j] * psi_Hxy_1[i][j][k]
									   + ch_y_1[j] * (Ez(i,j,k) - Ez(i,j+1,k)) / dy;
					Hx(i,j,k) = Hx(i,j,k) + DB * psi_Hxy_1[i][j][k];

					psi_Hxy_2[i][jj][k] = bh_y_2[jj] * psi_Hxy_2[i][jj][k]
										+ ch_y_2[jj] * (Ez(i,j2,k) - Ez(i,j2+1,k)) / dy;
					Hx(i,j2,k) = Hx(i,j2,k) + DB * psi_Hxy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}

		for(i = is; i < Imax-1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				kk = nzPML_2 - 2;
				k2 = Kmax - nzPML_2;
				for(k = 1; k < nzPML_1; ++k)
				{
					psi_Hxz_1[i][j][k-1] = bh_z_1[k-1] * psi_Hxz_1[i][j][k-1]
										 + ch_z_1[k-1] * (Ey(i,j,k) - Ey(i,j,k-1)) / dz;
					Hx(i,j,k) = Hx(i,j,k) + DB * psi_Hxz_1[i][j][k-1];

					psi_Hxz_2[i][j][kk] = bh_z_2[kk] * psi_Hxz_2[i][j][kk]
										+ ch_z_2[kk] * (Ey(i,j,k2) - Ey(i,j,k2-1)) / dz;
					Hx(i,j,k2) = Hx(i,j,k2) + DB * psi_Hxz_2[i][j][kk];
					--kk;
					++k2;
				}
			}
		}

		ii = nxPML_2 - 2;
		i2 = Imax - nxPML_2;
		for(i = 0; i < nxPML_1-1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				for (k = 1; k < Kmax-1; ++k)
				{
					psi_Hyx_2[ii][j][k] = bh_x_2[ii] * psi_Hyx_2[ii][j][k]
							+ ch_x_2[ii] * (Ez(i2+1,j,k) - Ez(i2,j,k)) / dx;
					Hy(i2,j,k) = Hy(i2,j,k) + DB * psi_Hyx_2[ii][j][k];
				}
			}
			--ii;
			++i2;
		}

		for(i = is; i < Imax-1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				kk = nzPML_2 - 2;
				k2 = Kmax - nzPML_2;
				for(k = 1; k < nzPML_1; ++k)
				{
					psi_Hyz_1[i][j][k-1] = bh_z_1[k-1] * psi_Hyz_1[i][j][k-1]
										 + ch_z_1[k-1] * (Ex(i,j,k-1) - Ex(i,j,k)) / dz;
					Hy(i,j,k) = Hy(i,j,k) + DB * psi_Hyz_1[i][j][k-1];

					psi_Hyz_2[i][j][kk] = bh_z_2[kk] * psi_Hyz_2[i][j][kk]
										+ ch_z_2[kk] * (Ex(i,j,k2-1) - Ex(i,j,k2)) / dz;
					Hy(i,j,k2) = Hy(i,j,k2) + DB * psi_Hyz_2[i][j][kk];
					--kk;
					++k2;
				}
			}
		}

		ii = nxPML_2 - 2;
		i2 = Imax - nxPML_2;
		for(i = 0; i < nxPML_1-1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					psi_Hzx_2[ii][j][k] = bh_x_2[ii] * psi_Hzx_2[ii][j][k]
										+ ch_x_2[ii] * (Ey(i2,j,k) - Ey(i2+1,j,k))/dx;
					Hz(i2,j,k) = Hz(i2,j,k) + DB * psi_Hzx_2[ii][j][k];
				}
			}
			--ii;
			++i2;
		}

		for(i = is; i < Imax-1; ++i)
		{
			//........................................................
			//  PML for bottom Hz, y-direction
			//.........................................................
			jj = nyPML_2 - 2;
			j2 = Jmax - nyPML_2;
			for(j = 0; j < nyPML_1-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					psi_Hzy_1[i][j][k] = bh_y_1[j] * psi_Hzy_1[i][j][k]
									   + ch_y_1[j] * (Ex(i,j+1,k) - Ex(i,j,k)) / dy;
					Hz(i,j,k) = Hz(i,j,k) + DB*  psi_Hzy_1[i][j][k];
			//.........................................................
			//  PML for top Hz, y-direction
			//..........................................................
					psi_Hzy_2[i][jj][k] = bh_y_2[jj] * psi_Hzy_2[i][jj][k]
										+ ch_y_2[jj] * (Ex(i,j2+1,k) - Ex(i,j2,k)) / dy;
					Hz(i,j2,k) = Hz(i,j2,k) + DB * psi_Hzy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}
}

void computeFieldE_0_Dispersion()
{
	int i, j, k;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Ex
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i <= ie; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
#ifdef _POYNTING
					preEx(i,j,k) = Ex(i,j,k);
#endif
                    formulaEx_Dispersion(i, j, k);
				}
			}
		}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Ey
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is+1; i <= ie+1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
#ifdef _POYNTING
					preEy(i,j,k) = Ey(i,j,k);
#endif
                    formulaEy_Dispersion(i, j, k);
				}
			}
		}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Ez
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is+1; i <= ie+1; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				for(k = 1; k < Kmax-1; ++k)
				{
#ifdef _POYNTING
					preEz(i,j,k) = Ez(i,j,k);
#endif
                   /* 
                    if ( i == thicknessOfPml+1+planeWaveIndex && Ez(i, j, k) > maxEz(0, j, k) && Ez(i, j, k) > 0)
                    {
                        maxEz(0, j, k) = Ez(i, j, k);
                    }
                    if ( i == thicknessOfPml+1+planeWaveIndex && Ez(i, j, k) < minEz(0, j, k) && Ez(i, j, k) < 0)
                    {
                        minEz(0, j, k) = Ez(i, j, k);
                    }*/
                    formulaEz_Dispersion(i, j, k);
				}
			}
		}
/*
    for(j = 1; j < Jmax-1; ++j)
	{
		for(k = 1; k < Kmax-1; ++k)
		{
            if ( abs(Ez(is+thicknessOfPml+1+5+1, j, k)) > maxEz(j, k))//检查波源发生平面+1处的电场值
                maxEz(j, k) = abs(Ez(is+thicknessOfPml+1+5+1, j, k));
		}
	}*/
}

void computeFieldE_Dispersion()
{
	int i, j, k;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Ex
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i <= ie; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
#ifdef _POYNTING
					preEx(i,j,k) = Ex(i,j,k);
#endif
					formulaEx_Dispersion(i, j, k);
				}
			}
		}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Ey
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i <= ie; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
#ifdef _POYNTING
					preEy(i,j,k) = Ey(i,j,k);
#endif
					formulaEy_Dispersion(i, j, k);
				}
			}
		}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Ez
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i <= ie; ++i)
	{
		for(j = 1; j < Jmax-1; ++j)
		{
			for(k = 1; k < Kmax-1; ++k)
			{
#ifdef _POYNTING
				preEz(i,j,k) = Ez(i,j,k);
#endif
				formulaEz_Dispersion(i, j, k);
			}
		}
	}
}

void computeFieldE_right_Dispersion()
{
	int i, j, k;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Ey
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	i = ie+1;
			for(j = 0; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
#ifdef _POYNTING
					preEy(i,j,k) = Ey(i,j,k);
#endif
					formulaEy_Dispersion(i, j, k);
				}
			}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Ez
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	i = ie+1;
		for(j = 1; j < Jmax-1; ++j)
		{
			for(k = 1; k < Kmax-1; ++k)
			{
#ifdef _POYNTING
				preEz(i,j,k) = Ez(i,j,k);
#endif
				formulaEz_Dispersion(i, j, k);
			}
		}
}

void computeFieldE_nprocsSub1_Dispersion()
{
	int i, j, k;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Ex
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i < ie; ++i)//ie(max)= Imax -1
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
#ifdef _POYNTING
					preEx(i,j,k) = Ex(i,j,k);
#endif
					formulaEx_Dispersion(i, j, k);
				}
			}
		}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Ey
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i < ie; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
#ifdef _POYNTING
					preEy(i,j,k) = Ey(i,j,k);
#endif
					formulaEy_Dispersion(i, j, k);
				}
			}
		}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  UPDATE Ez
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(i = is; i < ie; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				for(k = 1; k < Kmax-1; ++k)
				{
#ifdef _POYNTING
					preEz(i,j,k) = Ez(i,j,k);
#endif
					formulaEz_Dispersion(i, j, k);
				}
			}
		}
}

void computePMLE_0_Dispersion()
{
	int i, j, k;
	int ii, i2, jj, j2, kk, k2;
	for(i = is; i <= ie; ++i)
		{
			jj = nyPML_2 - 1;
			j2 = Jmax - nyPML_2;
			for(j = 1; j < nyPML_1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					psi_Exy_1[i][j][k] = be_y_1[j] * psi_Exy_1[i][j][k]
									   + ce_y_1[j] * (Hz(i,j,k) - Hz(i,j-1,k))/dy;
					Ex(i,j,k) = Ex(i,j,k) + CB_PML[modelDataX[i][j][k]] * psi_Exy_1[i][j][k];
					psi_Exy_2[i][jj][k] = be_y_2[jj] * psi_Exy_2[i][jj][k]
										+ ce_y_2[jj] * (Hz(i,j2,k) - Hz(i,j2-1,k)) / dy;
					Ex(i,j2,k) = Ex(i,j2,k) + CB_PML[modelDataX[i][j2][k]] * psi_Exy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}

		for(i = is; i <= ie; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				kk = nzPML_2 - 1;
				k2 = Kmax - nzPML_2 - 1;
				for(k = 0; k < nzPML_1; ++k)
				{
					psi_Exz_1[i][j][k] = be_z_1[k] * psi_Exz_1[i][j][k]
									   + ce_z_1[k] * (Hy(i,j,k) - Hy(i,j,k+1)) / dz;
					Ex(i,j,k) = Ex(i,j,k) + CB_PML[modelDataX[i][j][k]] * psi_Exz_1[i][j][k];
					psi_Exz_2[i][j][kk] = be_z_2[kk] * psi_Exz_2[i][j][kk]
										+ ce_z_2[kk] * (Hy(i,j,k2) - Hy(i,j,k2+1)) / dz;
					Ex(i,j,k2) = Ex(i,j,k2) + CB_PML[modelDataX[i][j][k2]] * psi_Exz_2[i][j][kk];
					--kk;
					++k2;
				}
			}
		}

		ii = nxPML_2 - 1;
		i2 = Imax - nxPML_2;
		for(i = 1; i < nxPML_1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					psi_Eyx_1[i][j][k] = be_x_1[i] * psi_Eyx_1[i][j][k]
									   + ce_x_1[i] * (Hz(i-1,j,k) - Hz(i,j,k)) / dx;
					Ey(i,j,k) = Ey(i,j,k) + CB_PML[modelDataY[i][j][k]] * psi_Eyx_1[i][j][k];
				}
			}
			--ii;
			++i2;
		}

		for(i = is+1; i <= ie+1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				kk = nzPML_2 - 1;
				k2 = Kmax - nzPML_2 - 1;
				for(k = 0; k < nzPML_1; ++k)
				{
					psi_Eyz_1[i][j][k] = be_z_1[k] * psi_Eyz_1[i][j][k]
									   + ce_z_1[k] * (Hx(i,j,k+1) - Hx(i,j,k)) / dz;
					Ey(i,j,k) = Ey(i,j,k) + CB_PML[modelDataY[i][j][k]] * psi_Eyz_1[i][j][k];
					psi_Eyz_2[i][j][kk] = be_z_2[kk] * psi_Eyz_2[i][j][kk]
										+ ce_z_2[kk] * (Hx(i,j,k2+1) - Hx(i,j,k2)) / dz;
					Ey(i,j,k2) = Ey(i,j,k2) + CB_PML[modelDataY[i][j][k2]] * psi_Eyz_2[i][j][kk];
					--kk;
					++k2;
				}
			}
		}

		ii = nxPML_2 - 1;
		i2 = Imax - nxPML_2;
		for(i = 1; i < nxPML_1; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				for(k = 1; k < Kmax-1; ++k)
				{
					psi_Ezx_1[i][j][k] = be_x_1[i] * psi_Ezx_1[i][j][k]
									   + ce_x_1[i] * (Hy(i,j,k) - Hy(i-1,j,k)) / dx;
					Ez(i,j,k) = Ez(i,j,k) + CB_PML[modelDataZ[i][j][k]] * psi_Ezx_1[i][j][k];
				}
			}
			--ii;
			++i2;
		}

		for(i = is+1; i <= ie+1; ++i)
		{
			jj = nyPML_2 - 1;
			j2 = Jmax - nyPML_2;
			for(j = 1; j < nyPML_1; ++j)
			{
				for(k = 1; k < Kmax-1; ++k)
				{
					psi_Ezy_1[i][j][k] = be_y_1[j] * psi_Ezy_1[i][j][k]
									   + ce_y_1[j] * (Hx(i,j-1,k) - Hx(i,j,k)) / dy;
					Ez(i,j,k) = Ez(i,j,k) + CB_PML[modelDataZ[i][j][k]] * psi_Ezy_1[i][j][k];
					psi_Ezy_2[i][jj][k] = be_y_2[jj] * psi_Ezy_2[i][jj][k]
										+ ce_y_2[jj] * (Hx(i,j2-1,k) - Hx(i,j2,k)) / dy;
					Ez(i,j2,k) = Ez(i,j2,k) + CB_PML[modelDataZ[i][j2][k]] * psi_Ezy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}
}

void computePMLE_Dispersion()
{
	int i, j, k;
	int jj, j2, kk, k2;
	for(i = is; i <= ie; ++i)
		{
			jj = nyPML_2 - 1;
			j2 = Jmax - nyPML_2;
			for(j = 1; j < nyPML_1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					psi_Exy_1[i][j][k] = be_y_1[j] * psi_Exy_1[i][j][k]
									   + ce_y_1[j] * (Hz(i,j,k) - Hz(i,j-1,k))/dy;
					Ex(i,j,k) = Ex(i,j,k) + CB_PML[modelDataX[i][j][k]] * psi_Exy_1[i][j][k];
					psi_Exy_2[i][jj][k] = be_y_2[jj] * psi_Exy_2[i][jj][k]
										+ ce_y_2[jj] * (Hz(i,j2,k) - Hz(i,j2-1,k)) / dy;
					Ex(i,j2,k) = Ex(i,j2,k) + CB_PML[modelDataX[i][j2][k]] * psi_Exy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}

		for(i = is; i <= ie; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				kk = nzPML_2 - 1;
				k2 = Kmax - nzPML_2 - 1;
				for(k = 0; k < nzPML_1; ++k)
				{
					psi_Exz_1[i][j][k] = be_z_1[k] * psi_Exz_1[i][j][k]
									   + ce_z_1[k] * (Hy(i,j,k) - Hy(i,j,k+1)) / dz;
					Ex(i,j,k) = Ex(i,j,k) + CB_PML[modelDataX[i][j][k]] * psi_Exz_1[i][j][k];
					psi_Exz_2[i][j][kk] = be_z_2[kk] * psi_Exz_2[i][j][kk]
										+ ce_z_2[kk] * (Hy(i,j,k2) - Hy(i,j,k2+1)) / dz;
					Ex(i,j,k2) = Ex(i,j,k2) + CB_PML[modelDataX[i][j][k2]] * psi_Exz_2[i][j][kk];
					--kk;
					++k2;
				}
			}
		}
		for(i = is; i <= ie+1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				kk = nzPML_2 - 1;
				k2 = Kmax - nzPML_2 - 1;
				for(k = 0; k < nzPML_1; ++k)
				{
					psi_Eyz_1[i][j][k] = be_z_1[k] * psi_Eyz_1[i][j][k]
									   + ce_z_1[k] * (Hx(i,j,k+1) - Hx(i,j,k)) / dz;
					Ey(i,j,k) = Ey(i,j,k) + CB_PML[modelDataY[i][j][k]] * psi_Eyz_1[i][j][k];
					psi_Eyz_2[i][j][kk] = be_z_2[kk] * psi_Eyz_2[i][j][kk]
										+ ce_z_2[kk] * (Hx(i,j,k2+1) - Hx(i,j,k2)) / dz;
					Ey(i,j,k2) = Ey(i,j,k2) + CB_PML[modelDataY[i][j][k2]] * psi_Eyz_2[i][j][kk];
					--kk;
					++k2;
				}
			}
		}
		for(i = is; i <= ie+1; ++i)
		{
			jj = nyPML_2 - 1;
			j2 = Jmax - nyPML_2;
			for(j = 1; j < nyPML_1; ++j)
			{
				for(k = 1; k < Kmax-1; ++k)
				{
					psi_Ezy_1[i][j][k] = be_y_1[j] * psi_Ezy_1[i][j][k]
									   + ce_y_1[j] * (Hx(i,j-1,k) - Hx(i,j,k)) / dy;
					Ez(i,j,k) = Ez(i,j,k) + CB_PML[modelDataZ[i][j][k]] * psi_Ezy_1[i][j][k];
					psi_Ezy_2[i][jj][k] = be_y_2[jj] * psi_Ezy_2[i][jj][k]
										+ ce_y_2[jj] * (Hx(i,j2-1,k) - Hx(i,j2,k)) / dy;
					Ez(i,j2,k) = Ez(i,j2,k) + CB_PML[modelDataZ[i][j2][k]] * psi_Ezy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}
}

void computePMLE_nprocsSub1_Dispersion()
{
	int i, j, k;
	int ii, i2, jj, j2, kk, k2;
	for(i = is; i < Imax-1; ++i)
		{
			jj = nyPML_2 - 1;
			j2 = Jmax - nyPML_2;
			for(j = 1; j < nyPML_1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					psi_Exy_1[i][j][k] = be_y_1[j] * psi_Exy_1[i][j][k]
									   + ce_y_1[j] * (Hz(i,j,k) - Hz(i,j-1,k))/dy;
					Ex(i,j,k) = Ex(i,j,k) + CB_PML[modelDataX[i][j][k]] * psi_Exy_1[i][j][k];
					psi_Exy_2[i][jj][k] = be_y_2[jj] * psi_Exy_2[i][jj][k]
										+ ce_y_2[jj] * (Hz(i,j2,k) - Hz(i,j2-1,k)) / dy;
					Ex(i,j2,k) = Ex(i,j2,k) + CB_PML[modelDataX[i][j2][k]] * psi_Exy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}

		for(i = is; i < Imax-1; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				kk = nzPML_2 - 1;
				k2 = Kmax - nzPML_2 - 1;
				for(k = 0; k < nzPML_1; ++k)
				{
					psi_Exz_1[i][j][k] = be_z_1[k] * psi_Exz_1[i][j][k]
									   + ce_z_1[k] * (Hy(i,j,k) - Hy(i,j,k+1)) / dz;
					Ex(i,j,k) = Ex(i,j,k) + CB_PML[modelDataX[i][j][k]] * psi_Exz_1[i][j][k];
					psi_Exz_2[i][j][kk] = be_z_2[kk] * psi_Exz_2[i][j][kk]
										+ ce_z_2[kk] * (Hy(i,j,k2) - Hy(i,j,k2+1)) / dz;
					Ex(i,j,k2) = Ex(i,j,k2) + CB_PML[modelDataX[i][j][k2]] * psi_Exz_2[i][j][kk];
					--kk;
					++k2;
				}
			}
		}

		ii = nxPML_2 - 1;
		i2 = Imax - nxPML_2;
		for(i = 1; i < nxPML_1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					psi_Eyx_2[ii][j][k] = be_x_2[ii] * psi_Eyx_2[ii][j][k]
										+ ce_x_2[ii] * (Hz(i2-1,j,k) - Hz(i2,j,k)) / dx;
					Ey(i2,j,k) = Ey(i2,j,k) + CB_PML[modelDataY[i2][j][k]] * psi_Eyx_2[ii][j][k];
				}
			}
			--ii;
			++i2;
		}

		for(i = is; i < Imax-1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				kk = nzPML_2 - 1;
				k2 = Kmax - nzPML_2 - 1;
				for(k = 0; k < nzPML_1; ++k)
				{
					psi_Eyz_1[i][j][k] = be_z_1[k] * psi_Eyz_1[i][j][k]
									   + ce_z_1[k] * (Hx(i,j,k+1) - Hx(i,j,k)) / dz;
					Ey(i,j,k) = Ey(i,j,k) + CB_PML[modelDataY[i][j][k]] * psi_Eyz_1[i][j][k];
					psi_Eyz_2[i][j][kk] = be_z_2[kk] * psi_Eyz_2[i][j][kk]
										+ ce_z_2[kk] * (Hx(i,j,k2+1) - Hx(i,j,k2)) / dz;
					Ey(i,j,k2) = Ey(i,j,k2) + CB_PML[modelDataY[i][j][k2]] * psi_Eyz_2[i][j][kk];
					--kk;
					++k2;
				}
			}
		}

		ii = nxPML_2 - 1;
		i2 = Imax - nxPML_2;
		for(i = 1; i < nxPML_1; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				for(k = 1; k < Kmax-1; ++k)
				{
					psi_Ezx_2[ii][j][k] = be_x_2[ii] * psi_Ezx_2[ii][j][k]
										+ ce_x_2[ii] * (Hy(i2,j,k) - Hy(i2-1,j,k)) / dx;
					Ez(i2,j,k) = Ez(i2,j,k) + CB_PML[modelDataZ[i2][j][k]] * psi_Ezx_2[ii][j][k];
				}
			}
			--ii;
			++i2;
		}

		for(i = is; i < Imax-1; ++i)
		{
			jj = nyPML_2 - 1;
			j2 = Jmax - nyPML_2;
			for(j = 1; j < nyPML_1; ++j)
			{
				for(k = 1; k < Kmax-1; ++k)
				{
					psi_Ezy_1[i][j][k] = be_y_1[j] * psi_Ezy_1[i][j][k]
									   + ce_y_1[j] * (Hx(i,j-1,k) - Hx(i,j,k)) / dy;
					Ez(i,j,k) = Ez(i,j,k) + CB_PML[modelDataZ[i][j][k]] * psi_Ezy_1[i][j][k];
					psi_Ezy_2[i][jj][k] = be_y_2[jj] * psi_Ezy_2[i][jj][k]
										+ ce_y_2[jj] * (Hx(i,j2-1,k) - Hx(i,j2,k)) / dy;
					Ez(i,j2,k) = Ez(i,j2,k) + CB_PML[modelDataZ[i][j2][k]] * psi_Ezy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}
}

void initializeDispersion()
{
    if (nprocs == 1) 
    {
        fprintf(fp_log, "Dispersion calculation only run in parallel!\n");
        exit(SETTING_ERROR);
    }
    if (myrank == 0)
	{
		Dz = initArrayFloat((Imax+1) * Jmax * Kmax);
		Dy = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));
		Dx = initArrayFloat(Imax * Jmax * (Kmax-1));
		preDz = initArrayFloat((Imax+1) * Jmax * Kmax);
		preDy = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));
		preDx = initArrayFloat(Imax * Jmax * (Kmax-1));
		prepreDz = initArrayFloat((Imax+1) * Jmax * Kmax);
		prepreDy = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));
		prepreDx = initArrayFloat(Imax * Jmax * (Kmax-1));
		preEz = initArrayFloat((Imax+1) * Jmax * Kmax);
		preEy = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));
		preEx = initArrayFloat(Imax * Jmax * (Kmax-1));
		prepreEz = initArrayFloat((Imax+1) * Jmax * Kmax);
		prepreEy = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));
		prepreEx = initArrayFloat(Imax * Jmax * (Kmax-1));
	}
	else if (myrank == nprocs - 1)
	{
		Dz = initArrayFloat(Imax * Jmax * Kmax);
		Dy = initArrayFloat(Imax * (Jmax-1) * (Kmax-1));
		Dx = initArrayFloat(Imax * Jmax * (Kmax-1));
		preDz = initArrayFloat(Imax * Jmax * Kmax);
		preDy = initArrayFloat(Imax * (Jmax-1) * (Kmax-1));
		preDx = initArrayFloat(Imax * Jmax * (Kmax-1));
		prepreDz = initArrayFloat(Imax * Jmax * Kmax);
		prepreDy = initArrayFloat(Imax * (Jmax-1) * (Kmax-1));
		prepreDx = initArrayFloat(Imax * Jmax * (Kmax-1));
		preEz = initArrayFloat(Imax * Jmax * Kmax);
		preEy = initArrayFloat(Imax * (Jmax-1) * (Kmax-1));
		preEx = initArrayFloat(Imax * Jmax * (Kmax-1));
		prepreEz = initArrayFloat(Imax * Jmax * Kmax);
		prepreEy = initArrayFloat(Imax * (Jmax-1) * (Kmax-1));
		prepreEx = initArrayFloat(Imax * Jmax * (Kmax-1));
	}
	else
	{
		Dz = initArrayFloat((Imax+1) * Jmax * Kmax);
		Dy = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));
		Dx = initArrayFloat(Imax * Jmax * (Kmax-1));
		preDz = initArrayFloat((Imax+1) * Jmax * Kmax);
		preDy = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));
		preDx = initArrayFloat(Imax * Jmax * (Kmax-1));
		prepreDz = initArrayFloat((Imax+1) * Jmax * Kmax);
		prepreDy = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));
		prepreDx = initArrayFloat(Imax * Jmax * (Kmax-1));
		preEz = initArrayFloat((Imax+1) * Jmax * Kmax);
		preEy = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));
		preEx = initArrayFloat(Imax * Jmax * (Kmax-1));
		prepreEz = initArrayFloat((Imax+1) * Jmax * Kmax);
		prepreEy = initArrayFloat((Imax+1) * (Jmax-1) * (Kmax-1));
		prepreEx = initArrayFloat(Imax * Jmax * (Kmax-1));
	}
}

int loadMediaData_Dispersion(char* path, int size)
{
	int i;
	int id;
	double eps, sig, rho, spec_heat, K, B;
    double eps_inf, delay;
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
		fscanf(fp, "%d %lf %lf %lf %lf %lf", &id, &sig, &eps, &eps_inf, &delay, &rho);
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
		media[id].sigma= sig;
		media[id].epsilon = eps;
        media[id].epsilon_inf = eps_inf;
        media[id].delay = delay;
		media[id].rho= rho;
        if (need_compute_temperature_rise)
        {
			media[id].spec_heat = spec_heat;
            media[id].K = K;
            media[id].B = B;
        }
        double DisTemp1, DisTemp2;
        if (media[id].delay == 0)
        {
            CA[id] = ((1.0 - sig * dt / (2.0 * eps * epsO))
				    / (1.0 + sig * dt / (2.0 * eps * epsO)));
		    CB[id] = (dt / (eps * epsO))
		    		/ (1.0 + sig * dt / (2.0 * eps * epsO));
		    CB_PML[id] = (dt / (eps * epsO))
		    		/ (1.0 + sig * dt / (2.0 * eps * epsO));
        }
        else
        {
		    CA[id] = (1.0 - sig * dt / 2.0) / (1.0 + sig * dt / 2.0);
		    CB[id] = dt / (1.0 + sig * dt / 2.0);
		    CB_PML[id] = (dt / (eps * epsO))
		    		/ (1.0 + sig * dt / (2.0 * eps * epsO));
            DisDa[id] = 1+delay/dt;
            DisDb[id] = 1+2*delay/dt;
            DisDc[id] = delay/dt;
            DisTemp1 = (sig*delay + epsO*eps);
            DisTemp2 = delay*epsO*eps_inf/dt;
            DisEa[id] = sig*dt/2+DisTemp1+DisTemp2;
            DisEb[id] = sig*dt/2-DisTemp1-2*DisTemp2;
            DisEc[id] = DisTemp2;
            DisDa[id] /= DisEa[id];
            DisDb[id] /= DisEa[id];
            DisDc[id] /= DisEa[id];
            DisEb[id] /= DisEa[id];
            DisEc[id] /= DisEa[id];
        }
		if (id>=maxMedia) maxMedia = id+1;
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

int formulaEx_Dispersion(int i, int j, int k)
{
    if (media[modelDataX[i][j][k]].delay == 0)
    {
        Ex(i,j,k) = CA[modelDataX[i][j][k]] * Ex(i,j,k) + CB[modelDataX[i][j][k]] *
	    			  ((Hz(i,j,k) - Hz(i,j-1,k)) * den_ey[j] +
	    			  (Hy(i,j,k) - Hy(i,j,k+1)) * den_ez[k] );
    }
    else
    {
        prepreDx(i, j, k) = preDx(i, j, k);
        preDx(i, j, k) = Dx(i, j, k);
        prepreEx(i, j, k) = preEx(i, j, k);
        preEx(i, j, k) = Ex(i, j, k);
	    Dx(i,j,k) = Dx(i,j,k) + dt*((Hz(i,j,k) - Hz(i,j-1,k)) * den_ey[j] +
	    			  (Hy(i,j,k) - Hy(i,j,k+1)) * den_ez[k] );
        Ex(i, j, k) = (DisDa[modelDataX[i][j][k]]*Dx(i, j, k) - 
                       DisDb[modelDataX[i][j][k]]*preDx(i, j, k) + 
                       DisDc[modelDataX[i][j][k]]*prepreDx(i, j, k) -
                       DisEb[modelDataX[i][j][k]]*preEx(i, j, k) - 
                       DisEc[modelDataX[i][j][k]]*prepreEx(i, j, k));
    }

    return SUCCESS;
}

int formulaEy_Dispersion(int i, int j, int k)
{
    if (media[modelDataY[i][j][k]].delay == 0)
    {
        Ey(i,j,k) = CA[modelDataY[i][j][k]] * Ey(i,j,k) + CB[modelDataY[i][j][k]] *
								   ((Hz(i-1,j,k) - Hz(i,j,k)) * den_ex[i] +
								   (Hx(i,j,k+1) - Hx(i,j,k)) * den_ez[k] );
    }
    else
    {
        prepreDy(i, j, k) = preDy(i, j, k);
        preDy(i, j, k) = Dy(i, j, k);
        prepreEy(i, j, k) = preEy(i, j, k);
        preEy(i, j, k) = Ey(i, j, k);
	    Dy(i,j,k) = Dy(i,j,k) + dt*((Hz(i-1,j,k) - Hz(i,j,k)) * den_ex[i] +
				      (Hx(i,j,k+1) - Hx(i,j,k)) * den_ez[k] );
        Ey(i, j, k) = (DisDa[modelDataY[i][j][k]]*Dy(i, j, k) - 
                       DisDb[modelDataY[i][j][k]]*preDy(i, j, k) + 
                       DisDc[modelDataY[i][j][k]]*prepreDy(i, j, k) -
                       DisEb[modelDataY[i][j][k]]*preEy(i, j, k) - 
                       DisEc[modelDataY[i][j][k]]*prepreEy(i, j, k));
    }

    return SUCCESS;
}

int formulaEz_Dispersion(int i, int j, int k)
{
    if (media[modelDataZ[i][j][k]].delay == 0)
    {
        Ez(i,j,k) = CA[modelDataZ[i][j][k]] * Ez(i,j,k) + CB[modelDataZ[i][j][k]] * 
								((Hy(i,j,k) - Hy(i-1,j,k)) * den_ex[i] +
								(Hx(i,j-1,k) - Hx(i,j,k)) * den_ey[j]);
    }
    else
    {
        prepreDz(i, j, k) = preDz(i, j, k);
        preDz(i, j, k) = Dz(i, j, k);
        prepreEz(i, j, k) = preEz(i, j, k);
        preEz(i, j, k) = Ez(i, j, k);
	    Dz(i,j,k) = Dz(i,j,k) + dt*((Hy(i,j,k) - Hy(i-1,j,k)) * den_ex[i] +
					  (Hx(i,j-1,k) - Hx(i,j,k)) * den_ey[j]);
        Ez(i, j, k) = (DisDa[modelDataZ[i][j][k]]*Dz(i, j, k) - 
                       DisDb[modelDataZ[i][j][k]]*preDz(i, j, k) + 
                       DisDc[modelDataZ[i][j][k]]*prepreDz(i, j, k) -
                       DisEb[modelDataZ[i][j][k]]*preEz(i, j, k) - 
                       DisEc[modelDataZ[i][j][k]]*prepreEz(i, j, k));
    }

    return SUCCESS;
}
