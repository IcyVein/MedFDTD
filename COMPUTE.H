/*****************************************************************************************/
/*
 * Function name: compute
 * Description: The main compute function
 * Parameters: 
 * Return
 */
void compute()
{
	int i, n;
	double tStart, tCurrent;
	double tStepStart, tStepEnd;
	double tStep=0.0;
	double tNear20[20] = {0.0};
	int dn = 10;
    char path_check[MAX_SIZE_OF_PATH];
    sprintf(path_check, "%scheckRMS_%d.txt", path_save, myrank);
    checkData = 0, checkDataPast = 0;

#ifdef _DEBUG_L3
        cout<<"In computer @Cpu : "<<myrank<<endl;
#endif

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
#ifdef _DEBUG_L3
        cout<<"Cpu : "<<myrank<<endl;
#endif
#ifdef _DEBUG
        if (myrank == 0)
            cout<<"Current step: "<<n<<endl;
#endif
		tStepStart = MPI_Wtime();
#ifdef _DEBUG_L3
        cout<<"--Cpu : "<<myrank<<endl;
#endif
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

#ifdef _DEBUG
        if (myrank == 0)
            cout<<"Enter H"<<endl;
#endif
		/* H */
		t0 = time(NULL);
		if (0<myrank && myrank<nprocs-1)
			computeFieldH();
		else if (myrank == 0)
			computeFieldH_0();
		else /* myrank == nprocs - 1 */
			computeFieldH_nprocsSub1();
		t[0] += time(NULL)-t0;

#ifdef _DEBUG
        if (myrank == 0)
            cout<<"Enter H - PML"<<endl;
#endif
		/* H-PML */
		t0 = time(NULL);
		if (0<myrank && myrank<nprocs-1)
			computePMLH();
		else if (myrank == 0)
			computePMLH_0();
		else /* myrank == nprocs - 1 */
			computePMLH_nprocsSub1();
		t[1] += time(NULL)-t0;

        if( sourceType == 1/*If Plane Wave*/)
		{
			powerSourcePlaneWaveH(n);/* -1/2 - dx/2/C/dt */
		}

		//MPI_Barrier(MPI_COMM_WORLD);

//		/* Exchage Hy、Hz */
//		/* Send right, Recv left */
//		t0 = time(NULL);
//		if (0<myrank && myrank<nprocs-1)
//		{
//			MPI_Isend(&Hy(ie,0,0), Jmax * Kmax, MPI_FLOAT, myrank+1, 100+myrank+1, MPI_COMM_WORLD, &request_sendy_r);
//			MPI_Irecv(&Hy(is-1,0,0), Jmax * Kmax, MPI_FLOAT, myrank-1, 100+myrank, MPI_COMM_WORLD, &request_recvy_l);
//			MPI_Isend(&Hz(ie,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank+1, 200+myrank+1, MPI_COMM_WORLD, &request_sendz_r);
//			MPI_Irecv(&Hz(is-1,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank-1, 200+myrank, MPI_COMM_WORLD, &request_recvz_l);	
//		}
//		else if (myrank == 0)
//		{
//			MPI_Isend(&Hy(ie,0,0), Jmax * Kmax, MPI_FLOAT, myrank+1, 100+myrank+1, MPI_COMM_WORLD, &request_sendy_r);
//			MPI_Isend(&Hz(ie,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank+1, 200+myrank+1, MPI_COMM_WORLD, &request_sendz_r);
//		}
//		else /* myrank == nprocs - 1 */
//		{
//			MPI_Irecv(&Hy(is-1,0,0), Jmax * Kmax, MPI_FLOAT, myrank-1, 100+myrank, MPI_COMM_WORLD, &request_recvy_l);
//			MPI_Irecv(&Hz(is-1,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank-1, 200+myrank, MPI_COMM_WORLD, &request_recvz_l);
//		}
//		t[2] += time(NULL)-t0;
//
//		/* Send left, Recv right */
//		t0 = time(NULL);
//		if (0<myrank && myrank<nprocs-1)
//		{
//			MPI_Isend(&Hy(is,0,0), Jmax * Kmax, MPI_FLOAT, myrank-1, 300+myrank-1, MPI_COMM_WORLD, &request_sendy_l);
//			MPI_Irecv(&Hy(ie+1,0,0), Jmax * Kmax, MPI_FLOAT, myrank+1, 300+myrank, MPI_COMM_WORLD, &request_recvy_r);
//			MPI_Isend(&Hz(is,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank-1, 400+myrank-1, MPI_COMM_WORLD, &request_sendz_l);
//			MPI_Irecv(&Hz(ie+1,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank+1, 400+myrank, MPI_COMM_WORLD, &request_recvz_r);
//		}
//		else if (myrank == 0)
//		{
//			MPI_Irecv(&Hy(ie+1,0,0), Jmax * Kmax, MPI_FLOAT, myrank+1, 300+myrank, MPI_COMM_WORLD, &request_recvy_r);
//			MPI_Irecv(&Hz(ie+1,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank+1, 400+myrank, MPI_COMM_WORLD, &request_recvz_r);
//		}
//		else /* myrank == nprocs - 1 */
//		{
//			MPI_Isend(&Hy(is,0,0), Jmax * Kmax, MPI_FLOAT, myrank-1, 300+myrank-1, MPI_COMM_WORLD, &request_sendy_l);
//			MPI_Isend(&Hz(is,0,0), (Jmax-1) * (Kmax-1), MPI_FLOAT, myrank-1, 400+myrank-1, MPI_COMM_WORLD, &request_sendz_l);
//		}
//		t[3] += time(NULL)-t0;

#ifdef _DEBUG
        if (myrank == 0)
            cout<<"Trans data"<<endl;
#endif
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

#ifdef _DEBUG
        if (myrank == 0)
            cout<<"Enter E"<<endl;
#endif
		/* E */
		t0 = time(NULL);
		if (0<myrank && myrank<nprocs-1)
		{
			MPI_Wait(&request_recvy_l, &status_recvy_l);
			MPI_Wait(&request_sendy_l, &status_sendy_l);

			MPI_Wait(&request_recvz_l, &status_recvz_l);
			MPI_Wait(&request_sendz_l, &status_sendz_l);

			computeFieldE();

			MPI_Wait(&request_recvy_r, &status_recvy_r);
			MPI_Wait(&request_sendy_r, &status_sendy_r);
			MPI_Wait(&request_recvz_r, &status_recvz_r);
			MPI_Wait(&request_sendz_r, &status_sendz_r);
			computeFieldE_right();
		}
		else if (myrank == 0)
		{
			MPI_Wait(&request_recvy_r, &status_recvy_r);
			MPI_Wait(&request_sendy_r, &status_sendy_r);
			MPI_Wait(&request_recvz_r, &status_recvz_r);
			MPI_Wait(&request_sendz_r, &status_sendz_r);
			computeFieldE_0();
		}
		else /* myrank == nprocs - 1 */
		{
			MPI_Wait(&request_recvy_l, &status_recvy_l);
			MPI_Wait(&request_sendy_l, &status_sendy_l);
			MPI_Wait(&request_recvz_l, &status_recvz_l);
			MPI_Wait(&request_sendz_l, &status_sendz_l);

			computeFieldE_nprocsSub1();
		}
		t[4] += time(NULL)-t0;

#ifdef _DEBUG
        if (myrank == 0)
            cout<<"Enter E - PML"<<endl;
#endif
		/* E-PML */
		t0 = time(NULL);
		if (0<myrank && myrank<nprocs-1)
		{
			computePMLE();
		}
		else if (myrank == 0)
		{
			computePMLE_0();
		}
		else /* myrank == nprocs - 1 */
			computePMLE_nprocsSub1();
		t[5] += time(NULL)-t0;

		t0 = time(NULL);
		//MPI_Barrier(MPI_COMM_WORLD);
		t[6] += time(NULL)-t0;

	//-----------------------------------------------------------
	//   Apply source
	//-----------------------------------------------------------

#ifdef _DEBUG
        if (myrank == 0)
            cout<<"Source"<<endl;
#endif
		t0 = time(NULL);

		powerSource(n);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  WRITE TO OUTPUT FILES
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		saveData(n);

        saveMaxE();

		MPI_Barrier(MPI_COMM_WORLD);
#ifdef _DEBUG
        if (myrank == 0)
            cout<<"Sar"<<endl;
#endif
#ifdef _SAR

#ifdef _DEBUG_L2
        if (myrank == 0)
            cout<<"computeRMS"<<endl;
#endif

		for (i = 0; i<_spaceZ+1; ++i)
		{
			computeRMS(&pSAR[i], n);
		}
#ifdef _DEBUG_L2
        if (myrank == 0)
            cout<<"checkData"<<endl;
#endif
        if (n%(sourceT) == 0 && n+sourceT<=nMax)
        {
            checkDataPast = checkData;
#ifdef _DEBUG_L3
        if (myrank == 0)
            cout<<"checkRMS --> "<<myrank<<endl;
#endif
            checkData = checkRMS();
#ifdef _DEBUG_L3
        if (myrank == 0)
            cout<<"convergence --> "<<myrank<<endl;
#endif
            if (checkData != 0 && checkData != checkDataPast)
                convergence = 10.0*log10(fabs(checkData-checkDataPast)/checkData);
            else if (checkData == checkDataPast)
                convergence = convergenceTarget;
            else
                convergence = 0.0;
#ifdef _DEBUG_L3
        if (myrank == 0)
            cout<<"reset--> "<<myrank<<endl;
#endif
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
#ifdef _DEBUG_L2
        if (myrank == 0)
            cout<<"resetRMS"<<endl;
#endif
        if (n == (nMax-sourceT))
            resetRMS();
#endif

#ifdef _DEBUG
        if (myrank == 0)
            cout<<"End sar"<<endl;
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
    fprintf(fp_mem, "inC1:\t%d\n", mem_count);

    char fn_maxE[MAX_SIZE_OF_PATH];
    sprintf(fn_maxE, "%smaxE_%d.txt", path_save, myrank);
    FILE* fp_maxE = fopen(fn_maxE, "w+");
    writeField(fp_maxE, *maxE, outputRangeE);
    fclose(fp_maxE);

	freeFDTDData();
    fprintf(fp_mem, "inC2:\t%d\n", mem_count);

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
    fprintf(fp_mem, "inC3:\t%d\n", mem_count);

	freeArray3Char(modelDataX, Imax+1, Jmax, Kmax);
	freeArray3Char(modelDataY, Imax+1, Jmax, Kmax);
	freeArray3Char(modelDataZ, Imax+1, Jmax, Kmax);
    fprintf(fp_mem, "inC4:\t%d\n", mem_count);

#ifdef _SAR
	if (nXgSAR || whole_body_sar)
	{
		computeXgSAR(nXgSAR);
	}
#endif
    fprintf(fp_mem, "inC5:\t%d\n", mem_count);

#ifdef _POYNTING
	if (myrank == 0)
	{
		aveRadiationPower /= (int)(1/(freq*dt));
		printf("Averaged radiation power = %eW\n", aveRadiationPower);
		fclose(fp_POYNTING);
	}
#endif
        
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
void computeFieldH_0()
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

void computeFieldH()
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

void computeFieldH_nprocsSub1()
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

void computePMLH_0()
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

void computePMLH()
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

void computePMLH_nprocsSub1()
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

void computeFieldE_0()
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
					Ex(i,j,k) = CA[modelDataX[i][j][k]] * Ex(i,j,k) + CB[modelDataX[i][j][k]] *
								  ((Hz(i,j,k) - Hz(i,j-1,k)) * den_ey[j] +
								  (Hy(i,j,k) - Hy(i,j,k+1)) * den_ez[k] );
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
					Ey(i,j,k) = CA[modelDataY[i][j][k]] * Ey(i,j,k) + CB[modelDataY[i][j][k]] *
								   ((Hz(i-1,j,k) - Hz(i,j,k)) * den_ex[i] +
								   (Hx(i,j,k+1) - Hx(i,j,k)) * den_ez[k] );
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
					Ez(i,j,k) = CA[modelDataZ[i][j][k]] * Ez(i,j,k) + CB[modelDataZ[i][j][k]]
								* ((Hy(i,j,k) - Hy(i-1,j,k)) * den_ex[i] +
								(Hx(i,j-1,k) - Hx(i,j,k)) * den_ey[j]);
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

void computeFieldE()
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
					Ex(i,j,k) = CA[modelDataX[i][j][k]] * Ex(i,j,k) + CB[modelDataX[i][j][k]] *
								  ((Hz(i,j,k) - Hz(i,j-1,k)) * den_ey[j] +
								  (Hy(i,j,k) - Hy(i,j,k+1)) * den_ez[k] );
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
					Ey(i,j,k) = CA[modelDataY[i][j][k]] * Ey(i,j,k) + CB[modelDataY[i][j][k]] *
								   ((Hz(i-1,j,k) - Hz(i,j,k)) * den_ex[i] +
								   (Hx(i,j,k+1) - Hx(i,j,k)) * den_ez[k] );
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
				Ez(i,j,k) = CA[modelDataZ[i][j][k]] * Ez(i,j,k) + CB[modelDataZ[i][j][k]]
							* ((Hy(i,j,k) - Hy(i-1,j,k)) * den_ex[i] +
							(Hx(i,j-1,k) - Hx(i,j,k)) * den_ey[j]);
			}
		}
	}
}

void computeFieldE_right()
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
					Ey(i,j,k) = CA[modelDataY[i][j][k]] * Ey(i,j,k) + CB[modelDataY[i][j][k]] *
								   ((Hz(i-1,j,k) - Hz(i,j,k)) * den_ex[i] +
								   (Hx(i,j,k+1) - Hx(i,j,k)) * den_ez[k] );
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
				Ez(i,j,k) = CA[modelDataZ[i][j][k]] * Ez(i,j,k) + CB[modelDataZ[i][j][k]]
							* ((Hy(i,j,k) - Hy(i-1,j,k)) * den_ex[i] +
							(Hx(i,j-1,k) - Hx(i,j,k)) * den_ey[j]);
			}
		}
}

void computeFieldE_nprocsSub1()
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
					Ex(i,j,k) = CA[modelDataX[i][j][k]] * Ex(i,j,k) + CB[modelDataX[i][j][k]] *
								  ((Hz(i,j,k) - Hz(i,j-1,k)) * den_ey[j] +
								  (Hy(i,j,k) - Hy(i,j,k+1)) * den_ez[k] );
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
					Ey(i,j,k) = CA[modelDataY[i][j][k]] * Ey(i,j,k) + CB[modelDataY[i][j][k]] *
								   ((Hz(i-1,j,k) - Hz(i,j,k)) * den_ex[i] +
								   (Hx(i,j,k+1) - Hx(i,j,k)) * den_ez[k] );
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
					Ez(i,j,k) = CA[modelDataZ[i][j][k]] * Ez(i,j,k) + CB[modelDataZ[i][j][k]]
								* ((Hy(i,j,k) - Hy(i-1,j,k)) * den_ex[i] +
								(Hx(i,j-1,k) - Hx(i,j,k)) * den_ey[j]);
				}
			}
		}
}

void computePMLE_0()
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
					Ex(i,j,k) = Ex(i,j,k) + CB[modelDataX[i][j][k]] * psi_Exy_1[i][j][k];
					psi_Exy_2[i][jj][k] = be_y_2[jj] * psi_Exy_2[i][jj][k]
										+ ce_y_2[jj] * (Hz(i,j2,k) - Hz(i,j2-1,k)) / dy;
					Ex(i,j2,k) = Ex(i,j2,k) + CB[modelDataX[i][j2][k]] * psi_Exy_2[i][jj][k];
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
					Ex(i,j,k) = Ex(i,j,k) + CB[modelDataX[i][j][k]] * psi_Exz_1[i][j][k];
					psi_Exz_2[i][j][kk] = be_z_2[kk] * psi_Exz_2[i][j][kk]
										+ ce_z_2[kk] * (Hy(i,j,k2) - Hy(i,j,k2+1)) / dz;
					Ex(i,j,k2) = Ex(i,j,k2) + CB[modelDataX[i][j][k2]] * psi_Exz_2[i][j][kk];
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
					Ey(i,j,k) = Ey(i,j,k) + CB[modelDataY[i][j][k]] * psi_Eyx_1[i][j][k];
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
					Ey(i,j,k) = Ey(i,j,k) + CB[modelDataY[i][j][k]] * psi_Eyz_1[i][j][k];
					psi_Eyz_2[i][j][kk] = be_z_2[kk] * psi_Eyz_2[i][j][kk]
										+ ce_z_2[kk] * (Hx(i,j,k2+1) - Hx(i,j,k2)) / dz;
					Ey(i,j,k2) = Ey(i,j,k2) + CB[modelDataY[i][j][k2]] * psi_Eyz_2[i][j][kk];
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
					Ez(i,j,k) = Ez(i,j,k) + CB[modelDataZ[i][j][k]] * psi_Ezx_1[i][j][k];
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
					Ez(i,j,k) = Ez(i,j,k) + CB[modelDataZ[i][j][k]] * psi_Ezy_1[i][j][k];
					psi_Ezy_2[i][jj][k] = be_y_2[jj] * psi_Ezy_2[i][jj][k]
										+ ce_y_2[jj] * (Hx(i,j2-1,k) - Hx(i,j2,k)) / dy;
					Ez(i,j2,k) = Ez(i,j2,k) + CB[modelDataZ[i][j2][k]] * psi_Ezy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}
}

void computePMLE()
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
					Ex(i,j,k) = Ex(i,j,k) + CB[modelDataX[i][j][k]] * psi_Exy_1[i][j][k];
					psi_Exy_2[i][jj][k] = be_y_2[jj] * psi_Exy_2[i][jj][k]
										+ ce_y_2[jj] * (Hz(i,j2,k) - Hz(i,j2-1,k)) / dy;
					Ex(i,j2,k) = Ex(i,j2,k) + CB[modelDataX[i][j2][k]] * psi_Exy_2[i][jj][k];
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
					Ex(i,j,k) = Ex(i,j,k) + CB[modelDataX[i][j][k]] * psi_Exz_1[i][j][k];
					psi_Exz_2[i][j][kk] = be_z_2[kk] * psi_Exz_2[i][j][kk]
										+ ce_z_2[kk] * (Hy(i,j,k2) - Hy(i,j,k2+1)) / dz;
					Ex(i,j,k2) = Ex(i,j,k2) + CB[modelDataX[i][j][k2]] * psi_Exz_2[i][j][kk];
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
					Ey(i,j,k) = Ey(i,j,k) + CB[modelDataY[i][j][k]] * psi_Eyz_1[i][j][k];
					psi_Eyz_2[i][j][kk] = be_z_2[kk] * psi_Eyz_2[i][j][kk]
										+ ce_z_2[kk] * (Hx(i,j,k2+1) - Hx(i,j,k2)) / dz;
					Ey(i,j,k2) = Ey(i,j,k2) + CB[modelDataY[i][j][k2]] * psi_Eyz_2[i][j][kk];
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
					Ez(i,j,k) = Ez(i,j,k) + CB[modelDataZ[i][j][k]] * psi_Ezy_1[i][j][k];
					psi_Ezy_2[i][jj][k] = be_y_2[jj] * psi_Ezy_2[i][jj][k]
										+ ce_y_2[jj] * (Hx(i,j2-1,k) - Hx(i,j2,k)) / dy;
					Ez(i,j2,k) = Ez(i,j2,k) + CB[modelDataZ[i][j2][k]] * psi_Ezy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}
}

void computePMLE_nprocsSub1()
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
					Ex(i,j,k) = Ex(i,j,k) + CB[modelDataX[i][j][k]] * psi_Exy_1[i][j][k];
					psi_Exy_2[i][jj][k] = be_y_2[jj] * psi_Exy_2[i][jj][k]
										+ ce_y_2[jj] * (Hz(i,j2,k) - Hz(i,j2-1,k)) / dy;
					Ex(i,j2,k) = Ex(i,j2,k) + CB[modelDataX[i][j2][k]] * psi_Exy_2[i][jj][k];
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
					Ex(i,j,k) = Ex(i,j,k) + CB[modelDataX[i][j][k]] * psi_Exz_1[i][j][k];
					psi_Exz_2[i][j][kk] = be_z_2[kk] * psi_Exz_2[i][j][kk]
										+ ce_z_2[kk] * (Hy(i,j,k2) - Hy(i,j,k2+1)) / dz;
					Ex(i,j,k2) = Ex(i,j,k2) + CB[modelDataX[i][j][k2]] * psi_Exz_2[i][j][kk];
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
					Ey(i2,j,k) = Ey(i2,j,k) + CB[modelDataY[i2][j][k]] * psi_Eyx_2[ii][j][k];
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
					Ey(i,j,k) = Ey(i,j,k) + CB[modelDataY[i][j][k]] * psi_Eyz_1[i][j][k];
					psi_Eyz_2[i][j][kk] = be_z_2[kk] * psi_Eyz_2[i][j][kk]
										+ ce_z_2[kk] * (Hx(i,j,k2+1) - Hx(i,j,k2)) / dz;
					Ey(i,j,k2) = Ey(i,j,k2) + CB[modelDataY[i][j][k2]] * psi_Eyz_2[i][j][kk];
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
					Ez(i2,j,k) = Ez(i2,j,k) + CB[modelDataZ[i2][j][k]] * psi_Ezx_2[ii][j][k];
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
					Ez(i,j,k) = Ez(i,j,k) + CB[modelDataZ[i][j][k]] * psi_Ezy_1[i][j][k];
					psi_Ezy_2[i][jj][k] = be_y_2[jj] * psi_Ezy_2[i][jj][k]
										+ ce_y_2[jj] * (Hx(i,j2-1,k) - Hx(i,j2,k)) / dy;
					Ez(i,j2,k) = Ez(i,j2,k) + CB[modelDataZ[i][j2][k]] * psi_Ezy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}
}

void computeOneCPU()
{
	int i, j, k, ii, jj, kk, n;
	time_t tStart,tCurrent;
	double tStep=0.0;
	int dn = 10;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  BEGIN TIME STEP
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	printf("Begin time-stepping...\n");
	tStart = time(NULL);

	int i2, j2, k2;

	for(n = 1; n <= nMax; ++n)
	{
		if(n % dn == 0)
			printf("Step: %d, %fsec/step, %d/%d sec left/total.\n",
					n, tStep, (int)((nMax-n+1) * tStep), (int)(nMax * tStep) );

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  UPDATE Hx
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for(i = 0; i < Imax-1; ++i)
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

		for(i = 0; i < Imax-1; ++i)
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

		for(i = 0; i < Imax-1; ++i)
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

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  UPDATE Hy
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for(i = 0; i < Imax-1; ++i)
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

					psi_Hyx_2[ii][j][k] = bh_x_2[ii] * psi_Hyx_2[ii][j][k]
							+ ch_x_2[ii] * (Ez(i2+1,j,k) - Ez(i2,j,k)) / dx;
					Hy(i2,j,k) = Hy(i2,j,k) + DB * psi_Hyx_2[ii][j][k];
				}
			}
			--ii;
			++i2;
		}

		for(i = 0; i < Imax-1; ++i)
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

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  UPDATE Hz
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for(i = 0; i < Imax-1; ++i)
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
					psi_Hzx_2[ii][j][k] = bh_x_2[ii] * psi_Hzx_2[ii][j][k]
										+ ch_x_2[ii] * (Ey(i2,j,k) - Ey(i2+1,j,k))/dx;
					Hz(i2,j,k) = Hz(i2,j,k) + DB * psi_Hzx_2[ii][j][k];
				}
			}
			--ii;
			++i2;
		}

		for(i = 0; i < Imax-1; ++i)
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

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  UPDATE Ex
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for(i = 0; i < Imax-1; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					Ex(i,j,k) = CA[modelDataX[i][j][k]] * Ex(i,j,k) + CB[modelDataX[i][j][k]] *
								  ((Hz(i,j,k) - Hz(i,j-1,k)) * den_ey[j] +
								  (Hy(i,j,k) - Hy(i,j,k+1)) * den_ez[k] );
				}
			}
		}

		for(i = 0; i < Imax-1; ++i)
		{
			jj = nyPML_2 - 1;
			j2 = Jmax - nyPML_2;
			for(j = 1; j < nyPML_1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					psi_Exy_1[i][j][k] = be_y_1[j] * psi_Exy_1[i][j][k]
									   + ce_y_1[j] * (Hz(i,j,k) - Hz(i,j-1,k))/dy;
					Ex(i,j,k) = Ex(i,j,k) + CB[modelDataX[i][j][k]] * psi_Exy_1[i][j][k];
					psi_Exy_2[i][jj][k] = be_y_2[jj] * psi_Exy_2[i][jj][k]
										+ ce_y_2[jj] * (Hz(i,j2,k) - Hz(i,j2-1,k)) / dy;
					Ex(i,j2,k) = Ex(i,j2,k) + CB[modelDataX[i][j2][k]] * psi_Exy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}

		for(i = 0; i < Imax-1; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				kk = nzPML_2 - 1;
				k2 = Kmax - nzPML_2 - 1;
				for(k = 0; k < nzPML_1; ++k)
				{
					psi_Exz_1[i][j][k] = be_z_1[k] * psi_Exz_1[i][j][k]
									   + ce_z_1[k] * (Hy(i,j,k) - Hy(i,j,k+1)) / dz;
					Ex(i,j,k) = Ex(i,j,k) + CB[modelDataX[i][j][k]] * psi_Exz_1[i][j][k];
					psi_Exz_2[i][j][kk] = be_z_2[kk] * psi_Exz_2[i][j][kk]
										+ ce_z_2[kk] * (Hy(i,j,k2) - Hy(i,j,k2+1)) / dz;
					Ex(i,j,k2) = Ex(i,j,k2) + CB[modelDataX[i][j][k2]] * psi_Exz_2[i][j][kk];
					--kk;
					++k2;
				}
			}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  UPDATE Ey
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for(i = 1; i < Imax-1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					Ey(i,j,k) = CA[modelDataY[i][j][k]] * Ey(i,j,k) + CB[modelDataY[i][j][k]] *
								   ((Hz(i-1,j,k) - Hz(i,j,k)) * den_ex[i] +
								   (Hx(i,j,k+1) - Hx(i,j,k)) * den_ez[k] );
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
					Ey(i,j,k) = Ey(i,j,k) + CB[modelDataY[i][j][k]] * psi_Eyx_1[i][j][k];
					psi_Eyx_2[ii][j][k] = be_x_2[ii] * psi_Eyx_2[ii][j][k]
										+ ce_x_2[ii] * (Hz(i2-1,j,k) - Hz(i2,j,k)) / dx;
					Ey(i2,j,k) = Ey(i2,j,k) + CB[modelDataY[i2][j][k]] * psi_Eyx_2[ii][j][k];
				}
			}
			--ii;
			++i2;
		}

		for(i = 1; i < Imax-1; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				kk = nzPML_2 - 1;
				k2 = Kmax - nzPML_2 - 1;
				for(k = 0; k < nzPML_1; ++k)
				{
					psi_Eyz_1[i][j][k] = be_z_1[k] * psi_Eyz_1[i][j][k]
									   + ce_z_1[k] * (Hx(i,j,k+1) - Hx(i,j,k)) / dz;
					Ey(i,j,k) = Ey(i,j,k) + CB[modelDataY[i][j][k]] * psi_Eyz_1[i][j][k];
					psi_Eyz_2[i][j][kk] = be_z_2[kk] * psi_Eyz_2[i][j][kk]
										+ ce_z_2[kk] * (Hx(i,j,k2+1) - Hx(i,j,k2)) / dz;
					Ey(i,j,k2) = Ey(i,j,k2) + CB[modelDataY[i][j][k2]] * psi_Eyz_2[i][j][kk];
					--kk;
					++k2;
				}
			}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  UPDATE Ez
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for(i = 1; i < Imax-1; ++i)
		{
			for(j = 1; j < Jmax-1; ++j)
			{
				for(k = 1; k < Kmax-1; ++k)
				{
					Ez(i,j,k) = CA[modelDataZ[i][j][k]] * Ez(i,j,k) + CB[modelDataZ[i][j][k]]
								* ((Hy(i,j,k) - Hy(i-1,j,k)) * den_ex[i] +
								(Hx(i,j-1,k) - Hx(i,j,k)) * den_ey[j]);
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
					Ez(i,j,k) = Ez(i,j,k) + CB[modelDataZ[i][j][k]] * psi_Ezx_1[i][j][k];
					psi_Ezx_2[ii][j][k] = be_x_2[ii] * psi_Ezx_2[ii][j][k]
										+ ce_x_2[ii] * (Hy(i2,j,k) - Hy(i2-1,j,k)) / dx;
					Ez(i2,j,k) = Ez(i2,j,k) + CB[modelDataZ[i2][j][k]] * psi_Ezx_2[ii][j][k];
				}
			}
			--ii;
			++i2;
		}

		for(i = 1; i < Imax-1; ++i)
		{
			jj = nyPML_2 - 1;
			j2 = Jmax - nyPML_2;
			for(j = 1; j < nyPML_1; ++j)
			{
				for(k = 1; k < Kmax-1; ++k)
				{
					psi_Ezy_1[i][j][k] = be_y_1[j] * psi_Ezy_1[i][j][k]
									   + ce_y_1[j] * (Hx(i,j-1,k) - Hx(i,j,k)) / dy;
					Ez(i,j,k) = Ez(i,j,k) + CB[modelDataZ[i][j][k]] * psi_Ezy_1[i][j][k];
					psi_Ezy_2[i][jj][k] = be_y_2[jj] * psi_Ezy_2[i][jj][k]
										+ ce_y_2[jj] * (Hx(i,j2-1,k) - Hx(i,j2,k)) / dy;
					Ez(i,j2,k) = Ez(i,j2,k) + CB[modelDataZ[i][j2][k]] * psi_Ezy_2[i][jj][k];
				}
				--jj;
				++j2;
			}
		}

	//-----------------------------------------------------------
	//   Apply a point source (Hard)
	//-----------------------------------------------------------

		powerSource(n);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  WRITE TO OUTPUT FILES
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		saveData(n);
#ifdef _SAR
		for (i = 0; i<save_localSAR_amount+1; ++i)
		{
			computeRMS(&pSAR[i], n);
		}
#endif

		tCurrent = time(NULL);
		tStep = (tCurrent - tStart)/(double)(n + 1);
	}//  END TIME STEP
	free(Ex);
	free(Ey);
	free(Ez);
	free(Hx);
	free(Hy);
	free(Hz);

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

#ifdef _SAR
	if (nXgSAR)
	{
		computeXgSAR(nXgSAR);
	}
#endif
	tCurrent = time(NULL);
	printf("All computational tasks are completed.\n");
	printf("Total time consumption : %d second.\n", tCurrent - tStart);
	fprintf(fp_log, "Total time consumption : %d second.\n", tCurrent - tStart);
}

/*****************************************************************************************/
/*
 * Function name: radiationPower
 * Description: Compute radiation power
 * Parameters: 
 *			 rad_region		Location of the closed surface S
 * Return
 *			 Radiation power through the closed surface S
 */
#ifdef _POYNTING
double radiationPower(radiation_region rad_region)
{
	int xStart = rad_region.xStart;
	int xEnd = rad_region.xEnd;
	int yStart = rad_region.yStart;
	int yEnd = rad_region.yEnd;
	int zStart = rad_region.zStart;
	int zEnd = rad_region.zEnd;

	double ex1, ey1, ez1, hx1, hy1, hz1;
	double ex2, ey2, ez2, hx2, hy2, hz2;
	double radP_xy = 0, radP_yz = 0, radP_xz = 0;
	int i, j, k;

	if(xStart == -1 && xEnd == -1)
		return 0;

	for(i = xStart; i<=xEnd; ++i)
		for(j = yStart; j<=yEnd; ++j)
		{
			ex1 = ( ( (Ex(i,j,zStart-1)  +preEx(i,j,zStart-1)  )/2
					 +(Ex(i,j,zStart)    +preEx(i,j,zStart)    )/2 )/2
				  + ( (Ex(i-1,j,zStart-1)+preEx(i-1,j,zStart-1))/2
					 +(Ex(i-1,j,zStart)  +preEx(i-1,j,zStart)  )/2 )/2
				  )/2;

			ey1 = ( ( (Ey(i,j,zStart-1)  +preEy(i,j,zStart-1)  )/2
					 +(Ey(i,j,zStart)    +preEy(i,j,zStart)    )/2 )/2
				  + ( (Ey(i,j-1,zStart-1)+preEy(i,j-1,zStart-1))/2
					 +(Ey(i,j-1,zStart)  +preEy(i,j-1,zStart)  )/2 )/2
				  )/2;

			hx1 = ( Hx(i,j,zStart) + Hx(i,j-1,zStart) )/2;
			hy1 = ( Hy(i,j,zStart) + Hy(i-1,j,zStart) )/2;


			ex2 = ( ( (Ex(i,j,zEnd-1)  +preEx(i,j,zEnd-1)  )/2
					 +(Ex(i,j,zEnd)    +preEx(i,j,zEnd)    )/2 )/2
				  + ( (Ex(i-1,j,zEnd-1)+preEx(i-1,j,zEnd-1))/2
					 +(Ex(i-1,j,zEnd)  +preEx(i-1,j,zEnd)  )/2 )/2
				  )/2;
			ey2 = ( ( (Ey(i,j,zEnd-1)  +preEy(i,j,zEnd-1)  )/2
					 +(Ey(i,j,zEnd)    +preEy(i,j,zEnd)    )/2 )/2
				  + ( (Ey(i,j-1,zEnd-1)+preEy(i,j-1,zEnd-1))/2
					 +(Ey(i,j-1,zEnd)  +preEy(i,j-1,zEnd)  )/2 )/2
				  )/2;

			hx2 = ( Hx(i,j,zEnd) + Hx(i,j-1,zEnd) )/2;
			hy2 = ( Hy(i,j,zEnd) + Hy(i-1,j,zEnd) )/2;

			radP_xy += ex1*hy1-ey1*hx1 - (ex2*hy2-ey2*hx2);
		}
	radP_xy *= dx*dy;

	if (rad_region.computeX_1)
	{
		for(j = yStart; j<=yEnd; ++j)
			for(k = zStart; k<=zEnd; ++k)
			{
				ey1 = ( ( (Ey(xStart,j,k-1)  +preEy(xStart,j,k-1)  )/2
						 +(Ey(xStart,j,k)    +preEy(xStart,j,k)    )/2 )/2
					  + ( (Ey(xStart,j-1,k-1)+preEy(xStart,j-1,k-1))/2
						 +(Ey(xStart,j-1,k)  +preEy(xStart,j-1,k)  )/2 )/2
					  )/2;
				ez1 = ( Ez(xStart,j,k) + preEz(xStart,j,k) )/2;

				hy1 = ( Hy(xStart,j,k) + Hy(xStart-1,j,k) )/2;
				hz1 = ( (Hz(xStart,j,k)+Hz(xStart-1,j,k)+Hz(xStart,j-1,k)+Hz(xStart-1,j-1,k))/4
					  + (Hz(xStart,j,k-1)+Hz(xStart-1,j,k-1)+Hz(xStart,j-1,k-1)+Hz(xStart-1,j-1,k-1))/4
					  )/2;

				radP_yz += ey1*hz1-ez1*hy1;
			}
	}
	if (rad_region.computeX_2)
	{
		for(j = yStart; j<=yEnd; ++j)
			for(k = zStart; k<=zEnd; ++k)
			{
				ey2 = ( ( (Ey(xEnd,j,k-1)  +preEy(xEnd,j,k-1)  )/2
						 +(Ey(xEnd,j,k)    +preEy(xEnd,j,k)    )/2 )/2
					  + ( (Ey(xEnd,j-1,k-1)+preEy(xEnd,j-1,k-1))/2
						 +(Ey(xEnd,j-1,k)  +preEy(xEnd,j-1,k)  )/2 )/2
					  )/2;
				ez2 = ( Ez(xEnd,j,k) + preEz(xEnd,j,k) )/2;

				hy2 = ( Hy(xEnd,j,k) + Hy(xEnd-1,j,k) )/2;
				hz2 = ( (Hz(xEnd,j,k)+Hz(xEnd-1,j,k)+Hz(xEnd,j-1,k)+Hz(xEnd-1,j-1,k))/4
					  + (Hz(xEnd,j,k-1)+Hz(xEnd-1,j,k-1)+Hz(xEnd,j-1,k-1)+Hz(xEnd-1,j-1,k-1))/4
					  )/2;

				radP_yz += -(ey2*hz2-ez2*hy2);
			}
	}
	radP_yz *= dy*dz;

	for(i = xStart; i<=xEnd; ++i)
		for(k = zStart; k<=zEnd; ++k)
		{
			ex1 = ( ( (Ex(i,yStart,k-1)  +preEx(i,yStart,k-1)  )/2
					 +(Ex(i,yStart,k)    +preEx(i,yStart,k)    )/2 )/2
				  + ( (Ex(i-1,yStart,k-1)+preEx(i-1,yStart,k-1))/2
					 +(Ex(i-1,yStart,k)  +preEx(i-1,yStart,k)  )/2 )/2
				  )/2;
			ez1 = ( Ez(i,yStart,k) + preEz(i,yStart,k) )/2;

			hx1 = ( Hx(i,yStart,k) + Hx(i,yStart-1,k) )/2;
			hz1 = ( (Hz(i,yStart,k)+Hz(i-1,yStart,k)+Hz(i,yStart-1,k)+Hz(i-1,yStart-1,k))/4
				  + (Hz(i,yStart,k-1)+Hz(i-1,yStart,k-1)+Hz(i,yStart-1,k-1)+Hz(i-1,yStart-1,k-1))/4
				  )/2;

			ex2 = ( ( (Ex(i,yEnd,k-1)  +preEx(i,yEnd,k-1)  )/2
					 +(Ex(i,yEnd,k)    +preEx(i,yEnd,k)    )/2 )/2
				  + ( (Ex(i-1,yEnd,k-1)+preEx(i-1,yEnd,k-1))/2
					 +(Ex(i-1,yEnd,k)  +preEx(i-1,yEnd,k)  )/2 )/2
				  )/2;
			ez2 = ( Ez(i,yEnd,k) + preEz(i,yEnd,k) )/2;

			hx2 = ( Hx(i,yEnd,k) + Hx(i,yEnd-1,k) )/2;
			hz2 = ( (Hz(i,yEnd,k)+Hz(i-1,yEnd,k)+Hz(i,yEnd-1,k)+Hz(i-1,yEnd-1,k))/4
				  + (Hz(i,yEnd,k-1)+Hz(i-1,yEnd,k-1)+Hz(i,yEnd-1,k-1)+Hz(i-1,yEnd-1,k-1))/4
				  )/2;

			radP_xz += ez1*hx1-ex1*hz1 - (ez2*hx2-ex2*hz2);
		}
	radP_xz *= dx*dz;

	return (radP_xy+radP_xz+radP_yz);
}
#endif

void computeOneCPU_Mur2()
{
	int i, j, k, n;
	time_t tStart,tCurrent;
	double tStep=0.0;
	int dn = 10;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  BEGIN TIME STEP
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	printf("Begin time-stepping...\n");
	tStart = time(NULL);

	/* Mur2 */
	int is, ie, js, je, ks, ke;
	double C1, C2, C3;

	is = 0;
	ie = Imax - 1;
	js = 0;
	je = Jmax - 1;
	ks = 0;
	ke = Kmax - 1;

	C1 = (C*dt-dx)/(C*dt+dx);
	C2 = 2*(C*dt-dx)/dx;
	C3 = (C*C*dt*dt)/2/dx/(C*dt+dx);

	for(n = 1; n <= nMax; ++n)
	{
		if(n % dn == 0)
			printf("Step: %d, %fsec/step, %d/%d sec left/total.\n",
					n, tStep, (int)((nMax-n+1) * tStep), (int)(nMax * tStep) );

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  UPDATE Hx
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for(i = is; i <= ie; ++i)
		{
			for (j = js; j <= je-1; ++j)
			{
				for (k = ks+1; k <= ke-1; ++k)
				{
					Hx(i,j,k) = DA * Hx(i,j,k) + DB *
						((Ez(i,j,k) - Ez(i,j+1,k)) / dy  +
						(Ey(i,j,k) - Ey(i,j,k-1)) / dz );
				}
			}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  UPDATE Hy
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for(i = is; i <= ie-1; ++i)
		{
			for(j = js; j <= je; ++j)
			{
				for(k = ks+1; k <= ke-1; ++k)
				{
					Hy(i,j,k) = DA * Hy(i,j,k) + DB *
						((Ez(i+1,j,k) - Ez(i,j,k)) / dx +
						(Ex(i,j,k-1) - Ex(i,j,k)) / dz );
				}
			}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  UPDATE Hz
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for(i = is; i <= ie-1; ++i)
		{
			for(j = js; j<= je-1; ++j)
			{
				for(k = ks; k <= ke-1; ++k)
				{
					Hz(i,j,k) = DA * Hz(i,j,k) + DB
								* ((Ey(i,j,k) - Ey(i+1,j,k)) / dx +
								(Ex(i,j+1,k) - Ex(i,j,k)) / dy);
				}
			}
		}

		/* Mur2 save past value*/
		for(i = 0; i < Imax-1; ++i)
		{
			for(j = 0; j < Jmax; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					Ex_n_1(i,j,k) = Ex_n(i,j,k);
					Ex_n(i,j,k) = Ex(i,j,k);
				}
			}
		}
		for(i = 0; i < Imax; ++i)
		{
			for(j = 0; j < Jmax-1; ++j)
			{
				for(k = 0; k < Kmax-1; ++k)
				{
					Ey_n_1(i,j,k) = Ey_n(i,j,k);
					Ey_n(i,j,k) = Ey(i,j,k);
				}
			}
		}
		for(i = 0; i < Imax; ++i)
		{
			for(j = 0; j < Jmax; ++j)
			{
				for(k = 0; k < Kmax; ++k)
				{
					Ez_n_1(i,j,k) = Ez_n(i,j,k);
					Ez_n(i,j,k) = Ez(i,j,k);
				}
			}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  UPDATE Ex
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for(i = is; i <= ie-1; ++i)
		{
			for(j = js+1; j <= je-1; ++j)
			{
				for(k = ks; k <= ke-1; ++k)
				{
					Ex(i,j,k) = CA[modelDataX[i][j][k]] * Ex(i,j,k) + CB[modelDataX[i][j][k]] *
								  ((Hz(i,j,k) - Hz(i,j-1,k)) / dy +
								  (Hy(i,j,k) - Hy(i,j,k+1)) / dz );
				}
			}
		}

		/* Exy */
		j = js;
		for (i = is+1; i <= ie-2; ++i)
		{
			for (k = ks+2; k <= ke-3; ++k)
			{
				Ex(i,j,k) = -Ex_n_1(i,j+1,k) + C1*(Ex(i,j+1,k) + Ex_n_1(i,j,k)) - 
							C2*(Ex_n(i,j,k) + Ex_n(i,j+1,k)) + 
							C3*(Ex_n(i-1,j,k)+Ex_n(i+1,j,k)+Ex_n(i,j,k-1)+Ex_n(i,j,k+1) + 
								Ex_n(i-1,j+1,k)+Ex_n(i+1,j+1,k)+Ex_n(i,j+1,k-1)+Ex_n(i,j+1,k+1));
			}
		}
		k = ks + 1;
		for (i = is; i <= ie-1; ++i)
		{
			Ex(i,j,k) = Ex_n(i,j+1,k) + C1*(Ex(i,j+1,k) - Ex_n(i,j,k));
		}
		k = ke - 2;
		for (i = is; i <= ie-1; ++i)
		{
			Ex(i,j,k) = Ex_n(i,j+1,k) + C1*(Ex(i,j+1,k) - Ex_n(i,j,k));
		}
		i = is;
		for (k = ks+2; k <= ke-3; ++k)
		{
			Ex(i,j,k) = Ex_n(i,j+1,k) + C1*(Ex(i,j+1,k) - Ex_n(i,j,k));
		}
		i = ie - 1;
		for (k = ks+2; k <= ke-3; ++k)
		{
			Ex(i,j,k) = Ex_n(i,j+1,k) + C1*(Ex(i,j+1,k) - Ex_n(i,j,k));
		}
		
		j = je;
		for (i = is+1; i <= ie-2; ++i)
		{
			for (k = ks+2; k <= ke-3; ++k)
			{
				Ex(i,j,k) = -Ex_n_1(i,j-1,k) + C1*(Ex(i,j-1,k) + Ex_n_1(i,j,k)) - 
							C2*(Ex_n(i,j,k) + Ex_n(i,j-1,k)) + 
							C3*(Ex_n(i-1,j,k)+Ex_n(i+1,j,k)+Ex_n(i,j,k-1)+Ex_n(i,j,k+1) + 
								Ex_n(i-1,j-1,k)+Ex_n(i+1,j-1,k)+Ex_n(i,j-1,k-1)+Ex_n(i,j-1,k+1));
			}
		}
		k = ks + 1;
		for (i = is; i <= ie-1; ++i)
		{
			Ex(i,j,k) = Ex_n(i,j-1,k) + C1*(Ex(i,j-1,k) - Ex_n(i,j,k));
		}
		k = ke - 2;
		for (i = is; i <= ie-1; ++i)
		{
			Ex(i,j,k) = Ex_n(i,j-1,k) + C1*(Ex(i,j-1,k) - Ex_n(i,j,k));
		}
		i = is;
		for (k = ks+2; k <= ke-3; ++k)
		{
			Ex(i,j,k) = Ex_n(i,j-1,k) + C1*(Ex(i,j-1,k) - Ex_n(i,j,k));
		}
		i = ie - 1;
		for (k = ks+2; k <= ke-3; ++k)
		{
			Ex(i,j,k) = Ex_n(i,j-1,k) + C1*(Ex(i,j-1,k) - Ex_n(i,j,k));
		}

		/* Exz */
		k = ks;
		for (i = is+1; i <= ie-2; ++i)
		{
			for (j = js+2; j <= je-2; ++j)
			{
				Ex(i,j,k) = -Ex_n_1(i,j,k+1) + C1*(Ex(i,j,k+1) + Ex_n_1(i,j,k)) - 
							C2*(Ex_n(i,j,k) + Ex_n(i,j,k+1)) + 
							C3*(Ex_n(i-1,j,k)+Ex_n(i+1,j,k)+Ex_n(i,j-1,k)+Ex_n(i,j+1,k) + 
								Ex_n(i-1,j,k+1)+Ex_n(i+1,j,k+1)+Ex_n(i,j-1,k+1)+Ex_n(i,j+1,k+1));
			}
		}
		j = js + 1;
		for (i = is; i <= ie-1; ++i)
		{
			Ex(i,j,k) = Ex_n(i,j,k+1) + C1*(Ex(i,j,k+1) - Ex_n(i,j,k));
		}
		j = je - 1;
		for (i = is; i <= ie-1; ++i)
		{
			Ex(i,j,k) = Ex_n(i,j,k+1) + C1*(Ex(i,j,k+1) - Ex_n(i,j,k));
		}
		i = is;
		for (j = js+2; j <= je-2; ++j)
		{
			Ex(i,j,k) = Ex_n(i,j,k+1) + C1*(Ex(i,j,k+1) - Ex_n(i,j,k));
		}
		i = ie - 1;
		for (j = js+2; j <= je-2; ++j)
		{
			Ex(i,j,k) = Ex_n(i,j,k+1) + C1*(Ex(i,j,k+1) - Ex_n(i,j,k));
		}

		k = ke-1;
		for (i = is+1; i <= ie-2; ++i)
		{
			for (j = js+2; j <= je-2; ++j)
			{
				Ex(i,j,k) = -Ex_n_1(i,j,k-1) + C1*(Ex(i,j,k-1) + Ex_n_1(i,j,k)) - 
							C2*(Ex_n(i,j,k) + Ex_n(i,j,k-1)) + 
							C3*(Ex_n(i-1,j,k)+Ex_n(i+1,j,k)+Ex_n(i,j-1,k)+Ex_n(i,j+1,k) + 
								Ex_n(i-1,j,k-1)+Ex_n(i+1,j,k-1)+Ex_n(i,j-1,k-1)+Ex_n(i,j+1,k-1));
			}
		}
		j = js + 1;
		for (i = is; i <= ie-1; ++i)
		{
			Ex(i,j,k) = Ex_n(i,j,k-1) + C1*(Ex(i,j,k-1) - Ex_n(i,j,k));
		}
		j = je - 1;
		for (i = is; i <= ie-1; ++i)
		{
			Ex(i,j,k) = Ex_n(i,j,k-1) + C1*(Ex(i,j,k-1) - Ex_n(i,j,k));
		}
		i = is;
		for (j = js+2; j <= je-2; ++j)
		{
			Ex(i,j,k) = Ex_n(i,j,k-1) + C1*(Ex(i,j,k-1) - Ex_n(i,j,k));
		}
		i = ie - 1;
		for (j = js+2; j <= je-2; ++j)
		{
			Ex(i,j,k) = Ex_n(i,j,k-1) + C1*(Ex(i,j,k-1) - Ex_n(i,j,k));
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  UPDATE Ey
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for(i = is+1; i <= ie-1; ++i)
		{
			for(j = js; j <= je-1; ++j)
			{
				for(k = ks+1; k <= ke-2; ++k)
				{
					Ey(i,j,k) = CA[modelDataY[i][j][k]] * Ey(i,j,k) + CB[modelDataY[i][j][k]] *
								   ((Hz(i-1,j,k) - Hz(i,j,k)) / dx +
								   (Hx(i,j,k+1) - Hx(i,j,k)) / dz );
				}
			}
		}

		/* Eyx */
		i = is;
		for (j = js+1; j <= je-2; ++j)
		{
			for (k = ks+2; k <= ke-3; ++k)
			{
				Ey(i,j,k) = -Ey_n_1(i+1,j,k) + C1*(Ey(i+1,j,k) + Ey_n_1(i,j,k)) - 
							C2*(Ey_n(i,j,k) + Ey_n(i+1,j,k)) + 
							C3*(Ey_n(i,j-1,k)+Ey_n(i,j+1,k)+Ey_n(i,j,k-1)+Ey_n(i,j,k+1) + 
								Ey_n(i+1,j-1,k)+Ey_n(i+1,j+1,k)+Ey_n(i+1,j,k-1)+Ey_n(i+1,j,k+1));
			}
		}
		k = ks + 1;
		for (j = js; j <= je-1; ++j)
		{
			Ey(i,j,k) = Ey_n(i+1,j,k) + C1*(Ey(i+1,j,k) - Ey_n(i,j,k));
		}
		k = ke - 2;
		for (j = js; j <= je-1; ++j)
		{
			Ey(i,j,k) = Ey_n(i+1,j,k) + C1*(Ey(i+1,j,k) - Ey_n(i,j,k));
		}
		j = js;
		for (k = ks+2; k <= ke-3; ++k)
		{
			Ey(i,j,k) = Ey_n(i+1,j,k) + C1*(Ey(i+1,j,k) - Ey_n(i,j,k));
		}
		j = je - 1;
		for (k = ks+2; k <= ke-3; ++k)
		{
			Ey(i,j,k) = Ey_n(i+1,j,k) + C1*(Ey(i+1,j,k) - Ey_n(i,j,k));
		}

		i = ie;
		for (j = js+1; j <= je-2; ++j)
		{
			for (k = ks+2; k <= ke-3; ++k)
			{
				Ey(i,j,k) = -Ey_n_1(i-1,j,k) + C1*(Ey(i-1,j,k) + Ey_n_1(i,j,k)) - 
							C2*(Ey_n(i,j,k) + Ey_n(i-1,j,k)) + 
							C3*(Ey_n(i,j-1,k)+Ey_n(i,j+1,k)+Ey_n(i,j,k-1)+Ey_n(i,j,k+1) + 
								Ey_n(i-1,j-1,k)+Ey_n(i-1,j+1,k)+Ey_n(i-1,j,k-1)+Ey_n(i-1,j,k+1));
			}
		}
		k = ks + 1;
		for (j = js; j <= je-1; ++j)
		{
			Ey(i,j,k) = Ey_n(i-1,j,k) + C1*(Ey(i-1,j,k) - Ey_n(i,j,k));
		}
		k = ke - 2;
		for (j = js; j <= je-1; ++j)
		{
			Ey(i,j,k) = Ey_n(i-1,j,k) + C1*(Ey(i-1,j,k) - Ey_n(i,j,k));
		}
		j = js;
		for (k = ks+2; k <= ke-3; ++k)
		{
			Ey(i,j,k) = Ey_n(i-1,j,k) + C1*(Ey(i-1,j,k) - Ey_n(i,j,k));
		}
		j = je - 1;
		for (k = ks+2; k <= ke-3; ++k)
		{
			Ey(i,j,k) = Ey_n(i-1,j,k) + C1*(Ey(i-1,j,k) - Ey_n(i,j,k));
		}

		/* Eyz */
		k = ks;
		for (i = is+2; i <= ie-2; ++i)
		{
			for (j = js+1; j <= je-2; ++j)
			{
				Ey(i,j,k) = -Ey_n_1(i,j,k+1) + C1*(Ey(i,j,k+1) + Ey_n_1(i,j,k)) - 
							C2*(Ey_n(i,j,k) + Ey_n(i,j,k+1)) + 
							C3*(Ey_n(i,j-1,k)+Ey_n(i,j+1,k)+Ey_n(i-1,j,k)+Ey_n(i+1,j,k) + 
								Ey_n(i,j-1,k+1)+Ey_n(i,j+1,k+1)+Ey_n(i-1,j,k+1)+Ey_n(i+1,j,k+1));
			}
		}
		i = is + 1;
		for (j = js; j <= je-1; ++j)
		{
			Ey(i,j,k) = Ey_n(i,j,k+1) + C1*(Ey(i,j,k+1) - Ey_n(i,j,k));
		}
		i = ie - 1;
		for (j = js; j <= je-1; ++j)
		{
			Ey(i,j,k) = Ey_n(i,j,k+1) + C1*(Ey(i,j,k+1) - Ey_n(i,j,k));
		}
		j = js;
		for (i = is+2; i <= ie-2; ++i)
		{
			Ey(i,j,k) = Ey_n(i,j,k+1) + C1*(Ey(i,j,k+1) - Ey_n(i,j,k));
		}
		j = je - 1;
		for (i = is+2; i <= ie-2; ++i)
		{
			Ey(i,j,k) = Ey_n(i,j,k+1) + C1*(Ey(i,j,k+1) - Ey_n(i,j,k));
		}

		k = ke-1;
		for (i = is+2; i <= ie-2; ++i)
		{
			for (j = js+1; j <= je-2; ++j)
			{
				Ey(i,j,k) = -Ey_n_1(i,j,k-1) + C1*(Ey(i,j,k-1) + Ey_n_1(i,j,k)) - 
							C2*(Ey_n(i,j,k) + Ey_n(i,j,k-1)) + 
							C3*(Ey_n(i,j-1,k)+Ey_n(i,j+1,k)+Ey_n(i-1,j,k)+Ey_n(i+1,j,k) + 
								Ey_n(i,j-1,k-1)+Ey_n(i,j+1,k-1)+Ey_n(i-1,j,k-1)+Ey_n(i+1,j,k-1));
			}
		}
		i = is + 1;
		for (j = js; j <= je-1; ++j)
		{
			Ey(i,j,k) = Ey_n(i,j,k-1) + C1*(Ey(i,j,k-1) - Ey_n(i,j,k));
		}
		i = ie - 1;
		for (j = js; j <= je-1; ++j)
		{
			Ey(i,j,k) = Ey_n(i,j,k-1) + C1*(Ey(i,j,k-1) - Ey_n(i,j,k));
		}
		j = js;
		for (i = is+2; i <= ie-2; ++i)
		{
			Ey(i,j,k) = Ey_n(i,j,k-1) + C1*(Ey(i,j,k-1) - Ey_n(i,j,k));
		}
		j = je - 1;
		for (i = is+2; i <= ie-2; ++i)
		{
			Ey(i,j,k) = Ey_n(i,j,k-1) + C1*(Ey(i,j,k-1) - Ey_n(i,j,k));
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//  UPDATE Ez
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		for(i = is+1; i <= ie-1; ++i)
		{
			for(j = js+1; j <= je-1; ++j)
			{
				for(k = ks+1; k <= ke-1; ++k)
				{
					Ez(i,j,k) = CA[modelDataZ[i][j][k]] * Ez(i,j,k) + CB[modelDataZ[i][j][k]]
								* ((Hy(i,j,k) - Hy(i-1,j,k)) / dx +
								(Hx(i,j-1,k) - Hx(i,j,k)) / dy);
				}
			}
		}

		/* Ezx */
		i = is;
		for (j = js+2; j <= je-2; ++j)
		{
			for (k = ks+2; k <= ke-2; ++k)
			{
				Ez(i,j,k) = -Ez_n_1(i+1,j,k) + C1*(Ez(i+1,j,k) + Ez_n_1(i,j,k)) - 
							C2*(Ez_n(i,j,k) + Ez_n(i+1,j,k)) + 
							C3*(Ez_n(i,j-1,k)+Ez_n(i,j+1,k)+Ez_n(i,j,k-1)+Ez_n(i,j,k+1) + 
								Ez_n(i+1,j-1,k)+Ez_n(i+1,j+1,k)+Ez_n(i+1,j,k-1)+Ez_n(i+1,j,k+1));
			}
		}
		j = js + 1;
		for (k = ks+1; k <= ke-1; ++k)
		{
			Ez(i,j,k) = Ez_n(i+1,j,k) + C1*(Ez(i+1,j,k) - Ez_n(i,j,k));
		}
		j = je - 1;
		for (k = ks+1; k <= ke-1; ++k)
		{
			Ez(i,j,k) = Ez_n(i+1,j,k) + C1*(Ez(i+1,j,k) - Ez_n(i,j,k));
		}
		k = ks+1;
		for (j = js+2; j <= je-2; ++j)
		{
			Ez(i,j,k) = Ez_n(i+1,j,k) + C1*(Ez(i+1,j,k) - Ez_n(i,j,k));
		}
		k = ke-1;
		for (j = js+2; j <= je-2; ++j)
		{
			Ez(i,j,k) = Ez_n(i+1,j,k) + C1*(Ez(i+1,j,k) - Ez_n(i,j,k));
		}

		i = ie;
		for (j = js+2; j <= je-2; ++j)
		{
			for (k = ks+2; k <= ke-2; ++k)
			{
				Ez(i,j,k) = -Ez_n_1(i-1,j,k) + C1*(Ez(i-1,j,k) + Ez_n_1(i,j,k)) - 
							C2*(Ez_n(i,j,k) + Ez_n(i-1,j,k)) + 
							C3*(Ez_n(i,j-1,k)+Ez_n(i,j+1,k)+Ez_n(i,j,k-1)+Ez_n(i,j,k+1) + 
								Ez_n(i-1,j-1,k)+Ez_n(i-1,j+1,k)+Ez_n(i-1,j,k-1)+Ez_n(i-1,j,k+1));
			}
		}
		j = js + 1;
		for (k = ks+1; k <= ke-1; ++k)
		{
			Ez(i,j,k) = Ez_n(i-1,j,k) + C1*(Ez(i-1,j,k) - Ez_n(i,j,k));
		}
		j = je - 1;
		for (k = ks+1; k <= ke-1; ++k)
		{
			Ez(i,j,k) = Ez_n(i-1,j,k) + C1*(Ez(i-1,j,k) - Ez_n(i,j,k));
		}
		k = ks + 1;
		for (j = js+2; j <= je-2; ++j)
		{
			Ez(i,j,k) = Ez_n(i-1,j,k) + C1*(Ez(i-1,j,k) - Ez_n(i,j,k));
		}
		k = ke - 1;
		for (j = js+2; j <= je-2; ++j)
		{
			Ez(i,j,k) = Ez_n(i-1,j,k) + C1*(Ez(i-1,j,k) - Ez_n(i,j,k));
		}

		/* Ezy */
		j = js;
		for (i = is+2; i <= ie-2; ++i)
		{
			for (k = ks+2; k <= ke-2; ++k)
			{
				Ez(i,j,k) = -Ez_n_1(i,j+1,k) + C1*(Ez(i,j+1,k) + Ez_n_1(i,j,k)) - 
							C2*(Ez_n(i,j,k) + Ez_n(i,j+1,k)) + 
							C3*(Ez_n(i-1,j,k)+Ez_n(i+1,j,k)+Ez_n(i,j,k-1)+Ez_n(i,j,k+1) + 
								Ez_n(i-1,j+1,k)+Ez_n(i+1,j+1,k)+Ez_n(i,j+1,k-1)+Ez_n(i,j+1,k+1));
			}
		}
		i = is + 1;
		for (k = ks+1; k <= ke-1; ++k)
		{
			Ez(i,j,k) = Ez_n(i,j+1,k) + C1*(Ez(i,j+1,k) - Ez_n(i,j,k));
		}
		i = ie - 1;
		for (k = ks+1; k <= ke-1; ++k)
		{
			Ez(i,j,k) = Ez_n(i,j+1,k) + C1*(Ez(i,j+1,k) - Ez_n(i,j,k));
		}
		k = ks + 1;
		for (i = is+2; i <= ie-2; ++i)
		{
			Ez(i,j,k) = Ez_n(i,j+1,k) + C1*(Ez(i,j+1,k) - Ez_n(i,j,k));
		}
		k = ke - 1;
		for (i = is+2; i <= ie-2; ++i)
		{
			Ez(i,j,k) = Ez_n(i,j+1,k) + C1*(Ez(i,j+1,k) - Ez_n(i,j,k));
		}

		j = je;
		for (i = is+2; i <= ie-2; ++i)
		{
			for (k = ks+2; k <= ke-2; ++k)
			{
				Ez(i,j,k) = -Ez_n_1(i,j-1,k) + C1*(Ez(i,j-1,k) + Ez_n_1(i,j,k)) - 
							C2*(Ez_n(i,j,k) + Ez_n(i,j-1,k)) + 
							C3*(Ez_n(i-1,j,k)+Ez_n(i+1,j,k)+Ez_n(i,j,k-1)+Ez_n(i,j,k+1) + 
								Ez_n(i-1,j-1,k)+Ez_n(i+1,j-1,k)+Ez_n(i,j-1,k-1)+Ez_n(i,j-1,k+1));
			}
		}
		i = is + 1;
		for (k = ks+1; k <= ke-1; ++k)
		{
			Ez(i,j,k) = Ez_n(i,j-1,k) + C1*(Ez(i,j-1,k) - Ez_n(i,j,k));
		}
		i = ie - 1;
		for (k = ks+1; k <= ke-1; ++k)
		{
			Ez(i,j,k) = Ez_n(i,j-1,k) + C1*(Ez(i,j-1,k) - Ez_n(i,j,k));
		}
		k = ks + 1;
		for (i = is+2; i <= ie-2; ++i)
		{
			Ez(i,j,k) = Ez_n(i,j-1,k) + C1*(Ez(i,j-1,k) - Ez_n(i,j,k));
		}
		k = ke - 1;
		for (i = is+2; i <= ie-2; ++i)
		{
			Ez(i,j,k) = Ez_n(i,j-1,k) + C1*(Ez(i,j-1,k) - Ez_n(i,j,k));
		}

	//-----------------------------------------------------------
	//   Apply a point source (Hard)
	//-----------------------------------------------------------

		powerSource(n);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//  WRITE TO OUTPUT FILES
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		saveData(n);
#ifdef _SAR
		for (i = 0; i<save_localSAR_amount+1; ++i)
		{
			computeRMS(&pSAR[i], n);
		}
#endif

		tCurrent = time(NULL);
		tStep = (tCurrent - tStart)/(double)(n + 1);
	}//  END TIME STEP

	free(Ex);
	free(Ey);
	free(Ez);
	free(Hx);
	free(Hy);
	free(Hz);

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

#ifdef _SAR
	if (nXgSAR)
	{
		computeXgSAR(nXgSAR);
	}
#endif		
	tCurrent = time(NULL);
	printf("All computational tasks are completed.\n");
	printf("Total time consumption : %d second.\n", tCurrent - tStart);
	fprintf(fp_log, "Total time consumption : %d second.\n", tCurrent - tStart);
}

void saveMaxE() {
    int i, j, k;
    double currE = 0;
    for (i = outputRangeE[0]; i<outputRangeE[1]; ++i) {
        for (j = outputRangeE[2]; j<outputRangeE[3]; ++j) {
            for (k = outputRangeE[4]; k<outputRangeE[5]; ++k) {
                currE = sqrt( pow(Ex(i, j, k-1), 2) + pow(Ey(i, j, k-1), 2) + pow(Ez(i, j, k), 2) );
                if (currE > (*maxE)[i][j][k])
                    (*maxE)[i][j][k] = currE;
            }
        }
    }
}

