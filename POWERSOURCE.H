/*****************************************************************************************/
/*
 * Function name: powersource
 * Description: Exciation
 * Parameters: 
 *			 timestep
 * Return: 
 */
void powerSource(int timestep) 
{
	int i, j, k;
	int is_by_rank, ie_by_rank;
    double source;
    int T = (int)(1/(freq*dt));
    for (int nSrc = 0; nSrc<sourceNum; ++nSrc) {

    if (myrank<=1) //前方end-rings
    {
	    source = genSource(timestep + nSrc*T/sourceNum);
        //cout<<"rank = "<<myrank<<", src num = "<<sourceNum<<" @beg"<<endl;
    }
    else if (myrank >= nprocs-2) //后方end-rings，极性反向
    {
	    source = -genSource(timestep + nSrc*T/sourceNum);
    }
    else //Legs
    {
        source = genSource(timestep + nSrc*T/sourceNum/2.0 + T/4); //相位为方位角+pi/2
    }

	if (waveForm == 0)
	{
		source = timestep <= 2*T? source * (timestep/(double)T/2) : source;
	}

	switch (sourceType)
	{
	case 0 :
		if (waveForm == -1 || waveForm == 0)
		{
			if (port == 'x' || port == 'X')
				Ex(isource[nSrc],jsource[nSrc],ksource[nSrc]-1) = source / dx;
			else if (port == 'y' || port == 'Y')
				Ey(isource[nSrc],jsource[nSrc],ksource[nSrc]-1) = source / dy;
			else
				Ez(isource[nSrc],jsource[nSrc],ksource[nSrc]) = source / dz;
		}
		else /* Soft source when no sine source */
		{
			if (port == 'x' || port == 'X')
				Ex(isource[nSrc],jsource[nSrc],ksource[nSrc]-1) += source / dx;
			else if (port == 'y' || port == 'Y')
				Ey(isource[nSrc],jsource[nSrc],ksource[nSrc]-1) += source / dy;
			else
				Ez(isource[nSrc],jsource[nSrc],ksource[nSrc]) = source / dz;
		}
		break;
	case 1 :
		if (_isource[nSrc] && _jsource[nSrc] == 0 && _ksource[nSrc] == 0 && myrank == 0)
		{
			if (port == 'y')
			{
				for (j = thicknessOfPml+1+planeWaveIndex; j<Jmax-(thicknessOfPml+1+planeWaveIndex); ++j)
					for (k = thicknessOfPml+1+planeWaveIndex; k<Kmax-(thicknessOfPml+1+planeWaveIndex); ++k)
					{
						Ey(is+thicknessOfPml+1,j,k) += source;
					}
			}
			else if (port == 'z')
			{
				for (j = 0/*thicknessOfPml+1+5*/; j<Jmax/*(thicknessOfPml+1+5)*/; ++j)
					for (k = 1/*thicknessOfPml+1+5*/; k<Kmax-1/*(thicknessOfPml+1+5)*/; ++k)
					{
						Ez(is+thicknessOfPml+1+planeWaveIndex,j,k) += source;// PlaneWaveEz1
					}
			}
		}
		else if ( _isource[nSrc] == 0 && _jsource[nSrc] && _ksource[nSrc] == 0)
		{
			if (port == 'x')
			{
				is_by_rank = is;
				ie_by_rank = ie;

				for (i=is_by_rank; i<=ie_by_rank; ++i)
					for (k=0; k<Kmax-1; ++k)
					{
						Ex(i,thicknessOfPml+1,k) += source;
					}
			}
			else if (port == 'z')
			{
				if (myrank == 0)
				{
					is_by_rank = is;
					ie_by_rank = ie+1;
				}
				else if (myrank == nprocs - 1)
				{
					is_by_rank = is-1;
					ie_by_rank = ie;
				}
				else
				{
					is_by_rank = is-1;
					ie_by_rank = ie+1;
				}

				for (i=is_by_rank; i<=ie_by_rank; ++i)
					for (k=0; k<Kmax-1; ++k)
					{
						Ez(i,thicknessOfPml+1,k) += source;
					}
			}
		}
		else if ( _isource[nSrc] == 0 && _jsource[nSrc] == 0 && _ksource[nSrc])
		{
			if (myrank == 0)
			{
				is_by_rank = is;
				ie_by_rank = ie+1;
			}
			else if (myrank == nprocs - 1)
			{
				is_by_rank = is-1;
				ie_by_rank = ie;
			}
			else
			{
				is_by_rank = is-1;
				ie_by_rank = ie+1;
			}
			if (port == 'x')
			{
				for (i=is_by_rank; i<=ie_by_rank; ++i)
					for (j=0; j<Jmax-1; ++j)
					{
						Ex(i,j,thicknessOfPml+1) += source;
					}
			}
			else if (port == 'y')
			{
				for (i=is_by_rank; i<=ie_by_rank; ++i)
					for (j=0; j<Jmax-1; ++j)
					{
						Ey(i,j,thicknessOfPml+1) += source;
					}
			}
		}
		break;
	case 2 :
		source = (dt/epsO)*(1/dz/dz/dz)*(1e-10*exp(-pow((dt*timestep-t0*dt)/pulse_width/dt,2)))*(-2*(dt*timestep-t0*dt)/pulse_width/pulse_width/dt/dt);
		Ez(isource[nSrc],jsource[nSrc],ksource[nSrc]) -= source;
		break;
    case 3:
        for (i = isource[nSrc]; i<=isource[nSrc]; ++i)
            for (j = jsource[nSrc]; j<=jsource[nSrc]; ++j)
                for (k = ksource[nSrc]-3; k<=ksource[nSrc]; ++k)
                    Ez(i, j, k) = source;
        break;
    case 4: //线电流源
        for (k = ksource[nSrc] - pulse_width; k<=ksource[nSrc] + pulse_width; ++k)
        {
            Ez(isource[nSrc], jsource[nSrc], k) = source * (dt/epsO/dx/dy/dz);
        }
    case 5: //面电流源
        for (k = ksource[nSrc]-7; k<=ksource[nSrc]+7; ++k)
            for (j = jsource[nSrc]-7; j<=jsource[nSrc]+7; ++j)
            {
                if ((pow(abs(k-ksource[nSrc]), 2) + pow(abs(j-jsource[nSrc]), 2))<=49)
                    Ez(isource[nSrc], j, k) = source * (dt/epsO/dy/dz);
            }
        break;
	default :
		printf("Error : sourceTpye = %d\n", sourceType);
		fflush(stdout);
	}
    }
}

void powerSourcePlaneWaveH(int timestep)
{
	int i, j, k;
	int is_by_rank, ie_by_rank;
	double source;

/* TODO :
	switch (waveForm)
	{
	case 0 :
		source = amp * sin(2 * 4*atan(1.0) * freq * dt * (timestep - timeShift + (dx/C/dt)*spaceShift));
		break;
 * case 1 :
		source = amp * exp(-(4*pi*(timestep-t0)*(timestep-t0)/(double)pulse_width/(double)pulse_width));
		break;
	case 2 :
		if (timestep <= pulse_width)
			source = amp * (1-cos(2*pi*timestep*dt/pulse_width))/2;
		else
			source = 0;
		break;
	case 3 :
		source = amp * (timestep - t0)/pulse_width*exp(-(4*pi*(timestep-t0)*(timestep-t0)/pulse_width/pulse_width));
		break;
	case 4 :
		if (timestep < t0)
			source = amp * (10 - 15*cos(2*pi*timestep*dt)+6*cos(4*pi*timestep*dt)-cos(6*pi*timestep*dt));
		else
			source=0;
		break;
	case 5 :
		source = -amp * cos(2 * 4*atan(1.0) * freq * dt * timestep)*exp(-(4*pi*(timestep-t0)*(timestep-t0)/pulse_width/pulse_width));
		break;
	}
	if (waveForm == 0)
	{
		int T = (int)(1/(freq*dt));
		source = timestep <= 2*T? source * (timestep/(double)T/2) : source;
	}
*/

    for (int nSrc = 0; nSrc<sourceNum; ++nSrc) {
    if (_isource[nSrc] && _jsource[nSrc] == 0 && _ksource[nSrc] == 0)
    {
        if (myrank == 0)
        {
        	is_by_rank = is+thicknessOfPml+1+4;
        	ie_by_rank = ie+1;
        }
        else if (myrank == nprocs - 1)
        {
        	is_by_rank = is-1;
        	ie_by_rank = ie;
        }
        else
        {
        	is_by_rank = is-1;
        	ie_by_rank = ie+1;
        }
        if (port == 'z')
        {
            for (i = is_by_rank; i<=ie_by_rank; ++i)
                for (k = thicknessOfPml+1+planeWaveIndex; k<Kmax-(thicknessOfPml+1+planeWaveIndex); ++k)
                {
                    if (timestep > _global_is + i - (thicknessOfPml+1+planeWaveIndex-1))
                    {
                        //Hx(i, thicknessOfPml+1+4, k) -= getSource(timestep, 1/2, -1/2*(i-(thicknessOfPml+1+5)))/sqrt(muO/epsO)/sqrt(1);
                        //Hx(i, Jmax-1-(thicknessOfPml+1+5), k) -= getSource(timestep, 1/2, -1/2*(i-(thicknessOfPml+1+5)))/sqrt(muO/epsO)/sqrt(1);
                    }
                }
        }
    }

    if (_isource[nSrc] && _jsource[nSrc] == 0 && _ksource[nSrc] == 0 && myrank == 0)
    {
    	if (port == 'y')
    	{
    		for (j = 0; j<Jmax-1; ++j)
    			for (k = 0; k<Kmax-1; ++k)
    			{
    				//Ey(is+thicknessOfPml+1,j,k) += source;
                    Hz(is+thicknessOfPml+1+planeWaveIndex-1,j,k) -= getSource(timestep, 1/2, -1/2)/sqrt(muO/epsO)/sqrt(1);
    			}
    	}
    	else if (port == 'z')
    	{
    		for (j = 0/*thicknessOfPml+1+5*/; j<Jmax/*(thicknessOfPml+1+5)*/; ++j)
    			for (k = 1/*thicknessOfPml+1+5*/; k<Kmax-1/*(thicknessOfPml+1+5)*/; ++k)
    			{
    				//Ez(is+thicknessOfPml+1,j,k) += source;
                    Hy(is+thicknessOfPml+1+planeWaveIndex-1,j,k) -= getSource(timestep, 1/2, -1/2)/sqrt(muO/epsO)/sqrt(1);// PlaneWaveEz2
    			}
    	}
    }
    }
    /*
    else if ( _isource[nSrc] == 0 && _jsource[nSrc] && _ksource[nSrc] == 0)
    {
    	if (port == 'x')
    	{
    		is_by_rank = is;
    		ie_by_rank = ie;
    
    		for (i=is_by_rank; i<=ie_by_rank; ++i)
    			for (k=0; k<Kmax-1; ++k)
    			{
    				Ex(i,thicknessOfPml+1,k) += source;
    			}
    	}
    	else if (port == 'z')
    	{
    		if (myrank == 0)
    		{
    			is_by_rank = is;
    			ie_by_rank = ie+1;
    		}
    		else if (myrank == nprocs - 1)
    		{
    			is_by_rank = is-1;
    			ie_by_rank = ie;
    		}
    		else
    		{
    			is_by_rank = is-1;
    			ie_by_rank = ie+1;
    		}
    
    		for (i=is_by_rank; i<=ie_by_rank; ++i)
    			for (k=0; k<Kmax-1; ++k)
    			{
    				Ez(i,thicknessOfPml+1,k) += source;
    			}
    	}
    }
    else if ( _isource[nSrc] == 0 && _jsource[nSrc] == 0 && _ksource[nSrc])
    {
    	if (myrank == 0)
    	{
    		is_by_rank = is;
    		ie_by_rank = ie+1;
    	}
    	else if (myrank == nprocs - 1)
    	{
    		is_by_rank = is-1;
    		ie_by_rank = ie;
    	}
    	else
    	{
    		is_by_rank = is-1;
    		ie_by_rank = ie+1;
    	}
    	if (port == 'x')
    	{
    		for (i=is_by_rank; i<=ie_by_rank; ++i)
    			for (j=0; j<Jmax-1; ++j)
    			{
    				Ex(i,j,thicknessOfPml+1) += source;
    			}
    	}
    	else if (port == 'y')
    	{
    		for (i=is_by_rank; i<=ie_by_rank; ++i)
    			for (j=0; j<Jmax-1; ++j)
    			{
    				Ey(i,j,thicknessOfPml+1) += source;
    			}
    	}
    }*/
}

double getSource(int timeStep, double timeShift, double spaceShift)
{
    double source = amp * sin(2 * 4*atan(1.0) * freq * dt * (timeStep - timeShift + (dx/C/dt)*spaceShift));
    int T = (int)(1/(freq*dt));
    source = timeStep <= 2*T? source * (timeStep/(double)T/2) : source;
    return source;
}

double genSource(int timestep)
{
    double source = 0;
    switch (waveForm)
	{
	case -1 :
		if (timestep-1 < sizeSRC)
		{
			source = power[timestep-1];
		}
		else
			source = 0;
		break;
	case 0 :
		source = amp * sin(2 * pi * freq * dt * timestep);
		break;
	case 1 :
		source = amp * exp(-(4*pi*(timestep-t0)*(timestep-t0)/(double)pulse_width/(double)pulse_width));
		break;
	case 2 :
		if (timestep <= pulse_width)
			source = amp * (1-cos(2*pi*timestep*dt/pulse_width))/2;
		else
			source = 0;
		break;
	case 3 :
		source = (double)amp * (timestep - t0)/(double)pulse_width*exp(-(4.0*pi*(timestep-t0)*(timestep-t0)/(double)pulse_width/(double)pulse_width));
		break;
	case 4 :
		if (timestep < t0)
			source = amp * (10 - 15*cos(2*pi*timestep*dt)+6*cos(4*pi*timestep*dt)-cos(6*pi*timestep*dt));
		else
			source=0;
		break;
	case 5 :
		source = -amp * cos(2 * 4*atan(1.0) * freq * dt * timestep)*exp(-(4*pi*(timestep-t0)*(timestep-t0)/pulse_width/pulse_width));
		break;
    case 6 :
        source = amp * (exp(-0.0101*timestep*dt*1e9) - exp(-0.0349*timestep*dt*1e9)); // s --> ns
        break;
	}

    return source;
}
