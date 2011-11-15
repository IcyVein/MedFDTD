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
		source = amp * sin(2 * 4*atan(1.0) * freq * dt * timestep);
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

	switch (sourceType)
	{
	case 0 :
		if (waveForm == -1 || waveForm == 0)
		{
			if (port == 'x' || port == 'X')
				Ex(isource,jsource,ksource-1) = source / dx;
			else if (port == 'y' || port == 'Y')
				Ey(isource,jsource,ksource-1) = source / dy;
			else
				Ez(isource,jsource,ksource) = source / dz;
		}
		else /* Soft source when no sine source */
		{
			if (port == 'x' || port == 'X')
				Ex(isource,jsource,ksource-1) += source / dx;
			else if (port == 'y' || port == 'Y')
				Ey(isource,jsource,ksource-1) += source / dy;
			else
				Ez(isource,jsource,ksource) += source / dz;
		}
		break;
	case 1 :
		if (_isource && _jsource == 0 && _ksource == 0 && myrank == 0)
		{
			if (port == 'y')
			{
				for (j = 0; j<Jmax-1; ++j)
					for (k = 0; k<Kmax-1; ++k)
					{
						Ey(is+thicknessOfPml+paddingX_1+1,j,k) = source;
					}
			}
			else if (port == 'z')
			{
				for (j = 0; j<Jmax-1; ++j)
					for (k = 0; k<Kmax-1; ++k)
					{
						Ez(is+thicknessOfPml+paddingX_1+1,j,k) = source;
					}
			}
		}
		else if ( _isource == 0 && _jsource && _ksource == 0)
		{
			if (port == 'x')
			{
				is_by_rank = is;
				ie_by_rank = ie;

				for (i=is_by_rank; i<=ie_by_rank; ++i)
					for (k=0; k<Kmax-1; ++k)
					{
						Ex(i,thicknessOfPml+paddingY_1+1,k) = source;
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
						Ez(i,thicknessOfPml+paddingY_1+1,k) = source;
					}
			}
		}
		else if ( _isource == 0 && _jsource == 0 && _ksource)
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
						Ex(i,j,thicknessOfPml+paddingZ_1+1) = source;
					}
			}
			else if (port == 'y')
			{
				for (i=is_by_rank; i<=ie_by_rank; ++i)
					for (j=0; j<Jmax-1; ++j)
					{
						Ey(i,j,thicknessOfPml+paddingZ_1+1) = source;
					}
			}
		}
		break;
	case 2 :
		source = (dt/epsO)*(1/dz/dz/dz)*(1e-10*exp(-pow((dt*timestep-t0*dt)/pulse_width/dt,2)))*(-2*(dt*timestep-t0*dt)/pulse_width/pulse_width/dt/dt);
		Ez(isource,jsource,ksource) -= source;
		break;
	default :
		printf("Error : sourceTpye = %d\n", sourceType);
		fflush(stdout);
	}
}

