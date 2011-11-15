#ifdef _TF_SF
double* E_1D;
double* H_1D;
int spaceE, spaceH;
double ds;
double CA, CB;

int initialize1D();
int powersource1D(int);
int compute1D(int);

#endif

#ifdef _TF_SF

int initialize1D()
{
	CA = 1;
	CB = dt / epsO;
	if (_isource)
	{
		spaceE = spaceX+3;
		spaceH = spaceX+2;
		ds = dx;
	}
	else if (_jsource)
	{
		spaceE = spaceY+3;
		spaceH = spaceY+2;
		ds = dy;
	}
	else
	{
		spaceE = spaceZ+3;
		spaceH = spaceZ+2;
		ds = dz;
	}
	E_1D = initArray(spaceE);
	H_1D = initArray(spaceH);

	if (E_1D && H_1D)
		return 0;
	else
		return MEM_ERROR;
}

int powersource1D(int timestep)
{
	double source;
	t0=50;
	pulse_width = 40;

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
		source = exp(-(4*pi*(timestep-t0)*(timestep-t0)/pulse_width/pulse_width));
		break;
	case 2 :
		if (timestep <= pulse_width)
			source = (1-cos(2*pi*timestep*dt/pulse_width))/2;
		else
			source = 0;
		break;
	case 3 :
		source = (timestep - t0)/pulse_width*exp(-(4*pi*(timestep-t0)*(timestep-t0)/pulse_width/pulse_width));
		break;
	case 4 :
		if (timestep < t0)
			source = (10 - 15*cos(2*pi*timestep*dt)+6*cos(4*pi*timestep*dt)-cos(6*pi*timestep*dt));
		else
			source=0;
		break;
	case 5 :
		source = -cos(2 * 4*atan(1.0) * freq * dt * timestep)*exp(-(4*pi*(timestep-t0)*(timestep-t0)/pulse_width/pulse_width));
		break;
	}

	if (waveForm == 0)
	{
		int T = 1/(freq*dt);
		source = timestep <= 2*T? source * (timestep/(double)T/2) : source;
	}

	E_1D[padding+1] = -source;

	return 0;
}

int compute1D(int timestep)
{
	int i;
	double E_1D_n1, E_1D_nn;
	/* H */
	for (i = 0; i<spaceH; ++i)
	{
		H_1D[i] = DA*H_1D[i] + DB*(E_1D[i] - E_1D[i+1])/ds;
	}
	/* E */
	E_1D_n1 = E_1D[1];
	E_1D_nn = E_1D[spaceE-2];
	for (i = 1; i<spaceE-1; ++i)
	{
		E_1D[i] = CA*E_1D[i] + CB*(H_1D[i-1] - H_1D[i])/ds;
	}

	/* ABCs(Mur2) */
	E_1D[0] = E_1D_n1 + (C*dt - ds)/(C*dt + ds) * (E_1D[1] - E_1D[0]);
	E_1D[spaceE-1] = E_1D_nn + (C*dt - ds)/(C*dt + ds) * (E_1D[spaceE-2] - E_1D[spaceE-1]);

	return 0;
}

#endif

