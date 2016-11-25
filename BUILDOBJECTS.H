/*****************************************************************************************/
/*
 * Function name: buildCuboid
 * Description: Creat a cube
 * Parameters:
 *			 Iloc、Jloc、Kloc			: The nearest point from the origin
 *			 Ichange、Jchange、kchange	: Size
 *			 epsilon、sigma				: Relative dielectric constant, Conductivity
 * Return
 */
void buildCuboid(int Iloc, int Jloc, int Kloc,
				 int Ichange, int Jchange, int Kchange,
				 unsigned char mediaNo) 
{
	int i, j, k;
#ifdef _SAR
if (flag == 13)/* Will be executed before compute SAR */
{
	if (myrank == 0)
	{
		for(i = Iloc; i < Iloc + Ichange; ++i)
		{
			for(j = Jloc; j < Jloc + Jchange; ++j)
			{
				for(k = Kloc; k < Kloc + Kchange; ++k)
				{
					rho[i-1][j-1][k-1] = media[mediaNo].rho*dx*dy*dz;
				}
			}
		}
	}
}
#endif

	int change_by_is = 0;
	if (myrank == 0)
		change_by_is = -is;
	else
		change_by_is = -is + 1;
	Iloc = Iloc + paddingX_1 + thicknessOfPml;
	Jloc = Jloc + paddingY_1 + thicknessOfPml;
	Kloc = Kloc + paddingZ_1 + thicknessOfPml;

	if(myrank != 0)
		is = is-1;

	for(i = Iloc; i < Iloc + Ichange; ++i) 
	{
		if(is<=i && i<=ie+1)
		{
			for(j = Jloc; j < Jloc + Jchange; ++j) 
			{
				for(k = Kloc; k < Kloc + Kchange; ++k) 
				{
					yeeCube (i+change_by_is, j, k, mediaNo);
				}
			}
		}
	}

	if(myrank != 0)
		is = is+1;
}

/*
 * Function name: buildPlane
 * Description: Creat a Plane
 * Parameters:
 *			 (p1x,p1y,p1z)			: Starting Coordinate
 *			 (p2x,p2y,p2z)			: Ending coordinate
 *			 epsilon、sigma			: Relative dielectric constant, Conductivity
 * Return
 */
void buildPlane(int p1x, int p1y, int p1z,
				int p2x, int p2y, int p2z,
				unsigned char mediaNo) 
{
	int i, j, k;

	int change_by_is = 0;
	if (myrank == 0)
		change_by_is = -is;
	else
		change_by_is = -is + 1;
	p1x = p1x + paddingX_1 + thicknessOfPml;
	p1y = p1y + paddingX_1 + thicknessOfPml;
	p1z = p1z + paddingX_1 + thicknessOfPml;
	p2x = p2x + paddingX_1 + thicknessOfPml;
	p2y = p2y + paddingX_1 + thicknessOfPml;
	p2z = p2z + paddingX_1 + thicknessOfPml;

	if(myrank != 0)
		is = is-1;

	if (p1x == p2x && p1y != p2y && p1z != p2z && is<=p1x && p1x<=ie+1)
	{
		for (j = p1y; j <= p2y; ++j)
		{
			for (k = p1z; k <= p2z; ++k)
			{
				k = k-1;
				modelDataY[p1x+change_by_is][j][k] = mediaNo;
				modelDataY[p1x+change_by_is][j][k+1] = mediaNo;
				k = k+1;

				modelDataZ[p1x+change_by_is][j][k] = mediaNo;
				modelDataZ[p1x+change_by_is][j+1][k] = mediaNo;
			}
		}
	}

	if (p1x != p2x && p1y == p2y && p1z != p2z)
	{
		for (i = p1x; i <= p2x; ++i)
		{
			if (is<=i && i<=ie+1)
			{
				for (k = p1z; k <= p2z; ++k)
				{
					k = k-1;
					modelDataX[i+change_by_is][p1y][k] = mediaNo;
					modelDataX[i+change_by_is][p1y][k+1] = mediaNo;
					k = k+1;

					modelDataZ[i+change_by_is][p1y][k] = mediaNo;
					modelDataZ[i+change_by_is+1][p1y][k] = mediaNo;
				}
			}
		}
	}

	if (p1x != p2x && p1y != p2y && p1z == p2z)
	{
		for (i = p1x; i <= p2x; ++i)
		{
			if (is<=i && i<=ie+1)
			{
				for (j = p1y; j <= p2y; ++j)
				{
					p1z = p1z-1;
					modelDataX[i+change_by_is][j][p1z] = mediaNo;
					modelDataX[i+change_by_is][j+1][p1z] = mediaNo;

					modelDataY[i+change_by_is][j][p1z] = mediaNo;
					modelDataY[i+change_by_is+1][j][p1z] = mediaNo;
					p1z = p1z+1;
				}
			}
		}
	}

	if(myrank != 0)
		is = is+1;
}
/*
void buildPlaneCircle(int px, int py, int pz,
				      double radius, int direction,
				      unsigned char mediaNo) 
{
	int i, j, k;

	int change_by_is = 0;
	if (myrank == 0)
		change_by_is = -is;
	else
		change_by_is = -is + 1;
	px = px + paddingX_1 + thicknessOfPml;
	py = py + paddingX_1 + thicknessOfPml;
	pz = pz + paddingX_1 + thicknessOfPml;

	if(myrank != 0)
		is = is-1;

	if (p1x == p2x && p1y != p2y && p1z != p2z && is<=p1x && p1x<=ie+1)
	{
		for (j = p1y; j <= p2y; ++j)
		{
			for (k = p1z; k <= p2z; ++k)
			{
				k = k-1;
				modelDataY[p1x+change_by_is][j][k] = mediaNo;
				modelDataY[p1x+change_by_is][j][k+1] = mediaNo;
				k = k+1;

				modelDataZ[p1x+change_by_is][j][k] = mediaNo;
				modelDataZ[p1x+change_by_is][j+1][k] = mediaNo;
			}
		}
	}

	if(myrank != 0)
		is = is+1;
}*/

/*****************************************************************************************/
/*
 * Function name: buildLine
 * Description: Creat a Line 
 * Parameters:
 *			 (p1x,p1y,p1z)		: Starting Coordinate
 *			 (p2x,p2y,p2z)		: Ending coordinate
 * Return
 */
void buildLine(int p1x, int p1y, int p1z, int p2x, int p2y, int p2z, unsigned char mediaNo) 
{
	int i, j, k;
	p1x = p1x + paddingX_1 + thicknessOfPml;
	p1y = p1y + paddingY_1 + thicknessOfPml;
	p1z = p1z + paddingZ_1 + thicknessOfPml;
	p2x = p2x + paddingX_1 + thicknessOfPml;
	p2y = p2y + paddingY_1 + thicknessOfPml;
	p2z = p2z + paddingZ_1 + thicknessOfPml;
	int change_by_is = 0;
	if (myrank == 0)
		change_by_is = -is;
	else
		change_by_is = -is + 1;

	for(i = p1x; i <= p2x; ++i)
	{
		if(is<=i && i<=ie+1)
		{
			for(j = p1y; j <= p2y; ++j)
			{
				for(k = p1z; k <= p2z; ++k)
				{
					if (p1x != p2x)
					{
						k = k-1;
						modelDataX[i+change_by_is][j][k] = mediaNo;
						k = k+1;
					}
					else if (p1y != p2y)
					{
						k = k-1;
						modelDataY[i+change_by_is][j][k] = mediaNo;
						k = k+1;
					}
					else /* p1z != p2z */
					{
						modelDataZ[i+change_by_is][j][k] = mediaNo;
					}
				}
			}
		}
	}
}

void creatDipoleAntenna(_DIPOLE_ANTENNA _dipole_antenna)
{
    double length_of_cell = 0.0;
	switch(_dipole_antenna.direction)
	{
	case 1 :
		buildLine(_dipole_antenna.feed_x-_dipole_antenna.length_low, _dipole_antenna.feed_y, _dipole_antenna.feed_z,
				  _dipole_antenna.feed_x-1, _dipole_antenna.feed_y, _dipole_antenna.feed_z, 1);
		buildLine(_dipole_antenna.feed_x+1, _dipole_antenna.feed_y, _dipole_antenna.feed_z,
				  _dipole_antenna.feed_x+_dipole_antenna.length_high, _dipole_antenna.feed_y, _dipole_antenna.feed_z, 1);
        length_of_cell = dx;
		break;
	case 2 :
		buildLine(_dipole_antenna.feed_x, _dipole_antenna.feed_y-_dipole_antenna.length_low, _dipole_antenna.feed_z,
				  _dipole_antenna.feed_x, _dipole_antenna.feed_y-1, _dipole_antenna.feed_z, 1);
		buildLine(_dipole_antenna.feed_x, _dipole_antenna.feed_y+1, _dipole_antenna.feed_z,
				  _dipole_antenna.feed_x, _dipole_antenna.feed_y+_dipole_antenna.length_high, _dipole_antenna.feed_z, 1);
        length_of_cell = dy;
		break;
	case 3 :
		buildLine(_dipole_antenna.feed_x, _dipole_antenna.feed_y, _dipole_antenna.feed_z-_dipole_antenna.length_low,
				  _dipole_antenna.feed_x, _dipole_antenna.feed_y, _dipole_antenna.feed_z-1, 1);
		buildLine(_dipole_antenna.feed_x, _dipole_antenna.feed_y, _dipole_antenna.feed_z+1,
				  _dipole_antenna.feed_x, _dipole_antenna.feed_y, _dipole_antenna.feed_z+_dipole_antenna.length_high, 1);
        length_of_cell = dz;
		break;
	}

	if(_dipole_antenna.length_low == _dipole_antenna.length_high)/* Symmetric Oscillators */
	{
		double available_power = 0.0;
		double net_input_power = 0.0;
		double Rr = 0.0;
		double namuda = C/freq;
		double beta = 2.0*pi/namuda;

		if(_dipole_antenna.length_low*length_of_cell/namuda < 0.125)
			Rr = 20*pow((beta * _dipole_antenna.length_low*length_of_cell), 2);
		else if(_dipole_antenna.length_low*length_of_cell/namuda < 0.25)
			Rr = 24.7*pow((beta * _dipole_antenna.length_low*length_of_cell), 2.4);
		else if(_dipole_antenna.length_low*length_of_cell/namuda < 0.319)
			Rr = 11.14*pow((beta * _dipole_antenna.length_low*length_of_cell), 4.17);
		else
            Rr = 1;
        if (freq == 0)
            Rr = 1;

		if(_dipole_antenna.impedance && myrank == 0)
		{
			available_power = pow(amp/(_dipole_antenna.impedance+_dipole_antenna.impedance), 2)*_dipole_antenna.impedance/2;
			net_input_power = pow(amp/(_dipole_antenna.impedance+Rr), 2)*Rr/2;
			printf("Creat DipoleAntenna Success.\n");
			printf("Feed Point : (%4d, %4d, %4d), Available Power = %lf W.\n", _isource, _jsource, _ksource, available_power);
			printf("                                 Net input power = %lf W.\n", net_input_power);
			printf("                                 Radiation efficiency = %.2f%%\n", net_input_power/available_power*100.0);
            fprintf(fp_log, "Creat DipoleAntenna Success.\n");
			fprintf(fp_log, "Feed Point : (%4d, %4d, %4d), Available Power = %lf W.\n", _isource, _jsource, _ksource, available_power);
			fprintf(fp_log, "                                 Net input power = %lf W.\n", net_input_power);
			fprintf(fp_log, "                                 Radiation efficiency = %.2f%%\n", net_input_power/available_power*100.0);
            fflush(stdout);
		}

		amp = amp * Rr /(Rr + _dipole_antenna.impedance);

	}
}

unsigned char createNewMedia(int id, double sig, double eps, double rho, double spec_heat)
{
	media[id].sigma= sig;
	media[id].epsilon = eps;
	media[id].rho= rho;
	media[id].spec_heat = spec_heat;
	CA[id] = ((1.0 - sig * dt / (2.0 * eps * epsO))
				/ (1.0 + sig * dt / (2.0 * eps * epsO)));
	CB[id] = (dt / (eps * epsO))
				/ (1.0 + sig * dt / (2.0 * eps * epsO));
	if (id>=maxMedia) maxMedia = id+1;
	return (id);
}

