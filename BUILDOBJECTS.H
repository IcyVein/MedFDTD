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
				 double epsilon, double sigma, double _rho = 0.0) 
{
	int i, j, k;
#ifdef _SAR
	int Iloc0 = Iloc + paddingX_1;
	int Jloc0 = Jloc + paddingY_1;
	int Kloc0 = Kloc + paddingZ_1;
	if (myrank == 0)
	{
		for(i = Iloc0; i < Iloc0 + Ichange; ++i)
		{
			for(j = Jloc0; j < Jloc0 + Jchange; ++j)
			{
				for(k = Kloc0; k < Kloc0 + Kchange; ++k)
				{
					rho[i][j][k] = _rho;
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
					yeeCube (i+change_by_is, j, k, epsilon, sigma, _rho);
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
				double epsilon, double sigma, double _rho = 0) 
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
				ID2_EpsR[p1x+change_by_is][j][k] = epsilon;
				ID2_SigR[p1x+change_by_is][j][k] = sigma;
				ID2_EpsR[p1x+change_by_is][j][k+1] = epsilon;
				ID2_SigR[p1x+change_by_is][j][k+1] = sigma;
				k = k+1;

				ID3_EpsR[p1x+change_by_is][j][k] = epsilon;
				ID3_SigR[p1x+change_by_is][j][k] = sigma;
				ID3_EpsR[p1x+change_by_is][j+1][k] = epsilon;
				ID3_SigR[p1x+change_by_is][j+1][k] = sigma;
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
					ID1_EpsR[i+change_by_is][p1y][k] = epsilon;
					ID1_SigR[i+change_by_is][p1y][k] = sigma;
					ID1_EpsR[i+change_by_is][p1y][k+1] = epsilon;
					ID1_SigR[i+change_by_is][p1y][k+1] = sigma;
					k = k+1;

					ID3_EpsR[i+change_by_is][p1y][k] = epsilon;
					ID3_SigR[i+change_by_is][p1y][k] = sigma;
					ID3_EpsR[i+change_by_is+1][p1y][k] = epsilon;
					ID3_SigR[i+change_by_is+1][p1y][k] = sigma;
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
					ID1_EpsR[i+change_by_is][j][p1z] = epsilon;
					ID1_SigR[i+change_by_is][j][p1z] = sigma;
					ID1_EpsR[i+change_by_is][j+1][p1z] = epsilon;
					ID1_SigR[i+change_by_is][j+1][p1z] = sigma;

					ID2_EpsR[i+change_by_is][j][p1z] = epsilon;
					ID2_SigR[i+change_by_is][j][p1z] = sigma;
					ID2_EpsR[i+change_by_is+1][j][p1z] = epsilon;
					ID2_SigR[i+change_by_is+1][j][p1z] = sigma;
					p1z = p1z+1;
				}
			}
		}
	}

	if(myrank != 0)
		is = is+1;
}

/*****************************************************************************************/
/*
 * Function name: buildWire
 * Description: Creat a Wire
 * Parameters:
 *			 (p1x,p1y,p1z)		: Starting Coordinate
 *			 (p2x,p2y,p2z)		: Ending coordinate
 * Return
 */
void buildWire(int p1x, int p1y, int p1z, int p2x, int p2y, int p2z) 
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
						ID1_EpsR[i+change_by_is][j][k] = -1;
						ID1_SigR[i+change_by_is][j][k] = 0 ;
						k = k+1;
					}
					else if (p1y != p2y)
					{
						k = k-1;
						ID2_EpsR[i+change_by_is][j][k] = -1;
						ID2_SigR[i+change_by_is][j][k] = 0 ;
						k = k+1;
					}
					else /* p1z != p2z */
					{
						ID3_EpsR[i+change_by_is][j][k] = -1;
						ID3_SigR[i+change_by_is][j][k] = 0 ;
					}
				}
			}
		}
	}
}

void creatDipoleAntenna(_DIPOLE_ANTENNA _dipole_antenna)
{
	switch(_dipole_antenna.direction)
	{
	case 1 :
		buildWire(_dipole_antenna.feed_x-_dipole_antenna.length_low, _dipole_antenna.feed_y, _dipole_antenna.feed_z,
				  _dipole_antenna.feed_x-1, _dipole_antenna.feed_y, _dipole_antenna.feed_z);
		buildWire(_dipole_antenna.feed_x+1, _dipole_antenna.feed_y, _dipole_antenna.feed_z,
				  _dipole_antenna.feed_x+_dipole_antenna.length_high, _dipole_antenna.feed_y, _dipole_antenna.feed_z);
		break;
	case 2 :
		buildWire(_dipole_antenna.feed_x, _dipole_antenna.feed_y-_dipole_antenna.length_low, _dipole_antenna.feed_z,
				  _dipole_antenna.feed_x, _dipole_antenna.feed_y-1, _dipole_antenna.feed_z);
		buildWire(_dipole_antenna.feed_x, _dipole_antenna.feed_y+1, _dipole_antenna.feed_z,
				  _dipole_antenna.feed_x, _dipole_antenna.feed_y+_dipole_antenna.length_high, _dipole_antenna.feed_z);
		break;
	case 3 :
		buildWire(_dipole_antenna.feed_x, _dipole_antenna.feed_y, _dipole_antenna.feed_z-_dipole_antenna.length_low,
				  _dipole_antenna.feed_x, _dipole_antenna.feed_y, _dipole_antenna.feed_z-1);
		buildWire(_dipole_antenna.feed_x, _dipole_antenna.feed_y, _dipole_antenna.feed_z+1,
				  _dipole_antenna.feed_x, _dipole_antenna.feed_y, _dipole_antenna.feed_z+_dipole_antenna.length_high);
		break;
	}

	if(_dipole_antenna.length_low == _dipole_antenna.length_high)/* Symmetric Oscillators */
	{
		double power_antenna = 0;
		double Rr = 1;
		double namuda = C/freq;
		double beta = 2*pi/namuda;

		if(_dipole_antenna.length_low*LENGTH_OF_CELL/namuda < 0.125)
			Rr = 20*pow((beta * _dipole_antenna.length_low*LENGTH_OF_CELL), 2);
		else if(_dipole_antenna.length_low*LENGTH_OF_CELL/namuda < 0.25)
			Rr = 24.7*pow((beta * _dipole_antenna.length_low*LENGTH_OF_CELL), 2.4);
		else if(_dipole_antenna.length_low*LENGTH_OF_CELL/namuda < 0.319)
			Rr = 11.14*pow((beta * _dipole_antenna.length_low*LENGTH_OF_CELL), 4.17);
		else;

		if(Rr && myrank == 0)
		{
			power_antenna = amp * amp / Rr /8;
			printf("Creat DipoleAntenna Success.\n");
			printf("Feed Point : (%d, %d, %d), Power = %lf W.\n", _isource, _jsource, _ksource, power_antenna);
		}

		amp = amp * Rr /sqrt(Rr*Rr + _dipole_antenna.impedance*_dipole_antenna.impedance);

	}
}

