/*****************************************************************************************/
/*
 * Function name: saveData
 * Description:	Save The Field Data
 * Parameters: 
 *				timestep
 * Return:
 */
void saveData(int timestep)
{
	int i;
	for (i = 0; i < save_plane_amount; i++)
	{
		if (timestep >= fp_save_field_file[i].sp.start && 
			timestep <= fp_save_field_file[i].sp.end &&
			timestep % fp_save_field_file[i].sp.step ==0)
		{
			switch (fp_save_field_file[i].sp.plane_no)
			{
				case 1 :
					if (nprocs > 1)
						writeField_xy(fp_save_field_file[i].ex,
								  fp_save_field_file[i].ey,
								  fp_save_field_file[i].ez,
								  fp_save_field_file[i].hx,
								  fp_save_field_file[i].hy,
								  fp_save_field_file[i].hz,
								  fp_save_field_file[i].sp.slice);
					else
						writeField_xyOneCPU(fp_save_field_file[i].ex,
								  fp_save_field_file[i].ey,
								  fp_save_field_file[i].ez,
								  fp_save_field_file[i].hx,
								  fp_save_field_file[i].hy,
								  fp_save_field_file[i].hz,
								  fp_save_field_file[i].sp.slice);
					break;
				case 2 :
					if (nprocs > 1)
						writeField_yz(fp_save_field_file[i].ex,
								  fp_save_field_file[i].ey,
								  fp_save_field_file[i].ez,
								  fp_save_field_file[i].hx,
								  fp_save_field_file[i].hy,
								  fp_save_field_file[i].hz,
								  fp_save_field_file[i].sp.slice);
					else
						writeField_yzOneCPU(fp_save_field_file[i].ex,
								  fp_save_field_file[i].ey,
								  fp_save_field_file[i].ez,
								  fp_save_field_file[i].hx,
								  fp_save_field_file[i].hy,
								  fp_save_field_file[i].hz,
								  fp_save_field_file[i].sp.slice);
					break;
				case 3 :
					if (nprocs > 1)
						writeField_xz(fp_save_field_file[i].ex,
								  fp_save_field_file[i].ey,
								  fp_save_field_file[i].ez,
								  fp_save_field_file[i].hx,
								  fp_save_field_file[i].hy,
								  fp_save_field_file[i].hz,
								  fp_save_field_file[i].sp.slice);
					else
						writeField_xzOneCPU(fp_save_field_file[i].ex,
								  fp_save_field_file[i].ey,
								  fp_save_field_file[i].ez,
								  fp_save_field_file[i].hx,
								  fp_save_field_file[i].hy,
								  fp_save_field_file[i].hz,
								  fp_save_field_file[i].sp.slice);
					break;
			}
		}
	}
}

void writeField_xy (FILE* ex, FILE* ey, FILE* ez,
					FILE* hx, FILE* hy, FILE* hz,
					int slice)
{
	int i, j, k;

	k=slice;

	int saveis = is, saveie = ie;
	if (myrank == 0)
		saveis = nxPML_1;
	else if (myrank == nprocs - 1)
		saveie = Imax - nxPML_2 - 1;
	else;

	for(j = nyPML_1; j < Jmax - nyPML_2; j++)
	{
		for(i = saveis; i <= saveie; i++)
		{
			fprintf(ex,"%lf ", Ex(i,j,k-1));
			fprintf(ey,"%lf ", Ey(i,j,k-1));
			fprintf(ez,"%lf ", Ez(i,j,k));
			fprintf(hx,"%lf ", Hx(i,j,k));
			fprintf(hy,"%lf ", Hy(i,j,k));
			fprintf(hz,"%lf ", Hz(i,j,k-1));
		}
		fprintf(ex,"\n");
		fprintf(ey,"\n");
		fprintf(ez,"\n");
		fprintf(hx,"\n");
		fprintf(hy,"\n");
		fprintf(hz,"\n");
	}
}

void writeField_yz (FILE* ex, FILE* ey, FILE* ez,
					FILE* hx, FILE* hy, FILE* hz,
					int slice)
{
	int i, j, k;

	i=slice;

	int change_by_is = 0;
	if (myrank == 0)
		change_by_is = -_global_is;
	else
		change_by_is = -_global_is + 1;

	if(_global_is <= i &&  i<= _global_ie)
		i = i+change_by_is;
	else
		return;

	for(k = nzPML_1; k < Kmax - nzPML_2; k++)
	{
		for(j = nyPML_1; j < Jmax - nyPML_2; j++)
		{
			fprintf(ex,"%lf ", Ex(i,j,k-1));
			fprintf(ey,"%lf ", Ey(i,j,k-1));
			fprintf(ez,"%lf ", Ez(i,j,k));
			fprintf(hx,"%lf ", Hx(i,j,k));
			fprintf(hy,"%lf ", Hy(i,j,k));
			fprintf(hz,"%lf ", Hz(i,j,k-1));
		}
		fprintf(ex,"\n");
		fprintf(ey,"\n");
		fprintf(ez,"\n");
		fprintf(hx,"\n");
		fprintf(hy,"\n");
		fprintf(hz,"\n");
	}
}

void writeField_xz (FILE* ex, FILE* ey, FILE* ez,
					FILE* hx, FILE* hy, FILE* hz,
					int slice)
{
	int i, j, k;

	j=slice;

	int saveis = is, saveie = ie;
	if (myrank == 0)
		saveis = nxPML_1;
	else if (myrank == nprocs - 1)
		saveie = Imax - nxPML_2 - 1;
	else;

	for(k = nzPML_1; k < Kmax - nzPML_2; k++)
	{
		for(i = saveis; i <= saveie; i++)
		{
			fprintf(ex,"%lf ", Ex(i,j,k-1));
			fprintf(ey,"%lf ", Ey(i,j,k-1));
			fprintf(ez,"%lf ", Ez(i,j,k));
			fprintf(hx,"%lf ", Hx(i,j,k));
			fprintf(hy,"%lf ", Hy(i,j,k));
			fprintf(hz,"%lf ", Hz(i,j,k-1));
		}
		fprintf(ex,"\n");
		fprintf(ey,"\n");
		fprintf(ez,"\n");
		fprintf(hx,"\n");
		fprintf(hy,"\n");
		fprintf(hz,"\n");
	}
}

void writeField_xyOneCPU (FILE* ex, FILE* ey, FILE* ez,
					FILE* hx, FILE* hy, FILE* hz,
					int slice)
{
	int i, j, k;
	k=slice;
	for(j = nyPML_1; j < Jmax - nyPML_2; j++)
	{
		for(i = nxPML_1; i < Imax - nxPML_2; i++)
		{
			fprintf(ex,"%lf ", Ex(i,j,k-1));
			fprintf(ey,"%lf ", Ey(i,j,k-1));
			fprintf(ez,"%lf ", Ez(i,j,k));
			fprintf(hx,"%lf ", Hx(i,j,k));
			fprintf(hy,"%lf ", Hy(i,j,k));
			fprintf(hz,"%lf ", Hz(i,j,k-1));
		}
		fprintf(ex,"\n");
		fprintf(ey,"\n");
		fprintf(ez,"\n");
		fprintf(hx,"\n");
		fprintf(hy,"\n");
		fprintf(hz,"\n");
	}
}

void writeField_yzOneCPU (FILE* ex, FILE* ey, FILE* ez,
					FILE* hx, FILE* hy, FILE* hz,
					int slice)
{
	int i, j, k;
	i=slice;
	for(k = 1; k < Kmax; k++)
	{
		for(j = 0; j < Jmax-1; j++)
		{
			fprintf(ey,"%lf ", Ey(i,j,k-1));
			fprintf(hz,"%lf ", Hz(i,j,k-1));
		}
		fprintf(ey,"\n");
		fprintf(hz,"\n");
	}
	for(k = 1; k < Kmax-1; k++)
	{
		for(j = 0; j < Jmax; j++)
		{
			fprintf(ez,"%lf ", Ez(i,j,k));
			fprintf(hy,"%lf ", Hy(i,j,k));
		}
		fprintf(ez,"\n");
		fprintf(hy,"\n");
	}
	for(k = 1; k < Kmax; k++)
	{
		for(j = 0; j < Jmax; j++)
		{
			fprintf(ex,"%lf ", Ex(i,j,k-1));
		}
		fprintf(ex,"\n");
	}
	for(k = 1; k < Kmax-1; k++)
	{
		for(j = 0; j < Jmax - 1; j++)
		{
			fprintf(hx,"%lf ", Hx(i,j,k));
		}
		fprintf(hx,"\n");
	}
}

void writeField_xzOneCPU (FILE* ex, FILE* ey, FILE* ez,
					FILE* hx, FILE* hy, FILE* hz,
					int slice)
{
	int i, j, k;
	j=slice;
	for(k = 1; k < Kmax-1; k++)
	{
		for(i = 0; i < Imax; i++)
		{
			fprintf(ez,"%lf ", Ez(i,j,k));
			fprintf(hx,"%lf ", Hx(i,j,k));
		}
		fprintf(ez,"\n");
		fprintf(hx,"\n");
	}
	for(k = 1; k < Kmax; k++)
	{
		for(i = 0; i < Imax - 1; i++)
		{
			fprintf(ex,"%lf ", Ex(i,j,k-1));
			fprintf(hz,"%lf ", Hz(i,j,k-1));
		}
		fprintf(ex,"\n");
		fprintf(hz,"\n");
	}
	for(k = 1; k < Kmax; k++)
	{
		for(i = 0; i < Imax; i++)
		{
			fprintf(ey,"%lf ", Ey(i,j,k-1));
		}
		fprintf(ey,"\n");
	}
	for(k = 1; k < Kmax-1; k++)
	{
		for(i = 0; i < Imax - 1; i++)
		{
			fprintf(hy,"%lf ", Hy(i,j,k));
		}
		fprintf(hy,"\n");
	}
}

void saveData_Point(int timestep, int px, int py, int pz)
{
	int change_by_is = 0;
	if (myrank == 0)
		change_by_is = -_global_is;
	else
		change_by_is = -_global_is + 1;
	if (_global_is <= px && px <= _global_ie)
	{
			px += change_by_is;
			fprintf(fp_save_field_file[0].ex, "%lf\n", Ex(px,py,pz-1));
			fprintf(fp_save_field_file[0].ey, "%lf\n", Ey(px,py,pz-1));
			fprintf(fp_save_field_file[0].ez, "%lf\n", Ez(px,py,pz));
			fprintf(fp_save_field_file[0].hx, "%lf\n", Hx(px,py,pz));
			fprintf(fp_save_field_file[0].hy, "%lf\n", Hy(px,py,pz));
			fprintf(fp_save_field_file[0].hz, "%lf\n", Hz(px,py,pz-1));
	}
}

