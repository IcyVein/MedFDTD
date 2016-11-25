#include <stdio.h>

//float compute_whole_body_average_sar(float ***local_sar, float ***mass, const int *space_dim, FILE *fp_log)
float compute_whole_body_average_sar(float ***local_sar, unsigned char ***model, double *rho, const int *space_dim, FILE *fp_log)
{
    int x = 0, y = 0, z = 0;
    float whole_body_average_sar = 0, sum_mass = 0;
    int sum_whole_body_average_sar_voxel = 0;

    printf("-------- Compute whole body average SAR --------\n");
    fflush(stdout);/* %% */
    fprintf(fp_log, "-------- Compute whole body average SAR --------\n");

    for (z = 0; z < space_dim[2]; z++)
        for (y = 0; y < space_dim[1]; y++)
            for (x = 0; x < space_dim[0]; x++)
                if (local_sar[x][y][z] != 0 && /*mass[x][y][z]*/ rho[model[x][y][z]] != 0)
                {
                    whole_body_average_sar += local_sar[x][y][z]* /*mass[x][y][z]*/ rho[model[x][y][z]];
                    sum_mass += /*mass[x][y][z]*/ rho[model[x][y][z]];
                    sum_whole_body_average_sar_voxel += 1;
                }
	if (sum_mass != 0)
	{
    	whole_body_average_sar /= sum_mass;

    	printf("Total model voxel = %d\n", sum_whole_body_average_sar_voxel);
    	printf("Whole body average SAR = %lf W/kg\n", whole_body_average_sar);
    	fflush(stdout);/* %% */
    	fprintf(fp_log, "Total model voxel = %d\n", sum_whole_body_average_sar_voxel);
    	fprintf(fp_log, "Whole body average SAR = %lf W/kg\n", whole_body_average_sar);
		return whole_body_average_sar;
	}
	else
	{
		whole_body_average_sar = 0;
		printf("None body in computational space.\n");
    	printf("Whole body average SAR = %lf W/kg\n", whole_body_average_sar);
    	fflush(stdout);/* %% */
    	fprintf(fp_log, "None body in computational space.\n");
    	fprintf(fp_log, "Whole body average SAR = %lf W/kg\n", whole_body_average_sar);
	}
}

