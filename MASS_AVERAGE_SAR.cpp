#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "MASS_AVERAGE_SAR.H"

int compute_mass_average_sar(float ***local_sar, float ***mass, const int *space_dim, const float required_mass, const float air_voxel_rate_threshold, const char *output_path, FILE *fp_log)
{
    /**************************************************************
     * Compute mass average SAR according to IEEE Std C95.3-2002. *
     * The averaging volume is in the shape of a cube.            *
     **************************************************************/

    FILE *fp_sar;
    char mass_average_sar_file[MAX_SIZE_OF_PATH];/* %% */
    float ***mass_average_sar;
    int ***semi_side_length;
    int x = 0, y = 0, z = 0, tt = 0;
    int sum_mass_average_sar_voxel = 0;
    float max_mass_average_sar = 0;
    time_t timep;

    time(&timep);
    printf("------------------------\n");
    printf("%s", ctime(&timep));
    printf("------------------------\n");
    fflush(stdout); /* %% */
    fprintf(fp_log, "------------------------\n");
    fprintf(fp_log, "%s", ctime(&timep));
    fprintf(fp_log, "------------------------\n");

    printf("Required mass = %lf kg\n", required_mass);
    fflush(stdout); /* %% */
    fprintf(fp_log, "Required mass = %lf kg\n", required_mass);

    printf("-------- Initialize --------\n");
    fflush(stdout); /* %% */
    fprintf(fp_log, "-------- Initialize --------\n");
    fflush(fp_log); /* %% */

    sprintf(mass_average_sar_file, "%sSAR%.0fg.txt", output_path, 1000*required_mass);
    printf("<- %s\n", mass_average_sar_file);
    fflush(stdout); /* %% */
    fprintf(fp_log, "<- %s\n", mass_average_sar_file);
    if ((fp_sar = fopen(mass_average_sar_file, "w")) == NULL)
    {
        fprintf(stderr, "ERROR in openning %s.\n", mass_average_sar_file);
        fprintf(fp_log, "ERROR in openning %s.\n", mass_average_sar_file);
        exit(1);
    }

    mass_average_sar = (float ***)malloc(space_dim[0]*sizeof(float **));
    semi_side_length = (int ***)malloc(space_dim[0]*sizeof(int **));
    for(x = 0; x < space_dim[0]; x++)
    {
        mass_average_sar[x] = (float **)malloc(space_dim[1]*sizeof(float *));
        semi_side_length[x] = (int **)malloc(space_dim[1]*sizeof(int *));
        for(y = 0; y < space_dim[1]; y++)
        {
            mass_average_sar[x][y] = (float *)malloc(space_dim[2]*sizeof(float));
            semi_side_length[x][y] = (int *)malloc(space_dim[2]*sizeof(int));
            for(z = 0; z < space_dim[2]; z++)
            {
                mass_average_sar[x][y][z] = 0;
                semi_side_length[x][y][z] = 0;
            }
        }
    }

    printf("-------- Compute SAR of type 1 voxels --------\n");
    fflush(stdout); /* %% */
    fprintf(fp_log, "-------- Compute SAR of type 1 voxels --------\n");
    compute_type_1_voxel(local_sar, mass, semi_side_length, mass_average_sar, required_mass, air_voxel_rate_threshold, space_dim, fp_log);

    printf("-------- Compute SAR of type 2 voxels --------\n");
    fflush(stdout); /* %% */
    fprintf(fp_log, "-------- Compute SAR of type 2 voxels --------\n");
    compute_type_2_voxel(mass, mass_average_sar, semi_side_length, space_dim, fp_log);

    printf("-------- Compute SAR of type 3 voxels --------\n");
    fflush(stdout); /* %% */
    fprintf(fp_log, "-------- Compute SAR of type 3 voxels --------\n");
    compute_type_3_voxel(local_sar, mass, mass_average_sar, required_mass, space_dim, fp_log);

    printf("-------- Save data --------\n");
    fflush(stdout); /* %% */
    fprintf(fp_log, "-------- Save data --------\n");
    sum_mass_average_sar_voxel = 0;
    for (z = 0; z < space_dim[2]; z++)
    {
        printf("Process 1 / 3: |");
        fflush(stdout); /* %% */
        for (tt = 0; tt < (z+1)*20/space_dim[2]; tt++)
        {
            printf("=");
            fflush(stdout); /* %% */
        }
        for (tt = (z+1)*20/space_dim[2]; tt < 20; tt++)
        {
            printf(" ");
            fflush(stdout); /* %% */
        }
        printf("| %3.1f %%\r", (float)(z+1)/space_dim[2]*100);
        fflush(stdout); /* %% */

        for (y = 0; y < space_dim[1]; y++)
        {
            for (x = 0; x < space_dim[0]; x++)
            {
                fprintf(fp_sar, "%e ", mass_average_sar[x][y][z]);
                if (mass_average_sar[x][y][z] != 0)
                    sum_mass_average_sar_voxel += 1;
                if (mass_average_sar[x][y][z] > max_mass_average_sar)
                    max_mass_average_sar = mass_average_sar[x][y][z];
            }
            fprintf(fp_sar, "\n");
        }
    }

    printf("Total mass average SAR voxel = %d\n", sum_mass_average_sar_voxel);
    printf("Peak mass average SAR = %lf W/kg\n", max_mass_average_sar);
    fflush(stdout); /* %% */
    fprintf(fp_log, "Total mass average SAR voxel = %d\n", sum_mass_average_sar_voxel);
    fprintf(fp_log, "Peak mass average SAR = %lf W/kg\n", max_mass_average_sar);

    for (x = 0; x < space_dim[0]; x++)
        for (y = 0; y < space_dim[1]; y++)
        {
            free(mass_average_sar[x][y]);
            free(semi_side_length[x][y]);
        }
    for (x = 0; x < space_dim[0]; x++)
    {
        free(mass_average_sar[x]);
        free(semi_side_length[x]);
    }
    free(mass_average_sar);
    free(semi_side_length);

    if (fclose(fp_sar) != 0)
    {
        fprintf(stderr, "ERROR in closing %s.\n", mass_average_sar_file);
        fprintf(fp_log, "ERROR in closing %s.\n", mass_average_sar_file);
        exit(1);
    }

    time(&timep);
    printf("------------------------\n");
    printf("%s", ctime(&timep));
    printf("------------------------\n");
    printf("Compute %lf g average SAR. DONE!\n", 1000*required_mass);
    fflush(stdout); /* %% */

    fprintf(fp_log, "------------------------\n");
    fprintf(fp_log, "%s", ctime(&timep));
    fprintf(fp_log, "------------------------\n");
    fprintf(fp_log, "Compute %lf g average SAR. DONE!\n", required_mass);
    fflush(fp_log); /* %% */

    return 0;
}

int compute_type_1_voxel(float ***local_sar, float ***mass, int ***semi_side_length, float ***mass_average_sar, const float required_mass, const float air_voxel_rate_threshold, const int *space_dim, FILE *fp_log)
{
    /*******************************************************************
     * The "type 1" voxels are the voxels that locate in the center of *
     * "valid averaging volume" (IEEE Std C95.3-2002).                 *
     *******************************************************************/

    int x = 0, y = 0, z = 0, xx = 0, yy = 0, zz = 0, tt = 0;
    int semi_side_len = 0;
    float air_voxel_rate = 0, max_air_voxel_rate = 0;
    int part_mass_id = 0;
    float part_mass[4] = {0, 0, 0, 0};
    float fraction = 0, min_fraction = 1, max_fraction = 0;
    float cube_mass = 0;
    float current_cube_mass = 0, min_current_cube_mass = 10000, max_current_cube_mass = 0;
    float max_mass_average_sar = 0;
    int sum_type_1_voxel = 0;
    int max_mass_average_sar_x = 0, max_mass_average_sar_y = 0, max_mass_average_sar_z = 0;

    for (z = 0; z < space_dim[2]; z++)
    {
        printf("Process 1 / 3: |");
        fflush(stdout); /* %% */
        for (tt = 0; tt < (z+1)*20/space_dim[2]; tt++)
        {
            printf("=");
            fflush(stdout); /* %% */
        }
        for (tt = (z+1)*20/space_dim[2]; tt < 20; tt++)
        {
            printf(" ");
            fflush(stdout); /* %% */
        }
        printf("| %3.1f %%\r", (float)(z+1)/space_dim[2]*100);
        fflush(stdout); /* %% */

        for (y = 0; y < space_dim[1]; y++)
        {
            for (x = 0; x < space_dim[0]; x++)
            {
                if (mass[x][y][z] > 0) /* It's tissue. */
                {
                    semi_side_len = 1;
                    cube_mass = 0;
                    current_cube_mass = 0;
                    for (part_mass_id = 0; part_mass_id < 4; part_mass_id++)
                        part_mass[part_mass_id] = 0;
                    while (1)
                    {
                        if (find_empty_side(mass, semi_side_len, x, y, z, space_dim) == 0)
                        {
                            air_voxel_rate = compute_air_voxel_rate(mass, semi_side_len, x, y, z);
                            if (air_voxel_rate <= air_voxel_rate_threshold)
                            {
                                cube_mass = 0;
                                for (xx = x-semi_side_len; xx <= x+semi_side_len; xx++)
                                    for (yy = y-semi_side_len; yy <= y+semi_side_len; yy++)
                                        for (zz = z-semi_side_len; zz <= z+semi_side_len; zz++)
                                                cube_mass += mass[xx][yy][zz];
                                if (cube_mass < required_mass)
                                    semi_side_len += 1;
                                else
                                {
                                    compute_type_1_part_mass(mass, part_mass, semi_side_len, x, y, z);
                                    fraction = compute_fraction(part_mass, required_mass);

                                    current_cube_mass = part_mass[3]*fraction*fraction*fraction+part_mass[2]*fraction*fraction+part_mass[1]*fraction+part_mass[0];
                                    mass_average_sar[x][y][z] = compute_type_1_sar(local_sar, mass, semi_side_len, fraction, x, y, z)/current_cube_mass;

                                    semi_side_length[x][y][z] = semi_side_len;
                                    sum_type_1_voxel += 1;

                                    if (air_voxel_rate > max_air_voxel_rate)
                                        max_air_voxel_rate = air_voxel_rate;
                                    if (fraction < min_fraction)
                                        min_fraction = fraction;
                                    if (fraction > max_fraction)
                                        max_fraction = fraction;
                                    if (current_cube_mass < min_current_cube_mass)
                                        min_current_cube_mass = current_cube_mass;
                                    if (current_cube_mass > max_current_cube_mass)
                                        max_current_cube_mass = current_cube_mass;
                                    if (mass_average_sar[x][y][z] > max_mass_average_sar)
                                    {
                                        max_mass_average_sar = mass_average_sar[x][y][z];
                                        max_mass_average_sar_x = x;
                                        max_mass_average_sar_y = y;
                                        max_mass_average_sar_z = z;
                                    }
                                    break;
                                }
                            } /* if air voxel rate <= threshold */
                            else /* air voxel rate > threshold */
                                break;
                        } /* if there is not empty side */
                        else /* if there is empty side */
                            break;
                    } /* while */
                } /* if tissue */
            } /* for x */
        } /* for y */
    } /* for z */
    printf("\nMaximum air voxel rate = %lf\n", max_air_voxel_rate);
    printf("Fraction = [%lf, %lf]\n", min_fraction, max_fraction);
    printf("Averaging mass = [%lf, %lf] kg\n", min_current_cube_mass, max_current_cube_mass);
    printf("Total type 1 voxel = %d\n", sum_type_1_voxel);
    printf("Peak mass average SAR of type 1 voxel = %lf W/kg in (%d, %d, %d)\n", max_mass_average_sar, max_mass_average_sar_x, max_mass_average_sar_y, max_mass_average_sar_z);
    fflush(stdout); /* %% */

    fprintf(fp_log, "Maximum air voxel rate = %lf\n", max_air_voxel_rate);
    fprintf(fp_log, "Fraction = [%lf, %lf]\n", min_fraction, max_fraction);
    fprintf(fp_log, "Averaging mass = [%lf, %lf] kg\n", min_current_cube_mass, max_current_cube_mass);
    fprintf(fp_log, "Total type 1 voxel = %d\n", sum_type_1_voxel);
    fprintf(fp_log, "Peak mass average SAR of type 1 voxel = %lf W/kg in (%d, %d, %d)\n", max_mass_average_sar, max_mass_average_sar_x, max_mass_average_sar_y, max_mass_average_sar_z);
    fflush(fp_log); /* %% */

    return 0;
}

int find_empty_side(float ***mass, const int semi_side_len, const int px, const int py, const int pz, const int *space_dim)
{
    /*****************************************************************
     * return 0: There is not empty surface in the averaging volume. *
     * return 1: There is empty surface in the averaging volume.     *
     *****************************************************************/

    float cube_mass = 0;
    int x = 0, y = 0, z = 0;

    if (px-semi_side_len < 0 || px+semi_side_len >= space_dim[0] || py-semi_side_len < 0 || py+semi_side_len >= space_dim[1] || pz-semi_side_len < 0 || pz+semi_side_len >= space_dim[2])
        return 1;

    cube_mass = 0;
    for (z = pz-semi_side_len; z <= pz+semi_side_len; z++)
        for (y = py-semi_side_len; y <= py+semi_side_len; y++)
            cube_mass += mass[px-semi_side_len][y][z];
    if (cube_mass == 0)
        return 1;

    cube_mass = 0;
    for (z = pz-semi_side_len; z <= pz+semi_side_len; z++)
        for (y = py-semi_side_len; y <= py+semi_side_len; y++)
            cube_mass += mass[px+semi_side_len][y][z];
    if (cube_mass == 0)
        return 1;

    cube_mass = 0;
    for (x = px-semi_side_len; x <= px+semi_side_len; x++)
        for (z = pz-semi_side_len; z <= pz+semi_side_len; z++)
            cube_mass += mass[x][py-semi_side_len][z];
    if (cube_mass == 0)
        return 1;

    cube_mass = 0;
    for (x = px-semi_side_len; x <= px+semi_side_len; x++)
        for (z = pz-semi_side_len; z <= pz+semi_side_len; z++)
            cube_mass += mass[x][py+semi_side_len][z];
    if (cube_mass == 0)
        return 1;

    cube_mass = 0;
    for (y = py-semi_side_len; y <= py+semi_side_len; y++)
        for (x = px-semi_side_len; x <= px+semi_side_len; x++)
            cube_mass += mass[x][y][pz-semi_side_len];
    if (cube_mass == 0)
        return 1;

    cube_mass = 0;
    for (y = py-semi_side_len; y <= py+semi_side_len; y++)
        for (x = px-semi_side_len; x <= px+semi_side_len; x++)
            cube_mass += mass[x][y][pz+semi_side_len];
    if (cube_mass == 0)
        return 1;

    return 0;
}

float compute_air_voxel_rate(float ***mass, const int semi_side_len, const int px, const int py, const int pz)
{
    /**********************************************************
     * Compute the rate of air voxel in the averaging volume. *
     **********************************************************/

    int x = 0, y = 0, z = 0, sum_air_voxel = 0;

    for (z = pz-semi_side_len; z <= pz+semi_side_len; z++)
        for (y = py-semi_side_len; y <= py+semi_side_len; y++)
            for (x = px -semi_side_len; x <= px+semi_side_len; x++)
                if (mass[x][y][z] == 0)
                    sum_air_voxel += 1;

    return (float)sum_air_voxel/((2*semi_side_len+1)*(2*semi_side_len+1)*(2*semi_side_len+1));
}

int compute_type_1_part_mass(float ***mass, float *part_mass, const int semi_side_len, const int px, const int py, const int pz)
{
    /*******************************************
     * Compute mass of four parts of the cube. *
     *******************************************/

    int x = 0, y = 0, z = 0;

    /* core mass */
    for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
        for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                part_mass[0] += mass[x][y][z];
    /* 6 sides */
    for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
        for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
            part_mass[1] += mass[x][y][pz-semi_side_len]+mass[x][y][pz+semi_side_len];
    for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
        for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
            part_mass[1] += mass[px-semi_side_len][y][z]+mass[px+semi_side_len][y][z];
    for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
        for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
            part_mass[1] += mass[x][py-semi_side_len][z]+mass[x][py+semi_side_len][z];
    /* 12 edges */
    for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
    {
        part_mass[2] += mass[x][py-semi_side_len][pz-semi_side_len];
        part_mass[2] += mass[x][py+semi_side_len][pz-semi_side_len];
        part_mass[2] += mass[x][py-semi_side_len][pz+semi_side_len];
        part_mass[2] += mass[x][py+semi_side_len][pz+semi_side_len];
    }
    for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
    {
        part_mass[2] += mass[px-semi_side_len][y][pz-semi_side_len];
        part_mass[2] += mass[px+semi_side_len][y][pz-semi_side_len];
        part_mass[2] += mass[px-semi_side_len][y][pz+semi_side_len];
        part_mass[2] += mass[px+semi_side_len][y][pz+semi_side_len];
    }
    for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
    {
        part_mass[2] += mass[px-semi_side_len][py-semi_side_len][z];
        part_mass[2] += mass[px+semi_side_len][py-semi_side_len][z];
        part_mass[2] += mass[px-semi_side_len][py+semi_side_len][z];
        part_mass[2] += mass[px+semi_side_len][py+semi_side_len][z];
    }
    /* 8 corner voxels */
    part_mass[3] += mass[px-semi_side_len][py-semi_side_len][pz-semi_side_len];
    part_mass[3] += mass[px+semi_side_len][py-semi_side_len][pz-semi_side_len];
    part_mass[3] += mass[px-semi_side_len][py+semi_side_len][pz-semi_side_len];
    part_mass[3] += mass[px+semi_side_len][py+semi_side_len][pz-semi_side_len];
    part_mass[3] += mass[px-semi_side_len][py-semi_side_len][pz+semi_side_len];
    part_mass[3] += mass[px+semi_side_len][py-semi_side_len][pz+semi_side_len];
    part_mass[3] += mass[px-semi_side_len][py+semi_side_len][pz+semi_side_len];
    part_mass[3] += mass[px+semi_side_len][py+semi_side_len][pz+semi_side_len];

    return 0;
}

float compute_type_1_sar(float ***local_sar, float ***mass, const int semi_side_len, const float fraction, const int px, const int py, const int pz)
{
    /**************************************************
     * Compute mass average SAR of the type 1 voxels. *
     **************************************************/

    float core_energy = 0, corner_energy = 0, edge_energy = 0, side_energy = 0;
    int x = 0, y = 0, z = 0;

    corner_energy += local_sar[px-semi_side_len][py-semi_side_len][pz-semi_side_len]*mass[px-semi_side_len][py-semi_side_len][pz-semi_side_len];
    corner_energy += local_sar[px+semi_side_len][py-semi_side_len][pz-semi_side_len]*mass[px+semi_side_len][py-semi_side_len][pz-semi_side_len];
    corner_energy += local_sar[px-semi_side_len][py+semi_side_len][pz-semi_side_len]*mass[px-semi_side_len][py+semi_side_len][pz-semi_side_len];
    corner_energy += local_sar[px+semi_side_len][py+semi_side_len][pz-semi_side_len]*mass[px+semi_side_len][py+semi_side_len][pz-semi_side_len];
    corner_energy += local_sar[px-semi_side_len][py-semi_side_len][pz+semi_side_len]*mass[px-semi_side_len][py-semi_side_len][pz+semi_side_len];
    corner_energy += local_sar[px+semi_side_len][py-semi_side_len][pz+semi_side_len]*mass[px+semi_side_len][py-semi_side_len][pz+semi_side_len];
    corner_energy += local_sar[px-semi_side_len][py+semi_side_len][pz+semi_side_len]*mass[px-semi_side_len][py+semi_side_len][pz+semi_side_len];
    corner_energy += local_sar[px+semi_side_len][py+semi_side_len][pz+semi_side_len]*mass[px+semi_side_len][py+semi_side_len][pz+semi_side_len];

    for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
    {
        edge_energy += local_sar[x][py-semi_side_len][pz-semi_side_len]*mass[x][py-semi_side_len][pz-semi_side_len];
        edge_energy += local_sar[x][py+semi_side_len][pz-semi_side_len]*mass[x][py+semi_side_len][pz-semi_side_len];
        edge_energy += local_sar[x][py-semi_side_len][pz+semi_side_len]*mass[x][py-semi_side_len][pz+semi_side_len];
        edge_energy += local_sar[x][py+semi_side_len][pz+semi_side_len]*mass[x][py+semi_side_len][pz+semi_side_len];
    }
    for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
    {
        edge_energy += local_sar[px-semi_side_len][y][pz-semi_side_len]*mass[px-semi_side_len][y][pz-semi_side_len];
        edge_energy += local_sar[px+semi_side_len][y][pz-semi_side_len]*mass[px+semi_side_len][y][pz-semi_side_len];
        edge_energy += local_sar[px-semi_side_len][y][pz+semi_side_len]*mass[px-semi_side_len][y][pz+semi_side_len];
        edge_energy += local_sar[px+semi_side_len][y][pz+semi_side_len]*mass[px+semi_side_len][y][pz+semi_side_len];
    }
    for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
    {
        edge_energy += local_sar[px-semi_side_len][py-semi_side_len][z]*mass[px-semi_side_len][py-semi_side_len][z];
        edge_energy += local_sar[px+semi_side_len][py-semi_side_len][z]*mass[px+semi_side_len][py-semi_side_len][z];
        edge_energy += local_sar[px-semi_side_len][py+semi_side_len][z]*mass[px-semi_side_len][py+semi_side_len][z];
        edge_energy += local_sar[px+semi_side_len][py+semi_side_len][z]*mass[px+semi_side_len][py+semi_side_len][z];
    }

    for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
        for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
            side_energy += local_sar[x][y][pz-semi_side_len]*mass[x][y][pz-semi_side_len]+local_sar[x][y][pz+semi_side_len]*mass[x][y][pz+semi_side_len];
    for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
        for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
            side_energy += local_sar[px-semi_side_len][y][z]*mass[px-semi_side_len][y][z]+local_sar[px+semi_side_len][y][z]*mass[px+semi_side_len][y][z];
    for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
        for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
            side_energy += local_sar[x][py-semi_side_len][z]*mass[x][py-semi_side_len][z]+local_sar[x][py+semi_side_len][z]*mass[x][py+semi_side_len][z];

    for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
        for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                core_energy += local_sar[x][y][z]*mass[x][y][z];

    return corner_energy*fraction*fraction*fraction+edge_energy*fraction*fraction+side_energy*fraction+core_energy;

}

int compute_type_2_voxel(float ***mass, float ***mass_average_sar, int ***semi_side_length, const int *space_dim, FILE *fp_log)
{
    /****************************************************************
     * The "type 2" voxel is the voxel locates in the center of the *
     * "invalid averaging volume" (IEEE Std C95.3-2002)             *
     ****************************************************************/

    int x = 0, y = 0, z = 0, xx = 0, yy = 0, zz = 0, tt = 0;
    int max_semi_side_length = 0;
    int sum_type_2_voxel = 0;

    for (z = 0; z < space_dim[2]; z++)
        for (y = 0; y < space_dim[1]; y++)
            for (x = 0; x < space_dim[0]; x++)
                if (semi_side_length[x][y][z] > max_semi_side_length)
                    max_semi_side_length = semi_side_length[x][y][z];

    for (z = 0; z < space_dim[2]; z++)
    {
        printf("Process 2 / 3: |");
        fflush(stdout); /* %% */
        for (tt = 0; tt < (z+1)*20/space_dim[2]; tt++)
        {
            printf("=");
            fflush(stdout); /* %% */
        }
        for (tt = (z+1)*20/space_dim[2]; tt < 20; tt++)
        {
            printf(" ");
            fflush(stdout); /* %% */
        }
        printf("| %3.1f %%\r", (float)(z+1)/space_dim[2]*100);
        fflush(stdout); /* %% */

        for (y = 0; y < space_dim[1]; y++)
            for (x = 0; x < space_dim[0]; x++)
                if (mass[x][y][z] != 0 && semi_side_length[x][y][z] == 0)
                {
                    for (zz = z-max_semi_side_length; zz <= z+max_semi_side_length; zz++)
                        for (yy = y-max_semi_side_length; yy <= y+max_semi_side_length; yy++)
                            for (xx = x-max_semi_side_length; xx <= x+max_semi_side_length; xx++)
                                if (xx >= 0 && xx < space_dim[0] && yy >= 0 && yy < space_dim[1] && zz >= 0 && zz < space_dim[2] && semi_side_length[xx][yy][zz] != 0)
                                    if (x >= xx-semi_side_length[xx][yy][zz] && x <= xx+semi_side_length[xx][yy][zz] && y >= yy-semi_side_length[xx][yy][zz] && y <= yy+semi_side_length[xx][yy][zz] && z >= zz-semi_side_length[xx][yy][zz] && z <= zz+semi_side_length[xx][yy][zz])
                                        if (mass_average_sar[xx][yy][zz] > mass_average_sar[x][y][z])
                                            mass_average_sar[x][y][z] = mass_average_sar[xx][yy][zz];
                    if (mass_average_sar[x][y][z] != 0)
                        sum_type_2_voxel += 1;
                }
    }

    printf("\nTotal type 2 voxel = %d\n", sum_type_2_voxel);
    fflush(stdout);/* %% */
    fprintf(fp_log, "Total type 2 voxel = %d\n", sum_type_2_voxel);
    fflush(fp_log); /* %% */

    return 0;
}

int compute_type_3_voxel(float ***local_sar, float ***mass, float ***mass_average_sar, const float required_mass, const int *space_dim, FILE *fp_log)
{
    /****************************************************************************
     * Voxels which are not included in the types 1 and 2 are treated as type 3 *
     * (IEEE Std C95.3-2002).                                                   *
     ****************************************************************************/

    int x = 0, y = 0, z = 0, tt = 0, direction_id = 0;
    int layer_marker = 0, semi_side_len = 0, direction = 0;
    int sum_type_3_voxel = 0, sum_type_4_voxel = 0;
    float cube_mass[6] = {0, 0, 0, 0, 0, 0};
    float inner_layer_part_mass[4] = {0, 0, 0, 0}, outer_layer_part_mass[4] = {0, 0, 0, 0};
    float max_mass = 0, last_max_mass = 0, fraction = 0;
    float current_cube_mass = 0, min_current_cube_mass = 10000, max_current_cube_mass = 0;
    float max_type_3_mass_average_sar = 0;

    for (z = 0; z < space_dim[2]; z++)
    {
        printf("Process 3 / 3: |");
        fflush(stdout); /* %% */
        for (tt = 0; tt < (z+1)*20/space_dim[2]; tt++)
        {
            printf("=");
            fflush(stdout); /* %% */
        }
        for (tt = (z+1)*20/space_dim[2]; tt < 20; tt++)
        {
            printf(" ");
            fflush(stdout); /* %% */
        }
        printf("| %3.1f %%\r", (float)(z+1)/space_dim[2]*100);
        fflush(stdout); /* %% */

        for (y = 0; y < space_dim[1]; y++)
        {
            for (x = 0; x < space_dim[0]; x++)
            {
                if (mass[x][y][z] != 0 && mass_average_sar[x][y][z] == 0) /* it's type 3 voxel */
                {
                    semi_side_len = 0;
                    max_mass = 0;
                    while (1)
                    {
                        last_max_mass = max_mass;
                        if (max_mass == last_max_mass)
                        {
                            /**********************************************************
                             * The voxel would be treated as type 4 voxel if the mass *
                             * can't reach required mass and doesn't increase while   *
                             * the averaging volume expands.                          *
                             **********************************************************/
                            sum_type_4_voxel += 1;
                            break;
                        }

                        for (direction_id = 0; direction_id < 6; direction_id++)
                            cube_mass[direction_id] = 0;
                        max_mass = compute_type_3_cube_mass(mass, cube_mass, semi_side_len, x, y, z, space_dim);
                        if (max_mass < required_mass)
                            semi_side_len += 1;
                        else
                        {
                            for (direction_id = 0; direction_id < 6; direction_id++)
                                if (cube_mass[direction_id] == max_mass)
                                    direction = direction_id;
                            compute_type_3_part_mass(inner_layer_part_mass, outer_layer_part_mass, mass, semi_side_len, direction, x, y, z);
                            fraction = compute_type_3_fraction(inner_layer_part_mass, outer_layer_part_mass, required_mass);
                            if (fraction > 1) /* outer layer is divided */
                            {
                                fraction -= 1;
                                layer_marker = 1;
                                current_cube_mass = outer_layer_part_mass[0]
                                                   +outer_layer_part_mass[1]*fraction
                                                   +outer_layer_part_mass[2]*fraction*fraction
                                                   +outer_layer_part_mass[3]*fraction*fraction*fraction;
                            }
                            else /* inner layer is divided */
                            {
                                layer_marker = 0;
                                current_cube_mass = inner_layer_part_mass[0]
                                                   +inner_layer_part_mass[1]*fraction
                                                   +inner_layer_part_mass[2]*fraction*fraction
                                                   +inner_layer_part_mass[3]*fraction*fraction*fraction;
                            }
                            mass_average_sar[x][y][z] = compute_type_3_sar(local_sar, mass, direction, semi_side_len, fraction, layer_marker, max_mass, x, y, z)/current_cube_mass;
                            if (current_cube_mass < min_current_cube_mass)
                                min_current_cube_mass = current_cube_mass;
                            if (current_cube_mass > max_current_cube_mass)
                                max_current_cube_mass = current_cube_mass;
                            if (mass_average_sar[x][y][z] > max_type_3_mass_average_sar)
                                max_type_3_mass_average_sar = mass_average_sar[x][y][z];
                            sum_type_3_voxel += 1;
                            break;
                        }
                    }
                    /*
                    while (max_mass < required_mass)
                    {
                        last_max_mass = max_mass;
                        semi_side_len += 1;
                        for (direction_id = 0; direction_id < 6; direction_id++)
                            cube_mass[direction_id] = 0;
                        max_mass = compute_type_3_cube_mass(mass, cube_mass, semi_side_len, x, y, z, space_dim);
                        if (max_mass == last_max_mass)
                        {
                            sum_type_4_voxel += 1;
                            break;
                        }
                    }
                    if (max_mass >= required_mass)
                    {
                        for (direction_id = 0; direction_id < 6; direction_id++)
                            if (cube_mass[direction_id] == max_mass)
                                direction = direction_id;
                        compute_type_3_part_mass(inner_layer_part_mass, outer_layer_part_mass, mass, semi_side_len, direction, x, y, z);
                        fraction = compute_type_3_fraction(inner_layer_part_mass, outer_layer_part_mass, required_mass);
                        if (fraction > 1)
                        {
                            fraction -= 1;
                            layer_marker = 1;
                            current_cube_mass = outer_layer_part_mass[0]
                                               +outer_layer_part_mass[1]*fraction
                                               +outer_layer_part_mass[2]*fraction*fraction
                                               +outer_layer_part_mass[3]*fraction*fraction*fraction;
                        }
                        else
                        {
                            layer_marker = 0;
                            current_cube_mass = inner_layer_part_mass[0]
                                               +inner_layer_part_mass[1]*fraction
                                               +inner_layer_part_mass[2]*fraction*fraction
                                               +inner_layer_part_mass[3]*fraction*fraction*fraction;
                        }
                        mass_average_sar[x][y][z] = compute_type_3_sar(local_sar, mass, direction, semi_side_len, fraction, layer_marker, max_mass, x, y, z)/current_cube_mass;
                        if (current_cube_mass < min_current_cube_mass)
                            min_current_cube_mass = current_cube_mass;
                        if (current_cube_mass > max_current_cube_mass)
                            max_current_cube_mass = current_cube_mass;
                        if (mass_average_sar[x][y][z] > max_type_3_mass_average_sar)
                            max_type_3_mass_average_sar = mass_average_sar[x][y][z];
                        sum_type_3_voxel += 1;
                    }
                    */
                } /* it's type 3 voxel */
            } /* for x */
        } /* for y */
    } /* for z */

    printf("\n");
    printf("Averaging mass = [%lf, %lf] kg\n", min_current_cube_mass, max_current_cube_mass);
    printf("Total type 3 voxel = %d\n", sum_type_3_voxel);
    printf("Peak mass average SAR of type 3 voxel = %lf W/kg\n", max_type_3_mass_average_sar);
    fflush(stdout);/* %% */

    fprintf(fp_log, "Averaging mass = [%lf, %lf] kg\n", min_current_cube_mass, max_current_cube_mass);
    fprintf(fp_log, "Total type 3 voxel = %d\n", sum_type_3_voxel);
    fprintf(fp_log, "Peak mass average SAR of type 3 voxel = %lf W/kg\n", max_type_3_mass_average_sar);
    fflush(fp_log); /* %% */

    return 0;
}

float compute_type_3_cube_mass(float ***mass, float *cube_mass, const int semi_side_len, const int px, const int py, const int pz, const int *space_dim)
{
    /************************************************************
     * Enlarge cube in six directions.                          *
     * If the cube extends out of the computational region, the *
     * mass of the cube would be assigned 0.                    *
     ************************************************************/

    int x = 0, y = 0, z = 0, direction_id = 0;
    float max_mass = 0;

    if (px+2*semi_side_len < space_dim[0] && py-semi_side_len >= 0 && py+semi_side_len < space_dim[1] && pz-semi_side_len >= 0 && pz+semi_side_len < space_dim[2])
        for (y = py-semi_side_len; y <= py+semi_side_len; y++) /* x+ */
            for (z = pz-semi_side_len; z <= pz+semi_side_len; z++)
                for (x = px; x <= px+2*semi_side_len; x++)
                    cube_mass[0] += mass[x][y][z];
    else
        cube_mass[0] = 0;

    if (px-2*semi_side_len >= 0 && py-semi_side_len >= 0 && py+semi_side_len < space_dim[1] && pz-semi_side_len >= 0 && pz+semi_side_len < space_dim[2])
        for (y = py-semi_side_len; y <= py+semi_side_len; y++) /* x- */
            for (z = pz-semi_side_len; z <= pz+semi_side_len; z++)
                for (x = px-2*semi_side_len; x <= px; x++)
                    cube_mass[1] += mass[x][y][z];
    else
        cube_mass[1] = 0;

    if (py+2*semi_side_len < space_dim[1] && pz-semi_side_len >= 0 && pz+semi_side_len < space_dim[2] && px-semi_side_len >= 0 && px+semi_side_len < space_dim[0])
        for (z = pz-semi_side_len; z <= pz+semi_side_len; z++) /* y+ */
            for (x = px-semi_side_len; x <= px+semi_side_len; x++)
                for (y = py; y <= py+2*semi_side_len; y++)
                    cube_mass[2] += mass[x][y][z];
    else
        cube_mass[2] = 0;

    if (py-2*semi_side_len >= 0 && pz-semi_side_len >= 0 && pz+semi_side_len < space_dim[2] && px-semi_side_len >= 0 && px+semi_side_len < space_dim[0])
        for (z = pz-semi_side_len; z <= pz+semi_side_len; z++) /* y- */
            for (x = px-semi_side_len; x <= px+semi_side_len; x++)
                for (y = py-2*semi_side_len; y <= py; y++)
                    cube_mass[3] += mass[x][y][z];
    else
        cube_mass[3] = 0;

    if (pz+2*semi_side_len < space_dim[2] && px-semi_side_len >= 0 && px+semi_side_len < space_dim[0] && py-semi_side_len >= 0 && py+semi_side_len < space_dim[1])
        for (x = px-semi_side_len; x <= px+semi_side_len; x++) /* z+ */
            for (y = py-semi_side_len; y <= py+semi_side_len; y++)
                for (z = pz; z <= pz+2*semi_side_len; z++)
                    cube_mass[4] += mass[x][y][z];
    else
        cube_mass[4] = 0;

    if (pz-2*semi_side_len >= 0 && px-semi_side_len >= 0 && px+semi_side_len < space_dim[0] && py-semi_side_len >= 0 && py+semi_side_len < space_dim[1])
        for (x = px-semi_side_len; x <= px+semi_side_len; x++) /* z- */
            for (y = py-semi_side_len; y <= py+semi_side_len; y++)
                for (z = pz-2*semi_side_len; z <= pz; z++)
                    cube_mass[5] += mass[x][y][z];
    else
        cube_mass[5] = 0;

    for (direction_id = 0; direction_id < 6; direction_id++)
        if (cube_mass[direction_id] > max_mass)
            max_mass = cube_mass[direction_id];

    return max_mass;
}

int compute_type_3_part_mass(float *inner_layer_part_mass, float *outer_layer_part_mass, float ***mass, const int semi_side_len, const int direction, const int px, const int py, const int pz)
{
    /************************************************************
     * Compute mass of different parts in the averaging volume. *
     ************************************************************/

    int x = 0, y = 0, z = 0;
    float core_mass = 0;
    float corner_1_mass = 0, corner_2_mass = 0;
    float edge_1_mass = 0, edge_2_mass = 0, edge_3_mass = 0;
    float side_1_mass = 0, side_2_mass = 0, side_3_mass = 0;

    switch (direction)
    {
        case 0: /* x+ */
            /* core mass */
            for (x = px; x <= px+2*semi_side_len-2; x++)
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                    for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                        core_mass += mass[x][y][z];
            /* corner mass */
            corner_1_mass += mass[px+2*semi_side_len-1][py-semi_side_len][pz-semi_side_len];
            corner_1_mass += mass[px+2*semi_side_len-1][py+semi_side_len][pz-semi_side_len];
            corner_1_mass += mass[px+2*semi_side_len-1][py-semi_side_len][pz+semi_side_len];
            corner_1_mass += mass[px+2*semi_side_len-1][py+semi_side_len][pz+semi_side_len];
            corner_2_mass += mass[px+2*semi_side_len][py-semi_side_len][pz-semi_side_len];
            corner_2_mass += mass[px+2*semi_side_len][py+semi_side_len][pz-semi_side_len];
            corner_2_mass += mass[px+2*semi_side_len][py-semi_side_len][pz+semi_side_len];
            corner_2_mass += mass[px+2*semi_side_len][py+semi_side_len][pz+semi_side_len];
            /* edge mass */
            for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
            {
                edge_1_mass += mass[px+2*semi_side_len-1][y][pz-semi_side_len];
                edge_1_mass += mass[px+2*semi_side_len-1][y][pz+semi_side_len];
                edge_2_mass += mass[px+2*semi_side_len][y][pz-semi_side_len];
                edge_2_mass += mass[px+2*semi_side_len][y][pz+semi_side_len];
            }
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
            {
                edge_1_mass += mass[px+2*semi_side_len-1][py-semi_side_len][z];
                edge_1_mass += mass[px+2*semi_side_len-1][py+semi_side_len][z];
                edge_2_mass += mass[px+2*semi_side_len][py-semi_side_len][z];
                edge_2_mass += mass[px+2*semi_side_len][py+semi_side_len][z];
            }
            for (x = px; x <= px+2*semi_side_len-2; x++)
            {
                edge_3_mass += mass[x][py-semi_side_len][pz-semi_side_len];
                edge_3_mass += mass[x][py+semi_side_len][pz-semi_side_len];
                edge_3_mass += mass[x][py-semi_side_len][pz+semi_side_len];
                edge_3_mass += mass[x][py+semi_side_len][pz+semi_side_len];
            }
            /* side mass */
            for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                {
                    side_1_mass += mass[px+2*semi_side_len-1][y][z];
                    side_2_mass += mass[px+2*semi_side_len][y][z];
                }
            for (x = px; x <= px+2*semi_side_len-2; x++)
            {
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                {
                    side_3_mass += mass[x][y][pz-semi_side_len];
                    side_3_mass += mass[x][y][pz+semi_side_len];
                }
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                {
                    side_3_mass += mass[x][py-semi_side_len][z];
                    side_3_mass += mass[x][py+semi_side_len][z];
                }
            }
            break;
        case 1: /* x- */
            /* core mass */
            for (x = px; x >= px-2*semi_side_len+2; x--)
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                    for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                        core_mass += mass[x][y][z];
            /* corner mass */
            corner_1_mass += mass[px-2*semi_side_len+1][py-semi_side_len][pz-semi_side_len];
            corner_1_mass += mass[px-2*semi_side_len+1][py+semi_side_len][pz-semi_side_len];
            corner_1_mass += mass[px-2*semi_side_len+1][py-semi_side_len][pz+semi_side_len];
            corner_1_mass += mass[px-2*semi_side_len+1][py+semi_side_len][pz+semi_side_len];
            corner_2_mass += mass[px-2*semi_side_len][py-semi_side_len][pz-semi_side_len];
            corner_2_mass += mass[px-2*semi_side_len][py+semi_side_len][pz-semi_side_len];
            corner_2_mass += mass[px-2*semi_side_len][py-semi_side_len][pz+semi_side_len];
            corner_2_mass += mass[px-2*semi_side_len][py+semi_side_len][pz+semi_side_len];
            /* edge mass */
            for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
            {
                edge_1_mass += mass[px-2*semi_side_len+1][y][pz-semi_side_len];
                edge_1_mass += mass[px-2*semi_side_len+1][y][pz+semi_side_len];
                edge_2_mass += mass[px-2*semi_side_len][y][pz-semi_side_len];
                edge_2_mass += mass[px-2*semi_side_len][y][pz+semi_side_len];
            }
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
            {
                edge_1_mass += mass[px-2*semi_side_len+1][py-semi_side_len][z];
                edge_1_mass += mass[px-2*semi_side_len+1][py+semi_side_len][z];
                edge_2_mass += mass[px-2*semi_side_len][py-semi_side_len][z];
                edge_2_mass += mass[px-2*semi_side_len][py+semi_side_len][z];
            }
            for (x = px; x >= px-2*semi_side_len+2; x--)
            {
                edge_3_mass += mass[x][py-semi_side_len][pz-semi_side_len];
                edge_3_mass += mass[x][py+semi_side_len][pz-semi_side_len];
                edge_3_mass += mass[x][py-semi_side_len][pz+semi_side_len];
                edge_3_mass += mass[x][py+semi_side_len][pz+semi_side_len];
            }
            /* side mass */
            for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                {
                    side_1_mass += mass[px-2*semi_side_len+1][y][z];
                    side_2_mass += mass[px-2*semi_side_len][y][z];
                }
            for (x = px; x >= px-2*semi_side_len+2; x--)
            {
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                {
                    side_3_mass += mass[x][y][pz-semi_side_len];
                    side_3_mass += mass[x][y][pz+semi_side_len];
                }
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                {
                    side_3_mass += mass[x][py-semi_side_len][z];
                    side_3_mass += mass[x][py+semi_side_len][z];
                }
            }
            break;
        case 2: /* y+ */
            /* core mass */
            for (y = py; y <= py+2*semi_side_len-2; y++)
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                    for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                        core_mass += mass[x][y][z];
            /* corner mass */
            corner_1_mass += mass[px-semi_side_len][py+2*semi_side_len-1][pz-semi_side_len];
            corner_1_mass += mass[px+semi_side_len][py+2*semi_side_len-1][pz-semi_side_len];
            corner_1_mass += mass[px-semi_side_len][py+2*semi_side_len-1][pz+semi_side_len];
            corner_1_mass += mass[px+semi_side_len][py+2*semi_side_len-1][pz+semi_side_len];
            corner_2_mass += mass[px-semi_side_len][py+2*semi_side_len][pz-semi_side_len];
            corner_2_mass += mass[px+semi_side_len][py+2*semi_side_len][pz-semi_side_len];
            corner_2_mass += mass[px-semi_side_len][py+2*semi_side_len][pz+semi_side_len];
            corner_2_mass += mass[px+semi_side_len][py+2*semi_side_len][pz+semi_side_len];
            /* edge mass */
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
            {
                edge_1_mass += mass[px-semi_side_len][py+2*semi_side_len-1][z];
                edge_1_mass += mass[px+semi_side_len][py+2*semi_side_len-1][z];
                edge_2_mass += mass[px-semi_side_len][py+2*semi_side_len][z];
                edge_2_mass += mass[px+semi_side_len][py+2*semi_side_len][z];
            }
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
            {
                edge_1_mass += mass[x][py+2*semi_side_len-1][pz-semi_side_len];
                edge_1_mass += mass[x][py+2*semi_side_len-1][pz+semi_side_len];
                edge_2_mass += mass[x][py+2*semi_side_len][pz-semi_side_len];
                edge_2_mass += mass[x][py+2*semi_side_len][pz+semi_side_len];
            }
            for (y = py; y <= py+2*semi_side_len-2; y++)
            {
                edge_3_mass += mass[px-semi_side_len][y][pz-semi_side_len];
                edge_3_mass += mass[px+semi_side_len][y][pz-semi_side_len];
                edge_3_mass += mass[px-semi_side_len][y][pz+semi_side_len];
                edge_3_mass += mass[px+semi_side_len][y][pz+semi_side_len];
            }
            /* side mass */
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                {
                    side_1_mass += mass[x][py+2*semi_side_len-1][z];
                    side_2_mass += mass[x][py+2*semi_side_len][z];
                }
            for (y = py; y <= py+2*semi_side_len-2; y++)
            {
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                {
                    side_3_mass += mass[px-semi_side_len][y][z];
                    side_3_mass += mass[px+semi_side_len][y][z];
                }
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                {
                    side_3_mass += mass[x][y][pz-semi_side_len];
                    side_3_mass += mass[x][y][pz+semi_side_len];
                }
            }
            break;
        case 3: /* y- */
            /* core mass */
            for (y = py; y >= py-2*semi_side_len+2; y--)
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                    for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                        core_mass += mass[x][y][z];
            /* corner mass */
            corner_1_mass += mass[px-semi_side_len][py-2*semi_side_len+1][pz-semi_side_len];
            corner_1_mass += mass[px+semi_side_len][py-2*semi_side_len+1][pz-semi_side_len];
            corner_1_mass += mass[px-semi_side_len][py-2*semi_side_len+1][pz+semi_side_len];
            corner_1_mass += mass[px+semi_side_len][py-2*semi_side_len+1][pz+semi_side_len];
            corner_2_mass += mass[px-semi_side_len][py-2*semi_side_len][pz-semi_side_len];
            corner_2_mass += mass[px+semi_side_len][py-2*semi_side_len][pz-semi_side_len];
            corner_2_mass += mass[px-semi_side_len][py-2*semi_side_len][pz+semi_side_len];
            corner_2_mass += mass[px+semi_side_len][py-2*semi_side_len][pz+semi_side_len];
            /* edge mass */
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
            {
                edge_1_mass += mass[px-semi_side_len][py-2*semi_side_len+1][z];
                edge_1_mass += mass[px+semi_side_len][py-2*semi_side_len+1][z];
                edge_2_mass += mass[px-semi_side_len][py-2*semi_side_len][z];
                edge_2_mass += mass[px+semi_side_len][py-2*semi_side_len][z];
            }
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
            {
                edge_1_mass += mass[x][py-2*semi_side_len+1][pz-semi_side_len];
                edge_1_mass += mass[x][py-2*semi_side_len+1][pz+semi_side_len];
                edge_2_mass += mass[x][py-2*semi_side_len][pz-semi_side_len];
                edge_2_mass += mass[x][py-2*semi_side_len][pz+semi_side_len];
            }
            for (y = py; y >= py-2*semi_side_len+2; y--)
            {
                edge_3_mass += mass[px-semi_side_len][y][pz-semi_side_len];
                edge_3_mass += mass[px+semi_side_len][y][pz-semi_side_len];
                edge_3_mass += mass[px-semi_side_len][y][pz+semi_side_len];
                edge_3_mass += mass[px+semi_side_len][y][pz+semi_side_len];
            }
            /* side mass */
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                {
                    side_1_mass += mass[x][py-2*semi_side_len+1][z];
                    side_2_mass += mass[x][py-2*semi_side_len][z];
                }
            for (y = py; y >= py-2*semi_side_len+2; y--)
            {
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                {
                    side_3_mass += mass[px-semi_side_len][y][z];
                    side_3_mass += mass[px+semi_side_len][y][z];
                }
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                {
                    side_3_mass += mass[x][y][pz-semi_side_len];
                    side_3_mass += mass[x][y][pz+semi_side_len];
                }
            }
            break;
        case 4: /* z+ */
            /* core mass */
            for (z = pz; z <= pz+2*semi_side_len-2; z++)
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                    for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                        core_mass += mass[x][y][z];
            /* corner mass */
            corner_1_mass += mass[px-semi_side_len][py-semi_side_len][pz+2*semi_side_len-1];
            corner_1_mass += mass[px+semi_side_len][py-semi_side_len][pz+2*semi_side_len-1];
            corner_1_mass += mass[px-semi_side_len][py+semi_side_len][pz+2*semi_side_len-1];
            corner_1_mass += mass[px+semi_side_len][py+semi_side_len][pz+2*semi_side_len-1];
            corner_2_mass += mass[px-semi_side_len][py-semi_side_len][pz+2*semi_side_len];
            corner_2_mass += mass[px+semi_side_len][py-semi_side_len][pz+2*semi_side_len];
            corner_2_mass += mass[px-semi_side_len][py+semi_side_len][pz+2*semi_side_len];
            corner_2_mass += mass[px+semi_side_len][py+semi_side_len][pz+2*semi_side_len];
            /* edge mass */
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
            {
                edge_1_mass += mass[x][py-semi_side_len][pz+2*semi_side_len-1];
                edge_1_mass += mass[x][py+semi_side_len][pz+2*semi_side_len-1];
                edge_2_mass += mass[x][py-semi_side_len][pz+2*semi_side_len];
                edge_2_mass += mass[x][py+semi_side_len][pz+2*semi_side_len];
            }
            for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
            {
                edge_1_mass += mass[px-semi_side_len][y][pz+2*semi_side_len-1];
                edge_1_mass += mass[px+semi_side_len][y][pz+2*semi_side_len-1];
                edge_2_mass += mass[px-semi_side_len][y][pz+2*semi_side_len];
                edge_2_mass += mass[px+semi_side_len][y][pz+2*semi_side_len];
            }
            for (z = pz; z <= pz+2*semi_side_len-2; z++)
            {
                edge_3_mass += mass[px-semi_side_len][py-semi_side_len][z];
                edge_3_mass += mass[px+semi_side_len][py-semi_side_len][z];
                edge_3_mass += mass[px-semi_side_len][py+semi_side_len][z];
                edge_3_mass += mass[px+semi_side_len][py+semi_side_len][z];
            }
            /* side mass */
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                {
                    side_1_mass += mass[x][y][pz+2*semi_side_len-1];
                    side_2_mass += mass[x][y][pz+2*semi_side_len];
                }
            for (z = pz; z <= pz+2*semi_side_len-2; z++)
            {
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                {
                    side_3_mass += mass[x][py-semi_side_len][z];
                    side_3_mass += mass[x][py+semi_side_len][z];
                }
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                {
                    side_3_mass += mass[px-semi_side_len][y][z];
                    side_3_mass += mass[px+semi_side_len][y][z];
                }
            }
            break;
        case 5: /* z- */
            /* core mass */
            for (z = pz; z >= pz-2*semi_side_len+2; z--)
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                    for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                        core_mass += mass[x][y][z];
            /* corner mass */
            corner_1_mass += mass[px-semi_side_len][py-semi_side_len][pz-2*semi_side_len+1];
            corner_1_mass += mass[px+semi_side_len][py-semi_side_len][pz-2*semi_side_len+1];
            corner_1_mass += mass[px-semi_side_len][py+semi_side_len][pz-2*semi_side_len+1];
            corner_1_mass += mass[px+semi_side_len][py+semi_side_len][pz-2*semi_side_len+1];
            corner_2_mass += mass[px-semi_side_len][py-semi_side_len][pz-2*semi_side_len];
            corner_2_mass += mass[px+semi_side_len][py-semi_side_len][pz-2*semi_side_len];
            corner_2_mass += mass[px-semi_side_len][py+semi_side_len][pz-2*semi_side_len];
            corner_2_mass += mass[px+semi_side_len][py+semi_side_len][pz-2*semi_side_len];
            /* edge mass */
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
            {
                edge_1_mass += mass[x][py-semi_side_len][pz-2*semi_side_len+1];
                edge_1_mass += mass[x][py+semi_side_len][pz-2*semi_side_len+1];
                edge_2_mass += mass[x][py-semi_side_len][pz-2*semi_side_len];
                edge_2_mass += mass[x][py+semi_side_len][pz-2*semi_side_len];
            }
            for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
            {
                edge_1_mass += mass[px-semi_side_len][y][pz-2*semi_side_len+1];
                edge_1_mass += mass[px+semi_side_len][y][pz-2*semi_side_len+1];
                edge_2_mass += mass[px-semi_side_len][y][pz-2*semi_side_len];
                edge_2_mass += mass[px+semi_side_len][y][pz-2*semi_side_len];
            }
            for (z = pz; z >= pz-2*semi_side_len+2; z--)
            {
                edge_3_mass += mass[px-semi_side_len][py-semi_side_len][z];
                edge_3_mass += mass[px+semi_side_len][py-semi_side_len][z];
                edge_3_mass += mass[px-semi_side_len][py+semi_side_len][z];
                edge_3_mass += mass[px+semi_side_len][py+semi_side_len][z];
            }
            /* side mass */
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                {
                    side_1_mass += mass[x][y][pz-2*semi_side_len+1];
                    side_2_mass += mass[x][y][pz-2*semi_side_len];
                }
            for (z = pz; z >= pz-2*semi_side_len+2; z--)
            {
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                {
                    side_3_mass += mass[x][py-semi_side_len][z];
                    side_3_mass += mass[x][py+semi_side_len][z];
                }
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                {
                    side_3_mass += mass[px-semi_side_len][y][z];
                    side_3_mass += mass[px+semi_side_len][y][z];
                }
            }
            break;
        default:
            {
                printf("ERROR! Direction is wrong when compute part mass of the type 3 voxels!\n");
                fflush(stdout);/* %% */
                exit(2);
            }
    }

    inner_layer_part_mass[0] = core_mass;
    inner_layer_part_mass[1] = side_1_mass+0.50*side_3_mass;
    inner_layer_part_mass[2] = 0.50*edge_1_mass+0.25*edge_3_mass;
    inner_layer_part_mass[3] = 0.25*corner_1_mass;

    outer_layer_part_mass[0] = core_mass+side_1_mass+0.50*side_3_mass+0.50*edge_1_mass+0.25*edge_3_mass+0.25*corner_1_mass;
    outer_layer_part_mass[1] = side_2_mass+0.50*side_3_mass+0.50*edge_1_mass+0.50*edge_2_mass+0.50*edge_3_mass+0.50*corner_1_mass+0.25*corner_2_mass;
    outer_layer_part_mass[2] = 0.50*edge_2_mass+0.25*edge_3_mass+0.25*corner_1_mass+0.50*corner_2_mass;
    outer_layer_part_mass[3] = 0.25*corner_2_mass;

    return 0;
}

float compute_type_3_fraction(const float *inner_layer_part_mass, const float *outer_layer_part_mass, const float required_mass)
{
    /**************************************************************************************
     * Compute the fraction of the outmost layer which is added to the averaging voume.   *
     * The fraction which is in the range from 0 to 1 indicate the inner layer is divided *
     * while 1 to 2 indecate the outer layer is divided.                                  *
     **************************************************************************************/

    if (inner_layer_part_mass[0]+inner_layer_part_mass[1]+inner_layer_part_mass[2]+inner_layer_part_mass[3] > required_mass)
        return compute_fraction(inner_layer_part_mass, required_mass);
    else
        return 1+compute_fraction(outer_layer_part_mass, required_mass);
}

float compute_fraction(const float *part_mass, const float required_mass)
{
    /*************************************************************************************
     * Compute the fraction of the outmost layer which is added to the averaging volume. *
     *************************************************************************************/

    float mass_threshold = required_mass/100;
    float current_mass = 0;
    float fraction = 0.5;
    float fraction_temp = fraction;
    float fraction_step = fraction/2;
    float delta_mass = mass_threshold+1;

    while (delta_mass > mass_threshold)
    {
        fraction = fraction_temp;
        current_mass = part_mass[3]*fraction*fraction*fraction+part_mass[2]*fraction*fraction+part_mass[1]*fraction+part_mass[0];
        if (current_mass > required_mass)
            fraction_temp = fraction-fraction_step;
        else if (current_mass < required_mass)
            fraction_temp = fraction+fraction_step;
        else
            break;
        fraction_step /= 2;
        if (current_mass > required_mass)
            delta_mass = current_mass-required_mass;
        else if (current_mass < required_mass)
            delta_mass = required_mass-current_mass;
        else
            delta_mass = 0;
    }

    return fraction;
}

float compute_type_3_sar(float ***local_sar, float ***mass, const int direction, const int semi_side_len, const float fraction, const int layer_marker, const float cube_mass, const int px, const int py, const int pz)
{
    /******************************************************
     * Compute the mass average SAR of the type 3 voxels. *
     ******************************************************/

    float core_energy = 0;
    float corner_1_energy = 0, corner_2_energy = 0;
    float edge_1_energy = 0, edge_2_energy = 0, edge_3_energy = 0;
    float side_1_energy = 0, side_2_energy = 0, side_3_energy = 0;
    float inner_layer_part_energy[4] = {0, 0, 0, 0}, outer_layer_part_energy[4] = {0, 0, 0, 0};
    int x = 0, y = 0, z = 0;

    switch (direction)
    {
        case 0: /* x+ */
            /* core energy */
            for (x = px; x <= px+2*semi_side_len-2; x++)
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                    for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                        core_energy += local_sar[x][y][z]*mass[x][y][z];
            /* corner energy */
            corner_1_energy += local_sar[px+2*semi_side_len-1][py-semi_side_len][pz-semi_side_len]*mass[px+2*semi_side_len-1][py-semi_side_len][pz-semi_side_len];
            corner_1_energy += local_sar[px+2*semi_side_len-1][py+semi_side_len][pz-semi_side_len]*mass[px+2*semi_side_len-1][py+semi_side_len][pz-semi_side_len];
            corner_1_energy += local_sar[px+2*semi_side_len-1][py-semi_side_len][pz+semi_side_len]*mass[px+2*semi_side_len-1][py-semi_side_len][pz+semi_side_len];
            corner_1_energy += local_sar[px+2*semi_side_len-1][py+semi_side_len][pz+semi_side_len]*mass[px+2*semi_side_len-1][py+semi_side_len][pz+semi_side_len];
            corner_2_energy += local_sar[px+2*semi_side_len][py-semi_side_len][pz-semi_side_len]*mass[px+2*semi_side_len][py-semi_side_len][pz-semi_side_len];
            corner_2_energy += local_sar[px+2*semi_side_len][py+semi_side_len][pz-semi_side_len]*mass[px+2*semi_side_len][py+semi_side_len][pz-semi_side_len];
            corner_2_energy += local_sar[px+2*semi_side_len][py-semi_side_len][pz+semi_side_len]*mass[px+2*semi_side_len][py-semi_side_len][pz+semi_side_len];
            corner_2_energy += local_sar[px+2*semi_side_len][py+semi_side_len][pz+semi_side_len]*mass[px+2*semi_side_len][py+semi_side_len][pz+semi_side_len];
            /* edge energy */
            for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
            {
                edge_1_energy += local_sar[px+2*semi_side_len-1][y][pz-semi_side_len]*mass[px+2*semi_side_len-1][y][pz-semi_side_len];
                edge_1_energy += local_sar[px+2*semi_side_len-1][y][pz+semi_side_len]*mass[px+2*semi_side_len-1][y][pz+semi_side_len];
                edge_2_energy += local_sar[px+2*semi_side_len][y][pz-semi_side_len]*mass[px+2*semi_side_len][y][pz-semi_side_len];
                edge_2_energy += local_sar[px+2*semi_side_len][y][pz+semi_side_len]*mass[px+2*semi_side_len][y][pz+semi_side_len];
            }
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
            {
                edge_1_energy += local_sar[px+2*semi_side_len-1][py-semi_side_len][z]*mass[px+2*semi_side_len-1][py-semi_side_len][z];
                edge_1_energy += local_sar[px+2*semi_side_len-1][py+semi_side_len][z]*mass[px+2*semi_side_len-1][py+semi_side_len][z];
                edge_2_energy += local_sar[px+2*semi_side_len][py-semi_side_len][z]*mass[px+2*semi_side_len][py-semi_side_len][z];
                edge_2_energy += local_sar[px+2*semi_side_len][py+semi_side_len][z]*mass[px+2*semi_side_len][py+semi_side_len][z];
            }
            for (x = px; x <= px+2*semi_side_len-2; x++)
            {
                edge_3_energy += local_sar[x][py-semi_side_len][pz-semi_side_len]*mass[x][py-semi_side_len][pz-semi_side_len];
                edge_3_energy += local_sar[x][py+semi_side_len][pz-semi_side_len]*mass[x][py+semi_side_len][pz-semi_side_len];
                edge_3_energy += local_sar[x][py-semi_side_len][pz+semi_side_len]*mass[x][py-semi_side_len][pz+semi_side_len];
                edge_3_energy += local_sar[x][py+semi_side_len][pz+semi_side_len]*mass[x][py+semi_side_len][pz+semi_side_len];
            }
            /* side energy */
            for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                {
                    side_1_energy += local_sar[px+2*semi_side_len-1][y][z]*mass[px+2*semi_side_len-1][y][z];
                    side_2_energy += local_sar[px+2*semi_side_len][y][z]*mass[px+2*semi_side_len][y][z];
                }
            for (x = px; x <= px+2*semi_side_len-2; x++)
            {
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                {
                    side_3_energy += local_sar[x][y][pz-semi_side_len]*mass[x][y][pz-semi_side_len];
                    side_3_energy += local_sar[x][y][pz+semi_side_len]*mass[x][y][pz+semi_side_len];
                }
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                {
                    side_3_energy += local_sar[x][py-semi_side_len][z]*mass[x][py-semi_side_len][z];
                    side_3_energy += local_sar[x][py+semi_side_len][z]*mass[x][py+semi_side_len][z];
                }
            }
            break;
        case 1: /* x- */
            /* core energy */
            for (x = px; x >= px-2*semi_side_len+2; x--)
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                    for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                        core_energy += local_sar[x][y][z]*mass[x][y][z];
            /* corner energy */
            corner_1_energy += local_sar[px-2*semi_side_len+1][py-semi_side_len][pz-semi_side_len]*mass[px-2*semi_side_len+1][py-semi_side_len][pz-semi_side_len];
            corner_1_energy += local_sar[px-2*semi_side_len+1][py+semi_side_len][pz-semi_side_len]*mass[px-2*semi_side_len+1][py+semi_side_len][pz-semi_side_len];
            corner_1_energy += local_sar[px-2*semi_side_len+1][py-semi_side_len][pz+semi_side_len]*mass[px-2*semi_side_len+1][py-semi_side_len][pz+semi_side_len];
            corner_1_energy += local_sar[px-2*semi_side_len+1][py+semi_side_len][pz+semi_side_len]*mass[px-2*semi_side_len+1][py+semi_side_len][pz+semi_side_len];
            corner_2_energy += local_sar[px-2*semi_side_len][py-semi_side_len][pz-semi_side_len]*mass[px-2*semi_side_len][py-semi_side_len][pz-semi_side_len];
            corner_2_energy += local_sar[px-2*semi_side_len][py+semi_side_len][pz-semi_side_len]*mass[px-2*semi_side_len][py+semi_side_len][pz-semi_side_len];
            corner_2_energy += local_sar[px-2*semi_side_len][py-semi_side_len][pz+semi_side_len]*mass[px-2*semi_side_len][py-semi_side_len][pz+semi_side_len];
            corner_2_energy += local_sar[px-2*semi_side_len][py+semi_side_len][pz+semi_side_len]*mass[px-2*semi_side_len][py+semi_side_len][pz+semi_side_len];
            /* edge energy */
            for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
            {
                edge_1_energy += local_sar[px-2*semi_side_len+1][y][pz-semi_side_len]*mass[px-2*semi_side_len+1][y][pz-semi_side_len];
                edge_1_energy += local_sar[px-2*semi_side_len+1][y][pz+semi_side_len]*mass[px-2*semi_side_len+1][y][pz+semi_side_len];
                edge_2_energy += local_sar[px-2*semi_side_len][y][pz-semi_side_len]*mass[px-2*semi_side_len][y][pz-semi_side_len];
                edge_2_energy += local_sar[px-2*semi_side_len][y][pz+semi_side_len]*mass[px-2*semi_side_len][y][pz+semi_side_len];
            }
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
            {
                edge_1_energy += local_sar[px-2*semi_side_len+1][py-semi_side_len][z]*mass[px-2*semi_side_len+1][py-semi_side_len][z];
                edge_1_energy += local_sar[px-2*semi_side_len+1][py+semi_side_len][z]*mass[px-2*semi_side_len+1][py+semi_side_len][z];
                edge_2_energy += local_sar[px-2*semi_side_len][py-semi_side_len][z]*mass[px-2*semi_side_len][py-semi_side_len][z];
                edge_2_energy += local_sar[px-2*semi_side_len][py+semi_side_len][z]*mass[px-2*semi_side_len][py+semi_side_len][z];
            }
            for (x = px; x >= px-2*semi_side_len+2; x--)
            {
                edge_3_energy += local_sar[x][py-semi_side_len][pz-semi_side_len]*mass[x][py-semi_side_len][pz-semi_side_len];
                edge_3_energy += local_sar[x][py+semi_side_len][pz-semi_side_len]*mass[x][py+semi_side_len][pz-semi_side_len];
                edge_3_energy += local_sar[x][py-semi_side_len][pz+semi_side_len]*mass[x][py-semi_side_len][pz+semi_side_len];
                edge_3_energy += local_sar[x][py+semi_side_len][pz+semi_side_len]*mass[x][py+semi_side_len][pz+semi_side_len];
            }
            /* side energy */
            for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                {
                    side_1_energy += local_sar[px-2*semi_side_len+1][y][z]*mass[px-2*semi_side_len+1][y][z];
                    side_2_energy += local_sar[px-2*semi_side_len][y][z]*mass[px-2*semi_side_len][y][z];
                }
            for (x = px; x >= px-2*semi_side_len+2; x--)
            {
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                {
                    side_3_energy += local_sar[x][y][pz-semi_side_len]*mass[x][y][pz-semi_side_len];
                    side_3_energy += local_sar[x][y][pz+semi_side_len]*mass[x][y][pz+semi_side_len];
                }
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                {
                    side_3_energy += local_sar[x][py-semi_side_len][z]*mass[x][py-semi_side_len][z];
                    side_3_energy += local_sar[x][py+semi_side_len][z]*mass[x][py+semi_side_len][z];
                }
            }
            break;
        case 2: /* y+ */
            /* core energy */
            for (y = py; y <= py+2*semi_side_len-2; y++)
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                    for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                        core_energy += local_sar[x][y][z]*mass[x][y][z];
            /* corner energy */
            corner_1_energy += local_sar[px-semi_side_len][py+2*semi_side_len-1][pz-semi_side_len]*mass[px-semi_side_len][py+2*semi_side_len-1][pz-semi_side_len];
            corner_1_energy += local_sar[px+semi_side_len][py+2*semi_side_len-1][pz-semi_side_len]*mass[px+semi_side_len][py+2*semi_side_len-1][pz-semi_side_len];
            corner_1_energy += local_sar[px-semi_side_len][py+2*semi_side_len-1][pz+semi_side_len]*mass[px-semi_side_len][py+2*semi_side_len-1][pz+semi_side_len];
            corner_1_energy += local_sar[px+semi_side_len][py+2*semi_side_len-1][pz+semi_side_len]*mass[px+semi_side_len][py+2*semi_side_len-1][pz+semi_side_len];
            corner_2_energy += local_sar[px-semi_side_len][py+2*semi_side_len][pz-semi_side_len]*mass[px-semi_side_len][py+2*semi_side_len][pz-semi_side_len];
            corner_2_energy += local_sar[px+semi_side_len][py+2*semi_side_len][pz-semi_side_len]*mass[px+semi_side_len][py+2*semi_side_len][pz-semi_side_len];
            corner_2_energy += local_sar[px-semi_side_len][py+2*semi_side_len][pz+semi_side_len]*mass[px-semi_side_len][py+2*semi_side_len][pz+semi_side_len];
            corner_2_energy += local_sar[px+semi_side_len][py+2*semi_side_len][pz+semi_side_len]*mass[px+semi_side_len][py+2*semi_side_len][pz+semi_side_len];
            /* edge energy */
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
            {
                edge_1_energy += local_sar[px-semi_side_len][py+2*semi_side_len-1][z]*mass[px-semi_side_len][py+2*semi_side_len-1][z];
                edge_1_energy += local_sar[px+semi_side_len][py+2*semi_side_len-1][z]*mass[px+semi_side_len][py+2*semi_side_len-1][z];
                edge_2_energy += local_sar[px-semi_side_len][py+2*semi_side_len][z]*mass[px-semi_side_len][py+2*semi_side_len][z];
                edge_2_energy += local_sar[px+semi_side_len][py+2*semi_side_len][z]*mass[px+semi_side_len][py+2*semi_side_len][z];
            }
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
            {
                edge_1_energy += local_sar[x][py+2*semi_side_len-1][pz-semi_side_len]*mass[x][py+2*semi_side_len-1][pz-semi_side_len];
                edge_1_energy += local_sar[x][py+2*semi_side_len-1][pz+semi_side_len]*mass[x][py+2*semi_side_len-1][pz+semi_side_len];
                edge_2_energy += local_sar[x][py+2*semi_side_len][pz-semi_side_len]*mass[x][py+2*semi_side_len][pz-semi_side_len];
                edge_2_energy += local_sar[x][py+2*semi_side_len][pz+semi_side_len]*mass[x][py+2*semi_side_len][pz+semi_side_len];
            }
            for (y = py; y <= py+2*semi_side_len-2; y++)
            {
                edge_3_energy += local_sar[px-semi_side_len][y][pz-semi_side_len]*mass[px-semi_side_len][y][pz-semi_side_len];
                edge_3_energy += local_sar[px+semi_side_len][y][pz-semi_side_len]*mass[px+semi_side_len][y][pz-semi_side_len];
                edge_3_energy += local_sar[px-semi_side_len][y][pz+semi_side_len]*mass[px-semi_side_len][y][pz+semi_side_len];
                edge_3_energy += local_sar[px+semi_side_len][y][pz+semi_side_len]*mass[px+semi_side_len][y][pz+semi_side_len];
            }
            /* side energy */
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                {
                    side_1_energy += local_sar[x][py+2*semi_side_len-1][z]*mass[x][py+2*semi_side_len-1][z];
                    side_2_energy += local_sar[x][py+2*semi_side_len][z]*mass[x][py+2*semi_side_len][z];
                }
            for (y = py; y <= py+2*semi_side_len-2; y++)
            {
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                {
                    side_3_energy += local_sar[px-semi_side_len][y][z]*mass[px-semi_side_len][y][z];
                    side_3_energy += local_sar[px+semi_side_len][y][z]*mass[px+semi_side_len][y][z];
                }
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                {
                    side_3_energy += local_sar[x][y][pz-semi_side_len]*mass[x][y][pz-semi_side_len];
                    side_3_energy += local_sar[x][y][pz+semi_side_len]*mass[x][y][pz+semi_side_len];
                }
            }
            break;
        case 3: /* y- */
            /* core energy */
            for (y = py; y >= py-2*semi_side_len+2; y--)
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                    for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                        core_energy += local_sar[x][y][z]*mass[x][y][x];
            /* corner energy */
            corner_1_energy += local_sar[px-semi_side_len][py-2*semi_side_len+1][pz-semi_side_len]*mass[px-semi_side_len][py-2*semi_side_len+1][pz-semi_side_len];
            corner_1_energy += local_sar[px+semi_side_len][py-2*semi_side_len+1][pz-semi_side_len]*mass[px+semi_side_len][py-2*semi_side_len+1][pz-semi_side_len];
            corner_1_energy += local_sar[px-semi_side_len][py-2*semi_side_len+1][pz+semi_side_len]*mass[px-semi_side_len][py-2*semi_side_len+1][pz+semi_side_len];
            corner_1_energy += local_sar[px+semi_side_len][py-2*semi_side_len+1][pz+semi_side_len]*mass[px+semi_side_len][py-2*semi_side_len+1][pz+semi_side_len];
            corner_2_energy += local_sar[px-semi_side_len][py-2*semi_side_len][pz-semi_side_len]*mass[px-semi_side_len][py-2*semi_side_len][pz-semi_side_len];
            corner_2_energy += local_sar[px+semi_side_len][py-2*semi_side_len][pz-semi_side_len]*mass[px+semi_side_len][py-2*semi_side_len][pz-semi_side_len];
            corner_2_energy += local_sar[px-semi_side_len][py-2*semi_side_len][pz+semi_side_len]*mass[px-semi_side_len][py-2*semi_side_len][pz+semi_side_len];
            corner_2_energy += local_sar[px+semi_side_len][py-2*semi_side_len][pz+semi_side_len]*mass[px+semi_side_len][py-2*semi_side_len][pz+semi_side_len];
            /* edge energy */
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
            {
                edge_1_energy += local_sar[px-semi_side_len][py-2*semi_side_len+1][z]*mass[px-semi_side_len][py-2*semi_side_len+1][z];
                edge_1_energy += local_sar[px+semi_side_len][py-2*semi_side_len+1][z]*mass[px+semi_side_len][py-2*semi_side_len+1][z];
                edge_2_energy += local_sar[px-semi_side_len][py-2*semi_side_len][z]*mass[px-semi_side_len][py-2*semi_side_len][z];
                edge_2_energy += local_sar[px+semi_side_len][py-2*semi_side_len][z]*mass[px+semi_side_len][py-2*semi_side_len][z];
            }
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
            {
                edge_1_energy += local_sar[x][py-2*semi_side_len+1][pz-semi_side_len]*mass[x][py-2*semi_side_len+1][pz-semi_side_len];
                edge_1_energy += local_sar[x][py-2*semi_side_len+1][pz+semi_side_len]*mass[x][py-2*semi_side_len+1][pz+semi_side_len];
                edge_2_energy += local_sar[x][py-2*semi_side_len][pz-semi_side_len]*mass[x][py-2*semi_side_len][pz-semi_side_len];
                edge_2_energy += local_sar[x][py-2*semi_side_len][pz+semi_side_len]*mass[x][py-2*semi_side_len][pz+semi_side_len];
            }
            for (y = py; y >= py-2*semi_side_len+2; y--)
            {
                edge_3_energy += local_sar[px-semi_side_len][y][pz-semi_side_len]*mass[px-semi_side_len][y][pz-semi_side_len];
                edge_3_energy += local_sar[px+semi_side_len][y][pz-semi_side_len]*mass[px+semi_side_len][y][pz-semi_side_len];
                edge_3_energy += local_sar[px-semi_side_len][y][pz+semi_side_len]*mass[px-semi_side_len][y][pz+semi_side_len];
                edge_3_energy += local_sar[px+semi_side_len][y][pz+semi_side_len]*mass[px+semi_side_len][y][pz+semi_side_len];
            }
            /* side energy */
            for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                {
                    side_1_energy += local_sar[x][py-2*semi_side_len+1][z]*mass[x][py-2*semi_side_len+1][z];
                    side_2_energy += local_sar[x][py-2*semi_side_len][z]*mass[x][py-2*semi_side_len][z];
                }
            for (y = py; y >= py-2*semi_side_len+2; y--)
            {
                for (z = pz-semi_side_len+1; z <= pz+semi_side_len-1; z++)
                {
                    side_3_energy += local_sar[px-semi_side_len][y][z]*mass[px-semi_side_len][y][z];
                    side_3_energy += local_sar[px+semi_side_len][y][z]*mass[px+semi_side_len][y][z];
                }
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                {
                    side_3_energy += local_sar[x][y][pz-semi_side_len]*mass[x][y][pz-semi_side_len];
                    side_3_energy += local_sar[x][y][pz+semi_side_len]*mass[x][y][pz+semi_side_len];
                }
            }
            break;
        case 4: /* z+ */
            /* core energy */
            for (z = pz; z <= pz+2*semi_side_len-2; z++)
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                    for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                        core_energy += local_sar[x][y][z]*mass[x][y][z];
            /* corner energy */
            corner_1_energy += local_sar[px-semi_side_len][py-semi_side_len][pz+2*semi_side_len-1]*mass[px-semi_side_len][py-semi_side_len][pz+2*semi_side_len-1];
            corner_1_energy += local_sar[px+semi_side_len][py-semi_side_len][pz+2*semi_side_len-1]*mass[px+semi_side_len][py-semi_side_len][pz+2*semi_side_len-1];
            corner_1_energy += local_sar[px-semi_side_len][py+semi_side_len][pz+2*semi_side_len-1]*mass[px-semi_side_len][py+semi_side_len][pz+2*semi_side_len-1];
            corner_1_energy += local_sar[px+semi_side_len][py+semi_side_len][pz+2*semi_side_len-1]*mass[px+semi_side_len][py+semi_side_len][pz+2*semi_side_len-1];
            corner_2_energy += local_sar[px-semi_side_len][py-semi_side_len][pz+2*semi_side_len]*mass[px-semi_side_len][py-semi_side_len][pz+2*semi_side_len];
            corner_2_energy += local_sar[px+semi_side_len][py-semi_side_len][pz+2*semi_side_len]*mass[px+semi_side_len][py-semi_side_len][pz+2*semi_side_len];
            corner_2_energy += local_sar[px-semi_side_len][py+semi_side_len][pz+2*semi_side_len]*mass[px-semi_side_len][py+semi_side_len][pz+2*semi_side_len];
            corner_2_energy += local_sar[px+semi_side_len][py+semi_side_len][pz+2*semi_side_len]*mass[px+semi_side_len][py+semi_side_len][pz+2*semi_side_len];
            /* edge energy */
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
            {
                edge_1_energy += local_sar[x][py-semi_side_len][pz+2*semi_side_len-1]*mass[x][py-semi_side_len][pz+2*semi_side_len-1];
                edge_1_energy += local_sar[x][py+semi_side_len][pz+2*semi_side_len-1]*mass[x][py+semi_side_len][pz+2*semi_side_len-1];
                edge_2_energy += local_sar[x][py-semi_side_len][pz+2*semi_side_len]*mass[x][py-semi_side_len][pz+2*semi_side_len];
                edge_2_energy += local_sar[x][py+semi_side_len][pz+2*semi_side_len]*mass[x][py+semi_side_len][pz+2*semi_side_len];
            }
            for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
            {
                edge_1_energy += local_sar[px-semi_side_len][y][pz+2*semi_side_len-1]*mass[px-semi_side_len][y][pz+2*semi_side_len-1];
                edge_1_energy += local_sar[px+semi_side_len][y][pz+2*semi_side_len-1]*mass[px+semi_side_len][y][pz+2*semi_side_len-1];
                edge_2_energy += local_sar[px-semi_side_len][y][pz+2*semi_side_len]*mass[px-semi_side_len][y][pz+2*semi_side_len];
                edge_2_energy += local_sar[px+semi_side_len][y][pz+2*semi_side_len]*mass[px+semi_side_len][y][pz+2*semi_side_len];
            }
            for (z = pz; z <= pz+2*semi_side_len-2; z++)
            {
                edge_3_energy += local_sar[px-semi_side_len][py-semi_side_len][z]*mass[px-semi_side_len][py-semi_side_len][z];
                edge_3_energy += local_sar[px+semi_side_len][py-semi_side_len][z]*mass[px+semi_side_len][py-semi_side_len][z];
                edge_3_energy += local_sar[px-semi_side_len][py+semi_side_len][z]*mass[px-semi_side_len][py+semi_side_len][z];
                edge_3_energy += local_sar[px+semi_side_len][py+semi_side_len][z]*mass[px+semi_side_len][py+semi_side_len][z];
            }
            /* side energy */
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                {
                    side_1_energy += local_sar[x][y][pz+2*semi_side_len-1]*mass[x][y][pz+2*semi_side_len-1];
                    side_2_energy += local_sar[x][y][pz+2*semi_side_len]*mass[x][y][pz+2*semi_side_len];
                }
            for (z = pz; z <= pz+2*semi_side_len-2; z++)
            {
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                {
                    side_3_energy += local_sar[x][py-semi_side_len][z]*mass[x][py-semi_side_len][z];
                    side_3_energy += local_sar[x][py+semi_side_len][z]*mass[x][py+semi_side_len][z];
                }
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                {
                    side_3_energy += local_sar[px-semi_side_len][y][z]*mass[px-semi_side_len][y][z];
                    side_3_energy += local_sar[px+semi_side_len][y][z]*mass[px+semi_side_len][y][z];
                }
            }
            break;
        case 5: /* z- */
            /* core energy */
            for (z = pz; z >= pz-2*semi_side_len+2; z--)
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                    for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                        core_energy += local_sar[x][y][z]*mass[x][y][z];
            /* corner energy */
            corner_1_energy += local_sar[px-semi_side_len][py-semi_side_len][pz-2*semi_side_len+1]*mass[px-semi_side_len][py-semi_side_len][pz-2*semi_side_len+1];
            corner_1_energy += local_sar[px+semi_side_len][py-semi_side_len][pz-2*semi_side_len+1]*mass[px+semi_side_len][py-semi_side_len][pz-2*semi_side_len+1];
            corner_1_energy += local_sar[px-semi_side_len][py+semi_side_len][pz-2*semi_side_len+1]*mass[px-semi_side_len][py+semi_side_len][pz-2*semi_side_len+1];
            corner_1_energy += local_sar[px+semi_side_len][py+semi_side_len][pz-2*semi_side_len+1]*mass[px+semi_side_len][py+semi_side_len][pz-2*semi_side_len+1];
            corner_2_energy += local_sar[px-semi_side_len][py-semi_side_len][pz-2*semi_side_len]*mass[px-semi_side_len][py-semi_side_len][pz-2*semi_side_len];
            corner_2_energy += local_sar[px+semi_side_len][py-semi_side_len][pz-2*semi_side_len]*mass[px+semi_side_len][py-semi_side_len][pz-2*semi_side_len];
            corner_2_energy += local_sar[px-semi_side_len][py+semi_side_len][pz-2*semi_side_len]*mass[px-semi_side_len][py+semi_side_len][pz-2*semi_side_len];
            corner_2_energy += local_sar[px+semi_side_len][py+semi_side_len][pz-2*semi_side_len]*mass[px+semi_side_len][py+semi_side_len][pz-2*semi_side_len];
            /* edge energy */
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
            {
                edge_1_energy += local_sar[x][py-semi_side_len][pz-2*semi_side_len+1]*mass[x][py-semi_side_len][pz-2*semi_side_len+1];
                edge_1_energy += local_sar[x][py+semi_side_len][pz-2*semi_side_len+1]*mass[x][py+semi_side_len][pz-2*semi_side_len+1];
                edge_2_energy += local_sar[x][py-semi_side_len][pz-2*semi_side_len]*mass[x][py-semi_side_len][pz-2*semi_side_len];
                edge_2_energy += local_sar[x][py+semi_side_len][pz-2*semi_side_len]*mass[x][py+semi_side_len][pz-2*semi_side_len];
            }
            for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
            {
                edge_1_energy += local_sar[px-semi_side_len][y][pz-2*semi_side_len+1]*mass[px-semi_side_len][y][pz-2*semi_side_len+1];
                edge_1_energy += local_sar[px+semi_side_len][y][pz-2*semi_side_len+1]*mass[px+semi_side_len][y][pz-2*semi_side_len+1];
                edge_2_energy += local_sar[px-semi_side_len][y][pz-2*semi_side_len]*mass[px-semi_side_len][y][pz-2*semi_side_len];
                edge_2_energy += local_sar[px+semi_side_len][y][pz-2*semi_side_len]*mass[px+semi_side_len][y][pz-2*semi_side_len];
            }
            for (z = pz; z >= pz-2*semi_side_len+2; z--)
            {
                edge_3_energy += local_sar[px-semi_side_len][py-semi_side_len][z]*mass[px-semi_side_len][py-semi_side_len][z];
                edge_3_energy += local_sar[px+semi_side_len][py-semi_side_len][z]*mass[px+semi_side_len][py-semi_side_len][z];
                edge_3_energy += local_sar[px-semi_side_len][py+semi_side_len][z]*mass[px-semi_side_len][py+semi_side_len][z];
                edge_3_energy += local_sar[px+semi_side_len][py+semi_side_len][z]*mass[px+semi_side_len][py+semi_side_len][z];
            }
            /* side energy */
            for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                {
                    side_1_energy += local_sar[x][y][pz-2*semi_side_len+1]*mass[x][y][pz-2*semi_side_len+1];
                    side_2_energy += local_sar[x][y][pz-2*semi_side_len]*mass[x][y][pz-2*semi_side_len];
                }
            for (z = pz; z >= pz-2*semi_side_len+2; z--)
            {
                for (x = px-semi_side_len+1; x <= px+semi_side_len-1; x++)
                {
                    side_3_energy += local_sar[x][py-semi_side_len][z]*mass[x][py-semi_side_len][z];
                    side_3_energy += local_sar[x][py+semi_side_len][z]*mass[x][py+semi_side_len][z];
                }
                for (y = py-semi_side_len+1; y <= py+semi_side_len-1; y++)
                {
                    side_3_energy += local_sar[px-semi_side_len][y][z]*mass[px-semi_side_len][y][z];
                    side_3_energy += local_sar[px+semi_side_len][y][z]*mass[px+semi_side_len][y][z];
                }
            }
            break;
        default:
            {
                printf("Direction is wrong when compute SAR of the type 3 voxels!\n");
                fflush(stdout);/* %% */
                exit(2);
            }
    }

    if (layer_marker == 0)
    {
        inner_layer_part_energy[0] = core_energy;
        inner_layer_part_energy[1] = side_1_energy+0.50*side_3_energy;
        inner_layer_part_energy[2] = 0.50*edge_1_energy+0.25*edge_3_energy;
        inner_layer_part_energy[3] = 0.25*corner_1_energy;

        return (inner_layer_part_energy[0]+inner_layer_part_energy[1]*fraction+inner_layer_part_energy[2]*fraction*fraction+inner_layer_part_energy[3]*fraction*fraction*fraction);
    }
    else
    {
        outer_layer_part_energy[0] = core_energy+side_1_energy+0.50*side_3_energy+0.50*edge_1_energy+0.25*edge_3_energy+0.25*corner_1_energy;
        outer_layer_part_energy[1] = side_2_energy+0.50*side_3_energy+0.50*edge_1_energy+0.50*edge_2_energy+0.50*edge_3_energy+0.50*corner_1_energy+0.25*corner_2_energy;
        outer_layer_part_energy[2] = 0.50*edge_2_energy+0.25*edge_3_energy+0.25*corner_1_energy+0.50*corner_2_energy;
        outer_layer_part_energy[3] = 0.25*corner_2_energy;

        return (outer_layer_part_energy[0]+outer_layer_part_energy[1]*fraction+outer_layer_part_energy[2]*fraction*fraction+outer_layer_part_energy[3]*fraction*fraction*fraction);
    }
}

