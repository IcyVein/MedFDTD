#include <stdio.h>
#include "computeType2Voxel.h"

int computeType2Voxel(double ***MassAveragedSAR, int ***UsedMarker, int ***SemiSideLength, double ***Mass, const int *SpaceDim, FILE *fpLog)
{
    int sumType2Voxel = 0;
    double maxCubeMassAveragedSAR = 0;
    int x = 0, y = 0, z = 0, xx = 0, yy = 0, zz = 0;
    int maxSemiSideLength = 0;
    int sumLayerType2Voxel = 0;

    for (z = 0; z < SpaceDim[2]; z++)
        for (y = 0; y < SpaceDim[1]; y++)
            for (x = 0; x < SpaceDim[0]; x++)
                if (SemiSideLength[x][y][z] > maxSemiSideLength)
                    maxSemiSideLength = SemiSideLength[x][y][z];

    for (z = 0; z < SpaceDim[2]; z++)
    {
        if (z%20 == 0)
            printf("Type 2: z = %d / %d ", z+1, SpaceDim[2]);
		fflush(stdout);
        fprintf(fpLog, "Type 2: z = %d / %d ", z+1, SpaceDim[2]);
        sumLayerType2Voxel = 0;
        for (y = 0; y < SpaceDim[1]; y++)
        {
            for (x = 0; x < SpaceDim[0]; x++)
            {
                if (Mass[x][y][z] > 0 && MassAveragedSAR[x][y][z] == 0 && UsedMarker[x][y][z] > 0)
                {
                    maxCubeMassAveragedSAR = 0;
                    for (zz = z-maxSemiSideLength; zz <= z+maxSemiSideLength; zz++)
                    {
                        for (yy = y-maxSemiSideLength; yy <= y+maxSemiSideLength; yy++)
                        {
                            for (xx = x-maxSemiSideLength; xx <= x+maxSemiSideLength; xx++)
                            {
                                if (xx < 0 || xx >= SpaceDim[0] || yy < 0 || yy >= SpaceDim[1] || zz < 0 || zz >= SpaceDim[2])
                                    continue;
                                else
                                {
                                    if (SemiSideLength[xx][yy][zz] != 0)
                                    {
                                        if (x >= xx-SemiSideLength[xx][yy][zz] && x <= xx+SemiSideLength[xx][yy][zz] && y >= yy-SemiSideLength[xx][yy][zz] && y <= yy+SemiSideLength[xx][yy][zz] && z >= zz-SemiSideLength[xx][yy][zz] && z <= zz+SemiSideLength[xx][yy][zz])
                                        {
                                            if (MassAveragedSAR[xx][yy][zz] > maxCubeMassAveragedSAR)
                                                maxCubeMassAveragedSAR = MassAveragedSAR[xx][yy][zz];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (maxCubeMassAveragedSAR != 0)
                    {
                       MassAveragedSAR[x][y][z] = maxCubeMassAveragedSAR;
                       sumType2Voxel += 1;
                       sumLayerType2Voxel += 1;
                    }
                }
            }
        }
        if (z%20 == 0)
            printf("sumVoxel = %d / %d\n", sumLayerType2Voxel, sumType2Voxel);
		fflush(stdout);
        fprintf(fpLog, "sumVoxel = %d / %d\n", sumLayerType2Voxel, sumType2Voxel);
    }

    printf("sumType2Voxel = %d\n", sumType2Voxel);
	fflush(stdout);
    fprintf(fpLog, "sumType2Voxel = %d\n", sumType2Voxel);

    return 0;
}

