#include <stdio.h>
#include <stdlib.h>
#include "computeType1Voxel.h"
#include "computeFraction.h"

int computeType1Voxel(double ***LocalSAR, double ***Mass, int ***UsedMarker, int ***SemiSideLength, double ***MassAveragedSAR, const double requiredMass, const int *SpaceDim, FILE *fpLog)
{
    int sumType1Voxel = 0, sumTissueVoxel = 0;
    int emptySide = 0, outsideVoxel = 0, sumLayerType1Voxel = 0;
    int semiSideLength = 0, minSemiSideLength = 1000, maxSemiSideLength = 0;
    double minCubeMass = 1000, maxCubeMass = 0;
    double maxFraction = 0, minFraction = 2, maxMassAveragedSAR = 0;
    double mass = 0, fraction = 0;
    double PartMass[4] = {0, 0, 0, 0};
    int x = 0, y = 0, z = 0, xx = 0, yy = 0, zz = 0, t = 0, i = 0;

    for (z = 0; z < SpaceDim[2]; z++)
    {
        if (z%20 == 0)
            printf("Type 1: z = %d / %d ", z+1, SpaceDim[2]);
        fprintf(fpLog, "Type 1: z = %d / %d ", z+1, SpaceDim[2]);
		fflush(stdout);
        sumLayerType1Voxel = 0;
        for (y = 0; y < SpaceDim[1]; y++)
        {
            for (x = 0; x < SpaceDim[0]; x++)
            {
                if (Mass[x][y][z] > 0)
                {
                    sumTissueVoxel += 1;
                    semiSideLength = 1;
                    mass = 0;
                    for (t = 0; t < 4; t++)
                        PartMass[t] = 0;

                    while (mass < requiredMass)
                    {
                        emptySide = findEmptySide(Mass, semiSideLength, x, y, z, SpaceDim);
                        if (emptySide == 0)
                        {
                            mass = 0;
                            for (xx = x-semiSideLength; xx <= x+semiSideLength; xx++)
                                for (yy = y-semiSideLength; yy <= y+semiSideLength; yy++)
                                    for (zz = z-semiSideLength; zz <= z+semiSideLength; zz++)
                                        mass += Mass[xx][yy][zz];
                            if (mass < requiredMass)
                                semiSideLength = semiSideLength+1;
                            else if (mass > requiredMass)
                            {
                                computeMass(Mass, PartMass, semiSideLength, x, y, z);
                                fraction = computeFraction(PartMass, requiredMass);
                                break;
                            }
                            else
                            {
                                fraction = 1;
                                break;
                            }
                        }
                        else
                            break;
                    }
                    if (mass >= requiredMass)
                    {
                        if (fraction > maxFraction)
                            maxFraction = fraction;
                        if (fraction < minFraction)
                            minFraction = fraction;
                        if (mass != 0 && mass < minCubeMass)
                            minCubeMass = mass;
                        if (mass > maxCubeMass)
                            maxCubeMass = mass;
                        MassAveragedSAR[x][y][z] = computeSAR(LocalSAR, Mass, mass, semiSideLength, fraction, x, y, z);
                        sumType1Voxel += 1;
                        sumLayerType1Voxel += 1;
                        for (xx = x-semiSideLength; xx <= x+semiSideLength; xx++)
                            for (yy = y-semiSideLength; yy <= y+semiSideLength; yy++)
                                for (zz = z-semiSideLength; zz <= z+semiSideLength; zz++)
                                    if (Mass[xx][yy][zz] > 0)
                                        UsedMarker[xx][yy][zz] += 1;
                        SemiSideLength[x][y][z] = semiSideLength;
                        if (semiSideLength > maxSemiSideLength && semiSideLength != 0)
                            maxSemiSideLength = semiSideLength;
                        if (semiSideLength < minSemiSideLength && semiSideLength != 0)
                            minSemiSideLength = semiSideLength;
                    }
                }
            }
        }
        if (z%20 == 0)
            printf("sumVoxel = %d / %d\n", sumLayerType1Voxel, sumType1Voxel);
		fflush(stdout);
        fprintf(fpLog, "sumVoxel = %d / %d\n", sumLayerType1Voxel, sumType1Voxel);
    }

    printf("-------- Statistical data --------\n");
    fprintf(fpLog, "-------- Statistical data --------\n");
    printf("sumTissueVoxel = %d\n", sumTissueVoxel);
    fprintf(fpLog, "sumTissueVoxel = %d\n", sumTissueVoxel);
    printf("Fraction = [%15.12lf, %15.12lf]\n", minFraction, maxFraction);
    fprintf(fpLog, "Fraction = [%15.12lf, %15.12lf]\n", minFraction, maxFraction);
    printf("CubeMass = [%15.12lf, %15.12lf] g\n", minCubeMass, maxCubeMass);
    fprintf(fpLog, "CubeMass = [%15.12lf, %15.12lf] g\n", minCubeMass, maxCubeMass);
	fflush(stdout);
    int maxSARx = 0, maxSARy = 0, maxSARz = 0;
    int *SumSemiSideLength = (int *)malloc((maxSemiSideLength+1)*sizeof(int));
    for (i = 0; i <= maxSemiSideLength; i++)
        SumSemiSideLength[i] = 0;
    for (z = 0; z < SpaceDim[2]; z++)
    {
        for (y = 0; y < SpaceDim[1]; y++)
        {
            for (x = 0; x < SpaceDim[0]; x++)
            {
                SumSemiSideLength[SemiSideLength[x][y][z]] += 1;
                if (MassAveragedSAR[x][y][z] > maxMassAveragedSAR)
                {
                    maxMassAveragedSAR = MassAveragedSAR[x][y][z];
                    maxSARx = x;
                    maxSARy = y;
                    maxSARz = z;
                }
            }
        }
    }
    printf("SumSemiSideLength = ");
    for (i = 0; i <= maxSemiSideLength; i++)
	{
        printf("%d ", SumSemiSideLength[i]);
		fflush(stdout);
	}
    printf("\n");
    printf("semiSideLength = [%d, %d]\n", minSemiSideLength, maxSemiSideLength);
    fprintf(fpLog, "semiSideLength = [%d, %d]\n", minSemiSideLength, maxSemiSideLength);
    printf("maxMassAveragedSAR = %15.12lf W/kg at (%d, %d, %d)\n", maxMassAveragedSAR, maxSARx, maxSARy, maxSARz);
    fprintf(fpLog, "maxMassAveragedSAR = %15.12lf W/kg at (%d, %d, %d)\n", maxMassAveragedSAR, maxSARx, maxSARy, maxSARz);

    free(SumSemiSideLength);

    return 0;
}

int findEmptySide(double ***Mass, const int semiSideLen, const int px, const int py, const int pz, const int *SpaceDim)
{
    double mass = 0;
    int x = 0, y = 0, z = 0;

    if (px-semiSideLen < 0 || px+semiSideLen >= SpaceDim[0] || py-semiSideLen < 0 || py+semiSideLen >= SpaceDim[1] || pz-semiSideLen < 0 || pz+semiSideLen >= SpaceDim[2])
        return 1;

    mass = 0;
    for (z = pz-semiSideLen; z <= pz+semiSideLen; z++)
        for (y = py-semiSideLen; y <= py+semiSideLen; y++)
            mass += Mass[px-semiSideLen][y][z];
    if (mass == 0)
        return 1;

    mass = 0;
    for (z = pz-semiSideLen; z <= pz+semiSideLen; z++)
        for (y = py-semiSideLen; y <= py+semiSideLen; y++)
            mass += Mass[px+semiSideLen][y][z];
    if (mass == 0)
        return 1;

    mass = 0;
    for (x = px-semiSideLen; x <= px+semiSideLen; x++)
        for (z = pz-semiSideLen; z <= pz+semiSideLen; z++)
            mass += Mass[x][py-semiSideLen][z];
    if (mass == 0)
        return 1;

    mass = 0;
    for (x = px-semiSideLen; x <= px+semiSideLen; x++)
        for (z = pz-semiSideLen; z <= pz+semiSideLen; z++)
            mass += Mass[x][py+semiSideLen][z];
    if (mass == 0)
        return 1;

    mass = 0;
    for (y = py-semiSideLen; y <= py+semiSideLen; y++)
        for (x = px-semiSideLen; x <= px+semiSideLen; x++)
            mass += Mass[x][y][pz-semiSideLen];
    if (mass == 0)
        return 1;

    mass = 0;
    for (y = py-semiSideLen; y <= py+semiSideLen; y++)
        for (x = px-semiSideLen; x <= px+semiSideLen; x++)
            mass += Mass[x][y][pz+semiSideLen];
    if (mass == 0)
        return 1;

    return 0;
}

int computeMass(double ***Mass, double *PartMass, const int semiSideLen, const int px, const int py, const int pz)
{
    int x = 0, y = 0, z = 0;

    /* core mass */
    for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
        for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                PartMass[0] += Mass[x][y][z];

    /* 6 sides */
    for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
        for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
            PartMass[1] += Mass[x][y][pz-semiSideLen]+Mass[x][y][pz+semiSideLen];
    for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
        for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
            PartMass[1] += Mass[px-semiSideLen][y][z]+Mass[px+semiSideLen][y][z];
    for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
        for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
            PartMass[1] += Mass[x][py-semiSideLen][z]+Mass[x][py+semiSideLen][z];

    /* 12 edges */
    for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
    {
        PartMass[2] += Mass[x][py-semiSideLen][pz-semiSideLen];
        PartMass[2] += Mass[x][py+semiSideLen][pz-semiSideLen];
        PartMass[2] += Mass[x][py-semiSideLen][pz+semiSideLen];
        PartMass[2] += Mass[x][py+semiSideLen][pz+semiSideLen];
    }
    for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
    {
        PartMass[2] += Mass[px-semiSideLen][y][pz-semiSideLen];
        PartMass[2] += Mass[px+semiSideLen][y][pz-semiSideLen];
        PartMass[2] += Mass[px-semiSideLen][y][pz+semiSideLen];
        PartMass[2] += Mass[px+semiSideLen][y][pz+semiSideLen];
    }
    for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
    {
        PartMass[2] += Mass[px-semiSideLen][py-semiSideLen][z];
        PartMass[2] += Mass[px+semiSideLen][py-semiSideLen][z];
        PartMass[2] += Mass[px-semiSideLen][py+semiSideLen][z];
        PartMass[2] += Mass[px+semiSideLen][py+semiSideLen][z];
    }

    /* 8 corner voxels */
    PartMass[3] += Mass[px-semiSideLen][py-semiSideLen][pz-semiSideLen];
    PartMass[3] += Mass[px+semiSideLen][py-semiSideLen][pz-semiSideLen];
    PartMass[3] += Mass[px-semiSideLen][py+semiSideLen][pz-semiSideLen];
    PartMass[3] += Mass[px+semiSideLen][py+semiSideLen][pz-semiSideLen];
    PartMass[3] += Mass[px-semiSideLen][py-semiSideLen][pz+semiSideLen];
    PartMass[3] += Mass[px+semiSideLen][py-semiSideLen][pz+semiSideLen];
    PartMass[3] += Mass[px-semiSideLen][py+semiSideLen][pz+semiSideLen];
    PartMass[3] += Mass[px+semiSideLen][py+semiSideLen][pz+semiSideLen];

    return 0;
}

double computeSAR(double ***LocalSAR, double ***Mass, const double cubeMass, const int semiSideLen, const double fraction, const int px, const int py, const int pz)
{
    double massAveragedSAR = 0;
    double corePower = 0;
    int x = 0, y = 0, z = 0;
    for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
        for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                if (Mass[x][y][z] != 0)
                    corePower += LocalSAR[x][y][z]*Mass[x][y][z];

    double cornerPower = 0;
    cornerPower += LocalSAR[px-semiSideLen][py-semiSideLen][pz-semiSideLen]*Mass[px-semiSideLen][py-semiSideLen][pz-semiSideLen];
    cornerPower += LocalSAR[px+semiSideLen][py-semiSideLen][pz-semiSideLen]*Mass[px+semiSideLen][py-semiSideLen][pz-semiSideLen];
    cornerPower += LocalSAR[px-semiSideLen][py+semiSideLen][pz-semiSideLen]*Mass[px-semiSideLen][py+semiSideLen][pz-semiSideLen];
    cornerPower += LocalSAR[px+semiSideLen][py+semiSideLen][pz-semiSideLen]*Mass[px+semiSideLen][py+semiSideLen][pz-semiSideLen];
    cornerPower += LocalSAR[px-semiSideLen][py-semiSideLen][pz+semiSideLen]*Mass[px-semiSideLen][py-semiSideLen][pz+semiSideLen];
    cornerPower += LocalSAR[px+semiSideLen][py-semiSideLen][pz+semiSideLen]*Mass[px+semiSideLen][py-semiSideLen][pz+semiSideLen];
    cornerPower += LocalSAR[px-semiSideLen][py+semiSideLen][pz+semiSideLen]*Mass[px-semiSideLen][py+semiSideLen][pz+semiSideLen];
    cornerPower += LocalSAR[px+semiSideLen][py+semiSideLen][pz+semiSideLen]*Mass[px+semiSideLen][py+semiSideLen][pz+semiSideLen];

    double edgePower = 0;
    for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
    {
        edgePower += LocalSAR[x][py-semiSideLen][pz-semiSideLen]*Mass[x][py-semiSideLen][pz-semiSideLen];
        edgePower += LocalSAR[x][py+semiSideLen][pz-semiSideLen]*Mass[x][py+semiSideLen][pz-semiSideLen];
        edgePower += LocalSAR[x][py-semiSideLen][pz+semiSideLen]*Mass[x][py-semiSideLen][pz+semiSideLen];
        edgePower += LocalSAR[x][py+semiSideLen][pz+semiSideLen]*Mass[x][py+semiSideLen][pz+semiSideLen];
    }
    for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
    {
        edgePower += LocalSAR[px-semiSideLen][y][pz-semiSideLen]*Mass[px-semiSideLen][y][pz-semiSideLen];
        edgePower += LocalSAR[px+semiSideLen][y][pz-semiSideLen]*Mass[px+semiSideLen][y][pz-semiSideLen];
        edgePower += LocalSAR[px-semiSideLen][y][pz+semiSideLen]*Mass[px-semiSideLen][y][pz+semiSideLen];
        edgePower += LocalSAR[px+semiSideLen][y][pz+semiSideLen]*Mass[px+semiSideLen][y][pz+semiSideLen];
    }
    for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
    {
        edgePower += LocalSAR[px-semiSideLen][py-semiSideLen][z]*Mass[px-semiSideLen][py-semiSideLen][z];
        edgePower += LocalSAR[px+semiSideLen][py-semiSideLen][z]*Mass[px+semiSideLen][py-semiSideLen][z];
        edgePower += LocalSAR[px-semiSideLen][py+semiSideLen][z]*Mass[px-semiSideLen][py+semiSideLen][z];
        edgePower += LocalSAR[px+semiSideLen][py+semiSideLen][z]*Mass[px+semiSideLen][py+semiSideLen][z];
    }

    double sidePower = 0; /* 6 sides. */
    for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
        for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
            sidePower += LocalSAR[x][y][pz-semiSideLen]*Mass[x][y][pz-semiSideLen]+LocalSAR[x][y][pz+semiSideLen]*Mass[x][y][pz+semiSideLen];
    for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
        for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
            sidePower += LocalSAR[px-semiSideLen][y][z]*Mass[px-semiSideLen][y][z]+LocalSAR[px+semiSideLen][y][z]*Mass[px+semiSideLen][y][z];
    for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
        for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
            sidePower += LocalSAR[x][py-semiSideLen][z]*Mass[x][py-semiSideLen][z]+LocalSAR[x][py+semiSideLen][z]*Mass[x][py+semiSideLen][z];

    massAveragedSAR = (corePower+cornerPower*fraction*fraction*fraction+edgePower*fraction*fraction+sidePower*fraction)/cubeMass; /* W/requiredMass -> W/kg. */

    return massAveragedSAR;
}

