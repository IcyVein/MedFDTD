#include <stdio.h>
#include <stdlib.h>
#include "computeType3Voxel.h"
#include "computeFraction.h"

int computeType3Voxel(double ***LocalSAR, double ***Mass, int ***UsedMarker, double ***MassAveragedSAR, const double requiredMass, const int *SpaceDim, FILE *fpLog)
{
    int x = 0, y = 0, z = 0, i = 0;
    int layerMarker = 0, semiSideLen = 0, direction = 0;
    int sumCannotComputeVoxel = 0, sumType3Voxel = 0;
    double CubeMass[6] = {0, 0, 0, 0, 0, 0};
    double InnerLayerPartMass[4] = {0, 0, 0, 0}, OuterLayerPartMass[4] = {0, 0, 0, 0};
    double cubeMass = 0, minType3CubeMass = 1000, maxType3CubeMass = 0;
	double sumMass = 0, minSumMass = 1000, maxSumMass = 0;
	double minType3MassAveragedSAR = 1000, maxType3MassAveragedSAR = 0;
    double fraction = 0, minType3Fraction = 1000, maxType3Fraction = 0;

    for (z = 0; z < SpaceDim[2]; z++)
    {
        for (y = 0; y < SpaceDim[1]; y++)
        {
            for (x = 0; x < SpaceDim[0]; x++)
            {
                if (Mass[x][y][z] != 0 && UsedMarker[x][y][z] == 0) /* Type 3 voxel */
                {
                    semiSideLen = 0;
                    cubeMass = 0;
                    while (cubeMass < requiredMass)
                    {
                        semiSideLen += 1;
                        for (i = 0; i < 6; i++)
                            CubeMass[i] = 0;
                        cubeMass = computeType3CubeMass(Mass, CubeMass, semiSideLen, x, y, z, SpaceDim);
                        if (cubeMass == 0)
                        {
                            sumCannotComputeVoxel += 1;
                            fprintf(stderr, "Voxel(%d, %d, %d) can't be computed.\n", x, y, z);
                            fprintf(fpLog, "Voxel(%d, %d, %d) can't be computed.\n", x, y, z);
                            break;
                        }
                    }
                    if (cubeMass >= requiredMass)
                    {
                        for (i = 0; i < 6; i++)
                            if (CubeMass[i] == cubeMass)
                                direction = i;
                        for (i = 0; i < 4; i++)
                        {
                            InnerLayerPartMass[i] = 0;
                            OuterLayerPartMass[i] = 0;
                        }
                        sumMass = computeType3PartMass(InnerLayerPartMass, OuterLayerPartMass, Mass, semiSideLen, direction, x, y, z);
                        fraction = computeType3Fraction(InnerLayerPartMass, OuterLayerPartMass, requiredMass);
                        if (OuterLayerPartMass[0] < requiredMass)
                            layerMarker = 1; /* outer layer */
                        else
                            layerMarker = 0; /* inner layer */
                        MassAveragedSAR[x][y][z] = computeType3SAR(LocalSAR, Mass, direction, semiSideLen, fraction, layerMarker, cubeMass, x, y, z);
                        sumType3Voxel += 1;

						if (cubeMass < minType3CubeMass)
							minType3CubeMass = cubeMass;
						if (cubeMass > maxType3CubeMass)
							maxType3CubeMass = cubeMass;
						if (sumMass < minSumMass)
							minSumMass = sumMass;
						if (sumMass > maxSumMass)
							maxSumMass = sumMass;
                        if (fraction < minType3Fraction)
                            minType3Fraction = fraction;
                        if (fraction > maxType3Fraction)
                            maxType3Fraction = fraction;
                        if (MassAveragedSAR[x][y][z] < minType3MassAveragedSAR)
                            minType3MassAveragedSAR = MassAveragedSAR[x][y][z];
                        if (MassAveragedSAR[x][y][z] > maxType3MassAveragedSAR)
                            maxType3MassAveragedSAR = MassAveragedSAR[x][y][z];
                    }
                }
            }
        }
    }

    printf("Fraction = [%lf, %lf]\n", minType3Fraction, maxType3Fraction);
    fprintf(fpLog, "Fraction = [%lf, %lf]\n", minType3Fraction, maxType3Fraction);
    printf("Type3MassAveragedSAR = [%lf, %lf]\n", minType3MassAveragedSAR, maxType3MassAveragedSAR);
    fprintf(fpLog, "Type3MassAveragedSAR = [%lf, %lf]\n", minType3MassAveragedSAR, maxType3MassAveragedSAR);
    printf("sumCannotComputeVoxel = %d\n", sumCannotComputeVoxel);
    printf("sumType3Voxel = %d\n", sumType3Voxel);
    fprintf(fpLog, "sumType3Voxel = %d\n", sumType3Voxel);
	fflush(stdout);

    return 0;
}

double computeType3CubeMass(double ***Mass, double *CubeMass, const int semiSideLen, const int px, const int py, const int pz, const int *SpaceDim)
{
    int x = 0, y = 0, z = 0, i = 0;
    double cubeMass = 0;

    /* x+ */
    if (px+2*semiSideLen < SpaceDim[0] && py-semiSideLen >= 0 && py+semiSideLen < SpaceDim[1] && pz-semiSideLen >= 0 && pz+semiSideLen < SpaceDim[2])
        for (y = py-semiSideLen; y <= py+semiSideLen; y++)
            for (z = pz-semiSideLen; z <= pz+semiSideLen; z++)
                for (x = px; x <= px+2*semiSideLen; x++)
                    CubeMass[0] += Mass[x][y][z];
    else
        CubeMass[0] = 0;

    /* x- */
    if (px-2*semiSideLen >= 0 && py-semiSideLen >= 0 && py+semiSideLen < SpaceDim[1] && pz-semiSideLen >= 0 && pz+semiSideLen < SpaceDim[2])
        for (y = py-semiSideLen; y <= py+semiSideLen; y++)
            for (z = pz-semiSideLen; z <= pz+semiSideLen; z++)
                for (x = px-2*semiSideLen; x <= px; x++)
                    CubeMass[1] += Mass[x][y][z];
    else
        CubeMass[1] = 0;

    /* y+ */
    if (py+2*semiSideLen < SpaceDim[1] && pz-semiSideLen >= 0 && pz+semiSideLen < SpaceDim[2] && px-semiSideLen >= 0 && px+semiSideLen < SpaceDim[0])
        for (z = pz-semiSideLen; z <= pz+semiSideLen; z++)
            for (x = px-semiSideLen; x <= px+semiSideLen; x++)
                for (y = py; y <= py+2*semiSideLen; y++)
                    CubeMass[2] += Mass[x][y][z];
    else
        CubeMass[2] = 0;

    /* y- */
    if (py-2*semiSideLen >= 0 && pz-semiSideLen >= 0 && pz+semiSideLen < SpaceDim[2] && px-semiSideLen >= 0 && px+semiSideLen < SpaceDim[0])
        for (z = pz-semiSideLen; z <= pz+semiSideLen; z++)
            for (x = px-semiSideLen; x <= px+semiSideLen; x++)
                for (y = py-2*semiSideLen; y <= py; y++)
                    CubeMass[3] += Mass[x][y][z];
    else
        CubeMass[3] = 0;

    /* z+ */
    if (pz+2*semiSideLen < SpaceDim[2] && px-semiSideLen >= 0 && px+semiSideLen < SpaceDim[0] && py-semiSideLen >= 0 && py+semiSideLen < SpaceDim[1])
        for (x = px-semiSideLen; x <= px+semiSideLen; x++)
            for (y = py-semiSideLen; y <= py+semiSideLen; y++)
                for (z = pz; z <= pz+2*semiSideLen; z++)
                    CubeMass[4] += Mass[x][y][z];
    else
        CubeMass[4] = 0;

    /* z- */
    if (pz-2*semiSideLen >= 0 && px-semiSideLen >= 0 && px+semiSideLen < SpaceDim[0] && py-semiSideLen >= 0 && py+semiSideLen < SpaceDim[1])
        for (x = px-semiSideLen; x <= px+semiSideLen; x++)
            for (y = py-semiSideLen; y <= py+semiSideLen; y++)
                for (z = pz-2*semiSideLen; z <= pz; z++)
                    CubeMass[5] += Mass[x][y][z];
    else
        CubeMass[5] = 0;

    for (i = 0; i < 6; i++)
        if (CubeMass[i] > cubeMass)
            cubeMass = CubeMass[i];

    return cubeMass;
}

double computeType3PartMass(double *InnerLayerPartMass, double *OuterLayerPartMass, double ***Mass, const int semiSideLen, const int direction, const int px, const int py, const int pz)
{
    double coreMass = 0;
    double corner1Mass = 0, corner2Mass = 0;
    double edge1Mass = 0, edge2Mass = 0, edge3Mass = 0;
    double side1Mass = 0, side2Mass = 0, side3Mass = 0;
    int x = 0, y = 0, z = 0;

    switch (direction)
    {
        case 0: /* x+ */
            /* core mass */
            for (x = px; x <= px+2*semiSideLen-2; x++)
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                    for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                        coreMass += Mass[x][y][z];
            /* corner mass */
            corner1Mass += Mass[px+2*semiSideLen-1][py-semiSideLen][pz-semiSideLen];
            corner1Mass += Mass[px+2*semiSideLen-1][py+semiSideLen][pz-semiSideLen];
            corner1Mass += Mass[px+2*semiSideLen-1][py-semiSideLen][pz+semiSideLen];
            corner1Mass += Mass[px+2*semiSideLen-1][py+semiSideLen][pz+semiSideLen];
            corner2Mass += Mass[px+2*semiSideLen][py-semiSideLen][pz-semiSideLen];
            corner2Mass += Mass[px+2*semiSideLen][py+semiSideLen][pz-semiSideLen];
            corner2Mass += Mass[px+2*semiSideLen][py-semiSideLen][pz+semiSideLen];
            corner2Mass += Mass[px+2*semiSideLen][py+semiSideLen][pz+semiSideLen];
            /* edge mass */
            for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
            {
                edge1Mass += Mass[px+2*semiSideLen-1][y][pz-semiSideLen];
                edge1Mass += Mass[px+2*semiSideLen-1][y][pz+semiSideLen];
                edge2Mass += Mass[px+2*semiSideLen][y][pz-semiSideLen];
                edge2Mass += Mass[px+2*semiSideLen][y][pz+semiSideLen];
            }
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
            {
                edge1Mass += Mass[px+2*semiSideLen-1][py-semiSideLen][z];
                edge1Mass += Mass[px+2*semiSideLen-1][py+semiSideLen][z];
                edge2Mass += Mass[px+2*semiSideLen][py-semiSideLen][z];
                edge2Mass += Mass[px+2*semiSideLen][py+semiSideLen][z];
            }
            for (x = px; x <= px+2*semiSideLen-2; x++)
            {
                edge3Mass += Mass[x][py-semiSideLen][pz-semiSideLen];
                edge3Mass += Mass[x][py+semiSideLen][pz-semiSideLen];
                edge3Mass += Mass[x][py-semiSideLen][pz+semiSideLen];
                edge3Mass += Mass[x][py+semiSideLen][pz+semiSideLen];
            }
            /* side mass */
			for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
			{
				for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
				{
					side1Mass += Mass[px+2*semiSideLen-1][y][z];
					side2Mass += Mass[px+2*semiSideLen][y][z];
				}
            }
            for (x = px; x <= px+2*semiSideLen-2; x++)
            {
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                {
                    side3Mass += Mass[x][y][pz-semiSideLen];
                    side3Mass += Mass[x][y][pz+semiSideLen];
                }
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                {
                    side3Mass += Mass[x][py-semiSideLen][z];
                    side3Mass += Mass[x][py+semiSideLen][z];
                }
            }
            break;
        case 1: /* x- */
            /* core mass */
            for (x = px; x >= px-2*semiSideLen+2; x--)
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                    for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                        coreMass += Mass[x][y][z];
            /* corner mass */
            corner1Mass += Mass[px-2*semiSideLen+1][py-semiSideLen][pz-semiSideLen];
            corner1Mass += Mass[px-2*semiSideLen+1][py+semiSideLen][pz-semiSideLen];
            corner1Mass += Mass[px-2*semiSideLen+1][py-semiSideLen][pz+semiSideLen];
            corner1Mass += Mass[px-2*semiSideLen+1][py+semiSideLen][pz+semiSideLen];
            corner2Mass += Mass[px-2*semiSideLen][py-semiSideLen][pz-semiSideLen];
            corner2Mass += Mass[px-2*semiSideLen][py+semiSideLen][pz-semiSideLen];
            corner2Mass += Mass[px-2*semiSideLen][py-semiSideLen][pz+semiSideLen];
            corner2Mass += Mass[px-2*semiSideLen][py+semiSideLen][pz+semiSideLen];
            /* edge mass */
            for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
            {
                edge1Mass += Mass[px-2*semiSideLen+1][y][pz-semiSideLen];
                edge1Mass += Mass[px-2*semiSideLen+1][y][pz+semiSideLen];
                edge2Mass += Mass[px-2*semiSideLen][y][pz-semiSideLen];
                edge2Mass += Mass[px-2*semiSideLen][y][pz+semiSideLen];
            }
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
            {
                edge1Mass += Mass[px-2*semiSideLen+1][py-semiSideLen][z];
                edge1Mass += Mass[px-2*semiSideLen+1][py+semiSideLen][z];
                edge2Mass += Mass[px-2*semiSideLen][py-semiSideLen][z];
                edge2Mass += Mass[px-2*semiSideLen][py+semiSideLen][z];
            }
            for (x = px; x >= px-2*semiSideLen+2; x--)
            {
                edge3Mass += Mass[x][py-semiSideLen][pz-semiSideLen];
                edge3Mass += Mass[x][py+semiSideLen][pz-semiSideLen];
                edge3Mass += Mass[x][py-semiSideLen][pz+semiSideLen];
                edge3Mass += Mass[x][py+semiSideLen][pz+semiSideLen];
            }
            /* side mass */
            for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                {
                    side1Mass += Mass[px-2*semiSideLen+1][y][z];
                    side2Mass += Mass[px-2*semiSideLen][y][z];
                }
            for (x = px; x >= px-2*semiSideLen+2; x--)
            {
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                {
                    side3Mass += Mass[x][y][pz-semiSideLen];
                    side3Mass += Mass[x][y][pz+semiSideLen];
                }
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                {
                    side3Mass += Mass[x][py-semiSideLen][z];
                    side3Mass += Mass[x][py+semiSideLen][z];
                }
            }
            break;
        case 2: /* y+ */
            /* core mass */
            for (y = py; y <= py+2*semiSideLen-2; y++)
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                    for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                        coreMass += Mass[x][y][z];
            /* corner mass */
            corner1Mass += Mass[px-semiSideLen][py+2*semiSideLen-1][pz-semiSideLen];
            corner1Mass += Mass[px+semiSideLen][py+2*semiSideLen-1][pz-semiSideLen];
            corner1Mass += Mass[px-semiSideLen][py+2*semiSideLen-1][pz+semiSideLen];
            corner1Mass += Mass[px+semiSideLen][py+2*semiSideLen-1][pz+semiSideLen];
            corner2Mass += Mass[px-semiSideLen][py+2*semiSideLen][pz-semiSideLen];
            corner2Mass += Mass[px+semiSideLen][py+2*semiSideLen][pz-semiSideLen];
            corner2Mass += Mass[px-semiSideLen][py+2*semiSideLen][pz+semiSideLen];
            corner2Mass += Mass[px+semiSideLen][py+2*semiSideLen][pz+semiSideLen];
            /* edge mass */
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
            {
                edge1Mass += Mass[px-semiSideLen][py+2*semiSideLen-1][z];
                edge1Mass += Mass[px+semiSideLen][py+2*semiSideLen-1][z];
                edge2Mass += Mass[px-semiSideLen][py+2*semiSideLen][z];
                edge2Mass += Mass[px+semiSideLen][py+2*semiSideLen][z];
            }
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
            {
                edge1Mass += Mass[x][py+2*semiSideLen-1][pz-semiSideLen];
                edge1Mass += Mass[x][py+2*semiSideLen-1][pz+semiSideLen];
                edge2Mass += Mass[x][py+2*semiSideLen][pz-semiSideLen];
                edge2Mass += Mass[x][py+2*semiSideLen][pz+semiSideLen];
            }
            for (y = py; y <= py+2*semiSideLen-2; y++)
            {
                edge3Mass += Mass[px-semiSideLen][y][pz-semiSideLen];
                edge3Mass += Mass[px+semiSideLen][y][pz-semiSideLen];
                edge3Mass += Mass[px-semiSideLen][y][pz+semiSideLen];
                edge3Mass += Mass[px+semiSideLen][y][pz+semiSideLen];
            }
            /* side mass */
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
			{
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                {
                    side1Mass += Mass[x][py+2*semiSideLen-1][z];
                    side2Mass += Mass[x][py+2*semiSideLen][z];
                }
			}
            for (y = py; y <= py+2*semiSideLen-2; y++)
            {
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                {
                    side3Mass += Mass[px-semiSideLen][y][z];
                    side3Mass += Mass[px+semiSideLen][y][z];
                }
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                {
                    side3Mass += Mass[x][y][pz-semiSideLen];
                    side3Mass += Mass[x][y][pz+semiSideLen];
                }
            }
            break;
        case 3: /* y- */
            /* core mass */
            for (y = py; y >= py-2*semiSideLen+2; y--)
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                    for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                        coreMass += Mass[x][y][x];
            /* corner mass */
            corner1Mass += Mass[px-semiSideLen][py-2*semiSideLen+1][pz-semiSideLen];
            corner1Mass += Mass[px+semiSideLen][py-2*semiSideLen+1][pz-semiSideLen];
            corner1Mass += Mass[px-semiSideLen][py-2*semiSideLen+1][pz+semiSideLen];
            corner1Mass += Mass[px+semiSideLen][py-2*semiSideLen+1][pz+semiSideLen];
            corner2Mass += Mass[px-semiSideLen][py-2*semiSideLen][pz-semiSideLen];
            corner2Mass += Mass[px+semiSideLen][py-2*semiSideLen][pz-semiSideLen];
            corner2Mass += Mass[px-semiSideLen][py-2*semiSideLen][pz+semiSideLen];
            corner2Mass += Mass[px+semiSideLen][py-2*semiSideLen][pz+semiSideLen];
            /* edge mass */
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
            {
                edge1Mass += Mass[px-semiSideLen][py-2*semiSideLen+1][z];
                edge1Mass += Mass[px+semiSideLen][py-2*semiSideLen+1][z];
                edge2Mass += Mass[px-semiSideLen][py-2*semiSideLen][z];
                edge2Mass += Mass[px+semiSideLen][py-2*semiSideLen][z];
            }
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
            {
                edge1Mass += Mass[x][py-2*semiSideLen+1][pz-semiSideLen];
                edge1Mass += Mass[x][py-2*semiSideLen+1][pz+semiSideLen];
                edge2Mass += Mass[x][py-2*semiSideLen][pz-semiSideLen];
                edge2Mass += Mass[x][py-2*semiSideLen][pz+semiSideLen];
            }
            for (y = py; y >= py-2*semiSideLen+2; y--)
            {
                edge3Mass += Mass[px-semiSideLen][y][pz-semiSideLen];
                edge3Mass += Mass[px+semiSideLen][y][pz-semiSideLen];
                edge3Mass += Mass[px-semiSideLen][y][pz+semiSideLen];
                edge3Mass += Mass[px+semiSideLen][y][pz+semiSideLen];
            }
            /* side mass */
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                {
                    side1Mass += Mass[x][py-2*semiSideLen+1][z];
                    side2Mass += Mass[x][py-2*semiSideLen][z];
                }
            for (y = py; y >= py-2*semiSideLen+2; y--)
            {
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                {
                    side3Mass += Mass[px-semiSideLen][y][z];
                    side3Mass += Mass[px+semiSideLen][y][z];
                }
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                {
                    side3Mass += Mass[x][y][pz-semiSideLen];
                    side3Mass += Mass[x][y][pz+semiSideLen];
                }
            }
            break;
        case 4: /* z+ */
            /* core mass */
            for (z = pz; z <= pz+2*semiSideLen-2; z++)
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                    for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                        coreMass += Mass[x][y][z];
            /* corner mass */
            corner1Mass += Mass[px-semiSideLen][py-semiSideLen][pz+2*semiSideLen-1];
            corner1Mass += Mass[px+semiSideLen][py-semiSideLen][pz+2*semiSideLen-1];
            corner1Mass += Mass[px-semiSideLen][py+semiSideLen][pz+2*semiSideLen-1];
            corner1Mass += Mass[px+semiSideLen][py+semiSideLen][pz+2*semiSideLen-1];
            corner2Mass += Mass[px-semiSideLen][py-semiSideLen][pz+2*semiSideLen];
            corner2Mass += Mass[px+semiSideLen][py-semiSideLen][pz+2*semiSideLen];
            corner2Mass += Mass[px-semiSideLen][py+semiSideLen][pz+2*semiSideLen];
            corner2Mass += Mass[px+semiSideLen][py+semiSideLen][pz+2*semiSideLen];
            /* edge mass */
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
            {
                edge1Mass += Mass[x][py-semiSideLen][pz+2*semiSideLen-1];
                edge1Mass += Mass[x][py+semiSideLen][pz+2*semiSideLen-1];
                edge2Mass += Mass[x][py-semiSideLen][pz+2*semiSideLen];
                edge2Mass += Mass[x][py+semiSideLen][pz+2*semiSideLen];
            }
            for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
            {
                edge1Mass += Mass[px-semiSideLen][y][pz+2*semiSideLen-1];
                edge1Mass += Mass[px+semiSideLen][y][pz+2*semiSideLen-1];
                edge2Mass += Mass[px-semiSideLen][y][pz+2*semiSideLen];
                edge2Mass += Mass[px+semiSideLen][y][pz+2*semiSideLen];
            }
            for (z = pz; z <= pz+2*semiSideLen-2; z++)
            {
                edge3Mass += Mass[px-semiSideLen][py-semiSideLen][z];
                edge3Mass += Mass[px+semiSideLen][py-semiSideLen][z];
                edge3Mass += Mass[px-semiSideLen][py+semiSideLen][z];
                edge3Mass += Mass[px+semiSideLen][py+semiSideLen][z];
            }
            /* side mass */
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                {
                    side1Mass += Mass[x][y][pz+2*semiSideLen-1];
                    side2Mass += Mass[x][y][pz+2*semiSideLen];
                }
            for (z = pz; z <= pz+2*semiSideLen-2; z++)
            {
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                {
                    side3Mass += Mass[x][py-semiSideLen][z];
                    side3Mass += Mass[x][py+semiSideLen][z];
                }
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                {
                    side3Mass += Mass[px-semiSideLen][y][z];
                    side3Mass += Mass[px+semiSideLen][y][z];
                }
            }
            break;
        case 5: /* z- */
            /* core mass */
            for (z = pz; z >= pz-2*semiSideLen+2; z--)
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                    for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                        coreMass += Mass[x][y][z];
            /* corner mass */
            corner1Mass += Mass[px-semiSideLen][py-semiSideLen][pz-2*semiSideLen+1];
            corner1Mass += Mass[px+semiSideLen][py-semiSideLen][pz-2*semiSideLen+1];
            corner1Mass += Mass[px-semiSideLen][py+semiSideLen][pz-2*semiSideLen+1];
            corner1Mass += Mass[px+semiSideLen][py+semiSideLen][pz-2*semiSideLen+1];
            corner2Mass += Mass[px-semiSideLen][py-semiSideLen][pz-2*semiSideLen];
            corner2Mass += Mass[px+semiSideLen][py-semiSideLen][pz-2*semiSideLen];
            corner2Mass += Mass[px-semiSideLen][py+semiSideLen][pz-2*semiSideLen];
            corner2Mass += Mass[px+semiSideLen][py+semiSideLen][pz-2*semiSideLen];
            /* edge mass */
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
            {
                edge1Mass += Mass[x][py-semiSideLen][pz-2*semiSideLen+1];
                edge1Mass += Mass[x][py+semiSideLen][pz-2*semiSideLen+1];
                edge2Mass += Mass[x][py-semiSideLen][pz-2*semiSideLen];
                edge2Mass += Mass[x][py+semiSideLen][pz-2*semiSideLen];
            }
            for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
            {
                edge1Mass += Mass[px-semiSideLen][y][pz-2*semiSideLen+1];
                edge1Mass += Mass[px+semiSideLen][y][pz-2*semiSideLen+1];
                edge2Mass += Mass[px-semiSideLen][y][pz-2*semiSideLen];
                edge2Mass += Mass[px+semiSideLen][y][pz-2*semiSideLen];
            }
            for (z = pz; z >= pz-2*semiSideLen+2; z--)
            {
                edge3Mass += Mass[px-semiSideLen][py-semiSideLen][z];
                edge3Mass += Mass[px+semiSideLen][py-semiSideLen][z];
                edge3Mass += Mass[px-semiSideLen][py+semiSideLen][z];
                edge3Mass += Mass[px+semiSideLen][py+semiSideLen][z];
            }
            /* side mass */
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                {
                    side1Mass += Mass[x][y][pz-2*semiSideLen+1];
                    side2Mass += Mass[x][y][pz-2*semiSideLen];
                }
            for (z = pz; z >= pz-2*semiSideLen+2; z--)
            {
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                {
                    side3Mass += Mass[x][py-semiSideLen][z];
                    side3Mass += Mass[x][py+semiSideLen][z];
                }
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                {
                    side3Mass += Mass[px-semiSideLen][y][z];
                    side3Mass += Mass[px+semiSideLen][y][z];
                }
            }
            break;
        default:
            {
                printf("ERROR! Direction is wrong when compute part Mass of the type 3 voxels!\n");
				fflush(stdout);
                exit(2);
            }
    }

    InnerLayerPartMass[0] = coreMass;
    InnerLayerPartMass[1] = side1Mass+0.50*side3Mass;
    InnerLayerPartMass[2] = 0.50*edge1Mass+0.25*edge3Mass;
    InnerLayerPartMass[3] = 0.25*corner1Mass;

    OuterLayerPartMass[0] = coreMass+side1Mass+0.50*side3Mass+0.50*edge1Mass+0.25*edge3Mass+0.25*corner1Mass;
    OuterLayerPartMass[1] = side2Mass+0.50*side3Mass+0.50*edge1Mass+0.50*edge2Mass+0.5*edge3Mass+0.50*corner1Mass+0.25*corner2Mass;
    OuterLayerPartMass[2] = 0.50*edge2Mass+0.25*edge3Mass+0.25*corner1Mass+0.50*corner2Mass;
    OuterLayerPartMass[3] = 0.25*corner2Mass;

    return OuterLayerPartMass[0]+OuterLayerPartMass[1]+OuterLayerPartMass[2]+OuterLayerPartMass[3];
}

double computeType3Fraction(double *InnerLayerPartMass, double *OuterLayerPartMass, const double requiredMass)
{
    double fraction = 0;

    if (OuterLayerPartMass[0] < requiredMass)
        fraction = computeFraction(OuterLayerPartMass, requiredMass);
    else
        fraction = computeFraction(InnerLayerPartMass, requiredMass);

    return fraction;
}

double computeType3SAR(double ***LocalSAR, double ***Mass, const int direction, const int semiSideLen, const double fraction, const int layerMarker, const double cubeMass, const int px, const int py, const int pz)
{
    double massAveragedSAR = 0;
    int x = 0, y = 0, z = 0;
    double corner1Power = 0, corner2Power = 0;
    double edge1Power = 0, edge2Power = 0, edge3Power = 0;
    double side1Power = 0, side2Power = 0, side3Power = 0;
    double corePower = 0;
    double InnerLayerPartPower[4] = {0, 0, 0, 0}, OuterLayerPartPower[4] = {0, 0, 0, 0};

    switch (direction)
    {
        case 0: /* x+ */
            /* core power */
            for (x = px; x <= px+2*semiSideLen-2; x++)
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                    for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                        corePower += LocalSAR[x][y][z]*Mass[x][y][z];
            /* corner power */
            corner1Power += LocalSAR[px+2*semiSideLen-1][py-semiSideLen][pz-semiSideLen]*Mass[px+2*semiSideLen-1][py-semiSideLen][pz-semiSideLen];
            corner1Power += LocalSAR[px+2*semiSideLen-1][py+semiSideLen][pz-semiSideLen]*Mass[px+2*semiSideLen-1][py+semiSideLen][pz-semiSideLen];
            corner1Power += LocalSAR[px+2*semiSideLen-1][py-semiSideLen][pz+semiSideLen]*Mass[px+2*semiSideLen-1][py-semiSideLen][pz+semiSideLen];
            corner1Power += LocalSAR[px+2*semiSideLen-1][py+semiSideLen][pz+semiSideLen]*Mass[px+2*semiSideLen-1][py+semiSideLen][pz+semiSideLen];
            corner2Power += LocalSAR[px+2*semiSideLen][py-semiSideLen][pz-semiSideLen]*Mass[px+2*semiSideLen][py-semiSideLen][pz-semiSideLen];
            corner2Power += LocalSAR[px+2*semiSideLen][py+semiSideLen][pz-semiSideLen]*Mass[px+2*semiSideLen][py+semiSideLen][pz-semiSideLen];
            corner2Power += LocalSAR[px+2*semiSideLen][py-semiSideLen][pz+semiSideLen]*Mass[px+2*semiSideLen][py-semiSideLen][pz+semiSideLen];
            corner2Power += LocalSAR[px+2*semiSideLen][py+semiSideLen][pz+semiSideLen]*Mass[px+2*semiSideLen][py+semiSideLen][pz+semiSideLen];
            /* edge power */
            for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
            {
                edge1Power += LocalSAR[px+2*semiSideLen-1][y][pz-semiSideLen]*Mass[px+2*semiSideLen-1][y][pz-semiSideLen];
                edge1Power += LocalSAR[px+2*semiSideLen-1][y][pz+semiSideLen]*Mass[px+2*semiSideLen-1][y][pz+semiSideLen];
                edge2Power += LocalSAR[px+2*semiSideLen][y][pz-semiSideLen]*Mass[px+2*semiSideLen][y][pz-semiSideLen];
                edge2Power += LocalSAR[px+2*semiSideLen][y][pz+semiSideLen]*Mass[px+2*semiSideLen][y][pz+semiSideLen];
            }
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
            {
                edge1Power += LocalSAR[px+2*semiSideLen-1][py-semiSideLen][z]*Mass[px+2*semiSideLen-1][py-semiSideLen][z];
                edge1Power += LocalSAR[px+2*semiSideLen-1][py+semiSideLen][z]*Mass[px+2*semiSideLen-1][py+semiSideLen][z];
                edge2Power += LocalSAR[px+2*semiSideLen][py-semiSideLen][z]*Mass[px+2*semiSideLen][py-semiSideLen][z];
                edge2Power += LocalSAR[px+2*semiSideLen][py+semiSideLen][z]*Mass[px+2*semiSideLen][py+semiSideLen][z];
            }
            for (x = px; x <= px+2*semiSideLen-2; x++)
            {
                edge3Power += LocalSAR[x][py-semiSideLen][pz-semiSideLen]*Mass[x][py-semiSideLen][pz-semiSideLen];
                edge3Power += LocalSAR[x][py+semiSideLen][pz-semiSideLen]*Mass[x][py+semiSideLen][pz-semiSideLen];
                edge3Power += LocalSAR[x][py-semiSideLen][pz+semiSideLen]*Mass[x][py-semiSideLen][pz+semiSideLen];
                edge3Power += LocalSAR[x][py+semiSideLen][pz+semiSideLen]*Mass[x][py+semiSideLen][pz+semiSideLen];
            }
            /* side power */
			for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
				for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
				{
					side1Power += LocalSAR[px+2*semiSideLen-1][y][z]*Mass[px+2*semiSideLen-1][y][z];
					side2Power += LocalSAR[px+2*semiSideLen][y][z]*Mass[px+2*semiSideLen][y][z];
				}
            for (x = px; x <= px+2*semiSideLen-2; x++)
            {
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                {
                    side3Power += LocalSAR[x][y][pz-semiSideLen]*Mass[x][y][pz-semiSideLen];
                    side3Power += LocalSAR[x][y][pz+semiSideLen]*Mass[x][y][pz+semiSideLen];
                }
                for (z = pz-semiSideLen+1; z <= pz-semiSideLen-1; z++)
                {
                    side3Power += LocalSAR[x][py-semiSideLen][z]*Mass[x][py-semiSideLen][z];
                    side3Power += LocalSAR[x][py+semiSideLen][z]*Mass[x][py+semiSideLen][z];
                }
            }
            break;
        case 1: /* x- */
            /* core power */
            for (x = px; x >= px-2*semiSideLen+2; x--)
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                    for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                        corePower += LocalSAR[x][y][z]*Mass[x][y][z];
            /* corner power */
            corner1Power += LocalSAR[px-2*semiSideLen+1][py-semiSideLen][pz-semiSideLen]*Mass[px-2*semiSideLen+1][py-semiSideLen][pz-semiSideLen];
            corner1Power += LocalSAR[px-2*semiSideLen+1][py+semiSideLen][pz-semiSideLen]*Mass[px-2*semiSideLen+1][py+semiSideLen][pz-semiSideLen];
            corner1Power += LocalSAR[px-2*semiSideLen+1][py-semiSideLen][pz+semiSideLen]*Mass[px-2*semiSideLen+1][py-semiSideLen][pz+semiSideLen];
            corner1Power += LocalSAR[px-2*semiSideLen+1][py+semiSideLen][pz+semiSideLen]*Mass[px-2*semiSideLen+1][py+semiSideLen][pz+semiSideLen];
            corner2Power += LocalSAR[px-2*semiSideLen][py-semiSideLen][pz-semiSideLen]*Mass[px-2*semiSideLen][py-semiSideLen][pz-semiSideLen];
            corner2Power += LocalSAR[px-2*semiSideLen][py+semiSideLen][pz-semiSideLen]*Mass[px-2*semiSideLen][py+semiSideLen][pz-semiSideLen];
            corner2Power += LocalSAR[px-2*semiSideLen][py-semiSideLen][pz+semiSideLen]*Mass[px-2*semiSideLen][py-semiSideLen][pz+semiSideLen];
            corner2Power += LocalSAR[px-2*semiSideLen][py+semiSideLen][pz+semiSideLen]*Mass[px-2*semiSideLen][py+semiSideLen][pz+semiSideLen];
            /* edge power */
            for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
            {
                edge1Power += LocalSAR[px-2*semiSideLen+1][y][pz-semiSideLen]*Mass[px-2*semiSideLen+1][y][pz-semiSideLen];
                edge1Power += LocalSAR[px-2*semiSideLen+1][y][pz+semiSideLen]*Mass[px-2*semiSideLen+1][y][pz+semiSideLen];
                edge2Power += LocalSAR[px-2*semiSideLen][y][pz-semiSideLen]*Mass[px-2*semiSideLen][y][pz-semiSideLen];
                edge2Power += LocalSAR[px-2*semiSideLen][y][pz+semiSideLen]*Mass[px-2*semiSideLen][y][pz+semiSideLen];
            }
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
            {
                edge1Power += LocalSAR[px-2*semiSideLen+1][py-semiSideLen][z]*Mass[px-2*semiSideLen+1][py-semiSideLen][z];
                edge1Power += LocalSAR[px-2*semiSideLen+1][py+semiSideLen][z]*Mass[px-2*semiSideLen+1][py+semiSideLen][z];
                edge2Power += LocalSAR[px-2*semiSideLen][py-semiSideLen][z]*Mass[px-2*semiSideLen][py-semiSideLen][z];
                edge2Power += LocalSAR[px-2*semiSideLen][py+semiSideLen][z]*Mass[px-2*semiSideLen][py+semiSideLen][z];
            }
            for (x = px; x >= px-2*semiSideLen+2; x--)
            {
                edge3Power += LocalSAR[x][py-semiSideLen][pz-semiSideLen]*Mass[x][py-semiSideLen][pz-semiSideLen];
                edge3Power += LocalSAR[x][py+semiSideLen][pz-semiSideLen]*Mass[x][py+semiSideLen][pz-semiSideLen];
                edge3Power += LocalSAR[x][py-semiSideLen][pz+semiSideLen]*Mass[x][py-semiSideLen][pz+semiSideLen];
                edge3Power += LocalSAR[x][py+semiSideLen][pz+semiSideLen]*Mass[x][py+semiSideLen][pz+semiSideLen];
            }
            /* side power */
            for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                {
                    side1Power += LocalSAR[px-2*semiSideLen+1][y][z]*Mass[px-2*semiSideLen+1][y][z];
                    side2Power += LocalSAR[px-2*semiSideLen][y][z]*Mass[px-2*semiSideLen][y][z];
                }
            for (x = px; x >= px-2*semiSideLen+2; x--)
            {
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                {
                    side3Power += LocalSAR[x][y][pz-semiSideLen]*Mass[x][y][pz-semiSideLen];
                    side3Power += LocalSAR[x][y][pz+semiSideLen]*Mass[x][y][pz+semiSideLen];
                }
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                {
                    side3Power += LocalSAR[x][py-semiSideLen][z]*Mass[x][py-semiSideLen][z];
                    side3Power += LocalSAR[x][py+semiSideLen][z]*Mass[x][py+semiSideLen][z];
                }
            }
            break;
        case 2: /* y+ */
            /* core power */
            for (y = py; y <= py+2*semiSideLen-2; y++)
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                    for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                        corePower += LocalSAR[x][y][z]*Mass[x][y][z];
            /* corner power */
            corner1Power += LocalSAR[px-semiSideLen][py+2*semiSideLen-1][pz-semiSideLen]*Mass[px-semiSideLen][py+2*semiSideLen-1][pz-semiSideLen];
            corner1Power += LocalSAR[px+semiSideLen][py+2*semiSideLen-1][pz-semiSideLen]*Mass[px+semiSideLen][py+2*semiSideLen-1][pz-semiSideLen];
            corner1Power += LocalSAR[px-semiSideLen][py+2*semiSideLen-1][pz+semiSideLen]*Mass[px-semiSideLen][py+2*semiSideLen-1][pz+semiSideLen];
            corner1Power += LocalSAR[px+semiSideLen][py+2*semiSideLen-1][pz+semiSideLen]*Mass[px+semiSideLen][py+2*semiSideLen-1][pz+semiSideLen];
            corner2Power += LocalSAR[px-semiSideLen][py+2*semiSideLen][pz-semiSideLen]*Mass[px-semiSideLen][py+2*semiSideLen][pz-semiSideLen];
            corner2Power += LocalSAR[px+semiSideLen][py+2*semiSideLen][pz-semiSideLen]*Mass[px+semiSideLen][py+2*semiSideLen][pz-semiSideLen];
            corner2Power += LocalSAR[px-semiSideLen][py+2*semiSideLen][pz+semiSideLen]*Mass[px-semiSideLen][py+2*semiSideLen][pz+semiSideLen];
            corner2Power += LocalSAR[px+semiSideLen][py+2*semiSideLen][pz+semiSideLen]*Mass[px+semiSideLen][py+2*semiSideLen][pz+semiSideLen];
            /* edge power */
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
            {
                edge1Power += LocalSAR[px-semiSideLen][py+2*semiSideLen-1][z]*Mass[px-semiSideLen][py+2*semiSideLen-1][z];
                edge1Power += LocalSAR[px+semiSideLen][py+2*semiSideLen-1][z]*Mass[px+semiSideLen][py+2*semiSideLen-1][z];
                edge2Power += LocalSAR[px-semiSideLen][py+2*semiSideLen][z]*Mass[px-semiSideLen][py+2*semiSideLen][z];
                edge2Power += LocalSAR[px+semiSideLen][py+2*semiSideLen][z]*Mass[px+semiSideLen][py+2*semiSideLen][z];
            }
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
            {
                edge1Power += LocalSAR[x][py+2*semiSideLen-1][pz-semiSideLen]*Mass[x][py+2*semiSideLen-1][pz-semiSideLen];
                edge1Power += LocalSAR[x][py+2*semiSideLen-1][pz+semiSideLen]*Mass[x][py+2*semiSideLen-1][pz+semiSideLen];
                edge2Power += LocalSAR[x][py+2*semiSideLen][pz-semiSideLen]*Mass[x][py+2*semiSideLen][pz-semiSideLen];
                edge2Power += LocalSAR[x][py+2*semiSideLen][pz+semiSideLen]*Mass[x][py+2*semiSideLen][pz+semiSideLen];
            }
            for (y = py; y <= py+2*semiSideLen-2; y++)
            {
                edge3Power += LocalSAR[px-semiSideLen][y][pz-semiSideLen]*Mass[px-semiSideLen][y][pz-semiSideLen];
                edge3Power += LocalSAR[px+semiSideLen][y][pz-semiSideLen]*Mass[px+semiSideLen][y][pz-semiSideLen];
                edge3Power += LocalSAR[px-semiSideLen][y][pz+semiSideLen]*Mass[px-semiSideLen][y][pz+semiSideLen];
                edge3Power += LocalSAR[px+semiSideLen][y][pz+semiSideLen]*Mass[px+semiSideLen][y][pz+semiSideLen];
            }
            /* side power */
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                {
                    side1Power += LocalSAR[x][py+2*semiSideLen-1][z]*Mass[x][py+2*semiSideLen-1][z];
                    side2Power += LocalSAR[x][py+2*semiSideLen][z]*Mass[x][py+2*semiSideLen][z];
                }
            for (y = py; y <= py+2*semiSideLen-2; y++)
            {
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                {
                    side3Power += LocalSAR[px-semiSideLen][y][z]*Mass[px-semiSideLen][y][z];
                    side3Power += LocalSAR[px+semiSideLen][y][z]*Mass[px+semiSideLen][y][z];
                }
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                {
                    side3Power += LocalSAR[x][y][pz-semiSideLen]*Mass[x][y][pz-semiSideLen];
                    side3Power += LocalSAR[x][y][pz+semiSideLen]*Mass[x][y][pz+semiSideLen];
                }
            }
            break;
        case 3: /* y- */
            /* core power */
            for (y = py; y >= py-2*semiSideLen+2; y--)
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                    for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                        corePower += LocalSAR[x][y][z]*Mass[x][y][x];
            /* corner power */
            corner1Power += LocalSAR[px-semiSideLen][py-2*semiSideLen+1][pz-semiSideLen]*Mass[px-semiSideLen][py-2*semiSideLen+1][pz-semiSideLen];
            corner1Power += LocalSAR[px+semiSideLen][py-2*semiSideLen+1][pz-semiSideLen]*Mass[px+semiSideLen][py-2*semiSideLen+1][pz-semiSideLen];
            corner1Power += LocalSAR[px-semiSideLen][py-2*semiSideLen+1][pz+semiSideLen]*Mass[px-semiSideLen][py-2*semiSideLen+1][pz+semiSideLen];
            corner1Power += LocalSAR[px+semiSideLen][py-2*semiSideLen+1][pz+semiSideLen]*Mass[px+semiSideLen][py-2*semiSideLen+1][pz+semiSideLen];
            corner2Power += LocalSAR[px-semiSideLen][py-2*semiSideLen][pz-semiSideLen]*Mass[px-semiSideLen][py-2*semiSideLen][pz-semiSideLen];
            corner2Power += LocalSAR[px+semiSideLen][py-2*semiSideLen][pz-semiSideLen]*Mass[px+semiSideLen][py-2*semiSideLen][pz-semiSideLen];
            corner2Power += LocalSAR[px-semiSideLen][py-2*semiSideLen][pz+semiSideLen]*Mass[px-semiSideLen][py-2*semiSideLen][pz+semiSideLen];
            corner2Power += LocalSAR[px+semiSideLen][py-2*semiSideLen][pz+semiSideLen]*Mass[px+semiSideLen][py-2*semiSideLen][pz+semiSideLen];
            /* edge power */
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
            {
                edge1Power += LocalSAR[px-semiSideLen][py-2*semiSideLen+1][z]*Mass[px-semiSideLen][py-2*semiSideLen+1][z];
                edge1Power += LocalSAR[px+semiSideLen][py-2*semiSideLen+1][z]*Mass[px+semiSideLen][py-2*semiSideLen+1][z];
                edge2Power += LocalSAR[px-semiSideLen][py-2*semiSideLen][z]*Mass[px-semiSideLen][py-2*semiSideLen][z];
                edge2Power += LocalSAR[px+semiSideLen][py-2*semiSideLen][z]*Mass[px+semiSideLen][py-2*semiSideLen][z];
            }
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
            {
                edge1Power += LocalSAR[x][py-2*semiSideLen+1][pz-semiSideLen]*Mass[x][py-2*semiSideLen+1][pz-semiSideLen];
                edge1Power += LocalSAR[x][py-2*semiSideLen+1][pz+semiSideLen]*Mass[x][py-2*semiSideLen+1][pz+semiSideLen];
                edge2Power += LocalSAR[x][py-2*semiSideLen][pz-semiSideLen]*Mass[x][py-2*semiSideLen][pz-semiSideLen];
                edge2Power += LocalSAR[x][py-2*semiSideLen][pz+semiSideLen]*Mass[x][py-2*semiSideLen][pz+semiSideLen];
            }
            for (y = py; y >= py-2*semiSideLen+2; y--)
            {
                edge3Power += LocalSAR[px-semiSideLen][y][pz-semiSideLen]*Mass[px-semiSideLen][y][pz-semiSideLen];
                edge3Power += LocalSAR[px+semiSideLen][y][pz-semiSideLen]*Mass[px+semiSideLen][y][pz-semiSideLen];
                edge3Power += LocalSAR[px-semiSideLen][y][pz+semiSideLen]*Mass[px-semiSideLen][y][pz+semiSideLen];
                edge3Power += LocalSAR[px+semiSideLen][y][pz+semiSideLen]*Mass[px+semiSideLen][y][pz+semiSideLen];
            }
            /* side power */
            for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                {
                    side1Power += LocalSAR[x][py-2*semiSideLen+1][z]*Mass[x][py-2*semiSideLen+1][z];
                    side2Power += LocalSAR[x][py-2*semiSideLen][z]*Mass[x][py-2*semiSideLen][z];
                }
            for (y = py; y >= py-2*semiSideLen+2; y--)
            {
                for (z = pz-semiSideLen+1; z <= pz+semiSideLen-1; z++)
                {
                    side3Power += LocalSAR[px-semiSideLen][y][z]*Mass[px-semiSideLen][y][z];
                    side3Power += LocalSAR[px+semiSideLen][y][z]*Mass[px+semiSideLen][y][z];
                }
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                {
                    side3Power += LocalSAR[x][y][pz-semiSideLen]*Mass[x][y][pz-semiSideLen];
                    side3Power += LocalSAR[x][y][pz+semiSideLen]*Mass[x][y][pz+semiSideLen];
                }
            }
            break;
        case 4: /* z+ */
            /* core power */
            for (z = pz; z <= pz+2*semiSideLen-2; z++)
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                    for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                        corePower += LocalSAR[x][y][z]*Mass[x][y][z];
            /* corner power */
            corner1Power += LocalSAR[px-semiSideLen][py-semiSideLen][pz+2*semiSideLen-1]*Mass[px-semiSideLen][py-semiSideLen][pz+2*semiSideLen-1];
            corner1Power += LocalSAR[px+semiSideLen][py-semiSideLen][pz+2*semiSideLen-1]*Mass[px+semiSideLen][py-semiSideLen][pz+2*semiSideLen-1];
            corner1Power += LocalSAR[px-semiSideLen][py+semiSideLen][pz+2*semiSideLen-1]*Mass[px-semiSideLen][py+semiSideLen][pz+2*semiSideLen-1];
            corner1Power += LocalSAR[px+semiSideLen][py+semiSideLen][pz+2*semiSideLen-1]*Mass[px+semiSideLen][py+semiSideLen][pz+2*semiSideLen-1];
            corner2Power += LocalSAR[px-semiSideLen][py-semiSideLen][pz+2*semiSideLen]*Mass[px-semiSideLen][py-semiSideLen][pz+2*semiSideLen];
            corner2Power += LocalSAR[px+semiSideLen][py-semiSideLen][pz+2*semiSideLen]*Mass[px+semiSideLen][py-semiSideLen][pz+2*semiSideLen];
            corner2Power += LocalSAR[px-semiSideLen][py+semiSideLen][pz+2*semiSideLen]*Mass[px-semiSideLen][py+semiSideLen][pz+2*semiSideLen];
            corner2Power += LocalSAR[px+semiSideLen][py+semiSideLen][pz+2*semiSideLen]*Mass[px+semiSideLen][py+semiSideLen][pz+2*semiSideLen];
            /* edge power */
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
            {
                edge1Power += LocalSAR[x][py-semiSideLen][pz+2*semiSideLen-1]*Mass[x][py-semiSideLen][pz+2*semiSideLen-1];
                edge1Power += LocalSAR[x][py+semiSideLen][pz+2*semiSideLen-1]*Mass[x][py+semiSideLen][pz+2*semiSideLen-1];
                edge2Power += LocalSAR[x][py-semiSideLen][pz+2*semiSideLen]*Mass[x][py-semiSideLen][pz+2*semiSideLen];
                edge2Power += LocalSAR[x][py+semiSideLen][pz+2*semiSideLen]*Mass[x][py+semiSideLen][pz+2*semiSideLen];
            }
            for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
            {
                edge1Power += LocalSAR[px-semiSideLen][y][pz+2*semiSideLen-1]*Mass[px-semiSideLen][y][pz+2*semiSideLen-1];
                edge1Power += LocalSAR[px+semiSideLen][y][pz+2*semiSideLen-1]*Mass[px+semiSideLen][y][pz+2*semiSideLen-1];
                edge2Power += LocalSAR[px-semiSideLen][y][pz+2*semiSideLen]*Mass[px-semiSideLen][y][pz+2*semiSideLen];
                edge2Power += LocalSAR[px+semiSideLen][y][pz+2*semiSideLen]*Mass[px+semiSideLen][y][pz+2*semiSideLen];
            }
            for (z = pz; z <= pz+2*semiSideLen-2; z++)
            {
                edge3Power += LocalSAR[px-semiSideLen][py-semiSideLen][z]*Mass[px-semiSideLen][py-semiSideLen][z];
                edge3Power += LocalSAR[px+semiSideLen][py-semiSideLen][z]*Mass[px+semiSideLen][py-semiSideLen][z];
                edge3Power += LocalSAR[px-semiSideLen][py+semiSideLen][z]*Mass[px-semiSideLen][py+semiSideLen][z];
                edge3Power += LocalSAR[px+semiSideLen][py+semiSideLen][z]*Mass[px+semiSideLen][py+semiSideLen][z];
            }
            /* side power */
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                {
                    side1Power += LocalSAR[x][y][pz+2*semiSideLen-1]*Mass[x][y][pz+2*semiSideLen-1];
                    side2Power += LocalSAR[x][y][pz+2*semiSideLen]*Mass[x][y][pz+2*semiSideLen];
                }
            for (z = pz; z <= pz+2*semiSideLen-2; z++)
            {
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                {
                    side3Power += LocalSAR[x][py-semiSideLen][z]*Mass[x][py-semiSideLen][z];
                    side3Power += LocalSAR[x][py+semiSideLen][z]*Mass[x][py+semiSideLen][z];
                }
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                {
                    side3Power += LocalSAR[px-semiSideLen][y][z]*Mass[px-semiSideLen][y][z];
                    side3Power += LocalSAR[px+semiSideLen][y][z]*Mass[px+semiSideLen][y][z];
                }
            }
            break;
        case 5: /* z- */
            /* core power */
            for (z = pz; z >= pz-2*semiSideLen+2; z--)
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                    for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                        corePower += LocalSAR[x][y][z]*Mass[x][y][z];
            /* corner power */
            corner1Power += LocalSAR[px-semiSideLen][py-semiSideLen][pz-2*semiSideLen+1]*Mass[px-semiSideLen][py-semiSideLen][pz-2*semiSideLen+1];
            corner1Power += LocalSAR[px+semiSideLen][py-semiSideLen][pz-2*semiSideLen+1]*Mass[px+semiSideLen][py-semiSideLen][pz-2*semiSideLen+1];
            corner1Power += LocalSAR[px-semiSideLen][py+semiSideLen][pz-2*semiSideLen+1]*Mass[px-semiSideLen][py+semiSideLen][pz-2*semiSideLen+1];
            corner1Power += LocalSAR[px+semiSideLen][py+semiSideLen][pz-2*semiSideLen+1]*Mass[px+semiSideLen][py+semiSideLen][pz-2*semiSideLen+1];
            corner2Power += LocalSAR[px-semiSideLen][py-semiSideLen][pz-2*semiSideLen]*Mass[px-semiSideLen][py-semiSideLen][pz-2*semiSideLen];
            corner2Power += LocalSAR[px+semiSideLen][py-semiSideLen][pz-2*semiSideLen]*Mass[px+semiSideLen][py-semiSideLen][pz-2*semiSideLen];
            corner2Power += LocalSAR[px-semiSideLen][py+semiSideLen][pz-2*semiSideLen]*Mass[px-semiSideLen][py+semiSideLen][pz-2*semiSideLen];
            corner2Power += LocalSAR[px+semiSideLen][py+semiSideLen][pz-2*semiSideLen]*Mass[px+semiSideLen][py+semiSideLen][pz-2*semiSideLen];
            /* edge power */
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
            {
                edge1Power += LocalSAR[x][py-semiSideLen][pz-2*semiSideLen+1]*Mass[x][py-semiSideLen][pz-2*semiSideLen+1];
                edge1Power += LocalSAR[x][py+semiSideLen][pz-2*semiSideLen+1]*Mass[x][py+semiSideLen][pz-2*semiSideLen+1];
                edge2Power += LocalSAR[x][py-semiSideLen][pz-2*semiSideLen]*Mass[x][py-semiSideLen][pz-2*semiSideLen];
                edge2Power += LocalSAR[x][py+semiSideLen][pz-2*semiSideLen]*Mass[x][py+semiSideLen][pz-2*semiSideLen];
            }
            for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
            {
                edge1Power += LocalSAR[px-semiSideLen][y][pz-2*semiSideLen+1]*Mass[px-semiSideLen][y][pz-2*semiSideLen+1];
                edge1Power += LocalSAR[px+semiSideLen][y][pz-2*semiSideLen+1]*Mass[px+semiSideLen][y][pz-2*semiSideLen+1];
                edge2Power += LocalSAR[px-semiSideLen][y][pz-2*semiSideLen]*Mass[px-semiSideLen][y][pz-2*semiSideLen];
                edge2Power += LocalSAR[px+semiSideLen][y][pz-2*semiSideLen]*Mass[px+semiSideLen][y][pz-2*semiSideLen];
            }
            for (z = pz; z >= pz-2*semiSideLen+2; z--)
            {
                edge3Power += LocalSAR[px-semiSideLen][py-semiSideLen][z]*Mass[px-semiSideLen][py-semiSideLen][z];
                edge3Power += LocalSAR[px+semiSideLen][py-semiSideLen][z]*Mass[px+semiSideLen][py-semiSideLen][z];
                edge3Power += LocalSAR[px-semiSideLen][py+semiSideLen][z]*Mass[px-semiSideLen][py+semiSideLen][z];
                edge3Power += LocalSAR[px+semiSideLen][py+semiSideLen][z]*Mass[px+semiSideLen][py+semiSideLen][z];
            }
            /* side power */
            for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                {
                    side1Power += LocalSAR[x][y][pz-2*semiSideLen+1]*Mass[x][y][pz-2*semiSideLen+1];
                    side2Power += LocalSAR[x][y][pz-2*semiSideLen]*Mass[x][y][pz-2*semiSideLen];
                }
            for (z = pz; z >= pz-2*semiSideLen+2; z--)
            {
                for (x = px-semiSideLen+1; x <= px+semiSideLen-1; x++)
                {
                    side3Power += LocalSAR[x][py-semiSideLen][z]*Mass[x][py-semiSideLen][z];
                    side3Power += LocalSAR[x][py+semiSideLen][z]*Mass[x][py+semiSideLen][z];
                }
                for (y = py-semiSideLen+1; y <= py+semiSideLen-1; y++)
                {
                    side3Power += LocalSAR[px-semiSideLen][y][z]*Mass[px-semiSideLen][y][z];
                    side3Power += LocalSAR[px+semiSideLen][y][z]*Mass[px+semiSideLen][y][z];
                }
            }
            break;
        default:
            {
                printf("Direction is wrong when compute SAR of the type 3 voxels!\n");
				fflush(stdout);
                exit(2);
            }
    }

    if (layerMarker == 0)
    {
        InnerLayerPartPower[0] = corePower;
        InnerLayerPartPower[1] = side1Power+0.50*side3Power;
        InnerLayerPartPower[2] = 0.50*edge1Power+0.25*edge3Power;
        InnerLayerPartPower[3] = 0.25*corner1Power;

        return (InnerLayerPartPower[0]+InnerLayerPartPower[1]*fraction+InnerLayerPartPower[2]*fraction*fraction+InnerLayerPartPower[3]*fraction*fraction*fraction)/cubeMass;
    }
    else
    {
        OuterLayerPartPower[0] = corePower+side1Power+0.50*side3Power+0.50*edge1Power+0.25*edge3Power+0.25*corner1Power;
        OuterLayerPartPower[1] = side2Power+0.50*side3Power+0.50*edge1Power+0.50*edge2Power+0.5*edge3Power+0.50*corner1Power+0.25*corner2Power;
        OuterLayerPartPower[2] = 0.50*edge2Power+0.25*edge3Power+0.25*corner1Power+0.50*corner2Power;
        OuterLayerPartPower[3] = 0.25*corner2Power;

        return (OuterLayerPartPower[0]+OuterLayerPartPower[1]*fraction+OuterLayerPartPower[2]*fraction*fraction+OuterLayerPartPower[3]*fraction*fraction*fraction)/cubeMass;
    }
}

