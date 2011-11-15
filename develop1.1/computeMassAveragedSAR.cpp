#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "computeType1Voxel.h"
#include "computeType2Voxel.h"
#include "computeType3Voxel.h"

int computeMassAveragedSAR(double ***LocalSAR, double ***Rho, const int *SpaceDim, const int *Resolution, const double *requiredMass, char *path, FILE *fpLog)
{
    time_t timep;
    time(&timep);
    printf("------------------------\n");
    printf("%s", ctime(&timep));
    printf("------------------------\n");

    /*
    const int SpaceDim[3] = {222, 242, 238};
    const int Resolution[3] = {1, 1, 1};
    const double requiredMass = 10;
    */

    int x = 0, y = 0, z = 0;

    printf("path = %s\n", path);
    int pathLen = strlen(path);

    /*
    char *logFileName = "computeMassAveragedSAR.log";
    char *inFileLocalSARName = "LocalSARWithoutEar.txt";
    char *inFileRhoName = "RhoWithoutEar.txt";
    */
    char outFileMassAveragedSARName[256];
	sprintf(outFileMassAveragedSARName, "MassAveragedSAR%fKg.txt", *requiredMass);

    /*
    int logFileLen = strlen(logFileName);
    int inFileLocalSARLen = strlen(inFileLocalSARName);
    int inFileRhoLen = strlen(inFileRhoName);
    */
    int outFileMassAveragedSARLen = strlen(outFileMassAveragedSARName);

    /*
    char logFile[pathLen+logFileLen];
    char inFileLocalSAR[pathLen+inFileLocalSARLen];
    char inFileRho[pathLen+inFileRhoLen];
    */

	char* outFileMassAveragedSAR;
	outFileMassAveragedSAR = (char*) calloc (pathLen+outFileMassAveragedSARLen, sizeof(char));

    /*
    strcpy(logFile, path);
    strcpy(inFileLocalSAR, path);
    strcpy(inFileRho, path);
    */
    strcpy(outFileMassAveragedSAR, path);

    /*
    strcat(logFile, logFileName);
    strcat(inFileLocalSAR, inFileLocalSARName);
    strcat(inFileRho, inFileRhoName);
    */
    strcat(outFileMassAveragedSAR, outFileMassAveragedSARName);

    /*
    FILE *fpLog, *fpLocalSAR, *fpRho;
    */
    FILE *fpMassAveragedSAR;
    /*
    printf("<- %s\n", logFile);
    if ((fpLog = fopen(logFile, "a")) == NULL)
    {
        fprintf(stderr, "ERROR! Can't open %s.", logFile);
        exit(1);
    }
    printf("-> %s\n", inFileLocalSAR);
    fprintf(fpLog, "-> %s\n", inFileLocalSAR);
    if ((fpLocalSAR = fopen(inFileLocalSAR, "r")) == NULL)
    {
        fprintf(stderr, "ERROR! Can't open %s.", inFileLocalSAR);
        fprintf(fpLog, "ERROR! Can't open %s.", inFileLocalSAR);
        exit(2);
    }
    printf("-> %s\n", inFileRho);
    fprintf(fpLog, "-> %s\n", inFileRho);
    if ((fpRho = fopen(inFileRho, "r")) == NULL)
    {
        fprintf(stderr, "ERROR! Can't open %s.", inFileRho);
        fprintf(fpLog, "ERROR! Can't open %s.", inFileRho);
        exit(3);
    }
    */
    printf("<- %s\n", outFileMassAveragedSAR);
    fprintf(fpLog, "<- %s\n", outFileMassAveragedSAR);

    if ((fpMassAveragedSAR = fopen(outFileMassAveragedSAR, "w")) == NULL)
    {
        fprintf(stderr, "ERROR! Can't open %s.", outFileMassAveragedSAR);
        fprintf(fpLog, "ERROR! Can't open %s.", outFileMassAveragedSAR);
        exit(4);
    }

    /*
    double ***LocalSAR;
    */
    double ***Mass, ***MassAveragedSAR;
    int ***UsedMarker, ***SemiSideLength;
    /*
    LocalSAR = (double ***)malloc(SpaceDim[0]*sizeof(double **));
    */
    Mass = (double ***)malloc(SpaceDim[0]*sizeof(double **));
    MassAveragedSAR = (double ***)malloc(SpaceDim[0]*sizeof(double **));
    UsedMarker = (int ***)malloc(SpaceDim[0]*sizeof(int **));
    SemiSideLength = (int ***)malloc(SpaceDim[0]*sizeof(int **));
    for(x = 0; x < SpaceDim[0]; x++)
    {
        /*
        LocalSAR[x] = (double **)malloc(SpaceDim[1]*sizeof(double *));
        */
        Mass[x] = (double **)malloc(SpaceDim[1]*sizeof(double *));
        MassAveragedSAR[x] = (double **)malloc(SpaceDim[1]*sizeof(double *));
        UsedMarker[x] = (int **)malloc(SpaceDim[1]*sizeof(int *));
        SemiSideLength[x] = (int **)malloc(SpaceDim[1]*sizeof(int *));
        for(y = 0; y < SpaceDim[1]; y++)
        {
            /*
            LocalSAR[x][y] = (double *)malloc(SpaceDim[2]*sizeof(double));
            */
            Mass[x][y] = (double *)malloc(SpaceDim[2]*sizeof(double));
            MassAveragedSAR[x][y] = (double *)malloc(SpaceDim[2]*sizeof(double));
            UsedMarker[x][y] = (int *)malloc(SpaceDim[2]*sizeof(int));
            SemiSideLength[x][y] = (int *)malloc(SpaceDim[2]*sizeof(int));
            for(z = 0; z < SpaceDim[2]; z++)
            {
                /*
                LocalSAR[x][y][z] = 0;
                */
                Mass[x][y][z] = 0;
                MassAveragedSAR[x][y][z] = 0;
                UsedMarker[x][y][z] = 0;
                SemiSideLength[x][y][z] = 0;
            }
        }
    }

    /*
    printf("\n-------- Load data --------\n");
    fprintf(fpLog, "\n-------- Load data --------\n");
    printf("---- Load LocalSAR ----\n");
    fprintf(fpLog, "---- Load LocalSAR ----\n");
    double localSAR = 0, minLocalSAR = 100, maxLocalSAR = 0;
    int sumLocalSARVoxel = 0;
    while(!feof(fpLocalSAR))
    {
        fscanf(fpLocalSAR, "%d %d %d %lf\n", &x, &y, &z, &localSAR);
        LocalSAR[x-1][y-1][z-1] = localSAR;
        sumLocalSARVoxel += 1;
        if (LocalSAR[x-1][y-1][z-1] != 0 && LocalSAR[x-1][y-1][z-1] > maxLocalSAR)
            maxLocalSAR = LocalSAR[x-1][y-1][z-1];
        if (LocalSAR[x-1][y-1][z-1] != 0 && LocalSAR[x-1][y-1][z-1] < minLocalSAR)
            minLocalSAR = LocalSAR[x-1][y-1][z-1];
    }
    printf("sumLocalSARVoxel = %d, localSAR = [%15.12lf, %15.12lf] W/kg\n", sumLocalSARVoxel, minLocalSAR, maxLocalSAR);
    fprintf(fpLog, "sumLocalSARVoxel = %d, localSAR = [%15.12lf, %15.12lf] W/kg\n", sumLocalSARVoxel, minLocalSAR, maxLocalSAR);

    printf("---- Load density ----\n");
    fprintf(fpLog, "---- Load density ----\n");
    double rho = 0, minRho = 1000, maxRho = 0, minMass = 100, maxMass = 0;
    int sumRhoVoxel = 0;
    while(!feof(fpRho))
    {
        fscanf(fpRho, "%d %d %d %lf\n", &x, &y, &z, &rho);
        sumRhoVoxel += 1;
        rho = rho/1000;
        Mass[x-1][y-1][z-1] = rho*Resolution[0]*Resolution[1]*Resolution[2];
        if (rho > maxRho)
            maxRho = rho;
        if (rho != 0 && rho < minRho)
            minRho = rho;
        if (Mass[x-1][y-1][z-1] != 0 && Mass[x-1][y-1][z-1] > maxMass)
            maxMass = Mass[x-1][y-1][z-1];
        if (Mass[x-1][y-1][z-1] != 0 && Mass[x-1][y-1][z-1] < minMass)
            minMass = Mass[x-1][y-1][z-1];
    }
    printf("sumRhoVoxel = %d\n", sumRhoVoxel);
    printf("density = [%lf, %lf] g/mm3\n", minRho, maxRho);
    printf("voxel mass = [%lf, %lf] g/mm3\n", minMass, maxMass);
    fprintf(fpLog, "sumRhoVoxel = %d\n", sumRhoVoxel);
    fprintf(fpLog, "density = [%lf, %lf] g/mm3\n", minRho, maxRho);
    fprintf(fpLog, "voxel mass = [%lf, %lf] g/mm3\n", minMass, maxMass);
    */

    for (z = 0; z < SpaceDim[2]; z++)
        for (y = 0; y < SpaceDim[1]; y++)
            for (x = 0; x < SpaceDim[0]; x++)
                Mass[x][y][z] = Rho[x][y][z]*Resolution[0]*Resolution[1]*Resolution[2] * 1e-6;/* Need Modify */

    printf("\n-------- Compute type 1 voxels --------\n");
    fprintf(fpLog, "\n-------- Compute type 1 voxels --------\n");
	fflush(stdout);
    computeType1Voxel(LocalSAR, Mass, UsedMarker, SemiSideLength, MassAveragedSAR, *requiredMass, SpaceDim, fpLog);

    int sumType1Voxel = 0, sumType2Voxel = 0, sumType3Voxel = 0;
    for (z = 0; z < SpaceDim[2]; z++)
    {
        for (y = 0; y < SpaceDim[1]; y++)
        {
            for (x = 0; x < SpaceDim[0]; x++)
            {
                if (MassAveragedSAR[x][y][z] != 0)
                {
                    sumType1Voxel += 1;
                }
                if (MassAveragedSAR[x][y][z] == 0 && UsedMarker[x][y][z] != 0)
                {
                    sumType2Voxel += 1;
                }
                if (Mass[x][y][z] > 0 && UsedMarker[x][y][z] == 0)
                {
                    sumType3Voxel += 1;
                }
            }
        }
    }
    printf("sumType1Voxel + sumType2Voxel + sumType3Voxel = %d + %d + %d = %d\n", sumType1Voxel, sumType2Voxel, sumType3Voxel, sumType1Voxel+sumType2Voxel+sumType3Voxel);
    fprintf(fpLog, "sumType1Voxel + sumType2Voxel + sumType3Voxel = %d + %d + %d = %d\n", sumType1Voxel, sumType2Voxel, sumType3Voxel, sumType1Voxel+sumType2Voxel+sumType3Voxel);

    printf("\n-------- Compute type 2 voxels --------\n");
    fprintf(fpLog, "\n-------- Compute type 2 voxels --------\n");
	fflush(stdout);
    computeType2Voxel(MassAveragedSAR, UsedMarker, SemiSideLength, Mass, SpaceDim, fpLog);

    printf("\n-------- Save data --------\n");
    fprintf(fpLog, "\n-------- Save data --------\n");
    double minMassAveragedSAR = 1000, maxMassAveragedSAR = 0;
    int sumMassAveragedSARVoxel = 0;
    for (z = 0; z < SpaceDim[2]; z++)
    {
        for (y = 0; y < SpaceDim[1]; y++)
        {
            for (x = 0; x < SpaceDim[0]; x++)
            {
                if (MassAveragedSAR[x][y][z] != 0)
                {
                    fprintf(fpMassAveragedSAR, "%d %d %d %e\n", x, y, z, MassAveragedSAR[x][y][z]);
                    sumMassAveragedSARVoxel += 1;
                    if (MassAveragedSAR[x][y][z] > maxMassAveragedSAR)
                        maxMassAveragedSAR = MassAveragedSAR[x][y][z];
                    if (MassAveragedSAR[x][y][z] < minMassAveragedSAR)
                        minMassAveragedSAR = MassAveragedSAR[x][y][z];
                }
            }
        }
    }
    printf("sumMassAveragedSARVoxel = %d, MassAveragedSAR = [%15.12lf, %15.12lf] W/kg\n", sumMassAveragedSARVoxel, minMassAveragedSAR, maxMassAveragedSAR);
    fprintf(fpLog, "sumMassAveragedSARVoxel = %d, MassAveragedSAR = [%15.12lf, %15.12lf] W/kg\n", sumMassAveragedSARVoxel, minMassAveragedSAR, maxMassAveragedSAR);

    printf("\n-------- Release memory --------\n");
    fprintf(fpLog, "\n-------- Release memory --------\n");
    for (x = 0; x < SpaceDim[0]; x++)
    {
        for (y = 0; y < SpaceDim[1]; y++)
        {
            /*
            free(LocalSAR[x][y]);
            */
            free(Mass[x][y]);
            free(MassAveragedSAR[x][y]);
            free(UsedMarker[x][y]);
            free(SemiSideLength[x][y]);
        }
    }
    for (x = 0; x < SpaceDim[0]; x++)
    {
        /*
        free(LocalSAR[x]);
        */
        free(Mass[x]);
        free(MassAveragedSAR[x]);
        free(UsedMarker[x]);
        free(SemiSideLength[x]);
    }
    /*
    free(LocalSAR);
    */
    free(Mass);
    free(MassAveragedSAR);
    free(UsedMarker);
    free(SemiSideLength);

    time(&timep);
    printf("------------------------\n");
    printf("%s", ctime(&timep));
    printf("------------------------\n");
    printf("--------------------------------\n");

    fprintf(fpLog, "------------------------\n");
    fprintf(fpLog, "%s", ctime(&timep));
    fprintf(fpLog, "------------------------\n");
    fprintf(fpLog, "--------------------------------\n");

    /*
    fclose(fpLog);
    fclose(fpLocalSAR);
    fclose(fpRho);
    */
    fclose(fpMassAveragedSAR);
	fflush(stdout);/* Need Modify */
    return 0;
}

