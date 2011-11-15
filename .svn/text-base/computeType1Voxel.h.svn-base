#ifndef COMPUTETYPE1VOXEL_H
#define COMPUTETYPE1VOXEL_H

int computeType1Voxel(double ***LocalSAR, double ***Mass, int ***UsedMarker, int ***SemiSideLength, double ***MassAveragedSAR, const double requiredMass, const int *SpaceDim, FILE *fpLog);

int findEmptySide(double ***Mass, const int semiSideLen, const int px, const int py, const int pz, const int *SpaceDim);

int computeMass(double ***Mass, double *PartMass, const int semiSideLen, const int px, const int py, const int pz);

double computeSAR(double ***LocalSAR, double ***Mass, const double cubeMass, const int semiSideLen, const double fraction, const int px, const int py, const int pz);

#endif

