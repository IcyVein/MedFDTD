#ifndef COMPUTETYPE3VOXEL_H
#define COMPUTETYPE3VOXEL_H

int computeType3Voxel(double ***LocalSAR, double ***Mass, int ***UsedMarker, double ***MassAveragedSAR, const double requiredMass, const int *SpaceDim, FILE *fpLog);

double computeType3CubeMass(double ***Mass, double *CubeMass, const int semiSideLen, const int px, const int py, const int pz, const int *SpaceDim);

double computeType3PartMass(double *InnerLayerPartMass, double *OuterLayerPartMass, double ***Mass, const int semiSideLen, const int direction, const int px, const int py, const int pz);

double computeType3Fraction(double *InnerLayerPartMass, double *OuterLayerPartMass, const double requiredMass);

double computeType3SAR(double ***LocalSAR, double ***Mass, const int direction, const int semiSideLen, const double fraction, const int layerMarker, const double cubeMass, const int px, const int py, const int pz);

#endif

