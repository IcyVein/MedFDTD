#ifndef MAX_SIZE_OF_PATH
#define MAX_SIZE_OF_PATH 256
#endif
#ifndef COMPUTE_MASS_AVERAGE_SAR_H
#define COMPUTE_MASS_AVERAGE_SAR_H

//float compute_whole_body_average_sar(float ***local_sar, float ***mass, const int *space_dim, FILE *fp_log);
float compute_whole_body_average_sar(float ***local_sar, unsigned char ***model, double* rho, const int *space_dim, FILE *fp_log);
	
int compute_mass_average_sar(float ***local_SAR, float ***mass, const int *space_dim, const float required_mass, const float air_voxel_rate_threshold, const char *output_path, FILE *fp_log);

int compute_type_1_voxel(float ***local_sar, float ***mass, int ***semi_side_length, float ***mass_average_sar, const float required_mass, const float air_voxel_rate_threshold, const int *space_dim, FILE *fp_log);

int find_empty_side(float ***mass, const int semi_side_len, const int px, const int py, const int pz, const int *space_dim);

float compute_air_voxel_rate(float ***mass, const int semi_side_len, const int px, const int py, const int pz);

int compute_type_1_part_mass(float ***mass, float *part_mass, const int semi_side_len, const int px, const int py, const int pz);

float compute_type_1_sar(float ***local_sar, float ***mass, const int semi_side_len, const float fraction, const int px, const int py, const int pz);

int compute_type_2_voxel(float ***mass, float ***mass_average_sar, int ***semi_side_length, const int *space_dim, FILE *fp_log);

int compute_type_3_voxel(float ***local_sar, float ***mass, float ***mass_average_sar, const float required_mass, const int *space_dim, FILE *fp_log);

float compute_type_3_cube_mass(float ***mass, float *cube_mass, const int semi_side_len, const int px, const int py, const int pz, const int *space_dim);

int compute_type_3_part_mass(float *inner_layer_part_mass, float *outer_layer_part_mass, float ***mass, const int semi_side_len, const int direction, const int px, const int py, const int pz);

float compute_type_3_fraction(const float *inner_layer_part_mass, const float *outer_layer_part_mass, const float required_mass);

float compute_type_3_sar(float ***local_sar, float ***mass, const int direction, const int semi_side_len, const float fraction, const int layer_marker, const float cube_mass, const int px, const int py, const int pz);

float compute_fraction(const float *part_mass, const float required_mass);

#endif

