#define _SAR
//#define _POYNTING
//Project path
#ifndef MAX_SIZE_OF_PATH
#define MAX_SIZE_OF_PATH 256
#endif
#define MAX_NUM_OF_MEDIA 1000
char proj_name[MAX_SIZE_OF_PATH];
char path_log[MAX_SIZE_OF_PATH];

//Load power source path(if the waveForm = -1)
char pathSRC[MAX_SIZE_OF_PATH];/* User define source */

//Save data path
char path_save[MAX_SIZE_OF_PATH];

//Load model path, MKS units
char path_data[MAX_SIZE_OF_PATH];
char model_name[MAX_SIZE_OF_PATH];
char media_name[MAX_SIZE_OF_PATH];

double dx;
double dy;
double dz;

/* 
 * Project Data
 */
int nMax;/* Total number of time steps */
int flag = 0;/* Record the current state */

int need_compute_temperature_rise;
char path_sar[MAX_SIZE_OF_PATH];
double ds, dt_T;
int nMaxT;

//Size of Yee Grids¡¢Thickness of PML¡¢Thickness of padding
int _spaceX;
int _spaceY;
int _spaceZ;
int abcNo;
int thicknessOfPml;
int padding;

//Link boundary
#ifdef _TF_SF
int linkBC_x_1, linkBC_x_2;
int linkBC_y_1, linkBC_y_2;
int linkBC_z_1, linkBC_z_2;
int dis_from_padding = 2 + 1;//distance from padding

void computeLinkBC();
void fixLinkBC_E();
void fixLinkBC_H();
#endif

//Location of feed point£¨1¡ª¡ª_spaceX£©¡¢port£¨x¡¢y¡¢z£©
int _isource[MAX_NUM_OF_MEDIA];
int _jsource[MAX_NUM_OF_MEDIA];
int _ksource[MAX_NUM_OF_MEDIA];
int sourceNum;
char port;

//Save parameters
FILE* fp_log;
int save_plane_amount;
struct save_parameters
{
	int start;//start saving
	int end;//end saving
	int step;//save interval
	int plane_no;//1 : xy, 2 : yz, 3 : xz
	int slice;//location of the saving plane
};

struct field_file
{
	FILE *ex, *ey, *ez, *hx, *hy, *hz;
	save_parameters sp;
};

field_file* fp_save_field_file;
char** save_field_file_name;
/*
 * SAR files
 */
int whole_body_sar = 0;
int nXgSAR = 0;
double* XgSAR;
int save_localSAR_amount;
struct localSAR
{
	int length;
	int plane_no;
	int slice;
	float** RMS_ex;
	float** RMS_ey;
	float** RMS_ez;	
	float** localSARData;
	FILE* fp;
};
localSAR* pSAR;
float*** rho;
float checkData, checkDataPast;
float convergence;
float convergenceCurr, convergenceTarget;
FILE* fp_check;

//Power source parameters
int sourceType;//Source type: 0--Point source, 1--Plane source
int waveForm;//Wave form: -1--User, 0--Sine : freq; 1--Gauss : t0, pulse_width; 2--raised cosine : pulse_width;
			 //			   3--differential Gauss : t0, pulse_width; 4--3 truncated cosine : pulse_width; 5--modulation Gauss : freq, t0, pulse_width.
int sizeSRC;//Size of user define power source
double* power;//User define power source
double amp;//Amplitude
double freq;//Frequency
int sourceT;//Period
int t0;//Start exciation
int pulse_width;//Exciation length
int planeWaveIndex;

/*
 * Boundary parameters
 */
int paddingX_1;
int paddingX_2;
int paddingY_1;
int paddingY_2;
int paddingZ_1;
int paddingZ_2;

//Add on padding and PML
int spaceX;//Add on padding
int spaceY;
int spaceZ;
int Imax;//Add on padding and PML
int Jmax;
int Kmax;
int isource[MAX_NUM_OF_MEDIA];//Add on padding and PML
int jsource[MAX_NUM_OF_MEDIA];
int ksource[MAX_NUM_OF_MEDIA];

/*
 * Select a cube(compute radiation power)
 */
struct radiation_region
{
	int xStart, xEnd;
	int yStart, yEnd;
	int zStart, zEnd;
	int computeX_1, computeX_2;
};
radiation_region rad_region;

/*
 * Dipole antenna
 */
struct _DIPOLE_ANTENNA
{
	int feed_x, feed_y, feed_z;/* Location of Feed Point */
	int direction;/* Toward to x--1, y--2, z--3*/
	int length_low;//length of low oscillator
	int length_high;
	double impedance;/* Characteristic Impedance */
};
int antenna_amount = 0;
_DIPOLE_ANTENNA* _dipole_antenna;

/*
 * MPI
 */
int nprocs, myrank;
int is, ie;//start and end in X-axis
int _global_is, _global_ie;
int* nx_procs;//compute length on X-axis

/*
 * Physical parameters
 */
//Fundamental Constants (MKS units)
double pi = 3.14159265;
double C = 2.99792458E8;
double muO = 4.0 * pi * 1.0E-7;
double epsO = 1.0 / (C * C * muO);

//Specify Material Relative Permittivity and Conductivity
double epsR = 1.0;//free space

/*
 * Specify the CPML Thickness in Each Direction (Value of Zero
 * Corresponds to No PML, and the Grid is Terminated with a PEC)
 * PML thickness in each direction(With PEC).
 */
int nxPML_1, nxPML_2;
int nyPML_1, nyPML_2;
int nzPML_1, nzPML_2;

//time step increment
double dt;

//Specify the CPML Order and Other Parameters:
int m = 3, ma = 1;

double sig_x_max;
double sig_y_max;
double sig_z_max;
double alpha_x_max;
double alpha_y_max;
double alpha_z_max;
double kappa_x_max;
double kappa_y_max;
double kappa_z_max;

/*
 * H & E Field components
 */

float *Hx;
float *Hy;
float *Hz;
float *Ex;
float *Ey;
float *Ez;

#define Ex(i,j,k) *(Ex + (i)*Jmax*(Kmax-1) + (j)*(Kmax-1) +k)
#define Ey(i,j,k) *(Ey + (i)*(Jmax-1)*(Kmax-1) + (j)*(Kmax-1) +k)
#define Ez(i,j,k) *(Ez + (i)*Jmax*Kmax + (j)*Kmax +k)
#define Hx(i,j,k) *(Hx + (i)*(Jmax-1)*Kmax + (j)*Kmax +k)
#define Hy(i,j,k) *(Hy + (i)*Jmax*Kmax + (j)*Kmax +k)
#define Hz(i,j,k) *(Hz + (i)*(Jmax-1)*(Kmax-1) + (j)*(Kmax-1) +k)

//float *maxE;
//#define maxE(i,j,k) *(maxE + (i)*Jmax*Kmax + (j)*Kmax +k)
float *minEz;
#define minE(i,j,k) *(minE + (i)*Jmax*Kmax + (j)*Kmax +k)

typedef vector<double> vectorD1;
typedef vector<vectorD1> vectorD2;
typedef vector<vectorD2> vectorD3;

vectorD3* maxE;
void writeField(FILE*, vectorD3, vector< int >);
void setOutputRange();
void saveMaxE();
vector< int > outputRangeEx(6, 0);
vector< int > outputRangeEy(6, 0);
vector< int > outputRangeEz(6, 0);
vector< int > outputRangeHx(6, 0);
vector< int > outputRangeHy(6, 0);
vector< int > outputRangeHz(6, 0);
vector< int > outputRangeE(6, 0);

// Dispersion
int media_dispersion = 0;
float *Dx;
float *Dy;
float *Dz;
#define Dx(i,j,k) *(Dx + (i)*Jmax*(Kmax-1) + (j)*(Kmax-1) +k)
#define Dy(i,j,k) *(Dy + (i)*(Jmax-1)*(Kmax-1) + (j)*(Kmax-1) +k)
#define Dz(i,j,k) *(Dz + (i)*Jmax*Kmax + (j)*Kmax +k)
float *preDx;
float *preDy;
float *preDz;
#define preDx(i,j,k) *(preDx + (i)*Jmax*(Kmax-1) + (j)*(Kmax-1) +k)
#define preDy(i,j,k) *(preDy + (i)*(Jmax-1)*(Kmax-1) + (j)*(Kmax-1) +k)
#define preDz(i,j,k) *(preDz + (i)*Jmax*Kmax + (j)*Kmax +k)
float *prepreDx;
float *prepreDy;
float *prepreDz;
#define prepreDx(i,j,k) *(prepreDx + (i)*Jmax*(Kmax-1) + (j)*(Kmax-1) +k)
#define prepreDy(i,j,k) *(prepreDy + (i)*(Jmax-1)*(Kmax-1) + (j)*(Kmax-1) +k)
#define prepreDz(i,j,k) *(prepreDz + (i)*Jmax*Kmax + (j)*Kmax +k)
float *preEx;
float *preEy;
float *preEz;
#define preEx(i,j,k) *(preEx + (i)*Jmax*(Kmax-1) + (j)*(Kmax-1) +k)
#define preEy(i,j,k) *(preEy + (i)*(Jmax-1)*(Kmax-1) + (j)*(Kmax-1) +k)
#define preEz(i,j,k) *(preEz + (i)*Jmax*Kmax + (j)*Kmax +k)
float *prepreEx;
float *prepreEy;
float *prepreEz;
#define prepreEx(i,j,k) *(prepreEx + (i)*Jmax*(Kmax-1) + (j)*(Kmax-1) +k)
#define prepreEy(i,j,k) *(prepreEy + (i)*(Jmax-1)*(Kmax-1) + (j)*(Kmax-1) +k)
#define prepreEz(i,j,k) *(prepreEz + (i)*Jmax*Kmax + (j)*Kmax +k)
double CB_PML[MAX_NUM_OF_MEDIA];
double DisDa[MAX_NUM_OF_MEDIA];
double DisDb[MAX_NUM_OF_MEDIA];
double DisDc[MAX_NUM_OF_MEDIA];
double DisEa[MAX_NUM_OF_MEDIA];
double DisEb[MAX_NUM_OF_MEDIA];
double DisEc[MAX_NUM_OF_MEDIA];

void computeDispersion();
void initializeDispersion();
int loadMediaData_Dispersion(char*, int);
void computeFieldH_Dispersion();
void computeFieldH_0_Dispersion();
void computeFieldH_nprocsSub1_Dispersion();
void computePMLH_Dispersion();
void computePMLH_0_Dispersion();
void computePMLH_nprocsSub1_Dispersion();
void computeFieldE_Dispersion();
void computeFieldE_right_Dispersion();
void computeFieldE_0_Dispersion();
void computeFieldE_nprocsSub1_Dispersion();
void computePMLE_Dispersion();
void computePMLE_0_Dispersion();
void computePMLE_nprocsSub1_Dispersion();
int formulaEx_Dispersion(int, int, int);
int formulaEy_Dispersion(int, int, int);
int formulaEz_Dispersion(int, int, int);

/*
 * Model
 */
unsigned char*** modelData;
unsigned char*** modelDataX;
unsigned char*** modelDataY;
unsigned char*** modelDataZ;
int object_num = 0;
int **object_creat;
struct media_data
{
	double sigma;
	double epsilon;
    double epsilon_inf;
    double delay;
	double rho;
	double spec_heat;
    double K;
    double B;
};
media_data media[MAX_NUM_OF_MEDIA];
int mediaNum;
int maxMedia = 2; /* 0---vaccum, 1---PEC */
double CA[MAX_NUM_OF_MEDIA];
double CB[MAX_NUM_OF_MEDIA];

/*
 * CPML components (Taflove 3rd Edition, Chapter 7)
 */
float ***psi_Ezx_1;
float ***psi_Ezx_2;
float ***psi_Hyx_1;
float ***psi_Hyx_2;
float ***psi_Ezy_1;
float ***psi_Ezy_2;
float ***psi_Hxy_1;
float ***psi_Hxy_2;
float ***psi_Hxz_1;
float ***psi_Hxz_2;
float ***psi_Hyz_1;
float ***psi_Hyz_2;
float ***psi_Exz_1;
float ***psi_Exz_2;
float ***psi_Eyz_1;
float ***psi_Eyz_2;
float ***psi_Hzx_1;
float ***psi_Eyx_1;
float ***psi_Hzx_2;
float ***psi_Eyx_2;
float ***psi_Hzy_1;
float ***psi_Exy_1;
float ***psi_Hzy_2;
float ***psi_Exy_2;

double *be_x_1, *ce_x_1, *alphae_x_PML_1, *sige_x_PML_1, *kappae_x_PML_1;
double *bh_x_1, *ch_x_1, *alphah_x_PML_1, *sigh_x_PML_1, *kappah_x_PML_1;
double *be_x_2, *ce_x_2, *alphae_x_PML_2, *sige_x_PML_2, *kappae_x_PML_2;
double *bh_x_2, *ch_x_2, *alphah_x_PML_2, *sigh_x_PML_2, *kappah_x_PML_2;
double *be_y_1, *ce_y_1, *alphae_y_PML_1, *sige_y_PML_1, *kappae_y_PML_1;
double *bh_y_1, *ch_y_1, *alphah_y_PML_1, *sigh_y_PML_1, *kappah_y_PML_1;
double *be_y_2, *ce_y_2, *alphae_y_PML_2, *sige_y_PML_2, *kappae_y_PML_2;
double *bh_y_2, *ch_y_2, *alphah_y_PML_2, *sigh_y_PML_2, *kappah_y_PML_2;
double *be_z_1, *ce_z_1, *alphae_z_PML_1, *sige_z_PML_1, *kappae_z_PML_1;
double *bh_z_1, *ch_z_1, *alphah_z_PML_1, *sigh_z_PML_1, *kappah_z_PML_1;
double *be_z_2, *ce_z_2, *alphae_z_PML_2, *sige_z_PML_2, *kappae_z_PML_2;
double *bh_z_2, *ch_z_2, *alphah_z_PML_2, *sigh_z_PML_2, *kappah_z_PML_2;

//denominators for the update equations
double *den_ex;
double *den_hx;
double *den_ey;
double *den_hy;
double *den_ez;
double *den_hz;

//E&H field update coefficients
double DA;
double DB;

/*
 * Load data parameters
 */
int size[3];
char path_model[MAX_SIZE_OF_PATH];
char path_media[MAX_SIZE_OF_PATH];
long mem_count = 0;

/*
 * Functions about Project
 */
int openProject(FILE*);
int newProject();
/*
 * Functions about initialize
 */
//files and parameters
void initializeFile();//Disk files initialization
void initializePart1();//Memeory initialization
void initializePart2();//Memeory initialization
int setUp();//Coefficients, parameters etc will get computed
void setUpCPML();//CPML coefficient computation
//Arrays
float* initArrayFloat(int);
double* initArray(int);
double*** initArray3(int, int, int ,double);
float*** initArray3Float(int, int, int ,float);
unsigned char*** initArray3Char(int, int, int);
void freeArray3(double***, int, int, int);
void freeArray3Float(float***, int, int, int);
void freeArray3Char(unsigned char***, int, int, int);
void freeFDTDData(void);

/*
 * Functions about compute
 */
void compute();//E & H Field update equation
void computeOneCPU();
void computeOneCPU_Mur2();
void computeOneCPU_Mur2_test();
void computeFieldH();
void computeFieldH_0();
void computeFieldH_nprocsSub1();
void computePMLH();
void computePMLH_0();
void computePMLH_nprocsSub1();
void computeFieldE();
void computeFieldE_right();
void computeFieldE_0();
void computeFieldE_nprocsSub1();
void computePMLE();
void computePMLE_0();
void computePMLE_nprocsSub1();
void computePlane();//E & H Field update equation
/*
 * Functions about model
 */
void buildObject();//Creates the object geometry
void yeeCube (int, int, int, unsigned char);//Sets material properties to a cell
void yeeCubeCap (int, int, int, unsigned char);//Sets material properties to a cellµçÈÝÆ÷
void buildCuboid(int, int, int, int, int, int, unsigned char);
void buildPlane(int, int, int, int, int, int, unsigned char);
void buildLine(int, int, int, int, int, int, unsigned char);
void creatDipoleAntenna(_DIPOLE_ANTENNA);
unsigned char createNewMedia(int, double, double, double, double spec_heat = 0);
unsigned char*** loadData3(char*, int*);//Load model data
int loadMediaData(char*, int);
void loadModel();//Creat model
void loadModelOneCPU();
/*
 * Functions about power source
 */
void powerSource(int);//Generate source
void powerSourcePlaneWaveH(int);//Generate source
double getSource(int, double, double);
double genSource(int);
double* loadSRC(char*);//Load user define power source
/*
 * Functions about save data
 */
void saveData(int);
void writeField_xy(FILE*, FILE*, FILE*, FILE*, FILE*, FILE*, int);//Writes output
void writeField_yz(FILE*, FILE*, FILE*, FILE*, FILE*, FILE*, int);//Writes output
void writeField_xz(FILE*, FILE*, FILE*, FILE*, FILE*, FILE*, int);//Writes output
void writeField_xyOneCPU(FILE*, FILE*, FILE*, FILE*, FILE*, FILE*, int);
void writeField_yzOneCPU(FILE*, FILE*, FILE*, FILE*, FILE*, FILE*, int);
void writeField_xzOneCPU(FILE*, FILE*, FILE*, FILE*, FILE*, FILE*, int);

void saveData_Point(int, int, int, int);//Writes output in appoint point

void int2str(int, char*);

/*
 * Poynting vector
 */
#ifdef _POYNTING
double radiationPower(radiation_region);
double* preEx;
double* preEy;
double* preEz;
#define preEx(i,j,k) *(preEx + (i)*Jmax*(Kmax-1) + (j)*(Kmax-1) +k)
#define preEy(i,j,k) *(preEy + (i)*(Jmax-1)*(Kmax-1) + (j)*(Kmax-1) +k)
#define preEz(i,j,k) *(preEz + (i)*Jmax*Kmax + (j)*Kmax +k)
#endif

/*
 * Mur2
 */
double* Ex_n_1;
double* Ex_n;
double* Ey_n_1;
double* Ey_n;
double* Ez_n_1;
double* Ez_n;
#define Ex_n_1(i,j,k) *(Ex_n_1 + (i)*Jmax*(Kmax-1) + (j)*(Kmax-1) +k)
#define Ey_n_1(i,j,k) *(Ey_n_1 + (i)*(Jmax-1)*(Kmax-1) + (j)*(Kmax-1) +k)
#define Ez_n_1(i,j,k) *(Ez_n_1 + (i)*Jmax*Kmax + (j)*Kmax +k)
#define Ex_n(i,j,k) *(Ex_n + (i)*Jmax*(Kmax-1) + (j)*(Kmax-1) +k)
#define Ey_n(i,j,k) *(Ey_n + (i)*(Jmax-1)*(Kmax-1) + (j)*(Kmax-1) +k)
#define Ez_n(i,j,k) *(Ez_n + (i)*Jmax*Kmax + (j)*Kmax +k)

/*
 * Extension
 */
void addFunctions();

localSAR* initializeLocalSAR(localSAR*, int);
int freeLocalSARRMS(localSAR*, int);
int freeLocalSARData(localSAR*, int);
void computeRMS(localSAR*, int);
float checkRMS();
void resetRMS();
void computeLocalSAR(localSAR, localSAR, float***);
void writeLocalSAR(FILE* , float**);

int computeXgSAR(int);
float*** loadLocalSAR(int*);
float* trans3DTo1D(float***, int*);
float*** trans1DTo3D(float*, int*);

int computeTemperatureRise(char*);
int temperatureRise(char*, double***, unsigned char***, int*);
int loadBackUp(float*, int*, int);
unsigned char*** loadData3T(char*, int*);
double*** loadSar(char*, int*);

/* ERROE Codes */
#define SUCCESS 0
#define MEM_ERROR 1
#define FILE_ERROR 2
#define SETTING_ERROR 3

FILE* fp_mem;

//#define _DEBUG
//#define _DEBUG_L2
//#define _DEBUG_L3
//#define _DEBUG_L4
