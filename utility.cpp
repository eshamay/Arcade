#include "utility.h"


void coordinate_conversion::CartesianToSpherical (double *coord) {
	double rho = sqrt(coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2]); 
	double theta = acos (coord[2]/rho);
	//double phi = atan (coord[1]/coord[0]);
	double phi = atan2 (coord[1], coord[0]);	// atan2 returns the range [-pi, pi]

	coord[0] = rho;
	coord[1] = theta;
	coord[2] = phi;
	return;
}

// assumes an incoming double[9] for the dcm matrix in column-major format
void coordinate_conversion::DCM2EulerAngles_ZXZ (double * dcm, double * angles) {
	angles[0] = atan2(dcm[6], -dcm[7]);
	//angles[0] = atan(-dcm[6]/dcm[7]);
	angles[1] = acos(dcm[8]);
	angles[2] = atan2(dcm[2], dcm[5]);
	//angles[2] = atan(dcm[2]/dcm[5]);
	return;
}

// assumes an incoming double[9] for the dcm matrix in column-major format
void coordinate_conversion::DCM2EulerAngles_YZY (double * dcm, double * angles) {
	angles[0] = atan2(dcm[5],-dcm[3]);
	angles[1] = acos(dcm[4]);
	angles[2] = atan2(dcm[7],dcm[1]);
	return;
}

// assumes an incoming double[9] for the dcm matrix in column-major format
void coordinate_conversion::DCM2EulerAngles_ZXY (double * dcm, double * angles) {
	angles[0] = atan2(-dcm[3],dcm[4]);
	angles[1] = asin(dcm[5]);
	angles[2] = atan2(-dcm[2],dcm[8]);
	return;
}

// assumes an incoming double[9] for the dcm matrix in column-major format
void coordinate_conversion::DCM2EulerAngles_YZX (double * dcm, double * angles) {
	angles[0] = atan2(-dcm[2],dcm[0]);
	angles[1] = asin(dcm[1]);
	angles[2] = atan2(-dcm[7],dcm[4]);
	return;
}
