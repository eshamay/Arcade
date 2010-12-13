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
