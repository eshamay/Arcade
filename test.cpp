#include "test.h"

using namespace md_system;

int main () {

	md_analysis::Analyzer<AmberSystem> sys;
	
	sys.LoadAll();

	int num = 10;
	printf ("%d\n\n", num*4);
	for (int i = 30; i < 30+num; i++) {

		Water * wat = new Water(sys.sys_mols[i]);
		wat->SetAtoms();
		wat->SetOrderAxes ();
		MatR dcm = wat->DCMToLab();

		double angles[3];
		//coordinate_conversion::DCM2EulerAngles_ZXZ (&dcm(0,0), angles);
		coordinate_conversion::DCM2EulerAngles_ZXZ (&dcm(0,0), angles);
		//double conv = 180./M_PI;
		//printf ("% .3f % .3f % .3f\n", conv*angles[0], conv*angles[1], conv*angles[2]);

		VecR rot = Vector3d::UnitZ();

		MatR rotation;
		rotation = Eigen::AngleAxisd(angles[0], Vector3d::UnitZ())
			* Eigen::AngleAxisd(angles[1], Vector3d::UnitX())
			* Eigen::AngleAxisd(angles[2], Vector3d::UnitY());

		VecR o = wat->O()->Position();
		VecR h1 = wat->H1()->Position();
		VecR h2 = wat->H2()->Position();
		rot = (rotation * rot).normalized() + o;

		printf ("O  % .5f % .5f % .5f\n", o[x], o[y], o[z]);
		printf ("H  % .5f % .5f % .5f\n", h1[x], h1[y], h1[z]);
		printf ("H  % .5f % .5f % .5f\n", h2[x], h2[y], h2[z]);
		printf ("X  % .5f % .5f % .5f\n", rot[x], rot[y], rot[z]);

		delete wat;
	}

	/*
	printf ("2\n\n");
	printf ("X  0.0 0.0 0.0\n");
	VecR z = VecR::UnitZ();

	MatR rotation;
	rotation = Eigen::AngleAxisd(M_PI/2., Vector3d::UnitY())
		* Eigen::AngleAxisd(-M_PI/2., Vector3d::UnitZ())
		* Eigen::AngleAxisd(0., Vector3d::UnitY());
	
	z = rotation * z;
	printf ("O  %f %f %f\n", z(0), z(1), z(2));
	*/
	
	return 0;
}
