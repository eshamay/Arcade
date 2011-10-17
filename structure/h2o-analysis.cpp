#include "h2o-analysis.h"

namespace h2o_analysis {

	void WaterThetaPhiAnalysis::MoleculeCalculation () {
		// find the center of mass location of the succinic acid
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);
		this->mol->SetOrderAxes();

		v1 = axis;// the reference axis - perp to the surface
		if (!(this->position.first))
			v1 = -v1;

		v2 = this->mol->Z();
		v3 = this->mol->OH1();

		theta = acos(v2 < v1) * 180.0 / M_PI;
		//phi = fabs(acos(this->so2->Y() < axis)) * 180.0 / M_PI;
		phi = Dihedral::Angle(v1,v2,v3) * 180.0 / M_PI;
		phi = fabs(phi);
		if (phi > 90.0)
			phi = 180.0 - phi;

		histos (this->position.second, theta, phi);
	}




	void DistanceAngleAnalysis::Analysis () {
		h2os.FindWaterSurfaceLocation();

		this->PreAnalysis();
		for (Wat_it wat = h2os.begin(); wat != h2os.end(); wat++) {
			this->MainWaterCalculation (*wat);
		}
		this->PostAnalysis();
	}



	void BisectorAnalysis::MainWaterCalculation (WaterPtr wat) {
		// find the center of mass location of the succinic acid
		this->distance = this->h2os.TopOrBottom(wat->O()->Position()[WaterSystem::axis]);

		axis = VecR::UnitY();
		if (!this->distance.first)
			axis = -axis;

		this->angle = wat->Bisector() < axis;
		this->angle = acos(this->angle) * 180.0/M_PI;

		this->histo(distance.second, this->angle);
	}

	void TwistAnalysis::MainWaterCalculation (WaterPtr wat) {
		// find the center of mass location of the succinic acid
		this->distance = this->h2os.TopOrBottom(wat->O()->Position()[WaterSystem::axis]);

		axis = VecR::UnitY();
		if (!this->distance.first)
			axis = -axis;

		this->angle = Dihedral::Angle(axis,wat->Bisector(),wat->OH1()) * 180.0 / M_PI;

		this->histo(distance.second, fabs(this->angle));
	}

}
