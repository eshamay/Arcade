#include "h2o-analysis.h"

namespace h2o_analysis {

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
