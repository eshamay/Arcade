#include "h2o-analysis.h"

namespace h2o_analysis {


	void WaterDipoleZComponentAnalysis::MoleculeCalculation () {
		this->mol->SetOrderAxes();
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		if (this->position.second < max && this->position.second >= min) {
			v1 = axis;
			if (!this->position.first)
				v1 = -v1;

			bisector = this->mol->Z();
			cos_sq = acos(bisector < v1) * 180.0 / M_PI;
			

			int pos = int((this->position.second - min) / res);
			histo[pos] += cos_sq;
			++counts[pos];
		}
	}


	void WaterDipoleZComponentAnalysis::DataOutput () {

		rewind(this->output);

		for (int i = 0; i < histo.size(); i++) {
			dep = i*res + min;
			avg = histo[i] / counts[i];
			fprintf (this->output, "% 12.8f % 12.8f\n", dep, avg);
		}
	}

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




	void DistanceAngleAnalysis::MoleculeCalculation () {
		this->mol->SetOrderAxes();
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		v1 = axis;
		if (!this->position.first)
			v1 = -v1;

		angle = acos(this->mol->Z() < v1) * 180.0 / M_PI;

		histo(this->position.second, angle);
	}


}
