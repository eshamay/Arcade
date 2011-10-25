#include "diacid-analysis.h"

namespace diacid {

	void CarbonylThetaPhiAnalysis::MoleculeCalculation () {

		// find the center of mass location of the succinic acid
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		// get the bisector
		angles (this->mol->CarbonylBisector1(), this->mol->CO1(), this->position);
		angles (this->mol->CarbonylBisector2(), this->mol->CO2(), this->position);

		// and ref-bond of each carbonyl group
	}

}
