#include "diacid-analysis.h"

namespace diacid {

	void Test::MoleculeCalculation () {
		this->mol->LoadAtomGroups();

		for (alkane::Diacid::atom_group_list::iterator it = this->mol->methyls_begin(); 
				it != this->mol->methyls_end(); ++it) {
			printf ("\t%f\n", acos(it->Angle())*180./M_PI);
		}
	}

	void CarbonylThetaPhiAnalysis::MoleculeCalculation () {

		// find the center of mass location of the succinic acid
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);
		this->mol->LoadAtomGroups();

		// get the dihedral angle
		angles (this->mol->CarbonylBisector1(), this->mol->CO1(), this->position);
		angles (this->mol->CarbonylBisector2(), this->mol->CO2(), this->position);

	}

	void MethylThetaPhiAnalysis::MoleculeCalculation () {

		// find the center of mass location of the succinic acid
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		// run through each methyl group and grab the angles needed
		this->mol->LoadAtomGroups();
		for (alkane::Diacid::atom_group_list::iterator it = this->mol->methyls_begin(); 
				it != this->mol->methyls_end(); it++) {
			angles (it->Bisector(), it->Bond1(), this->position);
		}
	}



}
