#include "molecule-analysis.h"
// test


namespace molecule_analysis {

	void SuccinicAcidAnalysis::Analysis () {
		h2os.FindWaterSurfaceLocation();

		alkane::SuccinicAcid * succ;

		this->PreCalculation ();

		Mol_ptr_vec acids;
		algorithm_extra::copy_if (this->begin_mols(), this->end_mols(), std::back_inserter(acids), member_functional::mem_fun_eq (&Molecule::MolType, Molecule::SUCCINIC));

		//algorithm_extra::copy_if (this->begin_mols(), this->end_mols(), std::back_inserter(acids), member_functional::mem_fun_eq(&Molecule::MolType, Molecule::H2O));

		for (Mol_it it = acids.begin(); it != acids.end(); it++) {
			succ = static_cast<alkane::SuccinicAcid *>(*it);
			this->SuccinicAcidCalculation (succ);
		}

		this->PostCalculation ();

		return;
	}

}
