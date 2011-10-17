#include "malonic-analysis.h"

namespace malonic_analysis {


	void MalonicTest::Analysis () {
		this->LoadAll();

		int numalk = alkane::Alkane::numAlkanes; 
		int nummal = alkane::MalonicAcid::numMalonicAcid;
		int numate = alkane::MalonicAcid::numMalonate;
		int numdiate = alkane::MalonicAcid::numDimalonate;

		//std::cout << "num alkanes: " << numalk << std::endl;
		//std::cout << "num malonic: " << nummal << std::endl;
		//std::cout << "num malonate: " << numate << std::endl;
		//std::cout << "num dimalonate: " << numdiate << std::endl;

		//if (numate) {
			//alkane::MalonicAcid * mal = static_cast<alkane::MalonicAcid *>( Molecule::FindByType (this->begin_mols(), this->end_mols(), Molecule::MALONATE));

			//mal->Print();

		//}
	}

	// for each timestep we're just finding the distances between the given pairs of atoms and outputing them all
	// to a data file.
	void MalonicBondLengthAnalysis::Analysis () {
		this->LoadAll();

		for (std::vector< std::pair<int,int> >::iterator it = atom_ids.begin(); it != atom_ids.end(); it++) {
			fprintf (this->output, "%13.5f ", BondLength(it->first,it->second));
		}
		fprintf (this->output, "\n");
		fflush(this->output);
	}	

	double MalonicBondLengthAnalysis::BondLength (const int id1, const int id2) {

		atom1 = this->Atoms()[id1];
		atom2 = this->Atoms()[id2];

		distance = MDSystem::Distance(atom1,atom2).Magnitude();

		return distance;
	}


	void MethyleneTilt::Analysis () {
		this->LoadAll();

		for (Mol_it mol = this->begin_mols(); mol != this->end_mols(); mol++) {
			if ((*mol)->MolType() != Molecule::MALONIC 
					&& (*mol)->MolType() != Molecule::MALONATE 
					&& (*mol)->MolType() != Molecule::DIMALONATE) continue;

			mal = static_cast<alkane::MalonicAcid *>(*mol);

			// methylene carbon is index 3
			// methylene hydrogen IDs are 7 and 10
			// However, that corresponds to atom #:
			// carbon = 3, hyd1 = 6, hyd2 = 7
			bisector = Molecule::Bisector(mal->GetAtom(6), mal->GetAtom(3), mal->GetAtom(7));
			angle = bisector < VecR::UnitZ();
			angle = acos(angle)*180.0/M_PI;

			histo (angle);

		}

	}


}
