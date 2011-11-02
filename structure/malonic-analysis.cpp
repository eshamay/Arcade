#include "malonic-analysis.h"

namespace malonic {

	BondLengths::BondLengths (Analyzer * t) :
		MalonicAnalysis (t,
				std::string ("Malonic bondlengths"),
				std::string("MalonicBondLengths.dat")) { 
			int numsteps = Analyzer::timesteps;
			lengths[c1o1] = std::vector<double> (numsteps, 0.0);
			lengths[c1oh1] = std::vector<double> (numsteps, 0.0);
			lengths[c2o2] = std::vector<double> (numsteps, 0.0);
			lengths[c2oh2] = std::vector<double> (numsteps, 0.0);
			lengths[h1o2] = std::vector<double> (numsteps, 0.0);
			lengths[h2o1] = std::vector<double> (numsteps, 0.0);
			lengths[h1oh1] = std::vector<double> (numsteps, 0.0);
			lengths[h2oh2] = std::vector<double> (numsteps, 0.0);
		}


	void BondLengths::CalcDistance (AtomPtr atom1, AtomPtr atom2, bond_t bond) {
		lengths[bond].operator[](Analyzer::timestep) = MDSystem::Distance(atom1, atom2).norm();
	}

	void BondLengths::MoleculeCalculation () {
		double distance;
		this->mol->SetAtoms();

		CalcDistance(this->mol->C1(), this->mol->O1(), c1o1);
		CalcDistance(this->mol->C2(), this->mol->O2(), c2o2);
		CalcDistance(this->mol->C1(), this->mol->OH1(), c1oh1);
		CalcDistance(this->mol->C2(), this->mol->OH2(), c2oh2);

		if (this->mol->H1() != (AtomPtr)NULL) {
			CalcDistance(this->mol->H1(), this->mol->OH1(), h1oh1);
			CalcDistance(this->mol->H1(), this->mol->O2(), h1o2);
		}

		if (this->mol->H2() != (AtomPtr)NULL) {
			CalcDistance(this->mol->H2(), this->mol->OH2(), h2oh2);
			CalcDistance(this->mol->H2(), this->mol->O1(), h2o1);
		}
	}

	void BondLengths::OutputDataPoint (bond_t bond, int timestep) {
		fprintf (this->output, " %6.4f", lengths[bond].operator[](timestep));
	}

	void BondLengths::DataOutput () {
		rewind (this->output);

		fprintf (this->output, "c1o1 c2o2 c1oh1 c2oh2 h1oh1 h1o2 h2oh2 h2o1\n");
		for (int i = 0; i < Analyzer::timesteps; i++) {
			OutputDataPoint (c1o1, i);
			OutputDataPoint (c2o2, i);
			OutputDataPoint (c1oh1, i);
			OutputDataPoint (c2oh2, i);
			OutputDataPoint (h1oh1, i);
			OutputDataPoint (h1o2, i);
			OutputDataPoint (h2oh2, i);
			OutputDataPoint (h2o1, i);
			fprintf (this->output, "\n");
		}
	}




		void MalonicTest::MoleculeCalculation () {
			this->mol->SetAtoms();

			AtomPtr h;
			AtomPtr o;
			double distance;

			h = this->mol->H1();
			if (h != NULL) {
				o = this->mol->O2();
				distance = MDSystem::Distance(o,h).norm();
				printf ("distance = %f\n", distance);
			}

			h = this->mol->H2();
			if (h != NULL) {
				o = this->mol->O1();
				distance = MDSystem::Distance(o,h).norm();
				printf ("distance = %f\n", distance);
			}
		}



		/*
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
*/


}
