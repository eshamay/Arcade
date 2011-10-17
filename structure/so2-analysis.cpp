#include "so2-analysis.h"

namespace so2_analysis {

	void SO2RDFAnalysis::MoleculeCalculation () {
		this->mol->SetAtoms();
		AtomPtr s = this->mol->S();
		//AtomPtr o1 = this->mol->O1();
		//AtomPtr o2 = this->mol->O2();

		AtomPtr h1, h2;
		for (Mol_it wat = this->begin_wats(); wat != this->end_wats(); wat++) {
			distance = MDSystem::Distance(s, (*wat)->GetAtom("O")).norm();
			rdf(distance);
			/*
			h1 = (*wat)->GetAtom("H1");
			h2 = (*wat)->GetAtom("H2");
			distance = MDSystem::Distance(o1, h1).norm();
			rdf(distance);
			distance = MDSystem::Distance(o1, h2).norm();
			rdf(distance);
			distance = MDSystem::Distance(o2, h1).norm();
			rdf(distance);
			distance = MDSystem::Distance(o2, h2).norm();
			rdf(distance);
			*/
		}
	}

}	// namespace so2 analysis
