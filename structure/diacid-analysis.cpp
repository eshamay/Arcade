#include "diacid-analysis.h"

namespace diacid {

	void Dimers::MoleculeCalculation () {

		this->mol->LoadAtomGroups();
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		// run the calculation between the current molecule and all the other diacids
		for (Mol_it it = this->analysis_mols.begin(); it != this->analysis_mols.end(); it++) {
			if (*it == this->mol) continue;	// don't check against itself

			// calculate the number of H-bonds formed between the two diacids
		}

	}

	int Dimers::NumberOfHBonds (MolPtr) {
		// check all Hs on one molecule to all Os on the other to find 
		// possible H-bonds between the carboxy groups and methyl groups
		//

	}




	RDF::RDF (Analyzer * t) :
		molecule_analysis::DiacidAnalysis (t,
				std::string ("Malonic RDFs"),
				std::string("")),
		rdf (std::string("MalonicRDF.alcO-H.surface.dat"), 
		//rdf (std::string("MalonicRDF.carbO-H.surface.dat"), 
				0.5, 7.0, 0.05) { }

	void RDF::MoleculeCalculation () {
		this->mol->LoadAtomGroups();
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		if (this->position.second > -12.0) {

			WaterPtr wat;
			AtomPtr o, oh;
			for (alkane::Diacid::atom_group_list::iterator coo = this->mol->carbonyls_begin(); coo != this->mol->carbonyls_end(); coo++) {
				o = coo->Right();	// carbonyl
				oh = coo->Left();	// alcohol

				for (Wat_it it = this->h2os.begin(); it != this->h2os.end(); it++) {
					wat = *it;
					distance = MDSystem::Distance (oh, wat->H1()).norm();
					rdf(distance);
					distance = MDSystem::Distance (oh, wat->H2()).norm();
					rdf(distance);
				}
			}
		}
	}



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
