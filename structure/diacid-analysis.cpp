#include "diacid-analysis.h"

namespace diacid {

	void Dimers::PreCalculation () {
		if (!initialized) {
			molecule_analysis::SingleMoleculeAnalysis<alkane::Diacid>::PreCalculation ();

			std::vector<alkane::Diacid *> mols;
			for (Mol_it it = this->analysis_mols.begin(); it != this->analysis_mols.end(); it++)
				mols.push_back(static_cast<alkane::Diacid *>(*it));

			this->_graph.Initialize(mols.begin(), mols.end());

			std::for_each(mols.begin(), mols.end(), std::mem_fun(&alkane::Diacid::LoadAtomGroups));
			initialized = true;
		} 
		
		else {
			this->_graph.RecalculateBonds();

		}
	}

	void Dimers::MoleculeCalculation () {

		/*
			 this->mol->LoadAtomGroups();
			 this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
			 this->position = this->h2os.TopOrBottom(com);

			 alkane::Diacid * dia;
			 int n=0;
		// run the calculation between the current molecule and all the other diacids
		for (Mol_it it = this->analysis_mols.begin(); it != this->analysis_mols.end(); it++) {
		if (*it == this->mol) continue;	// don't check against itself

		dia = static_cast<alkane::Diacid *>(*it);
		//n += NumberOfHBonds(this->mol, dia);

		// calculate the number of H-bonds formed between the two diacids
		}

		fprintf (this->output, "%d\n", n);
		fflush(this->output);
		*/
	}

	void Dimers::PostCalculation () {

		DimerGraph::hbond_data_t nums = this->_graph.NumHBonds();
		fprintf (this->output, "%5d %5d %5d %5d\n", 
				nums.acids, nums.methyls, nums.carbonyls, this->_graph.InterAcidConnections());
		fflush(this->output);
	}





	RDF::RDF (Analyzer * t) :
		molecule_analysis::DiacidAnalysis (t,
				std::string ("Malonic RDFs"),
				std::string("")),
		rdf (std::string("MalonicRDF.alcO-H.surface.dat"), 
				0.5, 7.0, 0.05) { 
			//rdf (std::string("MalonicRDF.carbO-H.surface.dat"), 
		}

	void RDF::MoleculeCalculation () {
		this->mol->LoadAtomGroups();
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		if (this->position.second > -12.0) {

			WaterPtr wat;
			AtomPtr o, oh;
			for (alkane::Diacid::atom_group_list::const_iterator coo = this->mol->carbonyls_begin(); coo != this->mol->carbonyls_end(); coo++) {
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
		this->mol->Print();

		std::cout << "methyl Hs" << std::endl;
		Atom_ptr_vec Hs = this->mol->methyl_hydrogens();
		std::for_each(Hs.begin(), Hs.end(), std::mem_fun(&Atom::Print));

		std::cout << "carbonyl Os" << std::endl;
		Atom_ptr_vec Os = this->mol->carbonyl_oxygens();
		std::for_each(Os.begin(), Os.end(), std::mem_fun(&Atom::Print));

		std::cout << "carbonyl Hs" << std::endl;
		Hs = this->mol->carbonyl_hydrogens();
		std::for_each(Hs.begin(), Hs.end(), std::mem_fun(&Atom::Print));
	}




	void CarboxylicThetaPhiAnalysis::MoleculeCalculation () {

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
		for (alkane::Diacid::atom_group_list::const_iterator it = this->mol->methyls_begin(); 
				it != this->mol->methyls_end(); it++) {
			angles (it->Bisector(), it->Bond1(), this->position);
		}
	}

	void CarbonylThetaThetaAnalysis::MoleculeCalculation () {
		this->mol->LoadAtomGroups();
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		angles (this->mol->CO1(), this->mol->CO2(), this->position);

		return;
	}


} // namespace diacid
