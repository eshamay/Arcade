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

		DimerGraph::hbond_data_t nums = this->_graph.NumHBonds(this->mol);
			if (nums.acids != 0)
				std::cout << std::endl << "methyls = " << nums.methyls << " carbos =  " << nums.carbonyls << " mols = " << nums.acids << std::endl;
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

	/*
		 int Dimers::NumberOfHBonds (alkane::Diacid * left, alkane::Diacid * right) {
// check all Hs on one molecule to all Os on the other to find 
// possible H-bonds between the carboxy groups and methyl groups
//

int N = 0;
// first grab all the Hs on the left molecule, and the Os on the right
Atom_ptr_vec Hs;
algorithm_extra::copy_if (left->begin(), left->end(), std::back_inserter(Hs), member_functional::mem_fun_eq (&Atom::Element, Atom::H));

Atom_ptr_vec Os;
algorithm_extra::copy_if (right->begin(), right->end(), std::back_inserter(Os), member_functional::mem_fun_eq (&Atom::Element, Atom::O));

// check each combo of O-H to see if it's an H-bond.
double distance;
for (Atom_it H = Hs.begin(); H != Hs.end(); H++) {
for (Atom_it O = Os.begin(); O != Os.end(); O++) {
distance = MDSystem::Distance(*H,*O).norm();
if (distance < bondgraph::HBONDLENGTH)
++N;
}
}

// now repeat the process with the Hs on the right and the Os on the left molecules
Hs.clear(); Os.clear();
algorithm_extra::copy_if (right->begin(), right->end(), std::back_inserter(Hs), member_functional::mem_fun_eq (&Atom::Element, Atom::H));
algorithm_extra::copy_if (left->begin(), left->end(), std::back_inserter(Os), member_functional::mem_fun_eq (&Atom::Element, Atom::O));
for (Atom_it H = Hs.begin(); H != Hs.end(); H++) {
for (Atom_it O = Os.begin(); O != Os.end(); O++) {
distance = MDSystem::Distance(*H,*O).norm();
if (distance < bondgraph::HBONDLENGTH)
++N;
}
}

return N;
}
using namespace boost;
using namespace md_system;
*/

/*
	 int Dimers::OutputData () {
	 rewind (this->output);
	 fflush(this->output);
	 }
	 */



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
	for (alkane::Diacid::atom_group_list::const_iterator it = this->mol->methyls_begin(); 
			it != this->mol->methyls_end(); it++) {
		angles (it->Bisector(), it->Bond1(), this->position);
	}
}



}
