#include "alkane.h"

namespace alkane {
	using namespace md_system;

	int Alkane::numAlkanes = 0;
	int MalonicAcid::numMalonicAcid = 0;
	int MalonicAcid::numMalonate = 0;
	int MalonicAcid::numDimalonate = 0;
	int Formaldehyde::numFormaldehyde = 0;

	Alkane::Alkane ()
		: Molecule () {
			++numAlkanes;
			this->Rename("alkane");
			_moltype = Molecule::ALKANE;
		}

	Alkane::~Alkane () {
		numAlkanes--;
	}

	Alkane::Alkane (const Molecule& molecule) 
		: Molecule(molecule) {
			++numAlkanes;
		}


	MalonicAcid::MalonicAcid (Molecule_t moltype) : Alkane () {

		this->_moltype = moltype;

		switch (moltype) {

			case Molecule::MALONIC : 
				this->Rename("malonic acid");
				++numMalonicAcid;
				break;

			case Molecule::MALONATE : 
				this->Rename("malonate");
				++numMalonate;
				break;

			case Molecule::DIMALONATE : 
				this->Rename("dimalonate");
				++numDimalonate;
				break;

			default:
				std::cerr << "Creating a malonic acid type but not using a valid molecule type" << std::endl;
		}

	}

	MalonicAcid::~MalonicAcid () {
		switch (this->_moltype) {

			case Molecule::MALONIC : 
				--numMalonicAcid;
				break;

			case Molecule::MALONATE : 
				--numMalonate;
				break;

			case Molecule::DIMALONATE : 
				--numDimalonate;
				break;

			default:
				std::cerr << "Something funny when killing a malonic acid" << std::endl;
		}
		return;
	}

	void MalonicAcid::SetAtoms () {
		std::vector<AtomPtr> parsed;
		// clear everything out
		this->c1 = this->c2 = cm = this->oh1 = this->oh2 = this->o1 = this->o2 = this->h1 = this->h2 = hc1 = hc2 = (AtomPtr)NULL;
		bool first = true;

		// grab the carbons
		Atom_ptr_vec cs;
		algorithm_extra::copy_if (this->begin(), this->end(), std::back_inserter(cs), member_functional::mem_fun_eq (&Atom::Element, Atom::C));
		algorithm_extra::copy_if (this->begin(), this->end(), std::back_inserter(parsed), member_functional::mem_fun_eq (&Atom::Element, Atom::C));

		// also grab the oxygens
		Atom_ptr_vec os;
		algorithm_extra::copy_if (this->begin(), this->end(), std::back_inserter(os), member_functional::mem_fun_eq (&Atom::Element, Atom::O));

		Atom_ptr_vec hs;
		algorithm_extra::copy_if (this->begin(), this->end(), std::back_inserter(hs), member_functional::mem_fun_eq (&Atom::Element, Atom::H));

		AtomPtr oh, oc;
		double distance;
		// go through each C-O atom pair and find ones that are bound to each other
		for (Atom_it c = cs.begin(); c != cs.end(); c++) {
			oh = oc = (AtomPtr)NULL;

			for (Atom_it o = os.begin(); o != os.end(); o++) {
				distance = ((*c)->Position() - (*o)->Position()).norm();
				// establish the two oxygens in the carbonyl group
				if (distance < 1.5) {	
					parsed.push_back(*o);
					if (oh == NULL) {
						oh = *o;
						continue;
					}
					else {
						oc = *o;
						break;
					}
				}
			}

			// In the case that the given carbon has 2 oxygens connected to it
			// we set the triatomic carbonyl groups
			// also, find the hydrogen, if any, connected to the alcohol oxygen
			if (oh != NULL && oc != NULL) {

				// find the acid h, if it exists
				AtomPtr ht = (AtomPtr)NULL;
				for (Atom_it h = hs.begin(); h != hs.end(); h++) {
					distance = ((*h)->Position() - oh->Position()).norm();
					//printf ("distance = %f\n", distance);
					// if the oh has the h, then fine
					if (distance < 1.3) {
						parsed.push_back(*h);
						ht = *h;
						break;
					}
					// otherwise swap the oh & oc
					distance = ((*h)->Position() - oc->Position()).norm();
					if (distance < 1.3) {
						parsed.push_back(*h);
						ht = *h;
						AtomPtr t = oh;
						oh = oc;
						oc = t;
						break;
					}
				}

				if (first) {
					this->c1 = *c;
					this->oh1 = oh;
					this->o1 = oc;
					this->h1 = ht;
					this->carbonyl_1.SetAtoms(oh, *c, oc);
					first = false;
				}
				else {
					this->c2 = *c;
					this->oh2 = oh;
					this->o2 = oc;
					this->h2 = ht;
					this->carbonyl_2.SetAtoms(oh, *c, oc);
					break;
				}
			}
		}

		// grab the middle carbon and hydrogens set it
		for (Atom_it c = cs.begin(); c != cs.end(); c++) {
			if (*c != c1 && *c != c2) {
				cm = *c;
				break;
			}
		}

		hc1 = hc2 = (AtomPtr)NULL;
		for (Atom_it hc = hs.begin(); hc != hs.end(); hc++) {
			distance = ((*hc)->Position() - cm->Position()).norm();
			if (distance < 1.5) {
				parsed.push_back(*hc);
				if (hc1 == (AtomPtr)NULL) {
					hc1 = *hc;
				}
				else {
					hc2 = *hc;
					break;
				}
			}
		}


		// now remove from this molecule any atoms that are not in the list that was parsed earlier

		if (parsed.size() != this->size()) {
			std::sort (this->begin(), this->end(), Atom::id_cmp);
			std::sort (parsed.begin(), parsed.end(), Atom::id_cmp);
			std::vector<AtomPtr> difference;
			std::set_difference (this->begin(), this->end(), parsed.begin(), parsed.end(), std::back_inserter(difference));

			//printf("\ncurrent molecule = %d\n", this->size());
			//std::for_each (this->begin(), this->end(), std::mem_fun(&Atom::Print));

			//printf("parsed %zu\n", parsed.size());
			//std::for_each (parsed.begin(), parsed.end(), std::mem_fun(&Atom::Print));

			//printf ("difference = %zu\n", difference.size());
			//std::for_each (difference.begin(), difference.end(), std::mem_fun(&Atom::Print));

			this->_atoms.clear();
			std::copy (parsed.begin(), parsed.end(), std::back_inserter(this->_atoms));
			//printf("\ncurrent molecule = %d\n", this->size());
			//std::for_each (this->begin(), this->end(), std::mem_fun(&Atom::Print));
			//std::cout << std::endl;
		}

	}


	Formaldehyde::Formaldehyde ()
		: Alkane () {
			++numFormaldehyde;
			this->Rename("formaldehyde");
			_moltype = Molecule::FORMALDEHYDE;
		}

	Formaldehyde::~Formaldehyde () {
		numFormaldehyde--;
	}

	Formaldehyde::Formaldehyde (const Molecule& molecule) 
		: Alkane(molecule) {
			++numFormaldehyde;
			this->Rename("formaldehyde");
			_moltype = Molecule::FORMALDEHYDE;
		}

	Formaldehyde::Formaldehyde (const MolPtr& molecule) 
		: Alkane(*molecule) {
			++numFormaldehyde;
			this->Rename("formaldehyde");
			_moltype = Molecule::FORMALDEHYDE;
		}

	void Formaldehyde::SetAtoms () {

		this->_c = this->GetAtom(Atom::C);
		this->_o = this->GetAtom(Atom::O);

		this->_h1 = (AtomPtr)NULL; this->_h2 = (AtomPtr)NULL;
		for (Atom_it it = this->_atoms.begin(); it != this->_atoms.end(); it++) {
			if ((*it)->Element() == Atom::H) {
				if (_h1 == (AtomPtr)NULL)
					this->_h1 = *it;
				else
					this->_h2 = *it;
			}
		}

		this->SetBonds();
	}

	void Formaldehyde::SetBonds () {
		this->_ch1 = this->_h1->Position() - this->_c->Position();
		this->_ch2 = this->_h2->Position() - this->_c->Position();
		this->_co = this->_o->Position() - this->_c->Position();
	}

	SuccinicAcid::SuccinicAcid ()
	{
		this->Rename("sin");
		_moltype = Molecule::SUCCINIC;
	}

	void SuccinicAcid::SetDihedralAtoms () {
		this->dihedral_atoms[0] = this->GetAtom("C1");
		this->dihedral_atoms[1] = this->GetAtom("C2");
		this->dihedral_atoms[2] = this->GetAtom("C3");
		this->dihedral_atoms[3] = this->GetAtom("C4");
	}


	/*
		 void SuccinicAcid::SetMethyleneBisectors () {
// get the atoms of the two methylene/CH2 groups
AtomPtr c1 = this->GetAtom("C2");
AtomPtr h1_1 = this->GetAtom("H3");
AtomPtr h1_2 = this->GetAtom("H4");

ch2_1 = Molecule::Bisector(h1_1, c1, h1_2);

AtomPtr c2 = this->GetAtom("C3");
AtomPtr h2_1 = this->GetAtom("H5");
AtomPtr h2_2 = this->GetAtom("H6");

ch2_2 = Molecule::Bisector(h2_1, c2, h2_2);
}
*/




Diacid::Diacid () : Alkane () { this->Rename("Diacid"); _moltype = Molecule::DIACID; } 
Diacid::~Diacid () { } 
Diacid::Diacid (const Molecule& molecule) : Alkane(molecule) { }

AtomPtr Diacid::CarbonylCarbon2 () {
	// second carbonyl carbon is always the last carbon (e.g. "C5" if there are 5 carbons)
	// filter out only carbons, and then grab the last one
	std::list<AtomPtr> carbons;
	algorithm_extra::copy_if (this->begin(), this->end(), std::back_inserter(carbons), member_functional::mem_fun_eq(&Atom::Element, Atom::C));

	//std::for_each (carbons.begin(), carbons.end(), std::mem_fun(&Atom::Print));
	return carbons.back();
}

VecR Diacid::CarbonylBisector1 () { return carbonyl_groups.front().Bisector(); }
VecR Diacid::CarbonylBisector2 () { return carbonyl_groups.back().Bisector(); }

VecR Diacid::CO1 () { return VecR (this->GetAtom("O1")->Position() - this->CarbonylCarbon1()->Position()); }
VecR Diacid::CO2 () { return VecR (this->GetAtom("O2")->Position() - this->CarbonylCarbon2()->Position()); }

void Diacid::LoadCarbonylGroups () {
	carbonyl_groups.clear();
	carbonyl_groups.push_back (ThreeAtomGroup(this->GetAtom("O3"), this->CarbonylCarbon1(), this->GetAtom("O1")));
	carbonyl_groups.push_back (ThreeAtomGroup(this->GetAtom("O4"), this->CarbonylCarbon2(), this->GetAtom("O2")));
}

// loads the carbon-hydrogen map with the methyl group atoms
void Diacid::LoadMethylGroups () {
	// grab all the hydrogens
	hydrogens.clear();
	algorithm_extra::copy_if (this->begin(), this->end(), std::back_inserter(hydrogens), member_functional::mem_fun_eq (&Atom::Element, Atom::H));

	// now load all the methyl carbons (i.e. all but the carbonyl carbons)
	methyl_groups.clear();
	for (Atom_it it = this->begin(); it != this->end(); ++it) {
		if ((*it)->Element() != Atom::C || *it == CarbonylCarbon1() || *it == CarbonylCarbon2()) continue;
		// for each carbon, find the 2 closest hydrogens and store them in the hydrogen bond-map 
		std::pair<AtomPtr,AtomPtr> hyds = FindMethylHydrogens(*it);
		methyl_groups.push_back(ThreeAtomGroup(hyds.first, *it, hyds.second));
	}
}

std::pair<AtomPtr,AtomPtr> Diacid::FindMethylHydrogens (AtomPtr carbon) {
	// run through each hydrogen in the molecule to find the one it's bound to
	AtomPtr first, second;
	bool one = false;

	for (std::list<AtomPtr>::iterator it = hydrogens.begin(); it != hydrogens.end(); it++) {
		if ((carbon->Position() - (*it)->Position()).norm() < 1.6) {
			if (!one) {
				first = *it;
				one = true;
			}
			else {
				second = *it;
				break;
			}
		}
	}
	return std::make_pair(first, second);
}


Atom_ptr_vec Diacid::methyl_hydrogens () const {
	Atom_ptr_vec Hs;
	for (std::list<ThreeAtomGroup>::const_iterator it = this->methyl_groups.begin(); it != this->methyl_groups.end(); it++) {
		Hs.push_back(it->Left());
		Hs.push_back(it->Right());
	}

	return Hs;
}


Atom_ptr_vec Diacid::carbonyl_oxygens () const {
	Atom_ptr_vec Os;
	for (std::list<ThreeAtomGroup>::const_iterator it = this->carbonyl_groups.begin(); it != this->carbonyl_groups.end(); it++) {
		Os.push_back(it->Left());
		Os.push_back(it->Right());
	}

	return Os;
}

Atom_ptr_vec Diacid::carbonyl_hydrogens () const {
	Atom_ptr_vec Os = this->carbonyl_oxygens();
	Atom_ptr_vec methyl_Hs = this->methyl_hydrogens();
	std::sort(methyl_Hs.begin(), methyl_Hs.end());

	Atom_ptr_vec all_Hs;
	algorithm_extra::copy_if (this->begin(), this->end(), std::back_inserter(all_Hs), member_functional::mem_fun_eq (&Atom::Element, Atom::H));
	std::sort(all_Hs.begin(), all_Hs.end());
	
	Atom_ptr_vec carbonyl_Hs;
	std::set_difference(all_Hs.begin(), all_Hs.end(), methyl_Hs.begin(), methyl_Hs.end(), std::back_inserter(carbonyl_Hs));

	return carbonyl_Hs;
}

std::pair<double,double> Diacid::MalonicDihedralAngle (Diacid *acid) {

		// first vector is the O=>C bond
		VecR v1 = acid->CO1();
		// but it points from the O to the C, not C->O
		v1 = -v1;
		// 2nd vector is C1->C2
		VecR v2 = acid->GetAtom("C2")->Position() - acid->GetAtom("C1")->Position();
		// 3rd is C2->C3
		VecR v3 = acid->GetAtom("C3")->Position() - acid->GetAtom("C2")->Position();
		// the dihedral is calculated from these 3 vectors
		double psi1 = Dihedral::Angle(v1,v2,v3) * 180.0/M_PI;

		
		// repeat for the 2nd side of the molecule, but the atom chain is now O2=>C3->C2->C1
		v1 = acid->CO2();
		v1 = -v1;
		v2 = acid->GetAtom("C2")->Position() - acid->GetAtom("C3")->Position();
		// 3rd is C2->C3
		v3 = acid->GetAtom("C1")->Position() - acid->GetAtom("C2")->Position();
		// the dihedral is calculated from these 3 vectors
		double psi2 = Dihedral::Angle(v1,v2,v3) * 180.0/M_PI;

		return std::make_pair(psi1,psi2);
}

}	// namespace alkane
