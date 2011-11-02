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
		// clear everything out
		c1 = c2 = oh1 = oh2 = o1 = o2 = h1 = h2 = (AtomPtr)NULL;
		bool first = true;

		// grab the carbons
		Atom_ptr_vec cs;
		algorithm_extra::copy_if (this->begin(), this->end(), std::back_inserter(cs), member_functional::mem_fun_eq (&Atom::Element, Atom::C));

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
					// if the oh has the h, then fine
					if (distance < 1.2) {
						ht = *h;
						break;
					}
					// otherwise swap the oh & oc
					distance = ((*h)->Position() - oc->Position()).norm();
					if (distance < 1.2) {
						ht = *h;
						AtomPtr t = oh;
						oh = oc;
						oc = t;
						break;
					}
				}

				if (first) {
					c1 = *c;
					oh1 = oh;
					o1 = oc;
					h1 = ht;
					carbonyl_1.SetAtoms(oh, *c, oc);
					first = false;
				}
				else {
					c2 = *c;
					oh2 = oh;
					o2 = oc;
					h2 = ht;
					carbonyl_2.SetAtoms(oh, *c, oc);
					break;
				}
			}
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


}	// namespace alkane
