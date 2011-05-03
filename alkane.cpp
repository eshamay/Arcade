#include "alkane.h"

namespace alkane {
	using namespace md_system;

	int Alkane::numAlkanes = 0;
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

}	// namespace alkane
