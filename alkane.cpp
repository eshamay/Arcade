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


	MalonicAcid::MalonicAcid (Molecule_t moltype)
		: Alkane () {

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
		/*
		this->UnsetAtoms();

		distances.clear();
		distances.resize(this->size(), std::vector<double> (this->size(), 0.0));

				//		C-C = 1
				//		C-O = 2
				//		C-H = 3
				//		O-H = 4
				
		// generate the connectivity matrix between all the atoms
		// the connectivity matrix is the triangle of a matrix
		// the diagonals specify the total number of bonds
		AtomPtr a1, a2;
		for (int i = 0; i < this->size()-1; i++) {
			for (int j = i+1; j < this->size(); j++) {

				a1 = this->_atoms[i];
				a2 = this->_atoms[j];

				distances[i][j] = MDSystem::Distance(a1, a2);

			} 
		}
		*/
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

}	// namespace alkane
