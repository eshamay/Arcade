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
		: Molecule () {
			++numFormaldehyde;
			this->Rename("formaldehyde");
			_moltype = Molecule::FORMALDEHYDE;
		}

	Formaldehyde::~Formaldehyde () {
		numFormaldehyde--;
	}

	Formaldehyde::Formaldehyde (const Molecule& molecule) 
		: Molecule(molecule) {
			++numFormaldehyde;
			this->Rename("formaldehyde");
			_moltype = Molecule::FORMALDEHYDE;
		}

}	// namespace alkane
