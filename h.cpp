#include "h.h"

namespace md_system {
	int Proton::numProtons = 0;

	Proton::Proton () : Molecule()
	{
		this->Rename("h+");
		_moltype = Molecule::H;
		++numProtons;
	}

	Proton::~Proton () {
		--numProtons;
	}

	void Proton::SetAtoms () {
		_h = this->GetAtom("H");
		return;
	}



	int Chlorine::numChlorines = 0;

	Chlorine::Chlorine () : Molecule()
	{
		this->Rename("Cl-");
		_moltype = Molecule::CL;
		++numChlorines;
	}

	Chlorine::~Chlorine () {
		--numChlorines;
	}

	void Chlorine::SetAtoms () {
		_cl = this->GetAtom(Atom::Cl);
		return;
	}
}
