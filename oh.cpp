#include "oh.h"

namespace md_system {

	int Hydroxide::numHydroxides = 0;

	Hydroxide::Hydroxide () : Molecule()
	{
		this->Rename("oh");
		this->_moltype = Molecule::OH;
		++numHydroxides;
	}

	Hydroxide::~Hydroxide () {
		--numHydroxides;
	}

	void Hydroxide::SetAtoms () {
		_o = this->GetAtom("O");
		_h = this->GetAtom("H");

		_oh = _h->Position() - _o->Position();

		return;
	}

}
