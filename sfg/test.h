#include "moritah2o.h"

int main () {

	double o_pos[3] = {-2.0028348429,       -2.6097723874,        4.6540504062};
	double h1_pos[3] = {-1.8139819844,       -1.9832539814,        3.8948771625};
	double h2_pos[3] = {-1.0416101094,       -2.8776873607,        4.8320834121};
	double frc[3];

	AtomPtr o = new Atom ("O", o_pos, frc);
	AtomPtr h1 = new Atom ("H1", h1_pos, frc);
	AtomPtr h2 = new Atom ("H2", h2_pos, frc);

	MolPtr mol = new Molecule();
	mol->AddAtom(o);
	mol->AddAtom(h1);
	mol->AddAtom(h2);

	MoritaH2O * wat = new MoritaH2O(mol);
	wat->SetAtoms();

	wat->Print();
	return 0;
}
