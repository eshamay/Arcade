#include "molgraphfactory.h"

namespace molgraph {

	typedef std::map<Atom::Element_t, int>	atomcounter;


	// given a molgraph, find out what molecule it is, and then create a new one
	MolPtr MoleculeGraph2Molecule (MoleculeGraph& molgraph) {

		atomcounter atomcount;
		MolPtr newmol;

		Atom_ptr_vec atoms (molgraph.Atoms());
		for (Atom_it at = atoms.begin(); at != atoms.end(); at++) {
			atomcount[(*at)->Element()]++;
		}

		int C_count = AtomCount(atomcount, Atom::C);
		int N_count = AtomCount(atomcount, Atom::N);
		int O_count = AtomCount(atomcount, Atom::O);
		int S_count = AtomCount(atomcount, Atom::S);
		int H_count = AtomCount(atomcount, Atom::H);

		// check for organics/alkanes
		if (C_count == 3 && O_count == 4) {
			newmol = new Alkane ();
		}

		// check for nitrates
		else if (N_count == 1 && O_count == 3 && H_count == 1) {
			newmol = new NitricAcid ();
		}
		else if (N_count == 1 && O_count == 3 && H_count == 0) {
			newmol = new Nitrate ();
		}

		// check for non-organics with an oxygen
		else if (O_count == 1 && H_count == 1) {
			newmol = new Hydroxide ();
		}
		else if (O_count == 1 && H_count == 2) {
			newmol = new Water ();
		}
		else if (O_count == 1 && H_count == 3) {
			newmol = new Hydronium ();
		}

		// check for non-organics with an oxygen
		else if (S_count == 1 && O_count == 2) {
			newmol = new SulfurDioxide ();
		}

		//else if (H_count == 1) {
			//newmol = new Proton();
		//}

		else {
			std::cerr << "MoleculeGraphFactory :: Couldn't figure out the molecule based on the graph given" << std::endl;
			std::for_each (atoms.begin(), atoms.end(), std::mem_fun(&Atom::Print));
			exit(1);
		}
			
		// add all the atoms into the new molecule
		for (Atom_it at = atoms.begin(); at != atoms.end(); at++) {
			newmol.AddAtom(*at);
		}

		return newmol;
	}



	// see if the atomcount map has the given atom
	bool HasAtom (atomcounter& atomcount, Atom::Element_t elmt) {
		ret = false;
		if (atomcount.find(elmt) != std::map::end)
			ret = true;
		return ret;
	}

	int AtomCount (atomcounter& atomcount, Atom::Element_t elmt) {
		int count = 0;
		if (HasAtom(atomcount, elmt)) {
			count = atomcount[elmt];
		}
		return count;
	}

}	// namespace molgraph
