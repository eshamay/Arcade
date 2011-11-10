#include "molgraphfactory.h"

namespace molgraph {

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
		int Cl_count = AtomCount(atomcount, Atom::Cl);
		int Total_count = C_count + N_count + O_count + S_count + H_count + Cl_count;

		if (Cl_count == Total_count) {
			newmol = new Chlorine ();
		}

		// check for organics/alkanes
		else if (C_count == 3) {
			newmol = new alkane::MalonicAcid (Molecule::MALONIC);
			if (H_count > 4 || O_count != 4) {
				//std::cerr << "----------- Right here --------------" << std::endl;
			}
		}

		// parse out formaldehydes
		else if (C_count == 1 && O_count == 1 && H_count == 2 && Total_count == 4) {
			newmol = new alkane::Formaldehyde ();
		}

		// check for nitrates
		else if (N_count == 1 && O_count == 3 && H_count == 1) {
			newmol = new NitricAcid ();
		}
		else if (N_count == 1 && O_count == 3 && H_count == 0) {
			newmol = new Nitrate ();
		}

		else if (H_count == 1 && Total_count == 1) {
			newmol = new Proton ();
		}
		// check for non-organics with an oxygen
		else if (O_count == 1 && H_count == 1 && Total_count == 2) {
			newmol = new Hydroxide ();
		}
		else if (O_count == 1 && H_count == 2 && Total_count == 3) {
			newmol = new Water ();
		}
		else if (O_count == 1 && H_count == 3 && Total_count == 4) {
			newmol = new Hydronium ();
		}
		else if (O_count == 2 && H_count == 5 && Total_count == 7) {
			newmol = new Zundel ();
		}

		// check for non-organics with an oxygen (i.e. so2)
		else if (S_count == 1 && O_count == 2 && Total_count == 3) {
			newmol = new SulfurDioxide ();
		}

		// if some other glob of Os and Hs
		// we call these blobs
		//else if (Total_count == O_count + H_count) {
			//newmol = new BlobMolecule ();
		//}

		else {
			MolgraphIdentificationError (molgraph);
		}

		// add all the atoms into the new molecule
		for (Atom_it at = atoms.begin(); at != atoms.end(); at++) {
			newmol->AddAtom(*at);
		}

		newmol->SetAtoms();
		newmol->FixAtoms();

		return newmol;
	}




	void MolgraphIdentificationError (MoleculeGraph& molgraph) {
		std::cerr << "MoleculeGraphFactory :: Couldn't figure out the molecule based on the graph of the following atoms:" << std::endl;
		printf ("\n");
		Atom_ptr_vec atoms = molgraph.Atoms();
		for (Atom_it atom = atoms.begin(); atom != atoms.end(); atom++) {
			// print the atom
			printf ("%s(%d)", (*atom)->Name().c_str(), (*atom)->ID());

			// then print all the atoms it's attached to
			Atom_ptr_vec bonded = molgraph.BondedAtoms (*atom);
			double distance;
			for (Atom_it it = bonded.begin(); it != bonded.end(); it++) {
				distance = MDSystem::Distance(*atom, *it).norm();
				printf (" -%.2f-> %s(%d)", distance, (*it)->Name().c_str(), (*it)->ID());
			}
			printf ("\n");
		}
		exit(1);
	}


	// see if the atomcount map has the given atom
	bool HasAtom (atomcounter& atomcount, Atom::Element_t elmt) {
		bool ret = false;
		if (atomcount.find(elmt) != atomcount.end())
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

	// a means to determine the number of atoms in a malonic acid
	// for each C, add 1. For each O add 10, etc.
	// C = 1, O = 10, H = 100

	/*
	void SetMalonicAtoms (MolPtr mol) {
		alkane::MalonicAcid * malonic = static_cast<alkane::MalonicAcid *>(mol);
		malonic->UnsetAtoms();

		bondgraph::BondGraph graph;
		graph.UpdateGraph (malonic->begin(), malonic->end());

		// run through each vertex of the molgraph 
		for (Atom_it atom = malonic->begin(); atom != malonic->end(); atom++) {
			// find the atoms to which it is bound
			Atom_ptr_vec bonded = graph.BondedAtoms (*atom);
			// decide what type of atom it is based on it's own element, and also what is connected to it.
			int num = 0;
			for (Atom_it it = bonded.begin(); it != bonded.end(); it++) {
				if ((*it)->Element() == Atom::H) num += 1;
				else if ((*it)->Element() == Atom::C) num += 10;
				else if ((*it)->Element() == Atom::O) num += 100;
			}
			printf ("%s %d\n", (*atom)->Name().c_str(), num);
			// based on what it's connected to, we can say what type of atom it is.

			// a carbonyl C in malonic acid will be C-(C)-OO, so a value of 2xO + C == 210
			if (num == 210) malonic->carbonyl_c.push_back (*atom);
			// A carbonyl O is only connected to the carbonyl C and would have a value of 10
			else if (num == 10) malonic->carbonyl_o.push_back(*atom);
			// An aliphatic C will have 2xC and 2xH, so 22
			else if (num == 22) malonic->aliphatic_c.push_back(*atom);
			// the acid O will be connected to a C and an H
			else if (num == 11) malonic->acid_o.push_back(*atom);
		}
	}	// set malonic atoms
	*/

}	// namespace molgraph
