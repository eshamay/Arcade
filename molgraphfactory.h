#ifndef MOLGRAPHFACTORY_H_
#define MOLGRAPHFACTORY_H_
#include "molgraph.h"
#include "bondgraph.h"
#include <map>

namespace molgraph {

	typedef std::map<Atom::Element_t, int>	atomcounter;

	// given a molgraph, find out what molecule it is, and then create a new one
	MolPtr MoleculeGraph2Molecule (MoleculeGraph& molgraph);

	// error given when a molgraph can't be figured out
	void MolgraphIdentificationError (MoleculeGraph& molgraph);

	// predicate for checking if a given element is contained in the molgraph
	bool HasAtom (atomcounter& atomcount, Atom::Element_t elmt);

	// count the number of atoms of a given element
	int AtomCount (atomcounter& atomcount, Atom::Element_t elmt);

	// Set the atoms of a malonic acid molecule
	void SetMalonicAtoms (MolPtr);

};	// namespace mol graph

#endif
