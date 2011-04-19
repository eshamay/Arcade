#ifndef MOLGRAPHFACTORY_H_
#define MOLGRAPHFACTORY_H_
#include "molgraph.h"
#include <map>

namespace molgraph {

	// given a molgraph, find out what molecule it is, and then create a new one
	MolPtr MoleculeGraph2Molecule (MoleculeGraph& molgraph);

};	// namespace mol graph

#endif
