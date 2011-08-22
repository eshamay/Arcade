#include "moleculefactory.h"

namespace md_system {

	MolPtr MoleculeFactory (const std::string name) {

		MolPtr mol;

		if (name == "h2o" || name == "wat" || name == "SW" || name == "SM2")
			mol = new Water;
		else if (name == "oh")
			mol = new Hydroxide;
		else if (name == "h3o")
			mol = new Hydronium;
		else if (name == "h")
			mol = new Proton;
		else if (name == "so2" || name == "sog" || name == "soa" || name == "soq")
			mol = new SulfurDioxide;
		else if (name == "sin")
			mol = new alkane::SuccinicAcid;

		else {
			//mol = new Molecule;
			std::cerr << "Couldn't determine the molecule using the given name: " << name << std::endl;
			exit(1);
		}

		return mol;
	}	// MoleculeFactory

}	// namespace md system
