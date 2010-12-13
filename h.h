#ifndef H_H_
#define H_H_

#include "molecule.h"

namespace md_system {

	/******************************
	 * Proton (H+)
	 * ****************************/
	class Proton: public Molecule {

		AtomPtr _h;			// pointers to the atoms for easy access

		public:
		Proton ();
		~Proton ();
		Proton (const Molecule& molecule);	// copy constructor for casting from a molecule

		static int numProtons;			// total number of waters in the system

		void SetAtoms ();					// set the _oh bond vector
	};

}
#endif
