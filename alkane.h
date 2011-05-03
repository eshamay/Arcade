#ifndef ALKANE_H_
#define ALKANE_H_

#include "molecule.h"

namespace alkane {
	using namespace md_system;


	class Alkane : public Molecule {

		public:
			Alkane ();			// a default constructor
			virtual ~Alkane ();
			Alkane (const Molecule& molecule);		// copy constructor for casting from a molecule

			static int numAlkanes;			// total number of carbon chains in the system

			// Functions for analysis
			//virtual void SetAtoms () = 0;
	};

	class Formaldehyde : public Alkane {

		public:
			Formaldehyde ();			// a default constructor
			virtual ~Formaldehyde ();
			Formaldehyde (const MolPtr& molecule);
			Formaldehyde (const Molecule& molecule);		// copy constructor for casting from a molecule
			static int numFormaldehyde;

			void SetAtoms ();
			void SetBonds ();
			VecR ReferencePoint () const { return _c->Position(); }

		protected:
			AtomPtr _c, _o, _h1, _h2;
			VecR _co, _ch1, _ch2;
	};
}
#endif
