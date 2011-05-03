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

			VecR ReferencePoint () const { 
				return this->CenterOfMass(); 
			}
			// Functions for analysis
			//virtual void SetAtoms () = 0;
	
			void ClearAlkaneAtoms () {
				carbonyl_c.clear();
				carbonyl_o.clear();
				aliphatic_c.clear();
				acid_o.clear();
			}

			Atom_ptr_vec carbonyl_c;
			Atom_ptr_vec carbonyl_o;
			Atom_ptr_vec aliphatic_c;
			Atom_ptr_vec acid_o;
	};

	class Formaldehyde : public Alkane {

		public:
			Formaldehyde ();			// a default constructor
			virtual ~Formaldehyde ();
			Formaldehyde (const MolPtr& molecule);
			Formaldehyde (const Molecule& molecule);		// copy constructor for casting from a molecule
			static int numFormaldehyde;

			AtomPtr C () const { return _c; }
			AtomPtr O () const { return _o; }
			AtomPtr H1 () const { return _h1; }
			AtomPtr H2 () const { return _h2; }

			VecR CH1 () const { return _ch1; }
			VecR CH2 () const { return _ch2; }
			VecR CO () const { return _co; }

			void SetAtoms ();
			void SetBonds ();
			VecR ReferencePoint () const { return _c->Position(); }

		protected:
			AtomPtr _c, _o, _h1, _h2;
			VecR _co, _ch1, _ch2;
	};
}
#endif
