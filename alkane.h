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

			virtual VecR ReferencePoint () const { 
				return this->CenterOfMass(); 
			}
			// Functions for analysis
			//virtual void SetAtoms () = 0;
	};




	class MalonicAcid : public Alkane {

		public:
			MalonicAcid (Molecule_t moltype);
			virtual ~MalonicAcid ();

			static int numMalonicAcid;
			static int numMalonate;
			static int numDimalonate;

			void SetAtoms ();

			void UnsetAtoms () {
				carbonyl_c.clear();
				carbonyl_o.clear();
				aliphatic_c.clear();
				acid_o.clear();
			}

			Atom_ptr_vec carbonyl_c;
			Atom_ptr_vec carbonyl_o;
			Atom_ptr_vec aliphatic_c;
			Atom_ptr_vec acid_o;

			std::vector< std::vector<double> > distances;
	}; // class malonic




	class SuccinicAcid : public Alkane, public Dihedral {
		protected:
			VecR ch2_1, ch2_2;

		public:
			SuccinicAcid ();
			virtual void SetDihedralAtoms();
			VecR Bisector (AtomPtr left, AtomPtr center, AtomPtr right);

			//void SetMethyleneBisectors ();
			//VecR& CH2_1 () { return ch2_1; }
			//VecR& CH2_2 () { return ch2_2; }

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


	class Diacid : public Alkane {
		public:
			Diacid ();			// a default constructor
			virtual ~Diacid ();
			Diacid (const Molecule& molecule);		// copy constructor for casting from a molecule

			// refer to the carbonyl carbons
			AtomPtr CarbonylCarbon1 () { return this->GetAtom("C1"); }
			AtomPtr CarbonylCarbon2 ();

			// bisector of the O-C-O of the carbonyl group
			VecR CarbonylBisector1 ();
			VecR CarbonylBisector2 ();

			// bonds from the carbonyl carbon to the carbonyl oxygens
			VecR CO1 ();
			VecR CO2 ();

			void LoadAtomGroups () { LoadMethylGroups(); LoadCarbonylGroups(); }

			typedef std::list<ThreeAtomGroup> atom_group_list;
			atom_group_list::iterator methyls_begin () { return methyl_groups.begin(); }
			atom_group_list::iterator methyls_end () { return methyl_groups.end(); }

			atom_group_list::iterator carbonyls_begin () { return carbonyl_groups.begin(); }
			atom_group_list::iterator carbonyls_end () { return carbonyl_groups.end(); }

		protected:
			void LoadMethylGroups ();
			void LoadCarbonylGroups ();

			std::list<ThreeAtomGroup> methyl_groups;
			std::list<ThreeAtomGroup> carbonyl_groups;
			std::list<AtomPtr> hydrogens;

			std::pair<AtomPtr,AtomPtr> FindMethylHydrogens (AtomPtr carbon);
	};


}
#endif
