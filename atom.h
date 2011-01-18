#ifndef ATOM_H_
#define ATOM_H_

#include "vecr.h"
#include "utility.h"
#include <vector>
#include <string>

namespace md_system {

	class Molecule;
	typedef Molecule* MolPtr;

	class Atom {

		public:

			typedef enum {
				NO_ELEMENT = 0,
				H = 1, He = 2,
				B = 5, C = 6, N = 7, O = 8, F = 9, Ne = 10,
				Na = 11, Mg = 12, Al = 13, Si = 14, P = 15, S = 16, Cl = 17, Ar = 18,
				K = 19, Ca = 20, 
				I = 53, 
				Cs = 55, Ba = 56,
				D, DW, SW
			} Element_t;

			// constructors
			Atom (const std::string& name, const double * position, const double * force);
			Atom (const Atom& oldAtom);				// copy constructor for deep copies
			virtual ~Atom ();

			typedef Atom* AtomPtr;
			typedef std::vector<AtomPtr> Atom_ptr_vec;
			typedef Atom_ptr_vec::const_iterator Atom_it;
			typedef Atom_ptr_vec::iterator Atom_it_non_const;
			typedef Atom_ptr_vec::reverse_iterator Atom_rit;

			double operator- (const Atom& input) const;		// operator usage to determine the distance between two atoms
			double operator[] (const coord index) const;	// get the atom's position by coordinate
			bool operator< (const AtomPtr& rhs) const;
			bool operator< (const Atom& rhs) const;

			// Input

			void Name (const std::string& name) { _name = name; }

			void Position (const vector_base& position);
			void Position (const double X, const double Y, const double Z);

			void Position (coord const axis, double const value);

			void Force (const vector_base& force) { _force = force; }
			void Force (const double X, const double Y, const double Z) { _force.Set(X, Y, Z); }
			void Force (coord const axis, double const value) { _force.Set (axis, value); }

			void ID (const int id) { _ID = id; }
			//void Charge (double charge) { _charge = charge; }
			void SetAtomProperties ();
			void Residue (const std::string& residue) { _residue = residue; }

			// for setting the atom's position
			void X (double val) { _position[0] = val; }
			void Y (double val) { _position[1] = val; }
			void Z (double val) { _position[2] = val; }

			void MolID (const int mol) { _molid = mol; }	// sets the ID of the molecule containing this atom
			void ParentMolecule (const MolPtr mol) { _pmolecule = mol;  }	// sets a pointer to the molecule that contains the atom

			void Shift (vector_base& shift)			// shift the atom's position
			{ _position += shift; }

			// Output
			std::string Name () const 	{ return _name; }
			Element_t Element () const { return _element; }
			double Mass () const 	{ return _mass; }
			double Charge () const 	{ return _charge; }
			int ID () const 		{ return _ID; }
			const std::string& Residue () const { return _residue; }

			const vector_map& Position () const	{ return _position; }
			const vector_map& Force () const		{ return _force; }

			double X () const 		{ return _position.x(); }
			double Y () const		{ return _position.y(); }
			double Z () const 		{ return _position.z(); }
			int MolID () const		{ return _molid; }
			MolPtr ParentMolecule () const { return _pmolecule; }

			void Print () const;


			static AtomPtr FindByElement (const Atom_ptr_vec& apv, Element_t elmt);
			static AtomPtr FindByID (const Atom_ptr_vec& apv, int id);
			// predicate to determine if two atoms are equal (by element)
			static bool element_eq (const AtomPtr& first, const AtomPtr& second);
			// predicate compares two atom ids for sorting
			static bool id_cmp (const AtomPtr& first, const AtomPtr& second);
			// compares two atoms for equality (using the atom id)
			static bool id_eq (const AtomPtr& first, const AtomPtr& second);

			// tests if the combination of atoms supplied matches the element pair supplied
			static bool ElementCombo (const AtomPtr& ai, const AtomPtr& aj, const Element_t element_a, const Element_t element_b);
			static void KeepByElement (Atom_ptr_vec& u, const Element_t& elmt);

			// conversion routines
			static std::string Element2String (Element_t);
			static Element_t String2Element (const std::string&);

		protected:
			std::string			_name,				// human-readable identifier
				_residue;			// name of the parent-molecule 

			int    _ID;							// some numerical identifier in case the atom is in an ordered list
			int	   _molid;					// the molecule that contains this atom

			MolPtr _pmolecule;		// parent molecule in which the atom is attached

			double _mass, _charge;

			Element_t _element;			// the actual element based on the atom name - always upper-case and max length of two letters

			vector_map _position;				// Particle position
			vector_map _force;					// the external force on the atom at any given point in time

	};	// class Atom


	typedef Atom::AtomPtr AtomPtr;
	typedef Atom::Atom_ptr_vec Atom_ptr_vec;
	typedef Atom::Atom_it Atom_it;
	typedef Atom::Atom_it_non_const Atom_it_non_const;
	typedef Atom::Atom_rit Atom_rit;

}	// namespace md system
#endif
