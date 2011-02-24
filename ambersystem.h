#pragma once
#ifndef AMBERSYSTEM_H_
#define AMBERSYSTEM_H_

#include "mdsystem.h"
#include "crdfile.h"
//#include "forcefile.h"
#include "topfile.h"
//#include "graph.h"

namespace md_files {

	using namespace md_system;

	class AmberSystem : public md_system::MDSystem {

		protected:
			TOPFile		_topfile;
			CRDFile		_coords;
			bool _periodic;

			void _ParseAtomInformation ();
			void _ParseMolecules ();

			Atom_ptr_vec	_atoms;		// the atoms in the system
			Mol_ptr_vec		_mols;		// the molecules in the system

		public:
			// constructors
			AmberSystem (const std::string& prmtop, const std::string& mdcrd, const bool periodic=true);//, const std::string& mdvel = "");
			virtual ~AmberSystem ();

			// Controller & Calculation methods
			void LoadNext ();	 					// Update the system to the next timestep
			void LoadFirst ();
			void Rewind () { 
				_coords.Rewind(); 
			}

			bool eof () const { return _coords.eof(); }

			// Output
			const VecR&	Dimensions () 		const 	{ return _coords.Dimensions(); }		// returns the system size.

			//! The set of all molecules in a system
			Mol_ptr_vec& Molecules () { return _mols; }
			//! An iterator to the beginning of the set of molecules
			Mol_it begin_mols () const { return _mols.begin(); }
			//! An iterator to the end of the set of molecules
			Mol_it end_mols () const { return _mols.end(); }
			//! An indexing method for retrieving specific molecules in a system
			MolPtr Molecules (int index) const { return _mols[index]; }
			//! Returns the total number of molecules in a system
			int NumMols () const { return _mols.size(); }

			Atom_ptr_vec& Atoms () { return _atoms; }
			Atom_it begin () const { return _atoms.begin(); }
			Atom_it end () const { return _atoms.end(); }
			AtomPtr Atoms (const int index) const { return _atoms[index]; }
			AtomPtr operator[] (int index) const { return _atoms[index]; }
			int NumAtoms ()	const { return (int)_atoms.size(); }

			int size () const { return (int)_atoms.size(); }

			//int 	Current ()		const 	{ return _coords.Current(); }

			//void PrintCRDFile () const;						// to output a frame of the system in .crd format

			//bondgraph::BondGraph graph;

	};

}
#endif
