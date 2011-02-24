#ifndef XYZSYSTEM_H_
#define XYZSYSTEM_H_

#include "mdsystem.h"
#include "xyzfile.h"
#include "wannier.h"
#include "bondgraph.h"
#include <ext/algorithm>
#include <exception>

namespace md_files {

	using namespace md_system;

	class XYZSystem : public MDSystem {

		private:
			XYZFile					_xyzfile;		// Atomlist parsed from an xyz file
			WannierFile 		_wanniers;		// The wannier centers

			Mol_ptr_vec			_mols;

			// This is set to control how often the system molecules will be reparsed.
			// Note that the wanniers centers need to be loaded with every timestep,
			// but the molecules don't necessarily need to
			int _reparse_limit;					
			int _reparse_step;

			/* For debugging (and other useful things?) this will keep a list of all the atoms that have been processed into molecules. Any atoms left over at the end of the parsing routine are not included and ... can potentially cause problems */
			Atom_ptr_vec _unparsed;

			void _ParseMolecules ();		// take the atoms we have and stick them into molecules - general umbrella routine

			//! Parses a simple molecule that is composed of a central atom, and all other atoms are connected to it - i.e. h2o, no3, ccl4, etc
			template <typename T>
				void _ParseSimpleMolecule (const Atom::Element_t central_elmt, const Atom::Element_t outer_elmt, const unsigned numOuter);

			/* more specialized parsing routines that don't fall under "simple molecules" */

			void _ParseNitricAcids ();
			void _ParseProtons ();
			void _ParseWanniers ();


			void _UpdateUnparsedList (Atom_ptr_vec& parsed);	// fixes the list of unparsed atoms
			void _CheckForUnparsedAtoms () const;

		public:
			// constructors
			XYZSystem (const std::string& filepath, const VecR& size, const std::string& wannierpath = "");
			~XYZSystem ();

			bondgraph::BondGraph graph;

			void SetReparseLimit (const int limit) { _reparse_limit = limit; }

			Atom_ptr_vec CovalentBonds (const AtomPtr atom) const { return graph.BondedAtoms(atom, bondgraph::covalent); }
			Atom_ptr_vec BondedAtoms (const AtomPtr atom) const { return graph.BondedAtoms (atom); }

			VecR SystemDipole ();	// calculate the total system dipole and return it

			typedef std::exception xyzsysex;

			struct unaccountedex : public xyzsysex { 
				char const* what() const throw() { return "After parsing of the xyz system molecules an atom(s) was left unaccounted for"; }
			};


			//! Loads the next frame of an MD simulation data set
			virtual void LoadNext ();
			//! rewinds the coordinate files
			virtual void Rewind ();

			//! The set of all molecules in a system
			virtual Mol_ptr_vec& Molecules () { return _mols; }
			//! An iterator to the beginning of the set of molecules
			virtual Mol_it begin_mols () const { return _mols.begin(); }
			//! An iterator to the end of the set of molecules
			virtual Mol_it end_mols () const { return _mols.end(); }
			//! An indexing method for retrieving specific molecules in a system
			virtual MolPtr Molecules (int index) const { return _mols[index]; }
			//! Returns the total number of molecules in a system
			virtual int NumMols () const { return _mols.size(); }

			virtual Atom_ptr_vec& Atoms () { return _xyzfile.Atoms(); }
			virtual Atom_it begin () const { return _xyzfile.begin(); }
			virtual Atom_it end () const { return _xyzfile.end(); }
			virtual AtomPtr Atoms (const int index) const { return _xyzfile[index]; }
			virtual AtomPtr operator[] (int index) const { return _xyzfile[index]; }
			virtual int NumAtoms ()	const { return _xyzfile.size(); }

			virtual int size () const { return _xyzfile.size(); }


			vector_map_it begin_wanniers () const { return _wanniers.begin(); }
			vector_map_it end_wanniers () const { return _wanniers.end(); }


			template <typename U>
				class AtomPtr_In_List : public std::binary_function<AtomPtr,U,bool> {
					public:
						bool operator () (const AtomPtr ap, const U u) const {
							return find(u.begin(), u.end(), ap) != u.end();
						}
				};

	};	// xyz system class



	template <typename T>
		void XYZSystem::_ParseSimpleMolecule (const Atom::Element_t central_elmt, const Atom::Element_t outer_elmt, const unsigned numOuter) {

			for (Atom_it it = _xyzfile.begin(); it != _xyzfile.end(); it++) {

				if ((*it)->Element() != central_elmt) continue;

				// grab all the atoms connected to the central atom
				Atom_ptr_vec outers = graph.BondedAtoms (*it, bondgraph::covalent, outer_elmt);
				if (outers.size() != numOuter) continue;

				int molIndex = _mols.size();

				// save the central atom
				outers.push_back(*it);

				T * newmol = new T();
				_mols.push_back (newmol);

				newmol->MolID (molIndex);

				for (Atom_it jt = outers.begin(); jt != outers.end(); jt++) {
					newmol->AddAtom (*jt);
				}

				_UpdateUnparsedList(outers);
			}

			return;
		}	// parse simple molecule	

}	// namespace md files

#endif
