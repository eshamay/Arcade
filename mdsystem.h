#ifndef MDSYSTEM_H_
#define MDSYSTEM_H_

#include "vecr.h"
#include "atom.h"
#include "molecule.h"
#include "moleculefactory.h"
#include <string>
#include <vector>

namespace md_system {

	const double WANNIER_BOND = 0.7;

	class CoordinateFile {

		public:

			CoordinateFile (const std::string path, int const c_size); 
			CoordinateFile (const std::string path);
			CoordinateFile ();

			virtual ~CoordinateFile () { 
				if (_file != (FILE *)NULL) {
					fclose (_file); 
				}
			}

			char * ReadLine () { 
				fgets (_line, 1000, _file);
				return _line;
			}

			char * Line () { return _line; }

			virtual void LoadNext () = 0;

			virtual void Rewind () {
				rewind(this->_file);
				this->LoadNext();
				_frame = 1;
			}

			// retrieves coordinates as VecR (3-element vectors)
			//const coord_t& Coordinate (const int index) const { return _vectors[index]; }
			//const coord_t& operator() (const int index) const { return _vectors[index]; }
			const double * Coordinate (const int index) const { return &_coords[3*index]; }
			const double * operator() (const int index) const { return &_coords[3*index]; }

			//coord_it begin () const { return _vectors.begin(); }
			//coord_it end () const { return _vectors.end(); }

			// retrieves the coordinate array
			const std::vector<double>& GetArray () const { return _coords; }

			int size () 	const { return _size; }

			bool eof () 	const { return _eof; }
			bool Loaded () const { return !_eof; }
			int Frame () 	const { return _frame; }

		protected:
			FILE				*_file;				// the file listing all the atom coordinates
			std::string _path;

			int					_size;				// number of coordinates to parse in each frame (e.g. number of atoms in the system)

			std::vector<double>								_coords;				// array of atomic coordinates
			//coord_set_t												_vectors;				// set of vectors representing positions

			char _line[1000];
			int 			_frame;		// The current frame (number of timesteps processed)
			bool			_eof;		// end of file marker for the coord file

	};	// Coordinate file




	class MDSystem {

		protected:

			static VecR		_dimensions;		// system dimensions - size

			//! Parses out molecules from the set of atoms in an MD data set. This is typically done via topology files, or some other defined routine that determines connectivity between atoms to form molecules.
			virtual void _ParseMolecules () = 0;
			bool _parse_molecules;	// this gets set if the molecules are to be parsed to determine the specific types.

		public:

			virtual ~MDSystem();

			//! Loads the next frame of an MD simulation data set
			virtual void LoadNext () = 0;
			//! rewinds the coordinate files
			virtual void Rewind () = 0;

			//! The set of all molecules in a system
			virtual Mol_ptr_vec& Molecules () = 0;
			//! An iterator to the beginning of the set of molecules
			virtual Mol_it begin_mols () const = 0;
			//! An iterator to the end of the set of molecules
			virtual Mol_it end_mols () const = 0;
			//! An indexing method for retrieving specific molecules in a system
			virtual MolPtr Molecules (int index) const = 0;
			//! Returns the total number of molecules in a system
			virtual int NumMols () const = 0;

			virtual Atom_ptr_vec& Atoms () = 0;
			virtual Atom_it begin () const = 0;
			virtual Atom_it end () const = 0;
			virtual AtomPtr Atoms (const int index) const = 0;
			virtual AtomPtr operator[] (int index) const = 0;
			virtual int NumAtoms ()	const = 0;

			virtual int size () const = 0;

			static VecR Dimensions () { return MDSystem::_dimensions; }
			static void Dimensions (const vector_base& dimensions) { MDSystem::_dimensions = dimensions; }

			/* Beyond simple system stats, various computations are done routinely in a molecular dynamics system: */

			// Calculate the distance between two points within a system that has periodic boundaries
			static VecR Distance (const VecR v1, const VecR v2);

			// Calculate the distance between two atoms given the periodic boundaries of the system
			static VecR Distance (const AtomPtr atom1, const AtomPtr atom2);

			// Calculates the minimum distance between two molecules - i.e. the shortest inter-molecular atom-pair distance
			static double Distance (const MolPtr mol1, const MolPtr mol2);

			//! Calculates a molecular dipole moment using the "classical E&M" method - consider each atom in the molecule as a point-charge, and that the molecule has no net charge (calculation is independent of origin location). This returns a vector that is the sum of the position*charge (r*q) of each atom.
			static VecR CalcClassicDipole (MolPtr mol);

			//! Calculates the dipole of a molecule using the classic E&M method, but also account for the wannier localization centers (WLC). This assumes that the wannier localization centers have also been calculated and parsed into the molecule.
			static VecR CalcWannierDipole (MolPtr mol);
	};


}	// namespace md system

#endif
