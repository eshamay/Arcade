#ifndef XYZFILE_H_
#define XYZFILE_H_

#include "mdsystem.h"

namespace md_files {

	using namespace md_system;

	class XYZFile : public md_system::CoordinateFile {


		public:

			XYZFile (std::string path);
			~XYZFile ();

			static void WriteXYZ (Atom_ptr_vec& atoms);

			// Various control functions
			// see LoadFirst for the init arg
			void LoadNext ();
			void Rewind ();


			// output functions
			Atom_ptr_vec& Atoms () { return _atoms; }
			Atom_it begin () const { return _atoms.begin(); }
			Atom_it end () const { return _atoms.end(); }
			AtomPtr front () const { return _atoms.front(); }
			AtomPtr back () const { return _atoms.back(); }

			AtomPtr operator[] (int index) const { return _atoms[index]; }

		protected:

			Atom_ptr_vec  _atoms;		// The listing of the atoms in the file

			bool _initialized;				// To tell wether or not a file has been loaded
	};	 // class xyzfile
} // namespace md_files

#endif
