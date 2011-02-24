#include "xyzfile.h"

namespace md_files {

	XYZFile::XYZFile (std::string path) 
		: 
			md_system::CoordinateFile (path),
			_initialized(false) {

				_file = fopen (path.c_str(), "r");
				if (_file == (FILE *)NULL) {
					printf ("Couldn't load the XYZ coordinate file file:: %s\n", path.c_str());
					exit(1);
				}

				// Initialize the atoms
				//this->LoadNext();
			}

	XYZFile::~XYZFile () {
		for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++)
			delete *it;
	}

	void XYZFile::LoadNext () {

		// first process the frame's header
		fscanf (_file, " %d", &(this->_size));
		//fscanf (_file, " i = %d", &_currentstep);			// some output xyz files have timestep data added
		this->ReadLine();
		this->ReadLine();

		double X, Y, Z;
		char name[10];

		// if we haven't already done so, let's clear out the previous atoms and resize things
		if (!_initialized) {
			for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
				delete *it;
			}
			this->_coords.resize(3*_size, 0.0);
		}

		for (int i = 0; i < _size; i++) {
			// now parse each line's information into the atoms
			fscanf (_file, " %s %lf %lf %lf", name, &X, &Y, &Z);

			// if we haven't already done so, let's create all the atoms we'll need
			if (!_initialized) {
				double * frc = (double *)NULL;
				_atoms.push_back ( new Atom (std::string(name), &_coords[3*i], frc) );
				_atoms.back()->ID (i);
				_atoms.back()->SetAtomProperties();
			}

			// finally we set the position of each atom for the timestep
			this->_coords[3*i] = X;
			this->_coords[3*i+1] = Y;
			this->_coords[3*i+2] = Z;
		}

		_initialized = true;
		_frame++;

		return;
	}	// load next

	void XYZFile::Rewind () {
		rewind(this->_file);
		LoadNext();
		this->_frame = 1;
	} // rewind

	void XYZFile::WriteXYZ (Atom_ptr_vec& atoms) {
		printf ("%d\n", (int)atoms.size());
		printf ("\n");

		for (Atom_it it = atoms.begin(); it != atoms.end(); it++) {
			printf ("%s % 10.5f % 10.5f % 10.5f\n", (*it)->Name().c_str(), (*it)->Position()[x], (*it)->Position()[y], (*it)->Position()[z]);
		}
	}	// write xyz file

} // namespace md_files
