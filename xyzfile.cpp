#include "xyzfile.h"

namespace md_files {

	std::string XYZFile::system_energy;
	int XYZFile::timestep;

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
		//
		// this grabs the number of atoms in the timeframe
		fscanf (_file, " %d", &(this->_size));
		//fscanf (_file, " i = %d", &_currentstep);			// some output xyz files have timestep data added
		this->ReadLine();	// skip the rest of the first line
		// Process the 2nd header line
		this->ReadLine();	
		ParseXYZHeader (std::string(this->_line));

		double X, Y, Z;
		char name[10];

		// if we haven't already done so, let's clear out the previous atoms and resize things
		if (!_initialized) {
			for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
				delete *it;
			}
			this->_coords.resize(3*_size, 0.0);
		}

		//double frc[3];
		for (int i = 0; i < _size; i++) {
			// now parse each line's information into the atoms
			fscanf (_file, " %s %lf %lf %lf", name, &X, &Y, &Z);

			// if we haven't already done so, let's create all the atoms we'll need
			if (!_initialized) {
				AtomPtr new_atom = new Atom (std::string(name), &_coords[3*i]);//, frc);
				_atoms.push_back (new_atom); 
				new_atom->ID(i);
				new_atom->SetAtomProperties();
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


	// trim up the header line of each xyz frame and extract some info from it
	void XYZFile::ParseXYZHeader (std::string header) {
		//std::cout << header << std::endl;
		boost::erase_all(header, " ");
		boost::erase_all(header, "\n");
		//std::cout << header << std::endl;
		std::vector<std::string> strs;
		// split the header - comma delimited
		boost::split(strs, header, boost::is_any_of(","));
		// go through each piece of the header and grab out the relevant info
		std::vector<std::string> tokens;
		for (std::vector<std::string>::iterator it = strs.begin(); it != strs.end(); it++) {
			tokens.clear();
			boost::split(tokens, *it, boost::is_any_of("="));

			// parse the possible tokens
			if (tokens[0] == "E") {
				XYZFile::system_energy = tokens[1];
				//printf ("energy = %s\n", system_energy.c_str());
			}
			else if (tokens[0] == "i") {
				XYZFile::timestep = atoi(tokens[1].c_str());
			}
		}

	}

} // namespace md_files
