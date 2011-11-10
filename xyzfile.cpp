#include "xyzfile.h"

namespace md_files {

	std::string XYZFile::system_energy;
	int XYZFile::timestep;

	XYZFile::XYZFile (std::string path) 
		: 
			md_system::CoordinateFile (path),
			_initialized(false) {

				_file = fopen (path.c_str(), "rb");
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
		//fscanf (_file, " %d", &(this->_size));
		//fscanf (_file, " i = %d", &_currentstep);			// some output xyz files have timestep data added
		//this->ReadLine();	// skip the rest of the first line
		// Process the 2nd header line
		//this->ReadLine();	
		//ParseXYZHeader (std::string(this->_line));

		//double X, Y, Z;
		//char name[10];

		// if we haven't already done so, let's clear out the previous atoms and resize things
		if (!_initialized) {
			rewind (this->_file);
			fread (&(this->_size), sizeof(unsigned int), 1, this->_file);

			for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
				delete *it;
			}
			_atoms.clear();
			this->_coords.resize(3*_size, 0.0);

			char name[5];
			unsigned int len;
			for (int i = 0; i < this->_size; i++) {
				fread (&len, sizeof(unsigned int), 1, this->_file);
				fread (name, sizeof(char), len+1, this->_file);
				AtomPtr new_atom = new Atom (std::string(name), &(this->_coords[3*i]));
				_atoms.push_back (new_atom); 
				new_atom->ID(i);
				new_atom->SetAtomProperties();
			}
			_initialized = true;
		}

		for (int i = 0; i < this->_size; i++) {
			fread (&(this->_coords[3*i]), sizeof(double), 3, this->_file);
			//fread (&_coords[0], sizeof(double), 3*this->_size, _file);
			// set the position of each atom for the timestep
			//this->_coords[3*i] = pos[0];
			//this->_coords[3*i+1] = pos[1];
			//this->_coords[3*i+2] = pos[2];
		}

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
