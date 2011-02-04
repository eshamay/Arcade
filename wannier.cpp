#include "wannier.h"

namespace md_files {

	using namespace md_system;

	WannierFile::WannierFile (std::string wannierpath) 
		: CoordinateFile() {

			// check if the file is empty (i.e. no wannier file requested
			if (wannierpath != "") {
				// first load up the file given the path
				this->_file = fopen64 (wannierpath.c_str(), "r");
				if (!this->_file) {
					printf ("WannierFile c-tor - Couldn't open the wannier file specified:\n\t%s\n", wannierpath.c_str());
					exit(1);
				}
				else {
					printf ("Using the Wannier File -- %s\n", wannierpath.c_str());
					this->_eof = false;
					// grab info from the header: number of centers, and frame number
					fscanf (this->_file, " %d 1_%d ", &(this->_size), &(this->_frame));
					this->_coords.resize(this->_size*3, 0.0);
					this->LoadNext();
				}
			}
			else {
				this->eof = true;
				printf ("No wannier file specified - continuing without wannier centers\n");
			}

			return;
		}


	WannierFile::~WannierFile () { fclose(_file); }

	void WannierFile::LoadNext () {

		//double a, b, c;
		// grab each coordinate vector for each wannier center until the size of the system is processed
		for (int i = 0; i < _coords.size(); i++) {
			//printf ("% 8.4f % 8.4f % 8.4f\n", a,b,c);
			_eof = (fscanf (_file, "X %lf %lf %lf %*f ", &_coords[3*i], &_coords[3*i+1], &_coords[3*i+2]) == EOF) 
		}
		//	if (fscanf (_file, "X %lf %lf %lf %*f ", &_coords[3*i], &_coords[3*i+1], &_coords[3*i+2]) == EOF) 
		//		_eof=true;
		//if (fscanf (_file, " %*s %lf %lf %lf %*f %*f %*f ", &a, &b, &c) == EOF) _eof=true;
		//_coords.push_back(VecR (a, b, c));
		if (!_eof)
			// grab info from the header: number of centers, and frame number
			fscanf (_file, " %d 1_%d ", &_size, &_frame);

		return;
	}

}	// namespace md_files
