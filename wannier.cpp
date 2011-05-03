#include "wannier.h"

namespace md_files {

	using namespace md_system;

	std::map<Molecule::Molecule_t, int>	WannierFile::numWanniers;

	//std::map<Molecule::Molecule_t, int> WannierFile::numWanniers;
	
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
					this->_eof = false;
					// grab info from the header: number of centers, and frame number
					fscanf (this->_file, " %d 1_%d ", &(this->_size), &(this->_frame));
					this->_coords.resize(this->_size*3, 0.0);
					for (int i = 0; i < this->_size; i++) {
						_wanniers.push_back (Eigen::Map<VecR> (&(this->_coords[3*i])));
					}
				}
			}
			else {
				this->_eof = true;
				printf ("No wannier file specified - continuing without wannier centers\n");
			}

			this->SetNumWanniers();

			return;
		}

	void WannierFile::SetNumWanniers () {
		WannierFile::numWanniers.insert(std::pair<Molecule::Molecule_t, int> (Molecule::H2O, 4));
		WannierFile::numWanniers.insert(std::pair<Molecule::Molecule_t, int> (Molecule::SO2, 9));
		WannierFile::numWanniers.insert(std::pair<Molecule::Molecule_t, int> (Molecule::FORMALDEHYDE, 6));
	}

	WannierFile::~WannierFile () { 
		fclose(_file); 
		_file = (FILE *)NULL;
	}

	void WannierFile::LoadNext () {

		//double a, b, c;
		// grab each coordinate vector for each wannier center until the size of the system is processed
		for (int i = 0; i < this->_size; i++) {
			//printf ("% 8.4f % 8.4f % 8.4f\n", a,b,c);
			this->_eof = (fscanf (_file, "X %lf %lf %lf %*f ", &this->_coords[3*i], &this->_coords[3*i+1], &this->_coords[3*i+2]) == EOF);
			//this->_eof = (fscanf (_file, " X %lf %lf %lf %*f %*f %*f ", &this->_coords[3*i], &this->_coords[3*i+1], &this->_coords[3*i+2]) == EOF);
			//	if (fscanf (_file, "X %lf %lf %lf %*f ", &_coords[3*i], &_coords[3*i+1], &_coords[3*i+2]) == EOF) 
			//		_eof=true;
		}

		if (!_eof)
			// grab info from the header: number of centers, and frame number
			fscanf (this->_file, " %d 1_%d ", &this->_size, &this->_frame);

		return;
	}

}	// namespace md_files
