#include "crdfile.h"

namespace md_files {
	CRDFile::CRDFile (std::string const crdpath, int const c_size, const bool periodic) :
		CoordinateFile (crdpath, c_size),
		_periodic(periodic)
	{
		ReadLine (); // skip the first frame's header
		LoadNext ();	// load the first frame of the file
	}


	void CRDFile::LoadNext () {

		// grab each coordinate from the file
		int i = 0;
		for (std::vector<double>::iterator it = _coords.begin(); it != _coords.end(); it++) {
			_eof = (fscanf (_file, "%lf ", &(*it)) == EOF);
		}

		if (_periodic) {
			// process the next frame's header line (grab the box dimensions)
			fscanf (_file, " %lf %lf %lf ", &_dims[0], &_dims[1], &_dims[2]);
		}

		++_frame;

		return;
	}

}

