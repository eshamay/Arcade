#include "crdfile.h"

namespace md_files {
	CRDFile::CRDFile (std::string const crdpath, int const c_size, const bool periodic) :
		CoordinateFile (crdpath, c_size),
		_periodic(periodic) {
//			ReadLine (); // skip the first frame's header
			rewind (_file);
			LoadNext ();	// load the first frame of the file
		}


	void CRDFile::LoadNext () {

		// grab each coordinate from the file
		float t;
		for (std::vector<double>::iterator it = _coords.begin(); it != _coords.end(); it++) {
			_eof = feof(_file);

			fread (&t, sizeof(float), 1, _file);
			*it = t;
		}

		float dims[3];
		if (_periodic) {
			// process the next frame's header line (grab the box dimensions)
			fread (dims, sizeof(float), 3, _file);
		}
		_dimensions[0] = dims[0];
		_dimensions[1] = dims[1];
		_dimensions[2] = dims[2];

		++_frame;

		return;
	}

	void CRDFile::Rewind () {
		rewind(this->_file);
		//ReadLine();
		LoadNext();
		this->_frame = 1;
	}	// rewind
}

