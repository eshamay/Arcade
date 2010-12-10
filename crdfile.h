#ifndef CRDFILE_H_
#define CRDFILE_H_

#include "mdfiles.h"

// Class for parsing coordinate files from Amber molecular dynamics trajectories

namespace md_files {
	class CRDFile : public CoordinateFile {

		public:

			CRDFile (std::string const crdpath, int const c_size, const bool periodic = true);
			virtual ~CRDFile () { }

			// Various control functions
			void LoadNext ();

			const VecR& Dimensions () const { return _dims; }

		protected:
			VecR			_dims;		// Dimensions of the system (box size)
			bool			_periodic;	// are periodic boundaries being used
	};

}	// namespace md files
#endif
