#ifndef CRDFILE_H_
#define CRDFILE_H_

#include "mdsystem.h"

// Class for parsing coordinate files from Amber molecular dynamics trajectories

namespace md_files {
	class CRDFile : public md_system::CoordinateFile {

		public:

			CRDFile (std::string const crdpath, int const c_size, const bool periodic = true);
			virtual ~CRDFile () { }

			// Various control functions
			void LoadNext ();

			const VecR& Dimensions () const { return _dimensions; }

		protected:
			VecR			_dimensions;		// Dimensions of the system (box size)
			bool			_periodic;	// are periodic boundaries being used
	};

}	// namespace md files
#endif
