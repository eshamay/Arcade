#ifndef WANNIER_H_
#define WANNIER_H_

#include "mdsystem.h"


namespace md_files {

	using namespace md_system;

	class WannierFile : public CoordinateFile {

		protected:
			typedef std::vector<vector_map>	vec_vec_map;
			vec_vec_map _wanniers;

		public:
			WannierFile (std::string wannierpath);
			~WannierFile ();

			typedef vec_vec_map::const_iterator Wannier_it;
			Wannier_it begin () const { return _wanniers.begin(); }
			Wannier_it end () const { return _wanniers.end(); }

			// Various control functions
			void LoadNext ();
	};

}
#endif
