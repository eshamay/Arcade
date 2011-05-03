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

			typedef vector_map_it Wannier_it;
			Wannier_it begin () { return _wanniers.begin(); }
			Wannier_it end () { return _wanniers.end(); }

			// a mapping for the number of wannier centers for each type of molecule
			static std::map<Molecule::Molecule_t, int>	numWanniers;
			void SetNumWanniers ();

			vector_map& operator[] (const int index) { 
				if (index < 0 || index > (int)_wanniers.size()-1) {
					std::cerr << "bad wannier file index -> " << index << " of " << _wanniers.size() << std::endl;
					exit(1);
				}
				return _wanniers[index]; 
			}

			// Various control functions
			void LoadNext ();
	};

}
#endif
