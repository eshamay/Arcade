#ifndef XYZTOPOLOGYFILE_H_
#define XYZTOPOLOGYFILE_H_

#include <cstdio>
#include <string>
#include <vector>
#include <libconfig.h++>

namespace md_files {

	class MolecularTopologyFile	{

		// Format is:
		//		number of molecules
		//		molname_1 num_atoms_1 atom1 atom2 atom3 atom4 ...
		//		molname_2 num_atoms_2 atom1 atom2 atom3 ...
		//		...
		//		molname_x num_atoms_x atom1 atom2 atom3 ...

		public:

			typedef struct {
				char name[10];
				int size;
				std::vector<int> atoms;
			} mol_topology_t;
			
			typedef std::vector<mol_topology_t> mol_topology_vec;
			typedef mol_topology_vec::iterator mol_topology_it;

			MolecularTopologyFile (const std::string filepath) {
				this->LoadFile (filepath);
			}
			virtual ~MolecularTopologyFile () { }

			int size () const { return num_mols; }

			mol_topology_t  Topology (const int id) const { return mols[id]; }
			mol_topology_it begin () { return mols.begin(); }
			mol_topology_it end () { return mols.end(); }

		protected:
			FILE * _topologyFile;

			int num_mols;
			mol_topology_vec							mols;		// the names given to each molecule

			void LoadFile (const std::string filename);
	};	// class molecular topology file

}	// namespace md_files


#endif
