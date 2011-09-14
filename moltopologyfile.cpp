#include "moltopologyfile.h"
namespace md_files {

	// open the file for reading
	void MolecularTopologyFile::LoadFile (const std::string filepath) {

		_topologyFile = fopen (filepath.c_str(), "r");
		mols.clear();

		// grab the number of mols to be processed
		fscanf (_topologyFile, " %d ", &num_mols);

		int atom_id;

		// for each molecule (row) we parse the name, and number of atoms in the molecule, and then each atom id in turn
		for (int i = 0; i < num_mols; i++) {
			mol_topology_t mol;
			fscanf (_topologyFile, " %s %d ", mol.name, &mol.size);

			mol.atoms.clear();
			for (int atom = 0; atom < mol.size; atom++) {
				fscanf (_topologyFile, " %d ", &atom_id);
				mol.atoms.push_back(atom_id);
			}
			mols.push_back(mol);
		}

		fclose(_topologyFile);

		/*
		try {
			_topologyFile = new libconfig::Config();
			printf ("\nLoading the topology file: \"%s\"\n", filepath.c_str());

			// load the configuration file
			try {
				_topologyFile->readFile(filepath.c_str());
			} catch (const libconfig::FileIOException fioex) {
				std::cerr << "MolecularTopologyFile::LoadFile\n*** Couldn't find the topology file (" << filepath.c_str() << ")***" << std::endl;
				exit(EXIT_FAILURE);
			}
			catch(const libconfig::SettingTypeException &stex) {
				std::cerr << "MolecularTopologyFile::LoadFile\n*** Something is wrong with the topology file settings - check syntax\n" << std::endl;
				exit(EXIT_FAILURE);
			}
			catch(const libconfig::SettingNotFoundException &snfex) {
				std::cerr << "MolecularTopologyFile::LoadFile\n*** A setting is missing from the topology file!" << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		*/

	}

} // namespace md_files
