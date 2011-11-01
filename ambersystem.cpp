#include "ambersystem.h"

namespace md_files {

	AmberSystem::AmberSystem (const std::string& prmtop, const std::string& mdcrd, const bool periodic)//, const std::string& mdvel)
		// some initialization needs to happen here
		: 	
			MDSystem (),
			_topfile(prmtop),
			_coords(mdcrd, _topfile.NumAtoms()),
			_periodic(periodic)
			//_forces(mdvel, _topfile.NumAtoms())
	{
		_atoms = md_system::Atom_ptr_vec(_topfile.NumAtoms(), (md_system::AtomPtr)NULL);

		// A lot of functionality depends on knowing the system size - so we set it here
		MDSystem::Dimensions (_coords.Dimensions());

		// and parse all the info out of the topology file into the atoms
		this->_ParseAtomInformation ();
		// then with all the atoms setup, group them into molecules
		this->_ParseMolecules ();

		return;
	}

	AmberSystem::~AmberSystem () {
		for (md_system::Mol_ptr_vec::iterator it = _mols.begin(); it != _mols.end(); it++) {
			delete *it;
		}
		for (md_system::Atom_ptr_vec::iterator it = _atoms.begin(); it != _atoms.end(); it++) {
			delete *it;
		}
	}

	// While the crdfile holds spatial coordinate information, and the topology file holds atomic information, the data has to be processed into proper atoms in order to play around with them more effectively.
	// This function should only be used once when first loading the system up. This info doesn't change during the course of the MD run
	void AmberSystem::_ParseAtomInformation () {
		//double * frc;	// This isn't set yet!
		for (int i = 0; i < _topfile.NumAtoms(); i++) {
			_atoms[i] = new md_system::Atom(_topfile.AtomNames()[i], _coords(i));//, frc);
			_atoms[i]->ID(i);// set the atom's index number - because we may need to access ordered/list info elsewhere
		}

		return;
	}

	/* the mol pointers in the topology file list the beginning and end atoms of each molecule. This function will group all atoms in a molecule, form a molecule object, and add it to the _mols vector
	*/
	void AmberSystem::_ParseMolecules () {

		md_system::MolPtr newmol;
		// Here we run through each mol pointer, create a new molecule, and add in the appropriate atoms
		for (int mol = 0; mol < _topfile.NumMols(); mol++) {
			// At each new pointer we create a molecule and start adding in atoms
			// if the molecule is of a specific type for which a class has been created, then let's use that!
			std::string name = _topfile.MolNames()[mol];
			newmol = md_system::MoleculeFactory(name);
			newmol->Name(name);
			//_mols[mol]->Name (name);	// set the molecule's residue name
			newmol->MolID(mol);

			// Then add all the atoms between the indices of the molpointers in the topology file
			int molsize = _topfile.MolSizes()[mol];
			int molpointer = _topfile.MolPointers()[mol];
			for (int atomCount = 0; atomCount < molsize; atomCount++) {
				// sets the index of the current atom that's being added to the molecule
				int curAtom = molpointer + atomCount - 1;
				// Now we're going to bless this new atom with loads of information about itself and its molecule
				newmol->AddAtom( _atoms[curAtom] );			// add the atom into the molecule
			}

			// lastly, tack it into our running list
			_mols.push_back(newmol);

			/*
			// a topology file doesn't give us the type, so that will be determined when creating the new molecules
			if (name == "no3") {
			_mols.push_back (new Nitrate());
			}
			else if (name == "hno3") {
			_mols.push_back (new NitricAcid());
			}
			// we might also have water molecules!
			else if (name == "h2o") {
			_mols.push_back (new Water());
			}
			else if (name == "so2" || name == "sog" || name == "soq") {
			_mols.push_back (new SulfurDioxide());
			}
			else if (name == "dec" || name == "pds") {
			_mols.push_back (new Decane());
			}
			// otherwise add on a generic molecule
			else {
			printf ("AmberSystem::_ParseMolecules() -- Couldn't parse the molecule with name '%s'. Don't know what to do with it\n", name.c_str());
			exit(1);
			}
			*/
		}

		return;
	}

	void AmberSystem::LoadFirst () {
		//_coords.LoadFirst();
		//if (_forces.Loaded()) _forces.LoadFirst();
		//this->_ParseAtomVectors ();
		return;
	}

	void AmberSystem::LoadNext () {
		_coords.LoadNext ();							// load up coordinate information from the file
		//if (_forces.Loaded()) _forces.LoadNext ();		// also load the force information while we're at it
		//this->_ParseAtomVectors ();
		return;
	}

	/*
	void AmberSystem::PrintCRDFile () const {

		for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
			printf ("  % 4.7f  % 4.7f  % 4.7f", (*it)->Position()[x], (*it)->Position()[y], (*it)->Position()[z]);
		}
		printf ("  % 4.7f  % 4.7f  % 4.7f", this->Dims()[x], this->Dims()[y], this->Dims()[z]);

		return;
	}
	*/

}	// namespace md_files
