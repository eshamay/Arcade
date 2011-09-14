#include "watersystem.h"

namespace md_system { 

	libconfig::Config * WaterSystem::config_file;	
	double WaterSystem::posmin;
	double WaterSystem::posmax;
	double WaterSystem::pbcflip;
	coord WaterSystem::axis;
	VecR WaterSystem::ref_axis;

	Atom_ptr_vec WaterSystem::sys_atoms;
	Mol_ptr_vec WaterSystem::sys_mols;
	Mol_ptr_vec WaterSystem::int_wats;
	Mol_ptr_vec WaterSystem::int_mols;
	Atom_ptr_vec WaterSystem::int_atoms;


	WaterSystem::WaterSystem (const std::string configuration_filename) 
	{
		try {
			config_file = new libconfig::Config();
			printf ("\nLoading configuration file: \"%s\"\n", configuration_filename.c_str());

			// load the configuration file
			try {
				config_file->readFile(configuration_filename.c_str());
			} catch (const libconfig::FileIOException fioex) {
				std::cerr << "WaterSystem::c-tor\n*** Couldn't find the configuration file (" << configuration_filename.c_str() << ")***" << std::endl;
				exit(EXIT_FAILURE);
			}

			posmin = SystemParameterLookup("analysis.position-range")[0];
			posmax = SystemParameterLookup("analysis.position-range")[1];
			axis = (coord)((int)SystemParameterLookup("analysis.reference-axis"));
			ref_axis = VecR(
					SystemParameterLookup("analysis.reference-vector")[0], 
					SystemParameterLookup("analysis.reference-vector")[1], 
					SystemParameterLookup("analysis.reference-vector")[2]);
			pbcflip = SystemParameterLookup("analysis.PBC-flip");
		}
		catch(const libconfig::SettingTypeException &stex) {
			std::cerr << "Something is wrong with the configuration parameters or file - check syntax\n(watersystem.h)" << std::endl;
			exit(EXIT_FAILURE);
		}
		catch(const libconfig::SettingNotFoundException &snfex) {
			std::cerr << "A setting is missing from the configuration file!" << std::endl;
			exit(EXIT_FAILURE);
		}


		return;
	}

	WaterSystem::~WaterSystem () {
		delete config_file;
		return;
	}


	/*
		 void WaterSystem< gromacs::GMXSystem<gromacs::TRRFile> >::_InitializeSystem () {
		 std::string gro = this->SystemParameterLookup("system.files.gmx-grofile");
		 std::string trr = this->SystemParameterLookup("system.files.gmx-trrfile");
		 this->sys = new gromacs::GMXSystem< gromacs::TRRFile >(gro.c_str(), trr.c_str());
		 return;
		 }

		 void WaterSystem< gromacs::GMXSystem<gromacs::XTCFile> >::_InitializeSystem () {
		 std::string gro = this->SystemParameterLookup("system.files.gmx-grofile");
		 std::string xtc = this->SystemParameterLookup("system.files.gmx-xtcfile");
		 this->sys = new gromacs::GMXSystem< gromacs::XTCFile >(gro.c_str(), xtc.c_str());
		 return;
		 }
		 */

	void WaterSystem::LoadAll () {

		sys_mols.clear();
		sys_atoms.clear();
		int_mols.clear();
		int_atoms.clear();

		// copy the molecules and atoms into each container
		std::copy(sys->begin_mols(), sys->end_mols(), std::back_inserter(sys_mols));
		std::copy(sys->begin_mols(), sys->end_mols(), std::back_inserter(int_mols));

		std::copy(sys->begin(), sys->end(), std::back_inserter(sys_atoms));
		std::copy(sys->begin(), sys->end(), std::back_inserter(int_atoms));

		return;
	}

} // namespace watersystem
