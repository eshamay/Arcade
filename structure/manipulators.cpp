#include "manipulators.h"


namespace md_analysis {

	SystemManipulator::SystemManipulator (Analyzer * sys) : _system(sys) 
	{ 
		_system->LoadAll();
		std::copy(WaterSystem::sys_atoms.begin(), WaterSystem::sys_atoms.end(), std::back_inserter(all_atoms));
		std::copy(WaterSystem::sys_mols.begin(), WaterSystem::sys_mols.end(), std::back_inserter(all_mols));
		this->Reload();
	}

	void SystemManipulator::Reload () {
		analysis_atoms.clear();
		analysis_mols.clear();
		std::copy(all_atoms.begin(), all_atoms.end(), std::back_inserter(analysis_atoms));
		std::copy(all_mols.begin(), all_mols.end(), std::back_inserter(analysis_mols));
	}

} // namespace md_analysis







namespace h2o_analysis {

	H2OSystemManipulator::H2OSystemManipulator (system_t * t, const int number_of_waters_for_surface_calc) : 
		SystemManipulator(t), 
		reference_point(WaterSystem::SystemParameterLookup("analysis.reference-location")),
		top_surface(WaterSystem::SystemParameterLookup("analysis.top-surface")),
		number_surface_waters(number_of_waters_for_surface_calc) 
	{ 
		this->Reload();
	}

	void H2OSystemManipulator::Reload () {

		// first clear out the previous water set
		for (Wat_it wat = this->all_waters.begin(); wat != this->all_waters.end(); wat++) {
			delete *wat;
		}

		all_waters.clear();
		// then load in the new water set
		this->_system->LoadWaters();
		// gather all the system waters
		for (Mol_it it = WaterSystem::int_wats.begin(); it != WaterSystem::int_wats.end(); it++) {
			WaterPtr wat (new Water(*it));
			wat->SetAtoms();
			all_waters.push_back(wat);
		}

		// grab all the water atoms
		all_water_atoms.clear();
		std::copy(WaterSystem::int_atoms.begin(), WaterSystem::int_atoms.end(), std::back_inserter(all_water_atoms));

		// now update the analysis containers
		this->UpdateAnalysisWaters();
		return;
	}


	void H2OSystemManipulator::UpdateAnalysisWaters () {

		this->analysis_waters.clear();
		std::copy (all_waters.begin(), all_waters.end(), std::back_inserter(analysis_waters));
		//std::for_each(analysis_waters.begin(), analysis_waters.end(), std::mem_fun(&Water::SetAtoms));
		// now all_waters has... all the waters, and analysis wats is used to perform some analysis
		this->analysis_atoms.clear();
		std::copy (all_water_atoms.begin(), all_water_atoms.end(), std::back_inserter(this->analysis_atoms));
	}	// reload analysis wats


	void H2OSystemManipulator::FindWaterSurfaceLocation () {

		// get rid of everything above (or below) the reference point
		if (top_surface) {
			analysis_waters.erase(
					remove_if(analysis_waters.begin(), analysis_waters.end(), typename system_t::MoleculeAbovePosition(reference_point, WaterSystem::axis)), analysis_waters.end());
		}

		else if (!top_surface) {
			analysis_waters.erase(
					remove_if(analysis_waters.begin(), analysis_waters.end(), typename system_t::MoleculeBelowPosition(reference_point, WaterSystem::axis)), analysis_waters.end()); // bottom surface
		}

		// sort the waters by position along the reference axis - first waters are lowest, last are highest
		std::sort (analysis_waters.begin(), analysis_waters.end(), typename system_t::molecule_position_pred(Atom::O));

		if (top_surface) {
			// get the surface waters at the beginning of the list
			std::reverse(analysis_waters.begin(), analysis_waters.end());
		}	// top surface

		// resize the list to contain only the surface waters
		analysis_waters.resize(number_surface_waters);

		// get the position of the top-most waters
		std::vector<double> surface_water_positions;
		// grab all the locations
		std::transform (analysis_waters.begin(), analysis_waters.end(), std::back_inserter(surface_water_positions), WaterLocation());

		// calculate the statistics
		surface_location = gsl_stats_mean (&surface_water_positions[0], 1, number_surface_waters);
		surface_width = gsl_stats_sd (&surface_water_positions[0], 1, number_surface_waters);

		if (surface_width > 2.5) {
			std::cout << std::endl << "Check the pbc-flip or the reference point settings and decrease/increase is to fix this gigantic surface width" << std::endl;
			std::cout << "Here's the positions of the waters used to calculate the surface:" << std::endl;
			std::copy (surface_water_positions.begin(), surface_water_positions.end(), std::ostream_iterator<double>(std::cout, " "));
			printf ("\nSurface width = % 8.3f, Reference point = % 8.3f\n", surface_width, reference_point);
			fflush(stdout);
		}

	}	// find surface water location


	void H2OSystemManipulator::CalcCenterOfMass () {
		center_of_mass.Zero();
		double mass = 0.0;


		for (Atom_it it = this->begin_atoms(); it != this->end_atoms(); it++) {
			mass += (*it)->Mass();
			center_of_mass += (*it)->Position() * (*it)->Mass();
		}
		center_of_mass = center_of_mass / mass;
		return;
	}


	// sort all the waters in the system by distance to a given reference atom
	void H2OSystemManipulator::FindClosestWaters (const AtomPtr a) {
		this->UpdateAnalysisWaters();
		std::sort(analysis_waters.begin(), analysis_waters.end(), typename system_t::molecule_distance_cmp(a));
	} // find closest waters

}	// namespace h2o analysis






namespace so2_analysis {


	void SO2SystemManipulator::Initialize () {
		this->_system->LoadAll();
		this->FindSO2 ();
		this->so2->SetAtoms();
		//printf ("\nFound %zu Sulfur Dioxides in the system\n", so2s.size());
		this->FindAllSO2s();
	}


	void SO2SystemManipulator::UpdateSO2 () {
		delete this->so2;
		for (std::vector<SulfurDioxide *>::iterator it = so2s.begin(); it != so2s.end(); it++) {
			delete *it;
		}
		this->Initialize();
	}


	void SO2SystemManipulator::FindSO2 () {
		// find the so2 in the system and set some pointers up
		int id = WaterSystem::SystemParameterLookup("analysis.reference-molecule-id");
		MolPtr mol = WaterSystem::sys_mols[id];
		//MolPtr mol = Molecule::FindByType(this->_system->sys_mols, Molecule::SO2);
		this->so2 = new SulfurDioxide(mol);
	}

	void SO2SystemManipulator::FindAllSO2s () {
		so2s.clear();
		// load all the system so2s into the container - in case
		for (Mol_it it = WaterSystem::sys_mols.begin(); it != WaterSystem::sys_mols.end(); it++) {
			if ((*it)->Name() == "so2") {
				SulfurDioxide * mol = new SulfurDioxide(*it);
				mol->SetAtoms();
				so2s.push_back(mol);
			}
		}
	}


	vector_map_vec XYZSO2Manipulator::GetWanniers (const AtomPtr& atom, const int num) const {
		std::sort (this->so2->wanniers_begin(), this->so2->wanniers_end(), vecr_distance_cmp(atom->Position()));
		vector_map_vec ret;
		std::copy (this->so2->wanniers_begin(), this->so2->wanniers_begin()+num, std::back_inserter(ret));
		return ret;
	}


	void XYZSO2Manipulator::FindSO2 () {
		for (Mol_it it = WaterSystem::sys_mols.begin(); it != WaterSystem::sys_mols.end(); it++) {
			if ((*it)->MolType() == Molecule::SO2) {
				this->so2 = new SulfurDioxide(*it);
				this->so2->SetAtoms();
			}
		}
	}

} // namespace so2 analysis
