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
		upper_reference_point(WaterSystem::SystemParameterLookup("analysis.upper-reference-point")),
		lower_reference_point(WaterSystem::SystemParameterLookup("analysis.lower-reference-point")),
		top_surface(WaterSystem::SystemParameterLookup("analysis.top-surface")),
		number_surface_waters(number_of_waters_for_surface_calc) 
	{ 
		this->Reload();
		//this->upper_reference_point = MDSystem::Dimensions()[WaterSystem::axis];
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
					remove_if(analysis_waters.begin(), analysis_waters.end(), typename system_t::MoleculeAbovePosition(upper_reference_point, WaterSystem::axis)), analysis_waters.end());
		}

		else if (!top_surface) {
			analysis_waters.erase(
					remove_if(analysis_waters.begin(), analysis_waters.end(), typename system_t::MoleculeBelowPosition(lower_reference_point, WaterSystem::axis)), analysis_waters.end()); // bottom surface
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

		/*
		if (surface_width > 2.5) {
			std::cout << std::endl << "Check the pbc-flip or the reference point settings and decrease/increase is to fix this gigantic surface width" << std::endl;
			std::cout << "Here's the positions of the waters used to calculate the surface:" << std::endl;
			std::copy (surface_water_positions.begin(), surface_water_positions.end(), std::ostream_iterator<double>(std::cout, " "));
			printf ("\nSurface width = % 8.3f, Reference point = % 8.3f\n", surface_width, reference_point);
			fflush(stdout);
		}
		*/

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


	void H2ODoubleSurfaceManipulator::FindWaterSurfaceLocation () {


		VecR shift (0.0,0.0,0.0);
		shift[WaterSystem::axis] = MDSystem::Dimensions()[WaterSystem::axis];

		// wrap all waters below the pbc-flip boundary up a box
		for (Wat_it it = analysis_waters.begin(); it != analysis_waters.end(); it++) {
			if ((*it)->ReferencePoint()[WaterSystem::axis] < this->lower_reference_point) {
				(*it)->Shift(shift);
			}
		}

		/*
		printf ("shift vector = ");
		shift.Print();

		printf ("\nreference point = %f\n", this->reference_point);
		*/

		// wrap all waters above the top reference point down a box
		shift = -shift;
		for (Wat_it it = analysis_waters.begin(); it != analysis_waters.end(); it++) {
			if ((*it)->ReferencePoint()[WaterSystem::axis] > this->upper_reference_point) {
				(*it)->Shift(shift);
			}
		}

		// sort the waters by position along the reference axis - first waters are lowest, last are highest
		std::sort (analysis_waters.begin(), analysis_waters.end(), typename system_t::molecule_position_pred(Atom::O));
		// for both the bottom and top waters, grab a certain number of them and calculate the stats
		
		// get the position of the bottom-most waters
		std::vector<double> surface_water_positions;
		std::transform (analysis_waters.begin(), analysis_waters.begin()+number_surface_waters, 
				std::back_inserter(surface_water_positions), WaterLocation());

		// calculate the statistics for the bottom surface
		bottom_location = gsl_stats_mean (&surface_water_positions[0], 1, number_surface_waters);
		bottom_width = gsl_stats_sd (&surface_water_positions[0], 1, number_surface_waters);

		// find the top water statistics, similarly - starting from the other end of the water list
		surface_water_positions.clear();
		std::transform (analysis_waters.rbegin(), analysis_waters.rbegin()+number_surface_waters, 
				std::back_inserter(surface_water_positions), WaterLocation());

		// calculate the statistics
		top_location = gsl_stats_mean (&surface_water_positions[0], 1, number_surface_waters);
		top_width = gsl_stats_sd (&surface_water_positions[0], 1, number_surface_waters);

		if (top_width > 3.0) {
			std::for_each (analysis_waters.rbegin(), analysis_waters.rbegin()+number_surface_waters, std::mem_fun(&Molecule::Print));
		}

		//printf ("\n\ntop = %.2f  %.2f\n", top_location, top_width);
		//printf ("bottom = %.2f  %.2f\n", bottom_location, bottom_width);

	}	// find surface water location - double surface



	// find the shortest distance between the two points given wrapping the unit cell
	double H2ODoubleSurfaceManipulator::WrappedDistance (double a, double b, double dim) const {
		// take care of wrapping into the central unit cell
		while (fabs(a-b) > dim/2.0) {
			if (a < b) a += dim;
			else 		 a -= dim;
		}

		return b-a;
	}

	surface_distance_t H2ODoubleSurfaceManipulator::TopOrBottom (const double pos) const {

		double dim = MDSystem::Dimensions()[WaterSystem::axis];
		double top = WrappedDistance(top_location, pos, dim);
		double bottom = WrappedDistance(pos, bottom_location, dim);

		// if the distance from the position to the top is smaller than the distance from the position to the bottom surface, 
		// then it's closer to the top surface, and we return true.
		bool surface = fabs(top) < fabs(bottom);
		double distance = (surface) ? top : bottom;

		return std::make_pair(surface, distance);
	}


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


	void CycleManipulator::BuildGraph () { 

		// check that the reference atom was set
		this->ReferenceAtom();

		// order the atoms in the system in closest to furthest from a particular reference point
		nm.OrderAtomsByDistance(ref_atom);

		// grab only a certain number of the closest atoms for analysis
		this->analysis_atoms.clear();
		this->analysis_atoms.push_back(ref_atom);

		// check that the cycle size was set
		this->CycleSize();

		std::copy(nm.closest(), nm.closest() + cycle_size, std::back_inserter(this->analysis_atoms));

		graph.UpdateGraph(this->analysis_atoms); 

	}	// build graph


	void CycleManipulator::ParseCycles () { 

		// determine if within the graph there are any closed cycles
		std::list<AtomPtr> gray_sources;
		std::list<AtomPtr> gray_targets;
		bondgraph::bfs_atom_visitor vis (gray_sources, gray_targets);
		boost::breadth_first_search(graph.Graph(), *graph._FindVertex(ref_atom), visitor(vis));

		// grab the atoms that make up the cycle
		std::list<AtomPtr>::const_iterator gray_source, gray_target;
		gray_source = gray_sources.begin();
		gray_target = gray_targets.begin();

		cycle.clear();
		cycle_type.clear();

		//printf ("\n --> %d\n", (int)std::distance (gray_sources.begin(), gray_sources.end()));
		while (gray_source != gray_sources.end() && gray_target != gray_targets.end()) {
			cycle_t new_cycle_type = NO_CYCLE;
			cycle_list new_cycle;

			// fill the cycle's atom list here by first getting all the atoms from the target to the ref-atom, and then from the source so that both sides of the cycle are accounted for.
			// cycle list will look like S -- O1 -- etc. -- O2
			for (AtomPtr a = *gray_target; a != ref_atom; a = graph.Parent(a)) {
				new_cycle.push_front(a);
			}
			new_cycle.push_front(ref_atom);
			for (AtomPtr a = *gray_source; a != ref_atom; a = graph.Parent(a)) {
				new_cycle.push_back(a);
			}

			// the first atom in the cycle is the reference. Here we get the atom connected to the reference - the first non-ref atom - and the last atom in the cycle
			// based on the molecules these are connected to, we discover the type of cycle they're involved in
			cycle_it _first = new_cycle.begin(); _first++;	
			AtomPtr atom1 = *_first;
			AtomPtr atom2 = new_cycle.back();

			Molecule::Molecule_t mol1_t, mol2_t;
			mol1_t = atom1->ParentMolecule()->MolType(); 
			mol2_t = atom2->ParentMolecule()->MolType(); 

			if (atom1 != atom2) {

				if ((mol1_t != Molecule::SO2 && mol2_t == Molecule::SO2) ||
						(mol1_t == Molecule::SO2 && mol2_t != Molecule::SO2)) {
					new_cycle_type = HALFBRIDGE;
					//printf ("\nhalf-bridge\n");
				}
				else if (mol1_t != Molecule::SO2 && mol2_t != Molecule::SO2) {
					new_cycle_type = FULLCROWN; 
					//printf ("\nfull-crown\n");
				}
				else if (mol1_t == Molecule::SO2 && mol2_t == Molecule::SO2) {
					new_cycle_type = FULLBRIDGE;
					//printf ("\nfull-bridge\n");
				}
				else {
					new_cycle_type = UNKNOWN;
				}
			}

			else {
				if (mol1_t == Molecule::SO2 && mol2_t == Molecule::SO2) { 
					new_cycle_type = WATERLEG; 
					//printf ("\nwater-leg\n");
				}

				else if (mol1_t != Molecule::SO2 && mol2_t != Molecule::SO2) {
					new_cycle_type = HALFCROWN;
					//printf ("\nhalf-crown\n");
				}
				else {
					new_cycle_type = UNKNOWN;
					std::cerr << "found some other type of funky cycle\n" << std::endl;
				}
			}

			//std::for_each (cycle.begin(), cycle.end(), std::mem_fun(&Atom::Print));
			cycle.push_back(new_cycle);
			cycle_type.push_back(new_cycle_type);

			gray_source++; gray_target++;
		}	// while

	}// parse cycles


	void CycleManipulator::FindUniqueMembers (const cycle_list& _cycle) {

		// find all the unique atoms in the cycle
		unique_cycle_atoms.clear();
		std::copy (_cycle.begin(), _cycle.end(), std::back_inserter(unique_cycle_atoms));
		std::sort(unique_cycle_atoms.begin(), unique_cycle_atoms.end(), std::ptr_fun(&Atom::id_cmp));
		Atom_it unique_atoms_it = std::unique(unique_cycle_atoms.begin(), unique_cycle_atoms.end(), std::ptr_fun(&Atom::id_eq));
		unique_cycle_atoms.resize(unique_atoms_it - unique_cycle_atoms.begin());

		unique_cycle_mols.clear();
		std::transform(
				unique_cycle_atoms.begin(), 
				unique_cycle_atoms.end(), 
				std::back_inserter(unique_cycle_mols),
				std::mem_fun<MolPtr,Atom>(&Atom::ParentMolecule));

		std::sort(unique_cycle_mols.begin(), unique_cycle_mols.end(), std::ptr_fun(&Molecule::mol_cmp));
		Mol_it unique_mols_it = std::unique(unique_cycle_mols.begin(), unique_cycle_mols.end(), std::ptr_fun(&Molecule::mol_eq));
		unique_cycle_mols.resize(unique_mols_it - unique_cycle_mols.begin());

		//std::for_each(cycle.begin(), cycle.end(), std::mem_fun(&Atom::Print));
		//printf ("unique = %zu %zu\n", unique_cycle_atoms.size(), unique_cycle_mols.size());
	}

	int CycleManipulator::NumUniqueWaterAtoms () const {
		int num = 0;
		for (Atom_it it = unique_cycle_atoms.begin(); it != unique_cycle_atoms.end(); it++) {
			if ((*it)->ParentMolecule()->MolType() == Molecule::H2O) {
				++num;
			}
		}
		return num;
	}


	// go through the cycle and find the location of the atom that attaches the cycle to the bridge from the so2 leg. Return the first and second occurrence.
	typename CycleManipulator::cycle_pair_t CycleManipulator::CyclePointAtom (const cycle_list& _cycle) {
		cycle_it it = _cycle.begin(); it++;	// first non-ref atom
		cycle_it jt = _cycle.end(); jt--;		// last atom
		cycle_it i_point = it;
		cycle_it j_point = jt;
		while (true) {
			it++;
			jt--;
			if (*it != *jt)
				break;
			i_point = it;
			j_point = jt;
		} 
		return std::make_pair(i_point,j_point);
	}


	std::pair<int,int> CycleManipulator::WaterLegInformation (const cycle_list& _cycle) {
		cycle_pair_t point_atom = CyclePointAtom(_cycle);
		int leg_cycle_location = std::distance(_cycle.begin(), point_atom.first);
		int leg_cycle_size = std::distance(point_atom.first, point_atom.second);
		return std::make_pair(leg_cycle_location, leg_cycle_size); // don't count an atom that's part of the cycle in the point
	}

} // namespace md analysis
