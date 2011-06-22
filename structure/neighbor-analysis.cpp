#include "neighbor-analysis.h"

namespace neighbor_analysis {

	void CycleManipulator::BuildGraph () { 

		this->ReferenceAtom();

		// order the atoms in the system in closest to furthest from a particular reference point
		nm.OrderAtomsByDistance(ref_atom);

		// grab only a certain number of the closest atoms for analysis
		this->analysis_atoms.clear();
		this->analysis_atoms.push_back(ref_atom);

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

				if ((mol1_t == Molecule::H2O && mol2_t == Molecule::SO2) ||
						(mol1_t == Molecule::SO2 && mol2_t == Molecule::H2O)) {
					new_cycle_type = HALFBRIDGE;
					//printf ("\nhalf-bridge\n");
				}
				else if (mol1_t == Molecule::H2O && mol2_t == Molecule::H2O) {
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

				else if (mol1_t == Molecule::H2O && mol2_t == Molecule::H2O) {
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







	void SO2BondingCycleAnalysis::Analysis() {

		// find the locatin of the surface of the water slab
		h2os.Reload();
		h2os.FindWaterSurfaceLocation();
		//h2os.FindWaterSurfaceLocation(false);	 // bottom surface
		double distance = system_t::Position(so2s.S()) - h2os.SurfaceLocation(); // top surface
		//double distance = h2os.SurfaceLocation() - system_t::Position(so2s.S());	// bottom surface
		fprintf (this->output, " % 8.3f ", distance);

		// find the H-bonding on each of the so2 oxygens
		h2os.Reload();

		// select the reference atom to use for determining h-bonding to the so2.
		cm.SetReferenceAtom(so2s.O1());
		cm.SetCycleSize (20);
		cm.BuildGraph();
		// output the number of bonds on the first oxygen
		fprintf (this->output, " %5d ", cm.NumReferenceAtomHBonds());

		// do the same for the 2nd oxygen
		cm.SetReferenceAtom(so2s.O2());
		cm.BuildGraph();
		fprintf (this->output, " %5d ", cm.NumReferenceAtomHBonds());

		// select the sulfur atom, get the number of Os it's playing with, and find any cycles within the neighboring waters in the graph
		cm.SetReferenceAtom(so2s.S());
		cm.BuildGraph();
		fprintf (this->output, " %5d ", cm.NumReferenceInteractions());

		cm.ParseCycles();

		typename std::list<typename CycleManipulator::cycle_list>::const_iterator	cycle = cm.cycle_begin();	
		typename std::list<typename CycleManipulator::cycle_t>::const_iterator		cycle_type = cm.cycle_type_begin();
		// once a cycle is detected - do something
		while (cycle != cm.cycle_end() && cycle_type != cm.cycle_type_end()) {
			int number_unique_atoms_in_cycle = 0;
			int number_molecules_in_cycle = 0;

			cm.FindUniqueMembers(*cycle);

			// full point cycle
			if (*cycle_type == CycleManipulator::FULLBRIDGE) {
				number_unique_atoms_in_cycle = cm.NumUniqueAtomsInCycle() - 3;	// don't count the so2 atoms
				number_molecules_in_cycle = cm.NumUniqueMoleculesInCycle() - 1;		// don't count the SO2
			}
			// half point cycle
			else if (*cycle_type == CycleManipulator::HALFBRIDGE) {
				number_unique_atoms_in_cycle = cm.NumUniqueAtomsInCycle() - 2;	// don't count the so2 atoms
				number_molecules_in_cycle = cm.NumUniqueMoleculesInCycle() - 1;		// don't count the SO2
			}
			// full crown cycle
			else if (*cycle_type == CycleManipulator::FULLCROWN) {
				number_unique_atoms_in_cycle = cm.NumUniqueAtomsInCycle() - 1;	// don't count the so2 atoms
				number_molecules_in_cycle = cm.NumUniqueMoleculesInCycle() - 1;		// don't count the SO2
			}
			// water leg cycle
			else if (*cycle_type == CycleManipulator::WATERLEG || *cycle_type == CycleManipulator::HALFCROWN) {
				std::pair<int,int> water_leg = cm.WaterLegInformation(*cycle);	// find the atoms that delimit the cycle
				number_unique_atoms_in_cycle = water_leg.first;	// this is the position of the cycle
				number_molecules_in_cycle = water_leg.second;		// this is the size of the cycle (# of atoms)
			}

			fprintf (this->output, " %5d %5d %5d ", 
					(int)*cycle_type, number_unique_atoms_in_cycle,	number_molecules_in_cycle);

			cycle++; cycle_type++;
		}	// while
		fprintf (this->output, "\n");
		fflush(this->output);

	} // analysis



	void SO2HBondingAnalysis::Analysis () {
		// find the locatin of the surface of the water slab and the distance of the so2 to it
		h2os.Reload();
		//h2os.FindWaterSurfaceLocation(true);	// top surface
		h2os.FindWaterSurfaceLocation();	// bottom surface
		//double distance = system_t::Position(so2s.S()) - h2os.SurfaceLocation();	// top surface
		double distance;
		if (h2os.TopSurface())
			distance = system_t::Position(so2s.S()) - h2os.SurfaceLocation();	// bottom surface
		else
			distance = h2os.SurfaceLocation() - system_t::Position(so2s.S());	// top surface

		this->LoadAll();
		// sort the atoms in the system in order of distance from the SO2
		std::sort(
				this->begin(), 
				this->end(), 
				typename system_t::atomic_distance_cmp (so2s.S()));

		//nearest.clear();
		//std::copy (this->begin(), this->begin+20, std::back_inserter(nearest));
		// build the graph using a handful of the closest atoms
		graph.UpdateGraph(this->begin(), this->begin() + 20);

		// grab the number of atoms connected to the reference atom through bonding, and check if the bond looks right.
		bonded.clear();

		// check for when the so2 has a single bond on the S and no other bonds
		// O2S-OH2
		/*
			 bonded = graph.BondedAtoms (so2s.S(), bondgraph::interaction, Atom::O);
			 bool good = false;
			 if (bonded.size() == 1 && bonded[0]->ParentMolecule()->MolType() == Molecule::H2O) {
			 good = true;
			 }
			 bonded = graph.BondedAtoms (so2s.O1(), bondgraph::hbond, Atom::H);
			 if (bonded.size() >= 1)
			 good = false;
			 bonded = graph.BondedAtoms (so2s.O2(), bondgraph::hbond, Atom::H);
			 if (bonded.size() >= 1)
			 good = false;

			 if (good)
			 histo(distance);
			 */

		// checking for the OSO--HOH interaction on a single water - only a 1:1 interaction
		// One O bond, and no S bonds
		int bonds = 0;
		bonded = graph.BondedAtoms (so2s.O1(), bondgraph::hbond, Atom::H);
		if (bonded.size() == 1 && bonded[0]->ParentMolecule()->MolType() == Molecule::H2O)
			++bonds;
		bonded = graph.BondedAtoms (so2s.O2(), bondgraph::hbond, Atom::H);
		if (bonded.size() == 1 && bonded[0]->ParentMolecule()->MolType() == Molecule::H2O)
			++bonds;
		bonded = graph.BondedAtoms (so2s.S(), bondgraph::hbond, Atom::H);
		if (bonds == 1 && bonded.empty())
			histo(distance);

	}	// so2 h bonding analysis





	void SO2NearestNeighborAnalysis::Analysis () {

		h2os.FindWaterSurfaceLocation();

		// At first we iterate through and find the several atoms that bind to the SO2
		// To do this we check each atom (S, O1, O2) for any hbonds/interactions, and if we find them 
		// and add it into our list of the first several that bind. 
		if (first_pass) {
			Atom_ptr_vec graph_atoms;
			std::copy (h2os.begin_atoms(), h2os.end_atoms(), std::back_inserter(graph_atoms));
			graph_atoms.push_back (so2s.S());
			graph_atoms.push_back (so2s.O1());
			graph_atoms.push_back (so2s.O2());
			graph.UpdateGraph(graph_atoms.begin(), graph_atoms.end());

			// find any oxygens bound to the S
			Atom_ptr_vec bound_atoms = graph.BondedAtoms (so2s.S(), bondgraph::interaction);

			// add in any unaccounted-for atoms that have bound to the S
			if (!bound_atoms.empty()) {
				for (Atom_it it = bound_atoms.begin(); it != bound_atoms.end(); ++it) {
					if (std::find (first_os.begin(), first_os.end(), *it) != first_os.end()) {
						first_os.push_back(*it);
						std::cout << "found an O on S" << std::endl;
					}
				}
			}

			// do the same for the first and second oxygens
			bound_atoms = graph.BondedAtoms (so2s.O1(), bondgraph::hbond);

			// add in any unaccounted-for atoms that have bound
			if (!bound_atoms.empty()) {
				for (Atom_it it = bound_atoms.begin(); it != bound_atoms.end(); ++it) {
					if (std::find (first_hs_o1.begin(), first_hs_o1.end(), *it) != first_hs_o1.end()){
						first_hs_o1.push_back(*it);
						std::cout << "found an H on O1" << std::endl;
					}
				}
			}

			bound_atoms = graph.BondedAtoms (so2s.O2(), bondgraph::hbond);

			// add in any unaccounted-for atoms that have bound
			if (!bound_atoms.empty()) {
				for (Atom_it it = bound_atoms.begin(); it != bound_atoms.end(); ++it) {
					if (std::find (first_hs_o2.begin(), first_hs_o2.end(), *it) != first_hs_o2.end()) {
						first_hs_o2.push_back(*it);
						std::cout << "found an H on O2" << std::endl;
					}
				}
			}

			if (first_hs_o1.size() >= 2 && first_hs_o2.size() >= 2 && first_os.size() >= 2) {
				first_pass = false;
				printf ("\n***Rewinding***\n");
				this->_system->Rewind();
			}

		} else if (!first_pass) {

			// now after going back in time, grab the position of the bound atoms and find their positions relative to the ... Sulfur
			// first get the distance vector to the oxygen
			so2s.SO2()->SetOrderAxes();
			// we need the DCM from lab frame to SO2 frame
			MatR dcm = so2s.SO2()->DCMToLab();

			// print out the distance of the so2 to the surface
			double distance = system_t::Position(so2s.S()) - h2os.SurfaceLocation(); // top surface
			fprintf (this->output, "% 12.5f ", distance);

			AnglePrintout (first_os.begin(), first_os.end(), so2s.S(), dcm);
			AnglePrintout (first_hs_o1.begin(), first_hs_o1.end(), so2s.O1(), dcm);
			AnglePrintout (first_hs_o2.begin(), first_hs_o2.end(), so2s.O2(), dcm);

			fprintf (this->output, "\n");
		}

	} // analysis

	void SO2NearestNeighborAnalysis::AnglePrintout (Atom_it first, Atom_it last, AtomPtr ref, const MatR& dcm) const {
		VecR bound_atom;
		for (Atom_it it = first; it != last; ++it) {
			bound_atom = MDSystem::Distance(ref, *it);
			// then find the position of the atom in the so2 frame
			bound_atom = dcm.transpose() * bound_atom;

			// now get the spherical coordinates of the atom
			double coords[3] = {bound_atom[x],bound_atom[y],bound_atom[z]};
			coordinate_conversion::CartesianToSpherical (coords);
			// now print the 3 spherical coords
			fprintf (this->output, "% 12.5f % 12.5f % 12.5f ", coords[0], coords[1]*180.0/M_PI, coords[2]*180.0/M_PI);
		}
	}



	void SO2BondingAnalysis::Analysis () {
		FindSO2();
		BuildBondingGraph();

		bonds.clear();
		FindBonds(so2->S());
		FindBonds(so2->O1());
		FindBonds(so2->O2());

		PrintBondingInformation();
		return;
	}

	void SO2BondingAnalysis::FindSO2 () {
		this->LoadAll();
		MolPtr mol = Molecule::FindByType (this->begin_mols(), this->end_mols(), Molecule::SO2);
		so2 = static_cast<SulfurDioxide *>(mol);
		so2->SetAtoms();
	}

	void SO2BondingAnalysis::BuildBondingGraph () {
		graph.UpdateGraph(this->begin(), this->end());
	}

	void SO2BondingAnalysis::FindBonds (const AtomPtr& atom) {
		// get the atoms bound to the reference atom
		int atom_bonds;
		if (atom->Element() == Atom::S) {
			Atom_ptr_vec inters = graph.InteractingAtoms (atom);
			std::for_each (inters.begin(), inters.end(), std::mem_fun(&Atom::Print));
			atom_bonds = inters.size();
		}

		if (atom->Element() == Atom::O)
			atom_bonds = graph.NumHBonds (atom);

		bonds.push_back(atom_bonds);
	}

	void SO2BondingAnalysis::PrintBondingInformation () {
		for (int i = 0; i < (int)bonds.size(); i++) {
			printf ("%d  ", bonds[i]);
		}
		printf ("\n");
		//std::copy (bonds.begin(), bonds.end(), std::ostream_iterator<int>(std::cout));
	}

}
