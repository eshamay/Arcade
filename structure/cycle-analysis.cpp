#include "cycle-analysis.h"

namespace cycle_analysis {

	// given a cycle iterator (from a cycle manipulator), this tells if the cycle is an SO-half bridge type
	// with the cycle spanning from the so2-S to an oxygen.
	bool SO2CycleAnalyzer::IsSO_Cycle ( CycleManipulator::cycle_type_it cycle_type, CycleManipulator::cycle_list_it cycle) {

		bool ret = false;

		if (*cycle_type == CycleManipulator::HALFBRIDGE) {

			// figure out what atoms are at the start and end of the cycle.
			CycleManipulator::cycle_it _first = cycle->begin(); _first++;	
			AtomPtr atom1 = *_first;
			AtomPtr atom2 = cycle->back();

			// figure out what molecules the starting and ending atoms of the cycle are in. 
			// This tells us about the cycle type.
			Molecule::Molecule_t mol1_t, mol2_t;
			mol1_t = atom1->ParentMolecule()->MolType(); 
			mol2_t = atom2->ParentMolecule()->MolType(); 

			// the half-bridge starts on the so2-S, and ends with a water-O, thus we check here the molecule types
			if ((mol1_t == Molecule::SO2 && mol2_t != Molecule::SO2) 
					|| (mol1_t != Molecule::SO2 && mol2_t == Molecule::SO2)) 
			{ ret = true; }
		}

		return ret;
	}


	void SO2CycleAnalyzer::CheckCycles () {

		// this points to the cycle type of each cycle found.
		CycleManipulator::cycle_type_it cycle_type = cm.cycle_type_begin();

		// snoop through each cycle found (probably not too many)
		for (CycleManipulator::cycle_list_it cycle = cm.cycle_begin(); cycle != cm.cycle_end(); cycle++) {

			// check only the cycles that are of the half-bridge variety that Baer proposed
			if (this->CycleCheck(cycle_type, cycle)) {
				// do what needs to be done if a cycle is found
				this->CycleCheckAction (cycle);
			} 
			else {
				this->ACycleCheckAction (cycle);
			}

			cycle_type++;
		} // for each cycle

		return;
	}




	void SO2CycleCoordinationAnalyzer::CycleCheckAction (CycleManipulator::cycle_list_it cycle) {
		this->cm.FindUniqueMembers(*cycle);

		int num_mol = this->cm.NumUniqueMoleculesInCycle(); 
		int num_atoms = this->cm.NumUniqueAtomsInCycle(); 

		// single cycles (i.e. only 1 water involved) only have 2 molecules involved... the so2 and the water
		if (num_mol == 2 && num_atoms == 4)
			++single_cycles;

		else if (num_mol == 3)
			++double_cycles;

		// found a 3-water cycle
		else if (num_mol == 4) {
			++triple_cycles;
			// 3-waters + so2 should make for 8 atoms in the cycle (because the cycle has to be (starting on the S):
			//		S-O-H-O-H-O-H-O
			// This can be checked be doing triple_cycles - type_1 - type_2. If the answer isn't 0 then there is some other
			// type of triple cycle forming! So far, that isn't the case..
			if (num_atoms == 8) {
				std::pair<int,int> minmax = CountMoleculeAtoms(cycle);
				if (minmax.first == 1 && minmax.second == 3)
					++type_1_triple_cycles;
				else if (minmax.first == 2 && minmax.second == 2)
					++type_2_triple_cycles;
				else {
					std::cerr << "check it -- " << minmax.first << "," << minmax.second << "  !" << std::endl;
					std::for_each (cycle->begin(), cycle->end(), std::mem_fun(&Atom::Print));
					exit(1);
				}
			}
		} // 3-water cycles

		// some error check
		else if (num_mol <= 4) {
			std::cerr << "hey boo!" << std::endl;
			std::cerr << num_mol << " mols and " << num_atoms << " atoms" << std::endl;
			std::for_each (cycle->begin(), cycle->end(), std::mem_fun(&Atom::Print));
			exit(1);
		}
		return;
	}




	// goes through the atoms in the cycle and counts the number of atoms in each molecule (by molecule ID).
	// After sorting those out, we find the number of atoms contributing to the cycle from each molecule,
	// and the max/min number of atoms contributed. i.e. one water contributes 1 atom (an oxygen) but another
	// water contributes 3 atoms (2x H and 1x O) to the cycle. Thus we have a min/max of 1/3.
	std::pair<int,int> SO2CycleCoordinationAnalyzer::CountMoleculeAtoms (CycleManipulator::cycle_list_it& cycle) {

		typedef std::vector<int> int_vec;
		typedef int_vec::iterator int_it;

		int_vec mol_ids;
		std::transform (cycle->begin(), cycle->end(), std::back_inserter(mol_ids), std::mem_fun<int>(&Atom::MolID));
		std::sort(mol_ids.begin(), mol_ids.end());

		int_vec unique_ids;
		std::copy(mol_ids.begin(), mol_ids.end(), std::back_inserter(unique_ids));
		int_it it = std::unique(unique_ids.begin(), unique_ids.end());
		unique_ids.resize(it - unique_ids.begin());

		int_vec count;

		for (int_it jt = unique_ids.begin(); jt != unique_ids.end(); jt++) {
			count.push_back((int)std::count(mol_ids.begin(), mol_ids.end(), *jt));
		}

		int_it min = std::min_element(count.begin(), count.end());
		int_it max = std::max_element(count.begin(), count.end());
		return std::make_pair (*min,*max);
	}


	void SO2CycleCoordinationAnalyzer::Analysis () {
		this->FindCoordination ();

		cm.SetReferenceAtom(so2->S());
		cm.SetCycleSize (30);
		cm.BuildGraph();
		cm.ParseCycles();

		// check that we're in the right type of coordination for the so2
		if (this->coordination / 10 >= 1 && this->coordination % 10 >= 1) {
			++total;
			this->CheckCycles ();
		}
	}

	// prints out a single line of data with the following values:
	void SO2CycleCoordinationAnalyzer::DataOutput () {
		rewind (this->output);
		fprintf (this->output, "%d %d %d %d %d %d\n", 
				total, single_cycles, double_cycles, triple_cycles, type_1_triple_cycles, type_2_triple_cycles);
		fflush (this->output);
	}




	void SO2CycleLifespanAnalyzer::Analysis () {

		this->FindCoordination ();

		cm.SetReferenceAtom(so2->S());
		cm.SetCycleSize (30);
		cm.BuildGraph();
		cm.ParseCycles();

		// update the previous step's state, and then find the current step;
		prev_step_state = step_state;

		// if the so2 is minimally an SO coordination (at least 1 S and 1 O bond)
		//if (this->coordination / 10 >= 1 && this->coordination % 10 >= 1)
		this->CheckCycles ();

		// if we've formed or broken a cycle, reset the timeout
		if (prev_step_state != step_state) {
			timeout_counter = 0;
			inTimeout = true;
		}

		// during a timeout/debounce
		else if (inTimeout) {
			++timeout_counter;

			if (timeout_counter >= max_timeout) {
				// state didn't change after a timeout?
				if (step_state == cycle_state) {
					timeout_counter = 0;
				}

				else {
					// keep track of both the lifespans of cycles, and the lifespans of the breaks
					if (cycle_state)
						lifespans.push_back (lifespan_counter);
					else
						breakspans.push_back (lifespan_counter);
					lifespan_counter = 0;
					cycle_state = step_state;
				}
				inTimeout = false;
			}

		}

		++lifespan_counter;

		return;
	}

	// print out the list of lifespans and breakspans encountered
	void SO2CycleLifespanAnalyzer::DataOutput () {
		std::list<int>::iterator ls,ls_end, bs,bs_end;
		boost::tie(ls,ls_end) = std::make_pair(lifespans.begin(), lifespans.end());
		boost::tie(bs,bs_end) = std::make_pair(breakspans.begin(), breakspans.end());

		while (ls != ls_end) {
			//fprintf (this->output, "%8d\n", *ls);
			ls++;
		}

		fflush(this->output);
		lifespans.clear();
		breakspans.clear();
		return;
	}

} // namespace cycle analysis
