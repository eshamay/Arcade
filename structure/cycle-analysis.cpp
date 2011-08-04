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
		CycleManipulator::cycle_list_it cycle = cm.cycle_begin();

		// for now... just look at the smallest cycles
		
		// snoop through each cycle found (probably not too many)
		//for (CycleManipulator::cycle_list_it cycle = cm.cycle_begin(); cycle != cm.cycle_end(); cycle++) {

		// check only the cycles that are of the half-bridge variety that Baer proposed
		if (this->CycleCheck(cycle_type, cycle)) {
			// do what needs to be done if a cycle is found
			this->CycleCheckAction (cycle);
		} 
		else {
			this->ACycleCheckAction (cycle);
		}

		//cycle_type++;
		//} // for each cycle

		return;
	}




	void SO2CycleCoordinationAnalyzer::CycleCheckAction (CycleManipulator::cycle_list_it cycle) {
		this->cm.FindUniqueMembers(*cycle);

		int num_mol = this->cm.NumUniqueMoleculesInCycle(); 
		int num_atoms = this->cm.NumUniqueAtomsInCycle(); 

		if (num_mol > cycle_counts.size()) {
			printf ("num_mol = %d\n", num_mol);
			fflush(stdout);
			exit(1);
		}

		cycle_counts[num_mol]++;

		// found a 3-water cycle
		if (num_mol == 4) {
			// 3-waters + so2 should make for 8 atoms in the cycle (because the cycle has to be (starting on the S):
			//		S-O-H-O-H-O-H-O
			// This can be checked be doing triple_cycles - type_1 - type_2. If the answer isn't 0 then there is some other
			// type of triple cycle forming! So far, that isn't the case..
			if (num_atoms == 8) {
				std::pair<int,int> minmax = CountMoleculeAtoms(cycle);
				if (minmax.first == 1 && minmax.second == 3)
					++triple_type_B_cycles;
				else if (minmax.first == 2 && minmax.second == 2)
					++triple_type_A_cycles;
				else {
					std::cerr << "check it -- " << minmax.first << "," << minmax.second << "  !" << std::endl;
					std::for_each (cycle->begin(), cycle->end(), std::mem_fun(&Atom::Print));
					exit(1);
				}
			}
		} // 3-water cycles

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
		int num_o = this->coordination % 10;
		int num_s = this->coordination / 10;
		//printf ("%d %d %d\n", this->coordination, num_o, num_s);
		if (num_o >= 1 && num_s >= 1) {
			++this->total;
			this->CheckCycles ();
		}
	}

	// prints out a single line of data with the following values:
	void SO2CycleCoordinationAnalyzer::DataOutput () {
		rewind (this->output);
		// print out the total times we encountered an SO coordinated so2
		fprintf (this->output, "%d", total);

		// print out the total number of times a cycle was encountered (of any size)
		fprintf (this->output, " %d", std::accumulate(cycle_counts.begin(), cycle_counts.end(), 0));

		// print out each cycle-size count
		for (int i = 2; i < cycle_counts.size(); i++) {
			fprintf (this->output, " %d", cycle_counts[i]);
		}

		// lastly print out the breakdown of the 2 types of triple-size cycles
		fprintf (this->output, " %d %d", triple_type_A_cycles, triple_type_B_cycles);

		fflush (this->output);
	}




	void SO2CycleLifespanAnalyzer::Analysis () {

		this->FindCoordination ();

		cm.SetReferenceAtom(so2->S());
		cm.SetCycleSize (30);
		cm.BuildGraph();
		cm.ParseCycles();

		// check that we're in the right type of coordination for the so2
		int num_o = this->coordination % 10;
		int num_s = this->coordination / 10;
		//printf ("%d %d %d\n", this->coordination, num_o, num_s);
		if (num_o >= 1 && num_s >= 1) {
			// if the so2 is minimally an SO coordination (at least 1 S and 1 O bond)
			//if (this->coordination / 10 >= 1 && this->coordination % 10 >= 1)
			this->CheckCycles ();
		}

		/*
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
		*/

		return;
	}

	// print out the list of lifespans and breakspans encountered
	void SO2CycleLifespanAnalyzer::DataOutput () {
		//for (std::list<int>::iterator ls = lifespans.begin(); ls != lifespans.end(); ls++) {
		//fprintf (this->output, "%8d\n", *ls);
		//ls++;
		//}

		fflush(this->output);
		return;
	}

} // namespace cycle analysis
