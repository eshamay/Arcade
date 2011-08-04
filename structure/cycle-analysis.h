#ifndef CYCLE_ANALYSIS_H_
#define CYCLE_ANALYSIS_H_

#include "bond-analysis.h"
#include <boost/utility.hpp>

namespace cycle_analysis {

	using namespace md_analysis;

	class SO2CycleAnalyzer : public bond_analysis::SO2BondingAnalyzer {
		protected:

			// checks that the cycle is right, and updates the type of cycle encountered
			void CheckCycles ();

			// tests if the cycle is the type we're interested in
			virtual bool CycleCheck ( CycleManipulator::cycle_type_it cycle_type, CycleManipulator::cycle_list_it cycle) = 0;
			// does what we need to do once we have a cycle we're interested in
			virtual void CycleCheckAction (CycleManipulator::cycle_list_it cycle) = 0;
			virtual void ACycleCheckAction (CycleManipulator::cycle_list_it cycle) = 0;

			CycleManipulator cm;

		public:
			SO2CycleAnalyzer (Analyzer * t, std::string desc, std::string fn) :
				SO2BondingAnalyzer (t, desc, fn),
				cm(t) { }

			virtual void Analysis () = 0;
			virtual void DataOutput () = 0;


			// different cycle checks
			
			// Check if the cycle is minimally SO coordinated, and the cycles spans through one of the so2-SO bonds
			bool IsSO_Cycle ( CycleManipulator::cycle_type_it cycle_type, CycleManipulator::cycle_list_it cycle);
	};



	// Find the percentage of SO-coordination structures that are also cyclic from the S to the O through a single water, or through 2 waters.
	class SO2CycleCoordinationAnalyzer : public SO2CycleAnalyzer {
		protected:
			// the number of each type of cycle (ie cycles with a single water, two waters, and three waters)
			// and also the breakdown by type of the triple-water cycles:
			// Type 1 triple-water: each water contribues an O, an OH, and HOH to the cycle, respectively.
			// Type 2 triple-water: each water contributes an OH to the cycle.
			std::vector<int>	cycle_counts;
			int triple_type_B_cycles, triple_type_A_cycles;	// 2 different types of triple-cycles
			int total;	// total timesteps we hit an SO coordination

			// CountMoleculeAtoms tells us the min and max # of water atoms involved
			// so a min of 1 would mean, most likely, that one of the water only has its oxygen involved
			// whereas the max of 3 would mean that all 3 water atoms are involved in the cycle
			std::pair<int,int> CountMoleculeAtoms (CycleManipulator::cycle_list_it& cycle);

			bool CycleCheck ( CycleManipulator::cycle_type_it cycle_type, CycleManipulator::cycle_list_it cycle) 
			{ return this->IsSO_Cycle(cycle_type,cycle); }

			void CycleCheckAction (CycleManipulator::cycle_list_it cycle);
			void ACycleCheckAction (CycleManipulator::cycle_list_it cycle) { return; }

		public:
			SO2CycleCoordinationAnalyzer (Analyzer * t) :
				SO2CycleAnalyzer (t, 
						std::string ("so2 cyclic coordination analyzer"),
						std::string ("so2.cyclic-SO.dat")),
				cycle_counts (12,0),
				triple_type_B_cycles(0), triple_type_A_cycles(0), total(0)
		{ }

			virtual void Analysis ();
			virtual void DataOutput ();
	};




	class SO2CycleLifespanAnalyzer : public SO2CycleAnalyzer {

		private:

			bool CycleCheck ( CycleManipulator::cycle_type_it cycle_type, CycleManipulator::cycle_list_it cycle)
			{ return this->IsSO_Cycle(cycle_type,cycle); }

			// for when a cycle forms
			void CycleCheckAction (CycleManipulator::cycle_list_it cycle) { fprintf(this->output, "1\n"); }
			// for when a cycle is not found
			void ACycleCheckAction (CycleManipulator::cycle_list_it cycle) { fprintf(this->output, "0\n"); }

		public:
			SO2CycleLifespanAnalyzer (Analyzer * t) :
				SO2CycleAnalyzer (t, 
						std::string ("so2 cycle lifespan analyzer"),
						std::string ("so2.cyclic-SO-lifespan.dat")) { }

			void Analysis ();
			void DataOutput ();

	}; // class so2 cycle lifespan analyzer



} // namespace cycle analysis
#endif
