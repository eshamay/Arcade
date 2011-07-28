#ifndef CYCLE_ANALYSIS_H_
#define CYCLE_ANALYSIS_H_

#include "bond-analysis.h"

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
			SO2CycleAnalyzer (Analyzer * t, std::string& desc, std::string& fn) :
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
			int single_cycles, double_cycles, type_1_triple_cycles, type_2_triple_cycles, triple_cycles;
			int total;

			// CountMoleculeAtoms tells us the min and max # of water atoms involved
			// so a min of 1 would mean, most likely, that one of the water only has its oxygen involved
			// whereas the max of 3 would mean that all 3 water atoms are involved in the cycle
			std::pair<int,int> CountMoleculeAtoms (CycleManipulator::cycle_list_it& cycle);


			void CycleCheckAction (CycleManipulator::cycle_list_it cycle);
			void ACycleCheckAction (CycleManipulator::cycle_list_it cycle) { return; }

		public:
			SO2CycleCoordinationAnalyzer (Analyzer * t) :
				SO2CycleAnalyzer (t, 
						std::string ("so2 cyclic coordination analyzer"),
						std::string ("so2.cyclic-SO.dat")),
				single_cycles(0), double_cycles(0), 
				type_1_triple_cycles(0), type_2_triple_cycles(0), triple_cycles(0), total(0)
		{ }

			virtual void Analysis ();
			virtual void DataOutput ();
	};




	class SO2CycleLifespanAnalyzer : public SO2CycleAnalyzer {

		private:
			// cycle state takes into account the timeout/debouncing. Step state only looks at the current timestep
			bool cycle_state, prev_step_state, step_state;	// true == in a cycle

			// the maximum length of timesteps that a cycle has to be "broken" or "formed" for it to officially be considered "broken" or "formed"
			bool inTimeout;
			int timeout_counter;
			const int max_timeout;

			int lifespan_counter;						// how long the current cycle has been formed
			std::list<int> lifespans;			// the running list of cycle lifespans.
			std::list<int> breakspans;			// the running list of cycle breakspans.
			// once a cycle forms, the cycle_lifespan increments with each timestep, until
			// the cycle breaks. At that points the cycle_lifespan is stored in the lifespans
			// list.

			bool CycleCheck ( CycleManipulator::cycle_type_it cycle_type, CycleManipulator::cycle_list_it cycle)
			{ return this->IsSOCycle(cycle_type,cycle); }

			// for when a cycle forms
			void CycleCheckAction (CycleManipulator::cycle_list_it cycle) { step_state = true; }
			// for when a cycle is not found
			void ACycleCheckAction (CycleManipulator::cycle_list_it cycle) { step_state = false; }

		public:
			SO2CycleLifespanAnalyzer (Analyzer * t) :
				SO2CycleAnalyzer (t, 
						std::string ("so2 cycle lifespan analyzer"),
						std::string ("so2.cyclic-SO-lifespan.dat")),
				cycle_state (false), prev_step_state (false), step_state (false),
				inTimeout (false), timeout_counter (0), max_timeout (20),
				lifespan_counter (0) { }

			void Analysis ();
			void DataOutput ();

	}; // class so2 cycle lifespan analyzer



} // namespace cycle analysis
#endif
