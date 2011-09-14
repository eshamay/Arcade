#ifndef NEIGHBOR_ANALYSIS_H_
#define NEIGHBOR_ANALYSIS_H_

#include "analysis.h"
#include "manipulators.h"
#include "histogram-analysis.h"
//#include "h2o-analysis.h"
//#include "so2-system-analysis.h"
#include "bondgraph.h"

//#include <boost/iterator/iterator_facade.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

namespace neighbor_analysis {

	using namespace md_analysis;


	//************* finds the nearest neighbors to each of the SO2 atoms **************/
	//	This analysis outputs the hbonding coordination (e.g. single-acceptor, double-donor, AD, AAD, etc), as well as any information about the cycles formed within neighboring waters around the so2
	class SO2BondingCycleAnalysis : public AnalysisSet {
		public:
			typedef Analyzer system_t;

			SO2BondingCycleAnalysis (system_t * t) :
				AnalysisSet(t,
						std::string ("SO2 cycle bonding analysis"),
						std::string ("cycle-bonding.dat")),
				cm(t),
				so2s (t), h2os(t)
		{ h2os.ReferencePoint(WaterSystem::SystemParameterLookup("analysis.reference-location")); }

			~SO2BondingCycleAnalysis() { }

			void Analysis();

			typedef typename CycleManipulator::cycle_t cycle_t;

		protected:
			CycleManipulator								cm;
			so2_analysis::SO2SystemManipulator	so2s;
			h2o_analysis::H2OSystemManipulator	h2os;

	};	// cycle bonding analysis







	//************* Finds the amount of H-bonding going on with the SO2 molecule ****************/
	class SO2HBondingAnalysis : public AnalysisSet {
		public:
			typedef Analyzer system_t;

			SO2HBondingAnalysis (system_t * t) :
				AnalysisSet(t,
						std::string ("SO2 H-bonding analysis"),
						std::string ("")),
				so2s(t), h2os(t),
				histo (std::string("hbonding.dat"), WaterSystem::posmin, WaterSystem::posmax, Analyzer::posres) { 
					h2os.ReferencePoint(WaterSystem::SystemParameterLookup("analysis.reference-location")); 
				}

			~SO2HBondingAnalysis() { }

			void Analysis();
			void DataOutput () { histo.OutputData(); }

		protected:
			bondgraph::BondGraph									graph;
			Atom_ptr_vec													nearest, bonded;
			so2_analysis::SO2SystemManipulator	so2s;
			h2o_analysis::H2OSystemManipulator	h2os;
			Histogram1DAgent histo;


	};	// h-bonding analysis




	//******************** Find out some things about so2's nearest neighbor **************************/


	class SO2NearestNeighborAnalysis : public AnalysisSet {
		public:
			typedef Analyzer system_t;

			SO2NearestNeighborAnalysis (system_t * t) :
				AnalysisSet(t,
						std::string ("Track SO2's nearest neighbors"),
						std::string ("so2-neighbors.dat")),
				neighbors (t),
				so2s (t), h2os (t),
				first_pass (true) {
					h2os.ReferencePoint(WaterSystem::SystemParameterLookup("analysis.reference-location")); 
				}

			~SO2NearestNeighborAnalysis() { }
			void AnglePrintout (Atom_it first, Atom_it last, AtomPtr ref, const MatR& dcm) const;

			void Analysis();

		protected:
			NeighborManipulator	neighbors;
			so2_analysis::SO2SystemManipulator	so2s;
			h2o_analysis::H2OSystemManipulator	h2os;

			Atom_ptr_vec nearest_atoms;
			Atom_ptr_vec first_os;
			Atom_ptr_vec first_hs_o1;
			Atom_ptr_vec first_hs_o2;
			Mol_ptr_vec	nearest_neighbors;
			bondgraph::BondGraph	graph;

			bool first_pass;

	};	// nearest neighbor analysis




	/*
	 * Find what atoms are bound to the S, and to the Os.
	 */
	class SO2BondingAnalysis : public AnalysisSet {
		public:
			typedef Analyzer system_t;

			SO2BondingAnalysis (system_t * t) :
				AnalysisSet(t,
						std::string ("Find out what atoms are bound the the SO2"),
						std::string ("so2-bonding.dat")) {
				}

			void Analysis ();
			void FindSO2 ();
			void BuildBondingGraph ();
			void FindBonds (const AtomPtr&);
			void PrintBondingInformation ();

		private:
			SulfurDioxide * so2;	// reference molecule
			bondgraph::BondGraph		graph;	// connectivity/bonding graph
			std::vector<int>				bonds;
	}; // so2 bonding analysis


}	// namespace neighbor analysis

#endif
