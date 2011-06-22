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


	/******************* Neighbor Manipulator ******************/
	class NeighborManipulator : public SystemManipulator {
		public:

			typedef Analyzer system_t;

			NeighborManipulator (system_t * sys) : SystemManipulator(sys) { }
			~NeighborManipulator() { }

			// Order all the analysis atoms by distance starting with the closest to the given atom
			void OrderAtomsByDistance (AtomPtr ap) {
				this->Reload();
				reference_atom = ap;
				std::sort(
						this->analysis_atoms.begin(), 
						this->analysis_atoms.end(), 
						typename system_t::atomic_distance_cmp (ap));
			}

			// iterate over the atoms sorted by distance - but only over a particular element.
			Atom_it closest (const Atom::Element_t elmt) {
				Atom_it it = closest();
				if ((*it)->Element() != elmt)
					next_closest(it, elmt);
				return it;
			}

			Atom_it closest () {
				Atom_it it = this->analysis_atoms.begin();
				if (*it == reference_atom) it++;
				return it;
			}

			Atom_it end () { return this->analysis_atoms.end(); }

			// increment the iterator to the next occurrence of an atom with the specified element type
			void next_closest (Atom_it& it, const Atom::Element_t elmt) {
				it++;
				while (it != end() && (*it)->Element() != elmt) { it++; }
			}

			// retrieve the closest atoms to the given reference_atom
			Atom_range GetClosestAtoms (AtomPtr a, const int num) {
				this->OrderAtomsByDistance(a);
				return std::make_pair(this->analysis_atoms.begin(), this->analysis_atoms.begin() + num);
			}

			void next_closest (Atom_it& it) { it++; }

		private:
			AtomPtr									reference_atom;

	}; // neighbor manipulator





	/************* Cycle Manipulator *****************/

	class CycleManipulator : public SystemManipulator {
		public:
			typedef Analyzer system_t;

			CycleManipulator (system_t * t) : 
				SystemManipulator(t),
				nm(t),
				cycle_size(0), ref_atom((AtomPtr)NULL), cycle_type(NO_CYCLE) { }


			typedef enum {
				NO_CYCLE = 0,	// no cycle present in the local graph of the so2 neighbors
				FULLBRIDGE=1,		// A complete cycle links the two so2 oxygens
				WATERLEG=2,			// a cycle exists but does not link the two so2 oxygens - the cycle is in the waters off one so2-O
				HALFBRIDGE=3,			// a cycle runs from the S to a water-O and back to the so2-O
				FULLCROWN=4,			// S connected to two water-Os
				HALFCROWN=5,			// cycle happens in waters and is only peripherally connected to the S (like for a waterleg)
				UNKNOWN				// some other type that isn't figured out right now
			}	cycle_t;


			typedef std::list<AtomPtr> cycle_list;
			typedef cycle_list::const_iterator cycle_it;
			typedef std::pair<cycle_it,cycle_it> cycle_pair_t;


			void BuildGraph ();	// builds the graph of the given atoms closest to the reference atom
			void ParseCycles ();	// finds if a cycle exists in the topology, and parses it
			void FindUniqueMembers (const cycle_list&);

			void SetCycleSize (const int size) { cycle_size = size; }
			void SetReferenceAtom (const AtomPtr ref) { ref_atom = ref; }

			AtomPtr ReferenceAtom () const { 
				if (ref_atom == (AtomPtr)NULL) {
					std::cerr << "CycleManipulator -- Reference atom was never set before use!" << std::endl;
					exit(1);
				}
				return ref_atom; 
			}

			int CycleSize () const {
				if (cycle_size == 0) {
					std::cerr << "CycleManipulator -- Cycle size was never set before use!" << std::endl;
					exit(1);
				}
				return cycle_size; 
			}


			typename std::list<cycle_list>::const_iterator	cycle_begin() const { return cycle.begin(); }
			typename std::list<cycle_list>::const_iterator	cycle_end() const { return cycle.end(); }

			typename std::list<cycle_t>::const_iterator		cycle_type_begin() const { return cycle_type.begin(); }
			typename std::list<cycle_t>::const_iterator		cycle_type_end() const { return cycle_type.end(); }

			size_t NumUniqueAtomsInCycle () const { return unique_cycle_atoms.size(); }
			size_t NumUniqueMoleculesInCycle () const { return unique_cycle_mols.size(); }
			int NumUniqueWaterAtoms () const;

			int NumReferenceAtomHBonds () const { return graph.NumHBonds(ref_atom); }
			int NumReferenceInteractions () const { return graph.NumInteractions(ref_atom); }

			cycle_pair_t CyclePointAtom (const cycle_list&);		// used to find the first atom that appears both in the bridge to a leg-cycle, and in the cycle itself
			std::pair<int,int> WaterLegInformation (const cycle_list&);

		private:
			NeighborManipulator	nm;

			int							cycle_size;
			AtomPtr						ref_atom;
			std::list<cycle_t>					cycle_type;

			bondgraph::BondGraph		graph;
			std::list<cycle_list>							cycle;				// all the atoms comprising the cycle
			Atom_ptr_vec						unique_cycle_atoms;	// only unique atoms in the cycle
			Mol_ptr_vec							unique_cycle_mols;

	}; // cycle manipulator





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
