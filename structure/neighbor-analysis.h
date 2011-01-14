#ifndef NEIGHBOR_ANALYSIS_H_
#define NEIGHBOR_ANALYSIS_H_

#include "analysis.h"
#include "histogram-analysis.h"
#include "h2o-analysis.h"
#include "so2-system-analysis.h"
#include "bondgraph.h"

//#include <boost/iterator/iterator_facade.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

namespace neighbor_analysis {

	using namespace md_analysis;

	typedef enum {
		HBRIDGE,	// a single H that bridges the two SO2 oxygens
		HOHBRIDGE // a water molecule where the two Hs are bound to the two Os of the SO2
	} cycle_type;


	/******************* Neighbor Manipulator ******************/
	template <typename T>
		class NeighborManipulator : public SystemManipulator<T> {
			public:

				typedef Analyzer<T> system_t;

				NeighborManipulator (system_t * sys) : SystemManipulator<T>(sys) { }
				~NeighborManipulator() { }

				// Order all the analysis atoms by distance starting with the closest to the given atom
				void OrderAtomsByDistance (AtomPtr ap) {
					this->Reload();
					reference_atom = ap;
					std::sort(
							this->analysis_atoms.begin(), 
							this->analysis_atoms.end(), 
							system_t::atomic_reference_distance_pred (ap));
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

				void next_closest (Atom_it& it) { it++; }

			private:
				AtomPtr									reference_atom;

		}; // neighbor manipulator





	/************* Cycle Manipulator *****************/

	template <typename T>
		class CycleManipulator : public SystemManipulator<T> {
			public:
				typedef Analyzer<T> system_t;

				CycleManipulator (system_t * t) : 
					SystemManipulator<T>(t),
					nm(t),
					cycle_size(0), ref_atom((AtomPtr)NULL) { }


				void BuildGraph ();	// builds the graph of the given atoms closest to the reference atom
				void FindCycle ();	// finds if a cycle exists in the topology, and parses it
				void FindUniqueMembers ();

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


				bool HasCycle () const { return has_cycle; }
				AtomPtr Source () const { return gray_source; }
				AtomPtr Target () const { return gray_target; }

				typedef std::list<AtomPtr> cycle_t;
				typedef cycle_t::const_iterator cycle_it;

				cycle_it	begin() const { return cycle.begin(); }
				cycle_it	end() const { return cycle.end(); }
				AtomPtr	front() const { return cycle.front(); }
				AtomPtr	back() const { return cycle.back(); }

				typedef std::vector<int> int_vec;

				int_vec::const_iterator begin_mol_ids () const { return mol_ids.begin(); }
				int_vec::const_iterator end_mol_ids () const { return mol_ids.end(); }
				size_t NumUniqueMoleculesInCycle () const { return mol_ids.size(); }

				int_vec::const_iterator begin_atom_ids () const { return atom_ids.begin(); }
				int_vec::const_iterator end_atom_ids () const { return atom_ids.end(); }
				size_t NumUniqueAtomsInCycle () const { return atom_ids.size(); }
				int NumUniqueWaterAtoms () const;

				size_t size () const { return cycle.size(); }

				int NumReferenceAtomHBonds () const { return graph.NumHBonds(ref_atom); }


			private:
				NeighborManipulator<T>	nm;

				bondgraph::BondGraph		graph;
				AtomPtr									ref_atom;
				int_vec									mol_ids, atom_ids;
				int											cycle_size;
				bool										has_cycle;
				AtomPtr									gray_source, gray_target;

				cycle_t									cycle;

		}; // cycle manipulator

	template <typename T>
		void CycleManipulator<T>::BuildGraph () { 

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



	template <typename T>
		void CycleManipulator<T>::FindCycle () { 
			// determine if within the graph there are any closed cycles
			has_cycle = false;
			bondgraph::bfs_atom_visitor vis(has_cycle, gray_source, gray_target);
			boost::breadth_first_search(graph.Graph(), *graph._FindVertex(ref_atom), visitor(vis));

			// grab the atoms that make up the cycle
			if (has_cycle) {
				cycle.clear();

				// fill the cycle's atom list here by first getting all the atoms from the target to the ref-atom, and then from the source so that both sides of the cycle are accounted for.
				for (AtomPtr a = gray_target; a != ref_atom; a = graph.Parent(a)) {
					cycle.push_front(a);
				}
				cycle.push_front(ref_atom);
				for (AtomPtr a = gray_source; a != ref_atom; a = graph.Parent(a)) {
					cycle.push_back(a);
				}
			}
		}	// find cycle


	template <typename T>
		void CycleManipulator<T>::FindUniqueMembers () {

			mol_ids.clear();
			atom_ids.clear();
			// find the number of unique atoms and molecules involved in the cycle
			for (cycle_it it = begin(); it != end(); it++) {
				mol_ids.push_back((*it)->MolID());
				atom_ids.push_back((*it)->ID());
			}

			// get the unique atom and molecule ids 
			std::sort(mol_ids.begin(), mol_ids.end());
			std::sort(atom_ids.begin(), atom_ids.end());
			std::vector<int>::iterator mol_ids_it = std::unique(mol_ids.begin(), mol_ids.end());
			std::vector<int>::iterator atom_ids_it = std::unique(atom_ids.begin(), atom_ids.end());
			mol_ids.resize(mol_ids_it - mol_ids.begin());
			atom_ids.resize(atom_ids_it - atom_ids.begin());
		}

	template <typename T>
		int CycleManipulator<T>::NumUniqueWaterAtoms () const {
			int num = 0;
			for (cycle_it it = begin(); it != end(); it++) {
				if ((*it)->ParentMolecule()->MolType() == Molecule::H2O) {
					++num;
				}
			}
			return num;
		}







	//************* finds the nearest neighbors to each of the SO2 atoms **************/
	//	This analysis outputs the hbonding coordination (e.g. single-acceptor, double-donor, AD, AAD, etc), as well as any information about the cycles formed within neighboring waters around the so2
	template <typename T>
		class SO2BondingCycleAnalysis : public AnalysisSet<T> {
			public:
				typedef Analyzer<T> system_t;


				typedef enum {
					NO_CYCLE = 0,		// no cycle present in the local graph of the so2 neighbors
					COMPLETE = 1,		// A complete cycle links the two so2 oxygens
					PARTIAL = 2,			// a cycle exists but does not link the two so2 oxygens
					UNKNOWN = 3			// some other type that isn't figured out right now
				}	cycle_t;


				SO2BondingCycleAnalysis (system_t * t) :
					AnalysisSet<T>(t,
							std::string ("SO2 cycle bonding analysis"),
							std::string ("cycle-bonding.dat")),
					cm(t),
					so2s (t), h2os(t)
			{ h2os.ReferencePoint(WaterSystem<T>::SystemParameterLookup("analysis.reference-location")); }

				~SO2BondingCycleAnalysis() { }

				void Analysis();

			protected:
				CycleManipulator<T>								cm;
				h2o_analysis::H2OSystemManipulator<T>	h2os;
				so2_analysis::SO2SystemManipulator<T>	so2s;
				std::vector<int> mol_ids, h2o_ids;

		};	// cycle bonding analysis



	template <typename T>
		void SO2BondingCycleAnalysis<T>::Analysis() {

			// find the locatin of the surface of the water slab
			h2os.Reload();
			//h2os.FindWaterSurfaceLocation(true);
			h2os.FindWaterSurfaceLocation(false);	 // bottom surface
			//double distance = system_t::Position(so2s.S()) - h2os.SurfaceLocation();
			double distance = h2os.SurfaceLocation() - system_t::Position(so2s.S());	// bottom surface
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

			// select the sulfur atom to use and find any cycles within the neighboring waters in the graph
			cm.SetReferenceAtom(so2s.S());
			cm.BuildGraph();
			cm.FindCycle();
			cycle_t cycle_type = NO_CYCLE;
			int number_water_atoms_in_cycle = 0;
			int number_waters_in_cycle = 0;

			// once a cycle is detected - do something
			if (cm.HasCycle()) {

				cm.FindUniqueMembers();

				//printf ("\n");
				//for (typename CycleManipulator<T>::cycle_it it = cm.begin(); it != cm.end(); it++) {
					//(*it)->Print();
				//}

				typename CycleManipulator<T>::cycle_it begin = cm.begin();
				typename CycleManipulator<T>::cycle_it end = cm.end();
				begin++;	// second one
				end--;	// last one
				// cycle involving both so2-Os
				if (*begin != *end) {
					cycle_type = COMPLETE;
					number_water_atoms_in_cycle = cm.NumUniqueWaterAtoms();
					number_waters_in_cycle = cm.NumUniqueMoleculesInCycle() - 1;		// don't count the SO2
					//printf ("\n%d %d %d\n", cycle_type, number_water_atoms_in_cycle, number_waters_in_cycle);
				}
				else if (*begin == *end) {
					cycle_type = PARTIAL;
					number_water_atoms_in_cycle = cm.NumUniqueWaterAtoms();
					number_waters_in_cycle = cm.NumUniqueMoleculesInCycle();		// don't count the SO2
					//printf ("\n%d %d %d\n", cycle_type, number_water_atoms_in_cycle, number_waters_in_cycle);
				}
				else { cycle_type = UNKNOWN; }
			}

			// final print out of the cycle information
			fprintf (this->output, " %5d %5d %5d \n", 
					(int)cycle_type, number_water_atoms_in_cycle,	number_waters_in_cycle);
			fflush(this->output);
		} // analysis



	//************* Finds the amount of H-bonding going on with the SO2 molecule ****************/
	template <typename T>
		class SO2HBondingAnalysis : public AnalysisSet<T> {
			public:
				typedef Analyzer<T> system_t;

				SO2HBondingAnalysis (system_t * t) :
					AnalysisSet<T>(t,
							std::string ("SO2 H-bonding analysis"),
							std::string ("so2-h-bonding.dat")),
					cm(t),
					so2s (t), h2os(t)
			{ h2os.ReferencePoint(WaterSystem<T>::SystemParameterLookup("analysis.reference-location")); }

				~SO2HBondingAnalysis() { }

				void Analysis();

			protected:
				CycleManipulator<T>								cm;
				h2o_analysis::H2OSystemManipulator<T>	h2os;
				so2_analysis::SO2SystemManipulator<T>	so2s;
				std::vector<int> mol_ids, h2o_ids;

		};	// h-bonding analysis

	template<typename T>
		void SO2HBondingAnalysis<T>::Analysis () {
			// find the locatin of the surface of the water slab and the distance of the so2 to it
			h2os.Reload();
			//h2os.FindWaterSurfaceLocation(true);	// top surface
			h2os.FindWaterSurfaceLocation(false);	// bottom surface
			//double distance = system_t::Position(so2s.S()) - h2os.SurfaceLocation();	// top surface
			double distance = h2os.SurfaceLocation() - system_t::Position(so2s.S());	// bottom surface
			fprintf (this->output, " % 8.3f ", distance);

			// select the reference atom to use for determining h-bonding to the so2.
			cm.SetReferenceAtom(so2s.O1());
			cm.SetCycleSize (15);

			cm.BuildGraph();

			// output the number of bonds on the first oxygen
			fprintf (this->output, " %5d ", cm.NumReferenceAtomHBonds());

			cm.SetReferenceAtom(so2s.O2());
			cm.BuildGraph();

			// then for the 2nd oxygen
			fprintf (this->output, " %5d\n", cm.NumReferenceAtomHBonds());

		}	// so2 h bonding analysis

}	// namespace md analysis

#endif
