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

	/*
	typedef enum {
		HBRIDGE,	// a single H that points the two SO2 oxygens
		HOHBRIDGE // a water molecule where the two Hs are bound to the two Os of the SO2
	} cycle_type;
	*/


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
				NeighborManipulator<T>	nm;

				bondgraph::BondGraph		graph;
				AtomPtr						ref_atom;
				int							cycle_size;

				std::list<cycle_list>							cycle;				// all the atoms comprising the cycle
				Atom_ptr_vec						unique_cycle_atoms;	// only unique atoms in the cycle
				Mol_ptr_vec							unique_cycle_mols;
				std::list<cycle_t>					cycle_type;

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
		void CycleManipulator<T>::ParseCycles () { 

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


	template <typename T>
		void CycleManipulator<T>::FindUniqueMembers (const cycle_list& _cycle) {

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

	template <typename T>
		int CycleManipulator<T>::NumUniqueWaterAtoms () const {
			int num = 0;
			for (Atom_it it = unique_cycle_atoms.begin(); it != unique_cycle_atoms.end(); it++) {
				if ((*it)->ParentMolecule()->MolType() == Molecule::H2O) {
					++num;
				}
			}
			return num;
		}


	// go through the cycle and find the location of the atom that attaches the cycle to the bridge from the so2 leg. Return the first and second occurrence.
	template <typename T>
		CycleManipulator<T>::cycle_pair_t CycleManipulator<T>::CyclePointAtom (const cycle_list& _cycle) {
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

	template <typename T>
		std::pair<int,int> CycleManipulator<T>::WaterLegInformation (const cycle_list& _cycle) {
			cycle_pair_t point_atom = CyclePointAtom(_cycle);
			int leg_cycle_location = std::distance(_cycle.begin(), point_atom.first);
			int leg_cycle_size = std::distance(point_atom.first, point_atom.second);
			return std::make_pair(leg_cycle_location, leg_cycle_size); // don't count an atom that's part of the cycle in the point
		}




	//************* finds the nearest neighbors to each of the SO2 atoms **************/
	//	This analysis outputs the hbonding coordination (e.g. single-acceptor, double-donor, AD, AAD, etc), as well as any information about the cycles formed within neighboring waters around the so2
	template <typename T>
		class SO2BondingCycleAnalysis : public AnalysisSet<T> {
			public:
				typedef Analyzer<T> system_t;

				SO2BondingCycleAnalysis (system_t * t) :
					AnalysisSet<T>(t,
							std::string ("SO2 cycle bonding analysis"),
							std::string ("cycle-bonding.dat")),
					cm(t),
					so2s (t), h2os(t)
			{ h2os.ReferencePoint(WaterSystem<T>::SystemParameterLookup("analysis.reference-location")); }

				~SO2BondingCycleAnalysis() { }

				void Analysis();

				typedef typename CycleManipulator<T>::cycle_t cycle_t;

			protected:
				CycleManipulator<T>								cm;
				h2o_analysis::H2OSystemManipulator<T>	h2os;
				so2_analysis::SO2SystemManipulator<T>	so2s;

		};	// cycle bonding analysis



	template <typename T>
		void SO2BondingCycleAnalysis<T>::Analysis() {

			// find the locatin of the surface of the water slab
			h2os.Reload();
			h2os.FindWaterSurfaceLocation(true);
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

			typename std::list<typename CycleManipulator<T>::cycle_list>::const_iterator	cycle = cm.cycle_begin();	
			typename std::list<typename CycleManipulator<T>::cycle_t>::const_iterator		cycle_type = cm.cycle_type_begin();
			// once a cycle is detected - do something
			while (cycle != cm.cycle_end() && cycle_type != cm.cycle_type_end()) {
				int number_unique_atoms_in_cycle = 0;
				int number_molecules_in_cycle = 0;

				cm.FindUniqueMembers(*cycle);

				// full point cycle
				if (*cycle_type == CycleManipulator<T>::FULLBRIDGE) {
					number_unique_atoms_in_cycle = cm.NumUniqueAtomsInCycle() - 3;	// don't count the so2 atoms
					number_molecules_in_cycle = cm.NumUniqueMoleculesInCycle() - 1;		// don't count the SO2
				}
				// half point cycle
				else if (*cycle_type == CycleManipulator<T>::HALFBRIDGE) {
					number_unique_atoms_in_cycle = cm.NumUniqueAtomsInCycle() - 2;	// don't count the so2 atoms
					number_molecules_in_cycle = cm.NumUniqueMoleculesInCycle() - 1;		// don't count the SO2
				}
				// full crown cycle
				else if (*cycle_type == CycleManipulator<T>::FULLCROWN) {
					number_unique_atoms_in_cycle = cm.NumUniqueAtomsInCycle() - 1;	// don't count the so2 atoms
					number_molecules_in_cycle = cm.NumUniqueMoleculesInCycle() - 1;		// don't count the SO2
				}
				// water leg cycle
				else if (*cycle_type == CycleManipulator<T>::WATERLEG || *cycle_type == CycleManipulator<T>::HALFCROWN) {
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
