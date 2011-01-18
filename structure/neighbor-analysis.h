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
					cycle_size(0), ref_atom((AtomPtr)NULL), cycle_type(NO_CYCLE) { }


				typedef enum {
					NO_CYCLE = 0,		// no cycle present in the local graph of the so2 neighbors
					FULL = 1,		// A complete cycle links the two so2 oxygens
					PARTIAL = 2,			// a cycle exists but does not link the two so2 oxygens
					UNKNOWN = 3			// some other type that isn't figured out right now
				}	cycle_t;

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

				typedef std::list<AtomPtr> cycle_list;
				typedef cycle_list::const_iterator cycle_it;

				cycle_it	begin() const { return cycle.begin(); }
				cycle_it	end() const { return cycle.end(); }
				AtomPtr	front() const { return cycle.front(); }
				AtomPtr	back() const { return cycle.back(); }

				cycle_t CycleType() const { return cycle_type; }
				size_t size () const { return cycle.size(); }

				size_t NumUniqueAtomsInCycle () const { return unique_cycle_atoms.size(); }
				size_t NumUniqueMoleculesInCycle () const { return unique_cycle_mols.size(); }
				int NumUniqueWaterAtoms () const;

				int NumReferenceAtomHBonds () const { return graph.NumHBonds(ref_atom); }

			private:
				NeighborManipulator<T>	nm;

				bondgraph::BondGraph		graph;
				AtomPtr									ref_atom;
				int											cycle_size;
				bool										has_cycle;
				AtomPtr									gray_source, gray_target;

				cycle_list							cycle;				// all the atoms comprising the cycle
				Atom_ptr_vec						unique_cycle_atoms;	// only unique atoms in the cycle
				Mol_ptr_vec							unique_cycle_mols;
				cycle_t									cycle_type;

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
			cycle_type = NO_CYCLE;
			if (has_cycle) {
			  cycle.clear();

			  // fill the cycle's atom list here by first getting all the atoms from the target to the ref-atom, and then from the source so that both sides of the cycle are accounted for.
			  // cycle list will look like S -- O1 -- etc. -- O2
			  for (AtomPtr a = gray_target; a != ref_atom; a = graph.Parent(a)) {
				cycle.push_front(a);
			  }
			  cycle.push_front(ref_atom);
			  for (AtomPtr a = gray_source; a != ref_atom; a = graph.Parent(a)) {
				cycle.push_back(a);
			  }

			  // the first atom in the cycle is the reference. Here we get the atom connected to the reference - the first non-ref atom - and the last atom in the cycle
			  // based on the molecules these are connected to, we discover the type of cycle they're involved in
			  cycle_it _first = begin(); _first++;	
			  AtomPtr atom1 = *_first;
			  AtomPtr atom2 = back();
			  if (atom1 != atom2) {
				Molecule::Molecule_t atom1_t, atom2_t;
				atom1_t = atom1->ParentMolecule()->MolType(); 
				atom2_t = atom2->ParentMolecule()->MolType(); 
				if ((atom1_t == Molecule::H2O && atom2_t == Molecule::SO2) ||
					(atom1_t == Molecule::SO2 && atom2_t == Molecule::H2O)) {
				  cycle_type = HALF;
				}
				else if (atom1_t == Molecule::SO2 && atom1_t == Molecule::SO2) { 
				  cycle_type = FULL; 
				}
				else {
				  printf ("found some other type of funky cycle\n");
				  std::for_each (cycle.begin(), cycle.end(), std::mem_fun(&Atom::Print));
				}
			  }
			  else if (atom1 == atom2) { 
				cycle_type = PARTIAL; 
			  }
			  else { 
				cycle_type = UNKNOWN; 
			  }
			}
		}	// find cycle


	template <typename T>
	  void CycleManipulator<T>::FindUniqueMembers () {

		// find all the unique atoms in the cycle
		unique_cycle_atoms.clear();
		std::copy (cycle.begin(), cycle.end(), std::back_inserter(unique_cycle_atoms));
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
		int number_water_atoms_in_cycle = 0;
		int number_waters_in_cycle = 0;

		// once a cycle is detected - do something
		if (cm.HasCycle()) {

		  cm.FindUniqueMembers();

		  // cycle involving both so2-Os
		  if (cm.CycleType() == CycleManipulator<T>::FULL || cm.CycleType() == CycleManipulator<T>::PARTIAL) {
			number_water_atoms_in_cycle = cm.NumUniqueWaterAtoms();
			number_waters_in_cycle = cm.NumUniqueMoleculesInCycle() - 1;		// don't count the SO2
		  }

		  //for (typename CycleManipulator<T>::cycle_it it = cm.begin(); it != cm.end(); it++) {
		  //(*it)->Print();
		  //}
		  //printf ("%d %d %zu\n", (int)cm.CycleType(), cm.NumUniqueWaterAtoms(), cm.NumUniqueMoleculesInCycle() - 1);
		}

		// final print out of the cycle information
		fprintf (this->output, " %5d %5d %5d \n", 
			(int)cm.CycleType(), number_water_atoms_in_cycle,	number_waters_in_cycle);
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
