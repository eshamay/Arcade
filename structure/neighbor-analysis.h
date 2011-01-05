#ifndef NEIGHBOR_ANALYSIS_H_
#define NEIGHBOR_ANALYSIS_H_

#include "analysis.h"
#include "histogram-analysis.h"
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

	/*
		 class atom_element_iterator : 
		 public boost::iterator_facade <
		 Atom_it,
		 AtomPtr,
		 boost::forward_traversal_tag
		 >
		 {
		 public:

		 atom_element_iterator (const Atom_it atom, const Atom::Element_t elmt) : _atom(atom), _elmt(elmt) { 
		 if (_atom->Element() != _elmt)
		 this->increment();
		 }

		 void increment() { while ( (*(++_atom))->Element() != _elmt) { } }

		 bool equal(atom_element_iterator const& other) const {
		 return this->_atom == other._atom;
		 }

		 AtomPtr& dereference() const { return *_atom; }

		 protected:
		 friend class boost::iterator_core_access;
		 Atom_it					_atom;
		 Atom::Element_t	_elmt;
		 };
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

				void next_closest (Atom_it& it) {
					it++;
				}

			private:
				AtomPtr	reference_atom;

		}; // neighbor manipulator




	//************* finds the nearest neighbors to each of the SO2 atoms **************/
	template <typename T>
		class SO2BondingCycleAnalysis : public AnalysisSet<T> {
			public:
				typedef Analyzer<T> system_t;

				SO2BondingCycleAnalysis (system_t * t) :
					AnalysisSet<T>(t,
							std::string ("SO2 nearest neighbor analysis"),
							std::string ("")),
					nm(t),
					so2s (t)
		 	{ }

				//~SO2BondingCycleAnalysis() { }

				void Analysis() {

					// order the atoms in the system in closest to furthest from a particular reference point
					AtomPtr ref_atom = so2s.S();
					nm.OrderAtomsByDistance(ref_atom);

					// grab only a certain number of the closest atoms for analysis
					Atom_ptr_vec atoms;
					atoms.push_back(ref_atom);
					int num_close_atoms = 10;
					std::copy(nm.closest(), nm.closest() + num_close_atoms, std::back_inserter(atoms));

					// build a graph of the selected atoms and their bonding/connectivity
					graph.UpdateGraph(atoms);

					// determine if within the graph there are any closed cycles
					typedef std::vector<bondgraph::BondGraph::Vertex> pred_vec;
					pred_vec preds (num_close_atoms);
					bool has_cycle = false;
					//bondgraph::cycle_detector<pred_vec> vis(has_cycle, preds); 
					AtomPtr source, target;
					bondgraph::bfs_atom_visitor vis(atoms, has_cycle, source, target);
					boost::breadth_first_search(graph.Graph(), *graph._FindVertex(ref_atom), visitor(vis));


					// once a cycle is detected - do something
					if (has_cycle) {

						typedef std::list<AtomPtr> apl;
						apl cycle;
						for (AtomPtr a = target; a != so2s.S(); a = graph.Parent(a)) {
							cycle.push_front(a);
						}
						for (AtomPtr a = source; a != so2s.S(); a = graph.Parent(a)) {
							cycle.push_back(a);
						}
						// cycle involving both so2-Os
						if (cycle.front() != cycle.back() || cycle.size() == 3) {


							// find the number of unique waters involved in the cycle
							std::vector<int> mol_ids;
							std::vector<int> h2o_ids;
							printf ("\n");
							for (apl::const_iterator it = cycle.begin(); it != cycle.end(); it++) {
								(*it)->Print();
								mol_ids.push_back((*it)->MolID());
								h2o_ids.push_back((*it)->ID());
							}
							std::sort(mol_ids.begin(), mol_ids.end());
							std::sort(h2o_ids.begin(), h2o_ids.end());
							std::vector<int>::iterator mol_ids_it = std::unique(mol_ids.begin(), mol_ids.end());
							std::vector<int>::iterator h2o_ids_it = std::unique(h2o_ids.begin(), h2o_ids.end());
							mol_ids.resize(mol_ids_it - mol_ids.begin());
							h2o_ids.resize(h2o_ids_it - h2o_ids.begin());
							// print the number of h2o atoms involved
							printf ("\ncycle found : %zu  %zu\n", h2o_ids.size() - 2, mol_ids.size() - 1);
						}

					}
				} // analysis


			protected:
				NeighborManipulator<T>	nm;
				so2_analysis::SO2SystemManipulator<T>	so2s;
				bondgraph::BondGraph	graph;

		};	// test analysis

}	// namespace md analysis

#endif
