#ifndef NEIGHBOR_ANALYSIS_H_
#define NEIGHBOR_ANALYSIS_H_

#include "analysis.h"
#include "histogram-analysis.h"
#include "so2-system-analysis.h"
#include "bondgraph.h"

# include <boost/iterator/iterator_facade.hpp>

namespace neighbor_analysis {

	using namespace md_analysis;

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
		class SO2NearestNeighborAnalysis : public AnalysisSet<T> {
			public:
				typedef Analyzer<T> system_t;

				SO2NearestNeighborAnalysis (system_t * t) :
					AnalysisSet<T>(t,
							std::string ("SO2 nearest neighbor analysis"),
							std::string ("")),
					nm(t),
					s_bonds ("so2-nearest.S.dat", 2.5, 3.5, 0.05, 2.5, 3.7, 0.05),
					o_bonds ("so2-nearest.O.dat", 2.0, 3.75, 0.05, 2.5, 3.75, 0.05),
					so2s (t)
		 	{ }

				~SO2NearestNeighborAnalysis() { }

				void Analysis() {

					AtomPtr ref_atom = so2s.S();
					nm.OrderAtomsByDistance(ref_atom);
					int num_close_atoms = 5;

					Atom_ptr_vec atoms;
					atoms.push_back(ref_atom);
					std::copy(nm.closest(), nm.closest() + num_close_atoms, std::back_inserter(atoms));

					graph.UpdateGraph(atoms);

					typedef std::vector<bondgraph::BondGraph::Vertex> pred_vec;
					pred_vec preds (8);

					bool has_cycle = false;
					bondgraph::cycle_detector<pred_vec> vis(has_cycle, preds); 
					//boost::depth_first_search(graph.Graph(), visitor(vis));
					boost::breadth_first_search(graph.Graph(), *graph._FindVertex(ref_atom), visitor(vis));
					if (has_cycle) {
						for (Atom_it it = atoms.begin(); it != atoms.end(); it++) {
							printf ("%s(%d), ", (*it)->Name().c_str(), (*it)->ID());
						}
						printf ("\n");
					}

					// Internal compiler error if commented out
					//typedef property_map<bondgraph::BondGraph::graph_t,vertex_color_t>::type Color;

					/*
					bondgraph::BondGraph::Vertex_it vi, vi_end;
					for (boost::tie(vi, vi_end) = boost::vertices(graph.Graph()); vi != vi_end; ++vi)
						boost::get(boost::vertex_color, graph.Graph())[*vi] = boost::white_color;

					// Line B
					boost::depth_first_visit(graph.Graph(), boost::vertices(graph.Graph()).first[1], bondgraph::vertex_processor(),
							get(boost::vertex_color, graph.Graph()) ); 
							*/




					//bondgraph::bfs_atom_visitor vis;
					//boost::breadth_first_search(graph.Graph(), *graph._FindVertex(so2s.S()), visitor(vis));

					//bondgraph::BondGraph::dfs_atom_visitor vis;


					/*
						 Atom_it closest = nm.closest(Atom::O); 
						 nm.next_closest(closest,Atom::O); 
						 nm.next_closest(closest,Atom::O); // first non-covalent O
						 distance_1 = MDSystem::Distance(so2s.S(), *closest).norm();
						 nm.next_closest(closest,Atom::O); // second non-covalent O
						 distance_2 = MDSystem::Distance(so2s.S(), *closest).norm();
						 s_bonds(distance_1, distance_2);

						 nm.OrderAtomsByDistance(so2s.O1());
						 closest = nm.closest(Atom::H);
						 distance_1 = MDSystem::Distance(so2s.O1(), *closest).norm();
						 nm.next_closest(closest,Atom::H);
						 distance_2 = MDSystem::Distance(so2s.O1(), *closest).norm();
						 o_bonds(distance_1, distance_2);

						 nm.OrderAtomsByDistance(so2s.O2());
						 closest = nm.closest(Atom::H);
						 distance_1 = MDSystem::Distance(so2s.O2(), *closest).norm();
						 nm.next_closest(closest, Atom::H);
						 distance_2 = MDSystem::Distance(so2s.O2(), *closest).norm();
						 o_bonds(distance_1, distance_2);
						 */
				}

				void DataOutput() {
					/*
						 s_bonds.OutputData();
						 o_bonds.OutputData();
						 */
				}

			protected:
				NeighborManipulator<T>	nm;
				Histogram2DAgent			s_bonds;
				Histogram2DAgent			o_bonds;
				so2_analysis::SO2SystemManipulator<T>	so2s;
				bondgraph::BondGraph	graph;

				double distance_1, distance_2;
		};	// test analysis

}	// namespace md analysis

#endif
