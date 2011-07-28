#ifndef GRAPH_H_
#define GRAPH_H_

#include "mdsystem.h"
#include "utility.h"

#include <map>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/depth_first_search.hpp>

//#include <boost/property_map/property_map.hpp>
#include <boost/property_map/property_map.hpp>

#include <exception>

//#define ANGLE_CRITERIA

/* The idea of a connectivity matrix is the same as mapping the connections between different members of a system. In this case, we're dealing with atoms in a system, and the connections are defined by the distances between the atoms. We're going to employ a few tricks here in order to make data retrieval a bit more succinct, and also to store more information into one matrix.

	 Firstly, the matrix is represented by a 2D array of numbers. Each number will represent the distance between 2 atoms. The matrix is symmetric such that each row and each column represent the ordered list of atoms in the system. The upper triangle (all elements above the diagonal) hold the distances between atoms. Only using one half of a symmetric matrix ensures that we don't double our efforts and redo calculations between atoms.

	 Because we'll be working with water, a hydrogen bond is any bond less than 2.4 angstroms and greater than 1.3 (or so... this is to be defined below). If we want to find the number of hydrogen bonds an atom is involved in, then we look at that atom's column from the top to the diagonal, and from the diagonal to the end of the row, and count the number of bonds that fall within the distance of an H-bond.

	 As for the diagonal elements - instead of writing a routine that will count the number of h-bonds by looking over the rows and columns, as distances are calculated the diagonal elements for each atom will be updated to reflect the number of H-bonds formed. Thus at any time we can look at the diagonal and know immediately the number of H-bonds an atom is involved in.

	 This leaves the bottom-diagonal free to store more information. If two atoms are covalently bound, then the bottom diagonal element will mark this with a 1.0.
 */

namespace bondgraph {

	using namespace boost;
	using namespace md_system;

	// various bondlengths to be used
	const double OHBONDLENGTH = 1.15;				// used to be 1.1
	const double HBONDLENGTH  = 2.46;				// used to be 2.46
	//const double HBONDANGLECOS	= cos(30.0*M_PI/180.0);		// bonding angle has to be bigger than this cos (i.e. smaller than ~30 degrees
	const double NOBONDLENGTH = 2.0;
	const double NHBONDLENGTH = 1.3;		// uhmm... check this?
	const double SOBONDLENGTH = 1.64;
	const double SOINTERACTIONLENGTH = 3.25;		// S not covalently bound to O
	const double COBONDLENGTH = 1.5;
	const double CHBONDLENGTH = 1.5;
	const double CCBONDLENGTH = 2.0;



	// bond types
	typedef enum {null = 0, unbonded, hbond, interaction, covalent} bondtype;

	// useful for tracking distances between atom pairs
	typedef std::pair<double, AtomPtr>	distance_pair;
	typedef std::vector<distance_pair> 	distance_vec;
	typedef distance_vec::const_iterator	distance_it;

	/* Encoding of the different coordination types
	 * The numbering is based on each O having a value of 1, and each H haveing a value of 10 (i.e. add 1 for every O, and 10 for every H...). So a water in a state of OOHH bonding would have a coordination of 22, and a coordination of 13 would be OOOH, 12 = OOH, 11 = OH, 10 = H, etc.
	 */
	typedef enum {
		UNBOUND=0, O=1, OO=2, OOO=3, OOOO=4, 			// no H
		H=10, OH=11, OOH=12, OOOH=13, OOOOH=14,			// 1 H
		HH=20, OHH=21, OOHH=22, OOOHH=23, OOOOHH=24,		// 2 Hs
		HHH=30, OHHH=31, OOHHH=32, OOOHHH=33, OOOOHHH=34,	// 3 Hs
		HHHH=40, OHHHH=41, OOHHHH=42, OOOHHHH=43, OOOOHHHH=44
	} coordination;
	// And hopefully that covers all the bonding coordination types :)

	typedef std::map<coordination, std::string> coord_map;





	class BondGraph {

		public:

			// Vertices are atoms
			struct VertexProperties {
				AtomPtr atom;							// the atom itself
				VecR position;
				Atom::Element_t element;	// easier than probing the AtomPtr itself
				AtomPtr parent;						// for use when performing various graph-traversal routines
			};

			// edges are bonds between atoms
			struct EdgeProperties {
				EdgeProperties (const double b_length, const bondtype b_type) : distance(b_length), btype(b_type) { }
				double 		distance;	// bond length - distance between two atoms
				bondtype	btype;		// type of bond (i.e. hbond, covalent, etc.)
			};

			//typedef adjacency_list<listS, listS, undirectedS, VertexProperties, EdgeProperties> graph_t;
			typedef adjacency_list < vecS, vecS, undirectedS, VertexProperties, EdgeProperties > graph_t;
			typedef graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
			typedef vertex_descriptor Vertex;
			typedef graph_traits<graph_t>::vertex_iterator vertex_iterator;
			typedef vertex_iterator Vertex_it;
			typedef std::pair<Vertex_it, Vertex_it> Vertex_pair;

			typedef graph_traits<graph_t>::edge_descriptor edge_descriptor;
			typedef edge_descriptor Edge;
			typedef graph_traits<graph_t>::edge_iterator edge_iterator;
			typedef edge_iterator Edge_it;
			typedef graph_traits<graph_t>::adjacency_iterator adjacency_iterator;
			typedef adjacency_iterator Adj_it;

			/*
				 typedef graph_traits<graph_t>::out_edge_iterator out_edge_iterator;
				 typedef graph_traits<graph_t>::in_edge_iterator in_edge_iterator;
				 typedef graph_t::vertex_property_type vertex_property_type;
				 typedef graph_t::graph_tag graph_tag;

				 typedef graph_traits<graph_t>::directed_category      directed_category;
				 typedef graph_traits<graph_t>::edge_parallel_category edge_parallel_category;
				 typedef graph_traits<graph_t>::traversal_category     traversal_category;

				 typedef graph_traits<graph_t>::vertices_size_type     vertices_size_type;
				 typedef graph_traits<graph_t>::edges_size_type        edges_size_type;
				 typedef graph_traits<graph_t>::degree_size_type       degree_size_type;

				 typedef graph_t::stored_vertex stored_vertex;
			 */

			// generic property maps
			template <class T, class Property_T> 
				struct PropertyMap {
					typedef typename boost::property_map<graph_t, T Property_T::*>::type Type;
				};


			// edge properties
			static PropertyMap<double,EdgeProperties>::Type 			b_length;
			static PropertyMap<bondtype,EdgeProperties>::Type 		b_type;

			// vertex properties
			static PropertyMap<AtomPtr,VertexProperties>::Type 		v_atom;
			static PropertyMap<VecR,VertexProperties>::Type 			v_position;
			static PropertyMap<Atom::Element_t,VertexProperties>::Type	v_elmt;
			static PropertyMap<AtomPtr,VertexProperties>::Type		v_parent;

			void _ParseAtoms (Atom_it first, Atom_it last);
			void _ParseAtoms (const Atom_ptr_vec& atoms);
			void _ParseBonds ();
			void _ClearBonds ();
			void _ClearAtoms ();
			void _ResolveSharedHydrogens ();

			void _SetBond (const Vertex& vi, const Vertex& vj, const double bondlength, const bondtype btype);
			Edge _GetBond (const Vertex& vi, const Vertex& vj) const;
			Edge _GetBond (const AtomPtr a1, const AtomPtr a2) const;
			void _RemoveBond (const Vertex& vi, const Vertex& vj);
			void _RemoveBond (const AtomPtr a1, const AtomPtr a2);

			static graph_t _graph;


			// constructor builds the matrix based on number of atoms to analyze
			BondGraph ();
			BondGraph (const Atom_ptr_vec& atoms);
			~BondGraph ();

			graph_t& Graph() { return _graph; }

			void UpdateGraph (Atom_it, Atom_it);
			void UpdateGraph (const Atom_ptr_vec&);

			Atom_ptr_vec BondedAtoms (
					const AtomPtr ap,
					const bondtype btype = null,
					const Atom::Element_t elmt = Atom::NO_ELEMENT
					) const;

			Atom_ptr_vec InteractingAtoms ( const AtomPtr ap) const;

			Vertex_it _FindVertex (const AtomPtr atom) const;
			// Given a molecule, find the atom (of an optionally given name) that is closest to the molecule but not part of it.
			distance_pair ClosestAtom (const MolPtr&, const Atom::Element_t = Atom::NO_ELEMENT) const;
			// Given an atom, find the atom closest to it (of an optionally given name) that is not part of the same molecule
			// returns a distance pair - [distance, AtomPtr]
			distance_pair ClosestAtom (const AtomPtr, const Atom::Element_t = Atom::NO_ELEMENT, bool = false) const;

			// find the atoms that are closest to an atom
			distance_vec ClosestAtoms (const AtomPtr, const int = 1, const Atom::Element_t = Atom::NO_ELEMENT, bool = false) const;

			// find the atoms closest to a given molecule
			distance_vec ClosestAtoms (const MolPtr, const int = 1, const Atom::Element_t = Atom::NO_ELEMENT) const;

			int NumInteractions (const AtomPtr ap) const;
			int NumHBonds (const AtomPtr) const;
			int NumHBonds (const WaterPtr) const;
			coordination WaterCoordination (const WaterPtr) const;

			// returns the distance between two vertices in the graph
			double Distance (const Vertex&, const Vertex&) const;
			// returns the distance between two atoms
			double Distance (const AtomPtr, const AtomPtr) const;

			// the vertex's parent after running a graph-traversal such as breadth-first search
			AtomPtr Parent (const Vertex&) const;
			AtomPtr Parent (const AtomPtr) const;



			typedef std::exception graphex;

			struct unboundhex : public graphex {
				const char* what() const throw() { return "A hydrogen was found with no covalent bonds, and no hydrogen bonds"; }
			};

			struct multiplyboundhex : public graphex {
				const char* what() const throw() { return "A hydrogen was found with more than 1 covalent bonds"; }
			};


			class VertexIsAtom_pred : public std::binary_function<Vertex_it,AtomPtr,bool> {
				public:
					bool operator() (const Vertex_it& it, const AtomPtr& atom) const {
						return BondGraph::v_atom[*it] == atom;
					}
			};




			/*
			// returns the closest atoms of a given name to a given atom
			// input is the target atom's id, atomname is the name of the other atoms in the system we want returned,
			// and number is the number of nearest atoms
			// i.e. ClosestAtoms (5, O, 3) - returns the three closest O's to the atom with ID 5
			std::vector<Atom *> ClosestAtoms (const int input, const string atomname, const int number) const;
			 */
	};	// BondGraph





	template <class ParentMap>
		class cycle_detector : public bfs_visitor<> {
			public:
				cycle_detector(bool& has_cycle, ParentMap parent_map)
					: m_has_cycle(has_cycle), m_p(parent_map) { }

				template <class Edge, class Graph>
					void tree_edge(Edge e, Graph& g) {
						m_p[target(e, g)] = source(e, g);
					}

				template <class Edge, class Graph>
					void gray_target(Edge e, Graph& g) {
						if (m_p[source(e, g)] != target(e, g))
							m_has_cycle = true;
						AtomPtr i = BondGraph::v_atom[source(e,g)];
						AtomPtr j = BondGraph::v_atom[target(e,g)];
						printf ("\n%s(%d) <--> %s(%d)\n", i->Name().c_str(), i->ID(), j->Name().c_str(), j->ID());
					}
			protected:
				bool& m_has_cycle;
				ParentMap m_p;
		}; 

	/*
		 class cycle_detector : public dfs_visitor<> {
		 public:
		 cycle_detector(bool& has_cycle) : _has_cycle(has_cycle) { }

		 template <class Edge, class Graph>
		 void back_edge(Edge e, Graph& g) { _has_cycle = true; }

		 protected:
		 bool& _has_cycle;
		 };
	 */


	class vertex_processor : public dfs_visitor<> {
		public:
			template <class Vertex, class Graph>
				void discover_vertex(Vertex u, Graph& g)
				{ std::cout << "Discover " << u << "\n";}

			template <class Vertex, class Graph>
				void finish_vertex(Vertex u, Graph& g)
				{ std::cout << " <- " << u << "\n";}

			template <class Edge, class Graph>
				void examine_edge(Edge e, Graph& g)
				{ std::cout << " examine edge " << source(e,g) << " "
					<< target(e,g) << "\n"; }
	}; 





	class bfs_atom_visitor : public default_bfs_visitor {
		public:
			typedef std::list<AtomPtr>	Atom_ptr_list;
			// supply the list of atoms that was used to make the graph,
			// a bool for determining if a cycle was found
			// and a target atom to identify the gray target
			bfs_atom_visitor(Atom_ptr_list& gray_source, Atom_ptr_list& gray_target) : _gray_source(gray_source), _gray_target(gray_target) { }

			template < typename Vertex, typename Graph >
				void initialize_vertex (Vertex v, Graph & g) { 
					BondGraph::v_parent[v] = (AtomPtr)NULL;
				}	// initialize vertex

			// set the currently dequeued vertex as the parent
			template < typename Vertex, typename Graph >
				void examine_vertex (Vertex v, Graph & g) { 
					parent = BondGraph::v_atom[v]; 
				} // examine vertex

			// Mark each child's parent for reconstructing any traversals
			template < typename Edge, typename Graph >
				void tree_edge (Edge e, Graph & g) { 
					BondGraph::Vertex t_v = target(e,g);
					BondGraph::Vertex s_v = source(e,g);

					AtomPtr t = BondGraph::v_atom[t_v];
					AtomPtr s = BondGraph::v_atom[s_v];
					AtomPtr s_parent = BondGraph::v_parent[s_v];

					if (s_parent != t) {
						BondGraph::v_parent[t_v] = s;
						//printf ("%s(%d) --> %s(%d)\n", s->Name().c_str(), s->ID(), t->Name().c_str(), t->ID());
					}
				}	// tree edge


			template < typename Edge, typename Graph >
				void gray_target(Edge e, const Graph & g) {

					++_num_cycles;
					_gray_source.push_back(BondGraph::v_atom[source(e,g)]);
					_gray_target.push_back(BondGraph::v_atom[target(e,g)]);
					//printf ("\n%s(%d) <--> %s(%d)\n", _gray_source->Name().c_str(), _gray_source->ID(), _gray_target->Name().c_str(), _gray_target->ID());

					//for (Atom_it it = _atoms.begin(); it != _atoms.end(); it++) {
					//(*it)->Print();
					//}
					//fflush(stdout);
				}


			/*
				 template < typename Vertex, typename Graph >
				 void finish_vertex (Vertex v, Graph & g) { 
				 printf ("\n\n");
				 }
			 */


			int NumCycles () const { return _num_cycles; }

		private:
			int									_num_cycles;
			Atom_ptr_list&	_gray_source; // the running list of gray sources/targets for each cycle that's found
			Atom_ptr_list&	_gray_target;	
			AtomPtr							parent;
	};


	class WaterCoordination_pred : public std::unary_function <WaterPtr, coordination> {

		protected:
			coordination _c;
			BondGraph &_graph;

		public:
			WaterCoordination_pred (const coordination c, BondGraph& graph) : _c(c), _graph(graph) { }
			bool operator () (const WaterPtr wat) const {
				coordination coord = _graph.WaterCoordination(wat);
				bool ret = (coord == _c) ? true : false;
				return ret;
			}
	};



}	// namespace bondgraph

#endif
