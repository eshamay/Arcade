#ifndef CARBONCHAIN_H_
#define CARBONCHAIN_H_

#include "molecule.h"
#include "bondgraph.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map.hpp>

namespace bondgraph {
	class BondGraph;
}

namespace alkane {
	using namespace md_system;
	using namespace boost;


	// implement an alkane molecule as a graph
	// Nodes:
	//		atom pointer to the actual atom
	//		atom type to specify the specific type of atom (i.e. carbonyl C, or alcohol O, etc)
	//		links to all the attached atoms

	typedef enum {null, CARBONYL_C, CARBONYL_O, ALCOHOL_O, CH2_C, CH2_H, ALCOHOL_H} alkane_atom_t;

	// Vertices are atoms
	struct AtomProperties {
		AtomPtr atom;							// the atom itself
		alkane_atom_t	type;
	};

	typedef adjacency_list < listS, vecS, bidirectionalS, AtomProperties > graph_t;	// graph type for alkanes
	typedef graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
	typedef vertex_descriptor Vertex;
	typedef graph_traits<graph_t>::vertex_iterator vertex_iterator;
	typedef vertex_iterator Vertex_it;



	class Alkane : public Molecule {

		protected:
			graph_t	_graph;

			void BuildGraph (Vertex v, const bondgraph::BondGraph& graph);
			Vertex AddAtomToGraph (AtomPtr const atom);
			bool InGraph (AtomPtr const atom) const;

		public:
			Alkane ();			// a default constructor
			virtual ~Alkane ();
			Alkane (const Molecule& molecule);		// copy constructor for casting from a molecule

			static int numAlkanes;			// total number of carbon chains in the system

			void InitializeAlkane (AtomPtr const atom, const bondgraph::BondGraph& graph);
			// Functions for analysis
			//virtual void SetAtoms () = 0;
	};

}
#endif
