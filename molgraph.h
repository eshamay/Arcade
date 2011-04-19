#ifndef CARBONCHAIN_H_
#define CARBONCHAIN_H_

#include "molecule.h"
#include "bondgraph.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map.hpp>

namespace bondgraph {
	class BondGraph;
}

namespace molgraph {
	using namespace md_system;
	using namespace boost;


	// implement an molgraph molecule as a graph
	// Nodes:
	//		atom pointer to the actual atom
	//		atom type to specify the specific type of atom (i.e. carbonyl C, or alcohol O, etc)
	//		links to all the attached atoms

	// Vertices are atoms
	struct AtomProperties {
		AtomPtr atom;							// the atom itself
	};

	typedef adjacency_list < listS, vecS, bidirectionalS, AtomProperties > graph_t;	// graph type for molgraphs
	typedef graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
	typedef vertex_descriptor Vertex;
	typedef graph_traits<graph_t>::vertex_iterator vertex_iterator;
	typedef vertex_iterator Vertex_it;



	class MoleculeGraph : public Molecule {

		protected:
			graph_t	_graph;

			void BuildGraph (Vertex v, const bondgraph::BondGraph& graph);

		public:
			MoleculeGraph ();			// a default constructor
			virtual ~MoleculeGraph ();
			MoleculeGraph (const Molecule& molecule);		// copy constructor for casting from a molecule

			static int numMoleculeGraphs;			// total number of carbon chains in the system

			Vertex AddAtomToGraph (AtomPtr const atom);
			bool InGraph (AtomPtr const atom) const;

			void Initialize (AtomPtr const atom, const bondgraph::BondGraph& graph);

			Atom_ptr_vec Atoms () const;
	};

}	//namespace molgraph
#endif
