#include "molgraph.h"

namespace molgraph {
	using namespace md_system;
	using namespace boost;

	int MoleculeGraph::numMoleculeGraphs = 0;

	MoleculeGraph::MoleculeGraph ()
		: Molecule () {
			++numMoleculeGraphs;
		}

	MoleculeGraph::~MoleculeGraph () {
		numMoleculeGraphs--;
	}

	/*
	MoleculeGraph::MoleculeGraph (const Molecule& molecule) 
		: Molecule(molecule) {
			++numMoleculeGraphs;
		}
		*/





	// given one of the atoms of an molgraph (i.e. a carbon), this takes that atom and a graph structure (boost-graph) and builds the entire molgraph (Graph)
	void MoleculeGraph::Initialize (AtomPtr const atom, const bondgraph::BondGraph& graph) {
		//printf ("new molgraph with top atom = "); atom->Print();
		Vertex v = AddAtomToGraph(atom);
		BuildGraph (v, graph);
	}

	void MoleculeGraph::BuildGraph (Vertex v, const bondgraph::BondGraph& graph) {

		// build the queue of connected atoms
		Atom_ptr_vec bonded_atoms = graph.BondedAtoms (_graph[v].atom, bondgraph::covalent);	// only covalently-bound atoms
		std::list<Vertex> queue;
		//printf ("adding:\n");
		for (Atom_it it = bonded_atoms.begin(); it != bonded_atoms.end(); it++) {
			// only queue up atoms that are not already in the graph
			if (!InGraph(*it)) {
				//(*it)->Print();
				Vertex w = AddAtomToGraph(*it);
				queue.push_back(w);
			}
		}

		// connect to each atom in the queue, and then build the sub-graph of that atom
		for (std::list<Vertex>::iterator it = queue.begin(); it != queue.end(); it++) {
			add_edge (v, *it, _graph);
			this->BuildGraph(*it, graph);
		}
	}	// BuildGraph


	Vertex MoleculeGraph::AddAtomToGraph (AtomPtr const atom) { 
		Vertex v = add_vertex(_graph);
		_graph[v].atom = atom;
		//_graph[v].type = null;

		this->AddAtom(atom);
		return v;
	}


	bool MoleculeGraph::InGraph (AtomPtr const atom) const {
		bool ret = false;
		Vertex_it vi, vi_end, next;
		tie(vi, vi_end) = vertices(_graph);
		for (next = vi; vi != vi_end; vi = next) {
			++next;
			if (_graph[*vi].atom == atom) {
				//printf ("found the atom: "); _graph[*vi].atom->Print();
				ret = true;
				break;
			}
		}
		return ret;
	}

	Atom_ptr_vec MoleculeGraph::Atoms () const {

		// return a container of all the atoms in the graph
		Atom_ptr_vec atoms;

		Vertex_it vi, vi_end, next;
		tie(vi, vi_end) = vertices(_graph);
		for (next = vi; vi != vi_end; vi = next) {
			++next;
			atoms.push_back(_graph[*vi].atom);
		}

		return atoms;
	}

}	// namespace molgraph
