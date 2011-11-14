#include "dimergraph.h"

namespace diacid {

	DimerGraph::graph_t DimerGraph::_graph;

	DimerGraph::PropertyMap<alkane::Diacid *,DimerGraph::DiacidProperties>::Type DimerGraph::v_acid = get(&DiacidProperties::acid, _graph);

	DimerGraph::PropertyMap<h_bond_t,DimerGraph::BondProperties>::Type DimerGraph::b_type = get(&BondProperties::type, _graph);

	DimerGraph::PropertyMap<AtomPtr,DimerGraph::BondProperties>::Type DimerGraph::b_left_atom = get(&BondProperties::left_atom, _graph);
	DimerGraph::PropertyMap<AtomPtr,DimerGraph::BondProperties>::Type DimerGraph::b_right_atom = get(&BondProperties::right_atom, _graph);

	DimerGraph::PropertyMap<alkane::Diacid *,DimerGraph::BondProperties>::Type DimerGraph::b_left_acid = get(&BondProperties::left_acid, _graph);
	DimerGraph::PropertyMap<alkane::Diacid *,DimerGraph::BondProperties>::Type DimerGraph::b_right_acid = get(&BondProperties::right_acid, _graph);

	// clear out the graph and recalculate all the inter-acid connections
	void DimerGraph::Initialize (alkane::Diacid_it start, alkane::Diacid_it end) {

		this->_ClearBonds();
		this->_ClearAcids();
		
		for (alkane::Diacid_it it = start; it != end; it++) {
			_AddAcidToGraph(*it);
		}
	}

	void DimerGraph::RecalculateBonds () {
		this->_ClearBonds();

		Vertex_it vi, vi_end, next;
		tie(vi, vi_end) = vertices(_graph);
		for (next = vi; vi != vi_end; vi = next) {
			++next;
			this->_AddBondsToOtherAcids (*vi);
		}
	}

	void DimerGraph::_ClearBonds () {
		// Remove all the edges.
		Edge_it ei, e_end, next;
		tie(ei, e_end) = edges(_graph);
		for (next = ei; ei != e_end; ei = next) {
			next++;
			remove_edge(*ei, _graph);
		}
		return;
	}

	void DimerGraph::_ClearAcids () {
		// Remove all the vertices.
		Vertex_it vi, vi_end, next;
		tie(vi, vi_end) = vertices(_graph);
		for (next = vi; vi != vi_end; vi = next) {
			++next;
			remove_vertex(*vi, _graph);
		}
		return;
	}

	DimerGraph::Vertex DimerGraph::_AddAcidToGraph (alkane::Diacid * const acid) {
		// add the acid vertex into the data structure
		Vertex v = add_vertex(_graph);
		v_acid[v] = acid;
		// now find the connections to each of the other acids in the graph
		this->_AddBondsToOtherAcids (v);
		return v;
	}

	void DimerGraph::_AddBondsToOtherAcids (Vertex v) {

		// iterate through all the other acids in the graph finding the H-bonds formed between the incoming one
		// and the rest
		Vertex_it vi, vi_end, next;
		tie(vi, vi_end) = vertices(_graph);
		for (next = vi; vi != vi_end; vi = next) {
			++next;
			// skip the vertex matched to itself
			if (*vi == v) continue;
			// find all the bonds between hydrogens on the incoming acid and the oxygens on the other acids
			_FindBondsBetweenAcids(v, *vi);
			// then do the reverse
			_FindBondsBetweenAcids(*vi, v);
		}
	} // add dimer bonds


	void  DimerGraph::_FindBondsBetweenAcids (Vertex vleft, Vertex vright) {

		alkane::Diacid * left = v_acid[vleft];
		alkane::Diacid * right = v_acid[vright];

		// grab all the Hs on the left molecule, and the Os on the right
		left->LoadAtomGroups();
		Atom_ptr_vec methyl_Hs = left->methyl_hydrogens();
		Atom_ptr_vec carbonyl_Hs = left->carbonyl_hydrogens();
		Atom_ptr_vec carbonyl_Os = right->carbonyl_oxygens();

		bool b;
		Edge e;
		double distance;

		for (Atom_it O = carbonyl_Os.begin(); O != carbonyl_Os.end(); O++) {

			// find the carbonyl H-bonds
			for (Atom_it H = carbonyl_Hs.begin(); H != carbonyl_Hs.end(); H++) {
				distance = MDSystem::Distance(*O, *H).norm(); 
				//(*O)->Print(); (*H)->Print();
				//std::cout << distance << std::endl;
				if (distance <= bondgraph::HBONDLENGTH) {
					//std::cout << "found two" << std::endl;
					tie(e, b) = add_edge(vleft, vright, 
							BondProperties(CARBONYL_CARBONYL, *H, *O, left, right), _graph);
				}
			}

			// then find Hbonds from the oxygens to methyl hydrogens
			for (Atom_it H = methyl_Hs.begin(); H != methyl_Hs.end(); H++) {
				distance = MDSystem::Distance(*O, *H).norm();
				//(*O)->Print(); (*H)->Print();
				//std::cout << distance << std::endl;
				if (distance <= bondgraph::HBONDLENGTH) {
					//std::cout << "found one" << std::endl;
					tie(e, b) = add_edge(vleft, vright, 
							BondProperties(METHYL_CARBONYL, *H, *O, left, right), _graph);
				}
			}
		}
	}

	DimerGraph::Vertex_it DimerGraph::_FindVertex (const alkane::Diacid * acid) const {
		Vertex_it vi, vi_end, next;
		tie(vi, vi_end) = vertices(_graph);
		for (next = vi; vi != vi_end; vi = next) {
			++next;
			if (v_acid[*vi] == acid)
				break;
		}
		return (vi);
	}

	// given an acid, find the number of h-bonds it forms with another acid
	// returns 
	//		number of carbonyl-O -- methyl-H bonds
	//		number of carb-carb bonds
	//		number of other acids to which it's bound
	DimerGraph::hbond_data_t DimerGraph::NumHBonds (alkane::Diacid * acid) {
		Vertex_it v = _FindVertex(acid);

		std::set<alkane::Diacid *> acids;

		hbond_data_t dat;
		dat.methyls = 0;
		dat.carbonyls = 0;
		dat.acids = 0;

		out_edge_iterator out, end, next;
		tie(out, end) = out_edges(*v, _graph);
		for (next = out; out != end; out = next) {
			++next;
			if (b_type[*out] == METHYL_CARBONYL)
				++dat.methyls;
			if (b_type[*out] == CARBONYL_CARBONYL)
				++dat.carbonyls;
			acids.insert(v_acid[target(*out, _graph)]);
		}
		dat.acids = acids.size();

		return dat;
	}

}
