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
		// clear out the previous graph edges
		this->_ClearBonds();

		// go through each pair of acids and find the bonds between them
		Vertex_it vi, vj, vi_end;
		
		for (tie(vi, vi_end) = vertices(_graph); vi != vi_end-1; vi++) {
			for (vj = vi+1; vj != vi_end; vj++) {
				this->_AddBondsToOtherAcids (*vi, *vj);
			}
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
		return v;
	}

	void DimerGraph::_AddBondsToOtherAcids (Vertex left, Vertex right) {

		if (left == right) return;
		// find all the bonds between hydrogens on the incoming acid and the oxygens on the other acids
		_FindBondsBetweenAcids(left, right);
		// then do the reverse
		_FindBondsBetweenAcids(right, left);
	} // add dimer bonds


	void  DimerGraph::_FindBondsBetweenAcids (Vertex vleft, Vertex vright) {

		alkane::Diacid * left = v_acid[vleft];
		alkane::Diacid * right = v_acid[vright];

		// grab all the Hs on the left molecule, and the Os on the right
		left->SetAtoms();
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
		for (tie(vi, vi_end) = vertices(_graph); vi != vi_end; vi++) {
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
	DimerGraph::hbond_data_t DimerGraph::NumHBonds (alkane::Diacid * acid) const {
		Vertex_it v = _FindVertex(acid);

		std::set<alkane::Diacid *> acids;

		hbond_data_t dat;
		dat.methyls = 0;
		dat.carbonyls = 0;
		dat.acids = 0;

		//std::cout << "\nAcid #" << acid->MolID() << std::endl;
		out_edge_iterator out, end, next;
		for (tie(out, end) = out_edges(*v, _graph); out != end; out++) {
			if (b_type[*out] == METHYL_CARBONYL) {
				//std::cout << "\nmethyl carbonyl between:\n" << std::endl;
				//b_left_atom[*out]->Print();
				//b_right_atom[*out]->Print();
				++dat.methyls;
			}
			if (b_type[*out] == CARBONYL_CARBONYL) {
				//std::cout << "\ncarbonyl carbonyl between:\n" << std::endl;
				//b_left_atom[*out]->Print();
				//b_right_atom[*out]->Print();
				++dat.carbonyls;
			}
			acids.insert(v_acid[target(*out, _graph)]);
		}
		dat.acids = acids.size();

		return dat;
	}

	// returns the same hbonding data as for NumHBonds, but for the entire system, not just for a single acid
	DimerGraph::hbond_data_t DimerGraph::NumHBonds () const {

		std::set<alkane::Diacid *> acids;

		hbond_data_t dat;
		dat.methyls = 0; // total number of methyl-carbo hbonds
		dat.carbonyls = 0; // total number of carbo-carbo hbonds
		dat.acids = 0;	// number of acids in the system taking part in dimerization -- h-bonding

		//std::cout << "\nAcid #" << acid->MolID() << std::endl;
		Edge_it ei, end;
		for (tie(ei, end) = edges(_graph); ei != end; ei++) {
			if (b_type[*ei] == METHYL_CARBONYL) {
				//std::cout << "\nmethyl carbonyl between:\n" << std::endl;
				//b_left_atom[*out]->Print();
				//b_right_atom[*out]->Print();
				++dat.methyls;
				acids.insert(v_acid[target(*ei, _graph)]);
				acids.insert(v_acid[source(*ei, _graph)]);
			}
			if (b_type[*ei] == CARBONYL_CARBONYL) {
				//std::cout << "\ncarbonyl carbonyl between:\n" << std::endl;
				//b_left_atom[*out]->Print();
				//b_right_atom[*out]->Print();
				acids.insert(v_acid[target(*ei, _graph)]);
				acids.insert(v_acid[source(*ei, _graph)]);
				++dat.carbonyls;
			}
		}
		dat.acids = acids.size();

		return dat;
	}


	// check if the two given acids are connected in any way
	bool DimerGraph::_AcidsConnected (Vertex left, Vertex right) const {
		Edge e;
		bool b;
		tie(e,b) = edge(left, right, _graph);
		return b;
	}	

	// counts the number of bound acids. If any two acids are hbonded in some way, then that counts as a single connection. 
	// Returned is: 
	//		the total number of connections between acids
	//	e.g. 3 acids with 2 connections means a chain of 3 acids. 
	//	3 acids with 3 connections means a cyclic set of connections between the acids
	int DimerGraph::InterAcidConnections () const {

		Vertex_it vi, vj, end;
		int num_connections = 0;
		// iterate through every pair of acids
		for (tie (vi,end) = vertices(_graph); vi != end-1; vi++) {
			for (vj = vi+1; vj != end; vj++) {
				if (_AcidsConnected(*vi,*vj)) 
					++num_connections;
			}
		}

		return num_connections;
	}


} // namespace diacid
