#include "bondgraph.h"

namespace bondgraph {
	BondGraph::graph_t BondGraph::_graph (0);

	// property map definitions
	BondGraph::PropertyMap<double,BondGraph::EdgeProperties>::Type BondGraph::b_length = get(&EdgeProperties::distance, _graph);
	BondGraph::PropertyMap<bondtype,BondGraph::EdgeProperties>::Type BondGraph::b_type = get(&EdgeProperties::btype, _graph);

	BondGraph::PropertyMap<AtomPtr,BondGraph::VertexProperties>::Type BondGraph::v_atom = get(&VertexProperties::atom, _graph);
	BondGraph::PropertyMap<VecR,BondGraph::VertexProperties>::Type BondGraph::v_position = get(&VertexProperties::position, _graph);
	BondGraph::PropertyMap<Atom::Element_t,BondGraph::VertexProperties>::Type BondGraph::v_elmt = get(&VertexProperties::element, _graph);
	BondGraph::PropertyMap<AtomPtr,BondGraph::VertexProperties>::Type BondGraph::v_parent = get(&VertexProperties::parent, _graph);


	/*
	BondGraph::BondGraph () :
		_sys_type("xyz") { 
			return; 
		}
		*/

	BondGraph::BondGraph () { }

	BondGraph::BondGraph (const Atom_ptr_vec& atoms) {
		this->UpdateGraph(atoms);
		return;
	}



	BondGraph::~BondGraph () {
		//this->_ClearBonds();
		//this->_ClearAtoms();
		return;
	}


	void BondGraph::_ParseAtoms (Atom_it first, Atom_it last) {
		// set all the vertices to contain the proper info
		//std::cout << "bondgraph::parseatoms " << last - first << std::endl;
		Vertex v;
		int i = 0;
		for (Atom_it it = first; it != last; it++) {
			v = boost::vertex(i, _graph);
			v_atom[v] = *it;
			v_position[v] = (*it)->Position();
			v_elmt[v] = (*it)->Element();
			i++;
		}
		//std::cout << "bondgraph::parseatoms " << last - first << std::endl;

		return;
	}

	void BondGraph::_ParseAtoms (const Atom_ptr_vec& atoms) {
		this->_ParseAtoms(atoms.begin(), atoms.end());
		return;
	}


	// Calculate the distance between each of the atoms and update the graph edge list with bonds between them
	// Edges/bonds should only exist between atoms that are bound in some fashion (h-hond, covalent, etc.)
	void BondGraph::_ParseBonds () {
		// first clear out all the bonds from before
		this->_ClearBonds();

		// now run through all atom combos and get their bondlengths in order to set their bondtypes
		AtomPtr ai, aj;

		// this little for-loop bit runs through all atom-pair combinations once and find the bond-types between them
		Vertex_it vi, vi_end, next_i;
		tie(vi, vi_end) = vertices(_graph);
		Vertex_it vj, vj_end, next_j;
		tie(vj, vj_end) = vertices(_graph);

		--vi_end;
		for (next_i = vi; vi != vi_end; vi = next_i) {
			++next_i;
			vj = next_i;
			for (next_j = vj; vj != vj_end; vj = next_j) {
				++next_j;

				ai = v_atom[*vi];	// first atom
				aj = v_atom[*vj];	// second atom
				// Don't connect oxygens to oxygens, and hydrogen to hydrogen...etc.
				//if (Atom::element_eq(ai,aj)) continue;

				// calculate the distance between the two atoms (taking into account the periodic boundaries)
				double bondlength = MDSystem::Distance (v_position[*vi], v_position[*vj]).Magnitude();
				if (bondlength > HBONDLENGTH && bondlength > SOINTERACTIONLENGTH) continue;
				// all bonds are considered unbound unless proven otherwise
				bondtype btype = unbonded;

				// first process O-H bonds
				if (Atom::ElementCombo(ai,aj,Atom::O,Atom::H))
				{
					// one type of bond is the O-H covalent
					if (bondlength <= OHBONDLENGTH) {
						btype = covalent;
					}

					// Or check if an H-bond is formed
					else if (bondlength <= HBONDLENGTH) {
#ifdef ANGLE_CRITERIA
						// additionally, let's check the angle-criteria for an H-bond.
						// This is done by looking at the angle formed from
						// o1 is covalently bound to h, and o2 is h-bound to h
						AtomPtr o1 = (AtomPtr)NULL, h = (AtomPtr)NULL, o2 = (AtomPtr)NULL;	

						if (ai->Element() == Atom::O) {		// ai is the O, and aj is the H
							o2 = ai;
							h = aj;
						}
						else if (aj->Element() == Atom::O) {
							o2 = aj;
							h = ai;
						}

						o1 = h->ParentMolecule()->GetAtom("O");

						if (h == (AtomPtr)NULL || o1 == (AtomPtr)NULL || o2 == (AtomPtr)NULL) {
							//throw (MALFORMED_H2O);
							//printf ("Something wrong in assigning the atoms O and H in forming an H-bond - graph.cpp\n");
							//exit(1);
						}

						VecR o1h = MDSystem::Distance (o1, h);	// the covalent bond
						VecR ho2 = MDSystem::Distance (h, o2);	// the H-bond

						double anglecos = o1h < ho2;		// cos(theta)

						if (anglecos > HBONDANGLECOS){
							//printf ("% 10.3f / %6.3f\n", acos(angle)*180.0/M_PI, HBONDANGLECOS);
#endif
							btype = hbond;
#ifdef ANGLE_CRITERIA
						}
#endif
					}	// check for h-bond

				}	// Check OH bond combos


				// process N-O bonds
				else if (Atom::ElementCombo (ai,aj, Atom::N, Atom::O) && (bondlength < NOBONDLENGTH)) {
					btype = covalent;
				}

				// now process SO2 molecules
				else if (Atom::ElementCombo (ai,aj, Atom::S, Atom::O)) { 
					if (bondlength < SOBONDLENGTH) {
						btype = covalent;
					} // process S-O covalent bonds
					else if (bondlength < SOINTERACTIONLENGTH) {
						btype = interaction;
					}	// process S-O interactions (sorta h-bonds)
				}

				else if (Atom::ElementCombo (ai,aj, Atom::C, Atom::O)) {
					if (bondlength < COBONDLENGTH) {
						btype = covalent;
					}
				}
				else if (Atom::ElementCombo (ai,aj, Atom::C, Atom::H)) {
					if (bondlength < CHBONDLENGTH) {
						btype = covalent;
					}
				}
				else if (Atom::ElementCombo (ai,aj, Atom::C, Atom::C)) {
					if (bondlength < CCBONDLENGTH) {
						btype = covalent;
					}
				}

				// add in the bond between two atoms
				if (btype != unbonded)
					this->_SetBond (*vi, *vj, bondlength, btype);
			}
		}

		/*
		// Now fix up any weird atom-sharing between molecules. At this point we have to consider if we want to divide the system into separate molecules, or if we're interested in other phenomena, such as contact-ion pairs, etc.
		if (_sys_type == "xyz")
		try {
		this->_ResolveSharedHydrogens ();
		} catch (unboundhex& ex) {
		std::cout << "Exception caught while resolving hydrogens shared between multiple molecules" << std::endl;
		throw;
		}
		*/
	}	// Parse Bonds

	void BondGraph::_ClearBonds () {
		// Remove all the edges.
		Edge_it ei, e_end, next;
		tie(ei, e_end) = edges(_graph);
		for (next = ei; ei != e_end; ei = next) {
			next++;
			remove_edge(*ei, _graph);
		}
		return;
	}

	void BondGraph::_ClearAtoms () {
		// Remove all the vertices.
		Vertex_it vi, vi_end, next;
		tie(vi, vi_end) = vertices(_graph);
		for (next = vi; vi != vi_end; vi = next) {
			++next;
			remove_vertex(*vi, _graph);
		}
		return;
	}

	void BondGraph::UpdateGraph (const Atom_it first, const Atom_it last) {
		// clear out the old graph info
		//std::cout << "1 bondgraph::updategraph " << last - first << std::endl;

		_graph.clear();
		//std::cout << "2 bondgraph::updategraph " << last - first << std::endl;

		// resize the graph with new vertices
		for (Atom_it it = first; it != last; it++) {
			boost::add_vertex(_graph);
		}
		//std::cout << "3 bondgraph::updategraph " << last - first << std::endl;

		// parse the atom info into the vertices
		this->_ParseAtoms(first, last);
		// then find all the needed bond information
		
		try {
			this->_ParseBonds ();
		} catch (graphex& ex) {
			std::cout << "An exception was caught while determining the bonds between atoms in the system" << std::endl;
		}

		//std::cout << "4 bondgraph::updategraph " << last - first << std::endl;
	}

	void BondGraph::UpdateGraph (const Atom_ptr_vec& atoms) {
		this->UpdateGraph (atoms.begin(), atoms.end());
		//this->_ParseAtoms(atoms.begin(), atoms.end());
		return;
	}

	void BondGraph::_SetBond (const Vertex& vi, const Vertex& vj, const double bondlength, const bondtype btype) {
		bool b;
		Edge e;

		tie(e, b) = add_edge(vi, vj, EdgeProperties(bondlength, btype), _graph);
		/*
			 if (!b) {
			 std::cout << "BondGraph::SetBond() - Tried to add a bond to an already-bonded atom pair" << std::endl;
			 v_atom[vi]->Print();
			 v_atom[vj]->Print();
			 exit(1);
			 }
			 */

		return;
	}

	BondGraph::Vertex_it BondGraph::_FindVertex (const AtomPtr atom) const {
		Vertex_it vi, vi_end, next;
		tie(vi, vi_end) = vertices(_graph);
		for (next = vi; vi != vi_end; vi = next) {
			++next;
			if (v_atom[*vi] == atom)
				break;
		}

		//Vertex_pair vp = vertices(_graph);
		//Vertex_it vi, vi_end, atom_vt;
		// find the vertex
		//Vertex_it atom_vt = std::find_if (vp.first, vp.second, std::bind2nd(VertexIsAtom_pred(), atom));
		//tie(vi, vi_end) = vertices(_graph);

		return (vi);
	}

	// returns a list of the atoms bonded to the given atom
	Atom_ptr_vec BondGraph::BondedAtoms (const AtomPtr ap, bondtype const btype, Atom::Element_t const elmt) const {
		Atom_ptr_vec atoms;
		Vertex_it va = this->_FindVertex(ap);
		Edge e;

		Adj_it vi, vi_end, next;
		tie(vi, vi_end) = adjacent_vertices(*va, _graph);
		for (next = vi; vi != vi_end; vi = next) {
			next++;

			// don't check an atom against itself
			if (v_atom[*va] == v_atom[*vi]) continue;

			// check the bondtype criteria - return only bonds that are specified by the bondtype argument, or if no argument is specified, return all hbond and covalent bonds.
			e = edge(*va, *vi, _graph).first;
			if (btype == b_type[e] || ((!btype) && ((b_type[e] == hbond) || (b_type[e] == covalent) || (b_type[e] == interaction)))) {
		//	if (btype == b_type[e] || (!btype && b_type[e] != unbonded)) {
				if (!elmt || (v_elmt[*vi] == elmt)) {
					atoms.push_back(v_atom[*vi]);
				}
			}
		}

		return (atoms);
	}	// Bonded atoms

	Atom_ptr_vec BondGraph::InteractingAtoms (const AtomPtr ap) const {
		return this->BondedAtoms(ap, interaction);
	}

	int BondGraph::NumInteractions (const AtomPtr ap) const {
		return (this->InteractingAtoms(ap).size());
	}

	int BondGraph::NumHBonds (const AtomPtr ap) const {
		return (this->BondedAtoms(ap, hbond).size());
	}

	int BondGraph::NumHBonds (const WaterPtr wat) const {
		int num = 0;
		for (Atom_it it = wat->begin(); it != wat->end(); it++) {
			num += this->NumHBonds(*it);
		}

		return (num);
	}

	// calculates the water bonding coordination of a given water molecule
	coordination BondGraph::WaterCoordination (const WaterPtr wat) const {
		int c = 0;
		for (Atom_it atom = wat->begin(); atom != wat->end(); atom++) {
			int bonds = NumHBonds (*atom);

			if ((*atom)->Element() == Atom::H) {
				c += 10 * bonds;
			}
			if ((*atom)->Element() == Atom::O) {
				c += bonds;
			}
		}
		return coordination(c);
	}	// water coordination


	void BondGraph::_ResolveSharedHydrogens () {
		/*
		// if any of the Hs are being shared between two molecules (i.e. a contact-ion pair) then we have to resolve which molecule gets the atom.
		// A simple solution is to decide that the molecule with the oxygen closer to the hydrogen is the one that wins out
		// This routine ***assumes*** that all the Hs will be bound to an Oxygen, and not some other atom in the system.
		// the closer oxygen becomes the covalent one, and the further oxygen becomes the h-bonded one

		// we'll go through each of the hydrogens and find if they are bound to one or two different molecules
		Vertex_it vi, vi_end, next;
		tie (vi, vi_end) = vertices(_graph);
		for (next = vi; vi != vi_end; vi = next) { next++;

		AtomPtr H = v_atom[*vi];

		if (H->Element() != Atom::H) continue;

		// here we check to see if it's bound to multiple atoms
		Atom_ptr_vec covalent_atoms = this->BondedAtoms (H, covalent);

		// great ... if the H is only covalently bound to one atom then don't worry about it
		if (covalent_atoms.size() == 1) continue;

		Atom_ptr_vec hbond_atoms = this->BondedAtoms (H, hbond);

		// The hydrogen isn't close enough to any oxygen to be covalent - it's probably shared between 2 as hbonds
		if (covalent_atoms.empty() && !hbond_atoms.empty()) {

		// find the closest atom to the hydrogen
		distance_pair dp = ClosestAtom (H);

		Vertex_it vO = this->_FindVertex(dp.second);

		// the closer bond is now changed from hbonded to covalent. The others are unchanged as they remain hbound
		this->_RemoveBond (*vi, *vO);								// cut the hbond to the nearest oxygen
		this->_SetBond(*vi, *vO, dp.first, covalent);	// and replace it with a covalent bond
		}

		// in the case that there are more than 1 covalent bond to the hydrogen (it's being shared, covalently, by more than two atoms!) then we have a physically very improbable scenario!
		else if (covalent_atoms.size() >= 2) {
		std::cout << "Found a hydrogen atom that is covalently bound to multiple oxygens. This is physically very improbable." << std::endl;
		H->Print();
		throw (multiplyboundhex());
		}

		else if (covalent_atoms.empty() && hbond_atoms.empty()) {
		std::cout << "Found a hydrogen with no covalent or hydrogen bonds to which it is attached (below)" << std::endl;
		H->Print();
		throw (unboundhex());
		}

		}

*/
		return;
	}	// Resolve shared hydrogens

	double BondGraph::Distance (const Vertex& vi, const Vertex& vj) const {
		bool b;
		Edge e;
		tie (e,b) = edge(vi, vj, _graph);
		if (!b) {
			printf ("The distance routine was run for the following two atoms, but they were never connected in the graph\n");
			v_atom[vi]->Print();
			v_atom[vj]->Print();
			exit(1);
		}
		return (b_length[e]);
	}

	double BondGraph::Distance (const AtomPtr a1, const AtomPtr a2) const {
		Vertex_it vi, vj;
		vi = this->_FindVertex (a1);
		vj = this->_FindVertex (a2);

		return Distance(*vi,*vj);
	}



	AtomPtr BondGraph::Parent (const Vertex& v) const { 
		return v_parent[v];
	}

	AtomPtr BondGraph::Parent (const AtomPtr a) const { 
		return Parent(*(this->_FindVertex(a)));
	}



	distance_pair BondGraph::ClosestAtom (const MolPtr& mol, const Atom::Element_t elmt) const {
		distance_vec distances;

		for (Atom_it it = mol->begin(); it != mol->end(); it++) {
			distances.push_back ((ClosestAtoms (*it, 1, elmt, false))[0]);
		}

		pair_utility::pair_sort_first(distances.begin(), distances.end());

		return distances[0];
	}


	distance_pair BondGraph::ClosestAtom (const AtomPtr atom, const Atom::Element_t elmt, bool SameMoleculeCheck) const {
		return ClosestAtoms (atom, 1, elmt, SameMoleculeCheck)[0];
	}	// Closest Atom

	distance_vec BondGraph::ClosestAtoms (const AtomPtr atom, const int num, const Atom::Element_t elmt, bool SameMoleculeCheck) const {
		distance_vec distances;

		// find the distance between the target atom and all the other atoms in the system
		Vertex_it vt, vi, vi_end, next;
		vt = _FindVertex (atom);

		tie(vi, vi_end) = vertices(_graph);
		for (next = vi; vi != vi_end; vi = next) {
			++next;

			// don't check the same atom against itself - that's not defined!
			if (vi == vt) continue;

			// don't consider atoms within the same molecule, or an atom compared to itself
			if (!SameMoleculeCheck && (v_atom[*vt]->ParentMolecule() == v_atom[*vi]->ParentMolecule() || vt == vi)) continue;

			// do an element check so that only atoms with the (optional) given element type are considered
			if (!elmt || elmt == v_elmt[*vi]) { 
				distances.push_back (std::make_pair (Distance(*vt, *vi), v_atom[*vi]));
			}
		}

		// sort all the distances to find the one closest
		pair_utility::pair_sort_first(distances.begin(), distances.end());
		distances.erase (distances.begin()+num, distances.end());

		return distances;
	}

	distance_vec BondGraph::ClosestAtoms (const MolPtr mol, const int num, const Atom::Element_t elmt) const {
		distance_vec distances;
		for (Atom_it atom = mol->begin(); atom != mol->end(); atom++) {
			distance_vec closest = ClosestAtoms (*atom, 10, elmt, false);
			std::copy (closest.begin(), closest.end(), std::back_inserter(distances));
		}

		pair_utility::pair_sort_first (distances.begin(), distances.end());
		distance_vec::const_iterator it = std::unique(distances.begin(), distances.end(), pair_utility::pair_equal_second_pred<distance_pair>());
		distances.resize(it - distances.begin());
		pair_utility::pair_sort_first (distances.begin(), distances.end());
		distances.resize(num);

		return distances;
	} // closest Atoms - molecular version



	BondGraph::Edge BondGraph::_GetBond (const Vertex& vi, const Vertex& vj) const {
		Edge e;
		bool b;
		tie (e,b) = edge(vi, vj, _graph);

		if (!b) {
			//throw(BOND_NOT_FOUND);
			/*
				 printf ("BondGraph::_GetBond - Couldn't find the requested Bond");
				 exit(1);
				 */
		}

		return (e);
	}	// _Get Bond

	BondGraph::Edge BondGraph::_GetBond (const AtomPtr a1, const AtomPtr a2) const {
		Vertex_it v1 = this->_FindVertex(a1);
		Vertex_it v2 = this->_FindVertex(a2);
		Edge e = this->_GetBond(*v1, *v2);

		return (e);
	}

	void BondGraph::_RemoveBond (const Vertex& vi, const Vertex& vj) {
		Edge e = this->_GetBond(vi, vj);
		remove_edge(e, _graph);

		return;
	}

	void BondGraph::_RemoveBond (const AtomPtr a1, const AtomPtr a2) {
		Edge e = this->_GetBond (a1, a2);
		remove_edge(e, _graph);

		return;
	}
}	// namespace bondgraph
