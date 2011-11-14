#ifndef DIMERGRAPH_H_
#define DIMERGRAPH_H_

#include "mdsystem.h"
#include "utility.h"
#include "molecule-analysis.h"
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>

namespace diacid {

	using namespace md_system;
	using namespace boost;

	typedef enum {
		METHYL_CARBONYL,
		CARBONYL_CARBONYL
	} h_bond_t;

	class DimerGraph {

		public:
			struct DiacidProperties {
				alkane::Diacid * acid;							// the diacid represented by the vertex
			};

			struct BondProperties {
				BondProperties (h_bond_t t, AtomPtr latom, AtomPtr ratom, alkane::Diacid * lacid, alkane::Diacid * racid) :
					type(t), left_atom(latom), right_atom(ratom), left_acid(lacid), right_acid(racid) { }

				h_bond_t type;
				AtomPtr left_atom, right_atom;		// atoms involved in the bond
				alkane::Diacid * left_acid, * right_acid;
			};

			typedef adjacency_list < listS, vecS, undirectedS, DiacidProperties, BondProperties > graph_t;	// graph type for molgraphs
			typedef graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
			typedef vertex_descriptor Vertex;
			typedef graph_traits<graph_t>::vertex_iterator vertex_iterator;
			typedef vertex_iterator Vertex_it;

			typedef graph_traits<graph_t>::edge_descriptor edge_descriptor;
			typedef edge_descriptor Edge;
			typedef graph_traits<graph_t>::edge_iterator edge_iterator;
			typedef edge_iterator Edge_it;
			typedef graph_traits<graph_t>::out_edge_iterator out_edge_iterator;
			typedef out_edge_iterator Out_it;

			typedef graph_traits<graph_t>::adjacency_iterator adjacency_iterator;
			typedef adjacency_iterator Adj_it;


			template <class T, class Property_T> 
				struct PropertyMap {
					typedef typename boost::property_map<graph_t, T Property_T::*>::type Type;
				};

			static PropertyMap<alkane::Diacid *,DiacidProperties>::Type 	v_acid;
			static PropertyMap<h_bond_t,BondProperties>::Type							b_type;
			static PropertyMap<AtomPtr,BondProperties>::Type							b_left_atom;;
			static PropertyMap<AtomPtr,BondProperties>::Type							b_right_atom;;
			static PropertyMap<alkane::Diacid *,BondProperties>::Type 		b_left_acid;;
			static PropertyMap<alkane::Diacid *,BondProperties>::Type 		b_right_acid;;

			void Initialize (alkane::Diacid_it start, alkane::Diacid_it end);
			void RecalculateBonds ();
			void _ClearBonds ();
			void _ClearAcids ();


			typedef struct { int methyls, carbonyls, acids; } hbond_data_t;
			hbond_data_t NumHBonds (alkane::Diacid * acid) const;
			hbond_data_t NumHBonds () const;
			int InterAcidConnections () const;

		private:
			static graph_t	_graph;

			Vertex _AddAcidToGraph (alkane::Diacid * const acid);
			void _AddBondsToOtherAcids (Vertex left, Vertex right);
			void  _FindBondsBetweenAcids (Vertex vleft, Vertex vright);
			Vertex_it _FindVertex (const alkane::Diacid * acid) const;
			bool _AcidsConnected (Vertex left, Vertex right) const;

	}; // dimer graph

} // namespace diacid


#endif
