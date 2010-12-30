#ifndef NEIGHBOR_ANALYSIS_H_
#define NEIGHBOR_ANALYSIS_H_

#include "analysis.h"
#include "h2o-analysis.h"
#include "so2-system-analysis.h"

# include <boost/iterator/iterator_facade.hpp>

namespace md_analysis {


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
					std::sort(
							this->analysis_atoms.begin(), 
							this->analysis_atoms.end(), 
							system_t::atomic_reference_distance_pred (reference_atom));
				}

				// iterate over the atoms sorted by distance - but only over a particular element.
				Atom_it begin (const Atom::Element_t elmt) {
					Atom_it it = this->analysis_atoms.begin();
					if ((*it)->Element() != elmt)
						this->Increment(it, elmt);
					return it;
				}

				Atom_it end () { return this->analysis_atoms.end(); }

				Atom_it Incremenent (const Atom_it it, const Atom::Element_t elmt) {
					while (it != this->analysis_atoms.end() && (*(++it))->Element() != elmt) { }
					return it;
				}

			private:

		}; // neighbor manipulator

	template <typename T>
		class TestAnalysis : public AnalysisSet<T> {
			typedef Analyzer<T> system_t;

			TestAnalysis (system_t * t) :
				AnalysisSet<T>(t,
						std::string ("test analysis"),
						std::string (""))
			{ 
				AtomPtr a = this->analysis_atoms[0];
				a->Print();
				OrderAtomsByDistance(a);
				for (Atom_it it = begin(Atom::O); it != end(); Increment(it,Atom::O)) {
					(*it)->Print();
				}
			}

			void Analysis() { }
		}

}	// namespace md analysis

#endif
