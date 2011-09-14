#ifndef MANIPULATORS_H_
#define MANIPULATORS_H_

#include "analysis.h"
#include <gsl/gsl_statistics.h>

namespace md_analysis {

	// manipulator superclass
	class SystemManipulator {
		public:

			SystemManipulator (Analyzer * sys);
			virtual ~SystemManipulator () { }

			// reload all the analysis atoms and molecules
			virtual void Reload ();

			Atom_it atoms_begin() { return analysis_atoms.begin(); }
			Atom_it atoms_end() { return analysis_atoms.end(); }

			Mol_it mols_begin() { return analysis_mols.begin(); }
			Mol_it mols_end() { return analysis_mols.end(); }

		protected:
			Analyzer * _system;

			Atom_ptr_vec		all_atoms;
			Mol_ptr_vec			all_mols;

			Atom_ptr_vec		analysis_atoms;
			Mol_ptr_vec			analysis_mols;
	};	// system manipulator





	/******************* Neighbor Manipulator ******************/
	class NeighborManipulator : public SystemManipulator {
		public:

			typedef Analyzer system_t;

			NeighborManipulator (system_t * sys) : SystemManipulator(sys) { }
			~NeighborManipulator() { }

			// Order all the analysis atoms by distance starting with the closest to the given atom
			void OrderAtomsByDistance (AtomPtr ap) {
				this->Reload();
				reference_atom = ap;
				std::sort(
						this->analysis_atoms.begin(), 
						this->analysis_atoms.end(), 
						typename system_t::atomic_distance_cmp (ap));
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

			// retrieve the closest atoms to the given reference_atom
			Atom_range GetClosestAtoms (AtomPtr a, const int num) {
				this->OrderAtomsByDistance(a);
				return std::make_pair(this->analysis_atoms.begin(), this->analysis_atoms.begin() + num);
			}

			void next_closest (Atom_it& it) { it++; }

		private:
			AtomPtr									reference_atom;

	}; // neighbor manipulator











	/************* Cycle Manipulator *****************/

	class CycleManipulator : public SystemManipulator {
		public:
			typedef Analyzer system_t;

			CycleManipulator (system_t * t) : 
				SystemManipulator(t),
				nm(t),
				cycle_size(0), ref_atom((AtomPtr)NULL), cycle_type(NO_CYCLE) { }


			typedef enum {
				NO_CYCLE = 0,	// no cycle present in the local graph of the so2 neighbors
				FULLBRIDGE=1,		// A complete cycle links the two so2 oxygens
				WATERLEG=2,			// a cycle exists but does not link the two so2 oxygens - the cycle is in the waters off one so2-O
				HALFBRIDGE=3,			// a cycle runs from the S to a water-O and back to the so2-O
				FULLCROWN=4,			// S connected to two water-Os
				HALFCROWN=5,			// cycle happens in waters and is only peripherally connected to the S (like for a waterleg)
				UNKNOWN				// some other type that isn't figured out right now
			}	cycle_t;


			typedef std::list<AtomPtr> cycle_list;
			typedef cycle_list::const_iterator cycle_it;

			typedef typename std::list<cycle_list>::const_iterator cycle_list_it;
			typedef typename std::list<cycle_t>::const_iterator cycle_type_it;

			typedef std::pair<cycle_it,cycle_it> cycle_pair_t;


			void BuildGraph ();	// builds the graph of the given atoms closest to the reference atom
			void ParseCycles ();	// finds if a cycle exists in the topology, and parses it
			void FindUniqueMembers (const cycle_list&);

			void SetCycleSize (const int size) { cycle_size = size; }
			void SetReferenceAtom (const AtomPtr ref) { ref_atom = ref; }

			AtomPtr ReferenceAtom () const { 
				if (ref_atom == (AtomPtr)NULL) {
					std::cerr << "CycleManipulator -- Reference atom was never set before use!" << std::endl;
					exit(1);
				}
				return ref_atom; 
			}

			int CycleSize () const {
				if (cycle_size == 0) {
					std::cerr << "CycleManipulator -- Cycle size was never set before use!" << std::endl;
					exit(1);
				}
				return cycle_size; 
			}

			cycle_list_it	cycle_begin() const { return cycle.begin(); }
			cycle_list_it	cycle_end() const { return cycle.end(); }

			cycle_type_it		cycle_type_begin() const { return cycle_type.begin(); }
			cycle_type_it		cycle_type_end() const { return cycle_type.end(); }

			size_t NumUniqueAtomsInCycle () const { return unique_cycle_atoms.size(); }
			size_t NumUniqueMoleculesInCycle () const { return unique_cycle_mols.size(); }
			int NumUniqueWaterAtoms () const;

			int NumReferenceAtomHBonds () const { return graph.NumHBonds(ref_atom); }
			int NumReferenceInteractions () const { return graph.NumInteractions(ref_atom); }

			cycle_pair_t CyclePointAtom (const cycle_list&);		// used to find the first atom that appears both in the bridge to a leg-cycle, and in the cycle itself
			std::pair<int,int> WaterLegInformation (const cycle_list&);

		private:
			NeighborManipulator	nm;

			int							cycle_size;
			AtomPtr						ref_atom;
			std::list<cycle_t>					cycle_type;

			bondgraph::BondGraph		graph;
			std::list<cycle_list>		cycle;				// all the atoms comprising the cycle
			Atom_ptr_vec						unique_cycle_atoms;	// only unique atoms in the cycle
			Mol_ptr_vec							unique_cycle_mols;

	}; // cycle manipulator

} // namespace md analysis





namespace h2o_analysis {

	using namespace md_analysis;

	// ****************** H2O Manipulator ********************** //

	class H2OSystemManipulator : public SystemManipulator {

		public:
			typedef Analyzer system_t;

			H2OSystemManipulator (system_t * t, const int number_of_waters_for_surface_calc = 70);

			virtual ~H2OSystemManipulator () { 
				for (Wat_it it = all_waters.begin(); it != all_waters.end(); it++)
					delete *it;
			}

			// find all the system waters and load them into the useable containers
			void Reload ();

			void UpdateAnalysisWaters ();

			virtual void FindWaterSurfaceLocation ();		
			virtual void FindClosestWaters (const AtomPtr);

			double ReferencePoint() const { return reference_point; }
			void ReferencePoint (const double point) { reference_point = point; }
			double SurfaceLocation () const { return surface_location; }	// mean location
			double SurfaceWidth () const { return surface_width; }	// standard deviation

			bool TopSurface () const { return top_surface; }

			Atom_it begin_atoms() const { return this->analysis_atoms.begin(); }
			Atom_it end_atoms() const { return this->analysis_atoms.end(); }

			Wat_it begin() { return analysis_waters.begin(); }
			Wat_it end() { return analysis_waters.end(); }
			Wat_rit rbegin() { return analysis_waters.rbegin(); }
			Wat_rit rend() { return analysis_waters.rend(); }


			// center of mass of all the waters
			VecR CenterOfMass () { 
				CalcCenterOfMass();
				return center_of_mass; 
			}

			void CalcCenterOfMass ();

			// function object that returns the distance of a given molecule to the given reference point
			class ReferenceDistance : public std::unary_function<MolPtr, double> {
				private:
					double reference;
					bool top_surface;

				public:
					ReferenceDistance (const double ref_position, bool top) : reference (ref_position), top_surface(top) { }
					double operator () (const MolPtr mol) const {
						double result;

						if (top_surface)
							result = system_t::Position(mol->ReferencePoint()) - this->reference; // top surface

						else
							result = reference - system_t::Position(mol->ReferencePoint()); // bottom surface

						return result;
					}
			};

		protected:
			Water_ptr_vec all_waters, analysis_waters;
			Atom_ptr_vec all_water_atoms;

			double reference_point;	// the original location of the so2 along the reference axis
			int number_surface_waters;
			bool	top_surface;

			VecR center_of_mass;

			double surface_location;	// location of the water surface along the reference axis
			double surface_width;			// standard deviation of the positions of waters used to calculate the surface_location

			// functor to grab a water's location based on the oxygen position
			class WaterLocation : public std::unary_function <WaterPtr, double> {
				public:
					double operator() (const WaterPtr wat) const {
						return system_t::Position(wat);
					}
			}; // water location

	};	// class H2OSystemManipulator


	typedef std::pair<bool,double>	surface_distance_t;

	class H2ODoubleSurfaceManipulator : public H2OSystemManipulator {
		protected:
			double top_location, bottom_location;
			double top_width, bottom_width;

		public:

			H2ODoubleSurfaceManipulator (system_t * t, const int number_of_waters_for_surface_calc = 70) :
				H2OSystemManipulator (t, number_of_waters_for_surface_calc) 
		{ }

			virtual void FindWaterSurfaceLocation ();		

			double TopSurfaceLocation () const { return top_location; }
			double BottomSurfaceLocation () const { return bottom_location; }

			double TopSurfaceWidth () const { return top_width; }
			double BottomSurfaceWidth () const { return bottom_width; }

			double WrappedDistance (double, double, double) const;
			// is the given value closer to the top or the bottom surface?
			surface_distance_t TopOrBottom (const double) const;
	};


} // namespace h2o analysis





namespace so2_analysis {

	using namespace md_analysis;

	typedef SulfurDioxide*								so2_ptr;
	typedef std::vector<SulfurDioxide *>	so2_vec;
	typedef so2_vec::iterator							so2_it;
	typedef so2_vec::reverse_iterator			so2_rit;

	// a convenience class for working with systems comprised of at least 1 SO2 molecule and a whole bunch of waters
	class SO2SystemManipulator : public SystemManipulator {
		public:
			typedef Analyzer system_t;

			SO2SystemManipulator (system_t * t) : SystemManipulator(t) { }

			virtual void Initialize ();

			virtual ~SO2SystemManipulator () { 
				delete this->so2;
				for (std::vector<SulfurDioxide *>::iterator it = so2s.begin(); it != so2s.end(); it++) {
					delete *it;
				}
			}

			virtual void UpdateSO2 ();

			SulfurDioxide * SO2 () { return so2; }
			AtomPtr S () { return so2->S(); }
			AtomPtr O1 () { return so2->O1(); }
			AtomPtr O2 () { return so2->O2(); }

			so2_it begin () { return so2s.begin(); }
			so2_it end ()		{ return so2s.end(); }
			so2_rit rbegin () { return so2s.rbegin(); }
			so2_rit rend () { return so2s.rend(); }

		protected:
			SulfurDioxide * so2;	// the sulfur dioxide of interest
			std::vector<SulfurDioxide *> so2s;

			// The method by which the SO2 of interest is found in the system
			virtual void FindSO2 ();
			virtual void FindAllSO2s ();

	}; // so2 system manipulator


	class XYZSO2Manipulator : public SO2SystemManipulator {

		public:
			typedef Analyzer system_t;

			XYZSO2Manipulator (system_t * t) : SO2SystemManipulator(t) { }

			virtual void Initialize () {
				this->_system->LoadAll();
				this->FindSO2 ();
			}

			void UpdateSO2 () {
				delete this->so2;
				this->Initialize();
			}

			// get the wannier centers associated with a particular atom
			vector_map_vec GetWanniers (const AtomPtr& atom, const int num) const;

		private:

			virtual void FindSO2 ();
			virtual void FindAllSO2s() { }

	};	// xyz so2 manipulator

} // namespace so2 analysis

#endif
