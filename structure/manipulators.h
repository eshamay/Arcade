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
