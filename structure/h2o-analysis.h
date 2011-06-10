#ifndef H2O_ANALYSIS_H
#define H2O_ANALYSIS_H

#include "analysis.h"
#include <gsl/gsl_statistics.h>

namespace h2o_analysis {

	using namespace md_analysis;


	// ****************** H2O Manipulator ********************** //

	template <typename T>
		class H2OSystemManipulator : public SystemManipulator<T> {

			public:
				typedef Analyzer<T> system_t;

				H2OSystemManipulator (system_t * t, const int number_of_waters_for_surface_calc = 70) : 
					SystemManipulator<T>(t), 
					reference_point(WaterSystem<T>::SystemParameterLookup("analysis.reference-location")),
					top_surface(WaterSystem<T>::SystemParameterLookup("analysis.top-surface")),
					number_surface_waters(number_of_waters_for_surface_calc) { 
						this->Reload();
					}

				virtual ~H2OSystemManipulator () { 
					for (Wat_it it = all_waters.begin(); it != all_waters.end(); it++)
						delete *it;
				}


				// find all the system waters and load them into the useable containers
				void Reload () {
					// first clear out the previous water set
					for (Wat_it wat = this->all_waters.begin(); wat != this->all_waters.end(); wat++) {
						delete *wat;
					}
					all_waters.clear();
					// then load in the new water set
					this->_system->LoadWaters();
					// gather all the system waters
					for (Mol_it it = this->_system->int_wats.begin(); it != this->_system->int_wats.end(); it++) {
						WaterPtr wat (new Water(*it));
						wat->SetAtoms();
						all_waters.push_back(wat);
					}

					// grab all the water atoms
					all_water_atoms.clear();
					std::copy(this->_system->int_atoms.begin(), this->_system->int_atoms.end(), std::back_inserter(all_water_atoms));

					// now update the analysis containers
					this->UpdateAnalysisWaters();
					return;
				}


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


	template <typename T>
		void H2OSystemManipulator<T>::UpdateAnalysisWaters () {

			this->analysis_waters.clear();
			std::copy (all_waters.begin(), all_waters.end(), std::back_inserter(analysis_waters));
			//std::for_each(analysis_waters.begin(), analysis_waters.end(), std::mem_fun(&Water::SetAtoms));
			// now all_waters has... all the waters, and analysis wats is used to perform some analysis
			this->analysis_atoms.clear();
			std::copy (all_water_atoms.begin(), all_water_atoms.end(), std::back_inserter(this->analysis_atoms));
		}	// reload analysis wats


	template <typename T>
		void H2OSystemManipulator<T>::FindWaterSurfaceLocation () {

			// get rid of everything above (or below) the reference point
			if (top_surface) {
				analysis_waters.erase(
						remove_if(analysis_waters.begin(), analysis_waters.end(), typename system_t::MoleculeAbovePosition(reference_point, system_t::axis)), analysis_waters.end());
			}

			else if (!top_surface) {
				analysis_waters.erase(
						remove_if(analysis_waters.begin(), analysis_waters.end(), typename system_t::MoleculeBelowPosition(reference_point, system_t::axis)), analysis_waters.end()); // bottom surface
			}

			// sort the waters by position along the reference axis - first waters are lowest, last are highest
			std::sort (analysis_waters.begin(), analysis_waters.end(), typename system_t::molecule_position_pred(Atom::O));

			if (top_surface) {
				// get the surface waters at the beginning of the list
				std::reverse(analysis_waters.begin(), analysis_waters.end());
			}	// top surface

			// resize the list to contain only the surface waters
			analysis_waters.resize(number_surface_waters);

			// get the position of the top-most waters
			std::vector<double> surface_water_positions;
			// grab all the locations
			std::transform (analysis_waters.begin(), analysis_waters.end(), std::back_inserter(surface_water_positions), WaterLocation());

			// calculate the statistics
			surface_location = gsl_stats_mean (&surface_water_positions[0], 1, number_surface_waters);
			surface_width = gsl_stats_sd (&surface_water_positions[0], 1, number_surface_waters);

			if (surface_width > 2.5) {
				std::cout << std::endl << "Check the pbc-flip or the reference point settings and decrease/increase is to fix this gigantic surface width" << std::endl;
				std::cout << "Here's the positions of the waters used to calculate the surface:" << std::endl;
				std::copy (surface_water_positions.begin(), surface_water_positions.end(), std::ostream_iterator<double>(std::cout, " "));
				printf ("\nSurface width = % 8.3f, Reference point = % 8.3f\n", surface_width, reference_point);
				fflush(stdout);
			}

		}	// find surface water location


	template <typename T>
		void H2OSystemManipulator<T>::CalcCenterOfMass () {
			center_of_mass.Zero();
			double mass = 0.0;


			for (Atom_it it = this->begin_atoms(); it != this->end_atoms(); it++) {
				mass += (*it)->Mass();
				center_of_mass += (*it)->Position() * (*it)->Mass();
			}
			center_of_mass = center_of_mass / mass;
			return;
		}


	// sort all the waters in the system by distance to a given reference atom
	template <typename T>
		void H2OSystemManipulator<T>::FindClosestWaters (const AtomPtr a) {
			this->UpdateAnalysisWaters();
			std::sort(analysis_waters.begin(), analysis_waters.end(), typename system_t::molecule_distance_cmp(a));
		} // find closest waters



	// functor takes a water and returns the values of the cos(angles) formed between the two oh-vectors. The first value of the pair is always the greater (magnitude) of the two values.
	class OHAngleCalculator : public std::unary_function <WaterPtr,std::pair<double,double> > {
		private:
			VecR axis;	// the reference axis to which the angles will be formed
		public:
			OHAngleCalculator (const VecR ax) : axis(ax) { }
			std::pair<double,double> operator() (const WaterPtr& wat) {
				double angle1 = wat->OH1() < axis;
				double angle2 = wat->OH2() < axis;
				std::pair<double,double> p = (fabs(angle1) > fabs(angle2)) 
					? std::make_pair(angle1,angle2) 
					: std::make_pair(angle2,angle1);
				return p;
			}
	};


}	// namespace h2o analysis



#endif
