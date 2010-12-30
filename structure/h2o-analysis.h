#ifndef H2O_ANALYSIS_H
#define H2O_ANALYSIS_H

#include "analysis.h"

namespace h2o_analysis {

	using namespace md_analysis;


	// ****************** H2O Manipulator ********************** //

	template <typename T>
		class H2OSystemManipulator : public SystemManipulator<T> {

			public:
				typedef Analyzer<T> system_t;

				H2OSystemManipulator (system_t * t) : 
					SystemManipulator<T>(t), 
					reference_point(WaterSystem<T>::SystemParameterLookup("analysis.reference-location")) { 

						this->_system->LoadWaters();
						// gather all the system waters into an analysis container
						for (Mol_it it = this->_system->int_wats.begin(); it != this->_system->int_wats.end(); it++) {
							WaterPtr wat (new Water(*(*it)));
							wat->SetAtoms();
							all_waters.push_back(wat);
						}

						all_water_atoms.clear();
						std::copy(this->_system->int_atoms.begin(), this->_system->int_atoms.end(), std::back_inserter(all_water_atoms));
						this->Reload();
					}

				virtual ~H2OSystemManipulator () { 
					for (Wat_it it = all_waters.begin(); it != all_waters.end(); it++)
						delete *it;
				}

				void Reload ();
				virtual void FindWaterSurfaceLocation ();		

				double ReferencePoint() const { return reference_point; }
				void ReferencePoint (const double point) { reference_point = point; }
				double SurfaceLocation () const { return surface_location; }

				Wat_it begin() { return analysis_waters.begin(); }
				Wat_it end() { return analysis_waters.end(); }

			protected:
				Water_ptr_vec all_waters, analysis_waters;
				Atom_ptr_vec all_water_atoms;

				double reference_point;	// the original location of the so2 along the reference axis
				double surface_location;	// location of the water surface along the reference axis

		};	// class H2OSystemManipulator


	template <typename T>
		void H2OSystemManipulator<T>::Reload () {
			analysis_waters.clear();
			std::copy (all_waters.begin(), all_waters.end(), std::back_inserter(analysis_waters));
			// now all_waters has... all the waters, and analysis wats is used to perform some analysis
			analysis_atoms.clear();
			std::copy (all_water_atoms.begin(), all_water_atoms.end(), std::back_inserter(analysis_atoms));
		}	// reload analysis wats


	template <typename T>
		void H2OSystemManipulator<T>::FindWaterSurfaceLocation () {
			// get rid of everything above the so2
			analysis_waters.erase(
					//remove_if(analysis_waters.begin(), analysis_waters.end(), system_t::MoleculeAbovePosition(reference_point, system_t::axis)), analysis_waters.end());
				remove_if(analysis_waters.begin(), analysis_waters.end(), system_t::MoleculeBelowPosition(reference_point, system_t::axis)), analysis_waters.end()); // bottom surface

			// sort the waters by position along the reference axis - first waters are lowest, last are highest
			std::sort (analysis_waters.begin(), analysis_waters.end(), system_t::molecule_position_pred(Atom::O));
			int numWats = 20;				// number of waters to use for calculating the location of the "top" of the water surface
			surface_location = 0.0;
			for (Wat_it it = analysis_waters.begin(); it != analysis_waters.begin() + numWats; it++) { // bottom surface
				//for (Wat_rit it = analysis_waters.rbegin(); it != analysis_waters.rbegin() + numWats; it++) {
				surface_location += system_t::Position((*it)->ReferencePoint());
				//printf ("% .3f\n", system_t::Position((*it)->ReferencePoint()));
			}
			surface_location /= numWats;
			}	// find surface water location


		}	// namespace md analysis



#endif
