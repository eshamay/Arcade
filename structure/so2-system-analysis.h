#ifndef SO2_ANALYSIS_H_
#define SO2_ANALYSIS_H_

#include "h2o-analysis.h"

namespace so2_analysis {

	using namespace md_analysis;


	// a convenience class for working with systems comprised of at least 1 SO2 molecule and a whole bunch of waters
	template <typename T>
		class SO2SystemAnalyzer : public h2o_analysis::H2OSystemAnalyzer<T> {
			public:
				typedef Analyzer<T> system_t;

				SO2SystemAnalyzer (std::string desc, std::string filename) :
					h2o_analysis::H2OSystemAnalyzer<T> (desc, filename) { }

				virtual ~SO2SystemAnalyzer () { 
					delete this->so2;
				}

				// The method by which the SO2 of interest is found in the system
				virtual void FindSO2 (system_t& t) {
					// find the so2 in the system and set some pointers up
					MolPtr mol = Molecule::FindByType(t.sys_mols, Molecule::SO2);
					this->so2 = new SulfurDioxide(mol);
				}

				// things to do after the standard setup
				virtual void PostSetup (system_t& t) { 
					this->FindSO2 (t);

					this->so2->SetAtoms();
					this->s = so2->S();
					this->o1 = so2->O1();
					this->o2 = so2->O2();

					// grab the first location of the so2 as a reference for later analyses
					this->reference_point = system_t::Position(this->s);
				}

			protected:
				SulfurDioxide * so2;	// the sulfur dioxide of interest
				AtomPtr s, o1, o2;	// sulfur dioxide's atoms
		};
		

	template <typename T>
		class SO2AngleAnalysis : public h2o_analysis::H2OSystemAngleAnalyzer<T> {
			protected:
				std::vector<SulfurDioxide *> so2s;

			public:

				typedef Analyzer<T> system_t;

				SO2AngleAnalysis () :
					h2o_analysis::H2OSystemAngleAnalyzer<T> (
							std::string("SO2 angle analysis")) { }

				~SO2AngleAnalysis () {
					for (std::vector<SulfurDioxide *>::iterator it = so2s.begin(); it != so2s.end(); it++) {
						delete *it;
					}
				}

				void Analysis (system_t& t);
				virtual void FindWaterSurfaceLocation ();		
				void BinAngles (MolPtr mol);

				void PostSetup (system_t& t) {
					this->reference_point = WaterSystem<T>::SystemParameterLookup("analysis.reference-location");

					t.LoadAll();
					so2s.clear();
					// load all the system so2s into the container - in case
					for (Mol_it it = t.sys_mols.begin(); it != t.sys_mols.end(); it++) {
						if ((*it)->Name() == "so2") {
							SulfurDioxide * mol = new SulfurDioxide(*it);
							mol->SetAtoms();
							so2s.push_back(mol);
						}
					}
					printf ("\nFound %zu Sulfur Dioxides in the system\n", so2s.size());
				}

		};	// so2 angle analysis


	template <typename T>
		void SO2AngleAnalysis<T>::FindWaterSurfaceLocation () {
			// get rid of everything above the so2
			this->analysis_wats.erase(
					remove_if(this->analysis_wats.begin(), this->analysis_wats.end(), system_t::MoleculeAbovePosition(this->reference_point, WaterSystem<T>::axis)), this->analysis_wats.end());

			// sort the waters by position along the reference axis - first waters are lowest, last are highest
			std::sort (this->analysis_wats.begin(), this->analysis_wats.end(), system_t::molecule_position_pred(Atom::O));
			int numWats = 20;				// number of waters to use for calculating the location of the "top" of the water surface
			this->surface_location = 0.0;
			for (Wat_it it = this->analysis_wats.begin(); it != this->analysis_wats.begin() + numWats; it++) { // bottom surface
			//for (Wat_it it = this->analysis_wats.end()-1; it != this->analysis_wats.end() - numWats; it--) {
				this->surface_location += system_t::Position((*it)->ReferencePoint());
				//printf ("% .3f\n", system_t::Position((*it)->ReferencePoint()));
			}
			this->surface_location /= numWats;
		}	// find surface water location

	template <typename T>
		void SO2AngleAnalysis<T>::BinAngles (MolPtr mol) {
			SulfurDioxide * so2 = new SulfurDioxide(mol);
			//SulfurDioxide * so2 = static_cast<SulfurDioxide *>(mol);
			so2->SetOrderAxes();
			double distance = this->surface_location - system_t::Position(so2->ReferencePoint());		// bottom surface
			//double distance = system_t::Position(so2->ReferencePoint()) - this->surface_location;

			so2->SetOrderAxes();
			this->r = so2->Bisector();

			this->_alpha(distance, -(this->r < VecR::UnitY()));		// bottom surface
			//this->_alpha(distance, (this->r < VecR::UnitY()));
			// get the projection of the molecule's x-axis onto the system's x-y plane
			this->_beta(distance, fabs(so2->Y() < VecR::UnitY()));

			delete so2;
			return;
		}

	template <typename T>
		void SO2AngleAnalysis<T>::Analysis (system_t& t) {

			this->FindWaterSurfaceLocation();

			for (std::vector<SulfurDioxide *>::iterator so2 = so2s.begin(); so2 != so2s.end(); so2++) {
				BinAngles(*so2);
			}
		}

}	// namespace so2_analysis


#endif
