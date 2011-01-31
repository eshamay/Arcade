#ifndef ANGLE_ANALYSIS_H
#define ANGLE_ANALYSIS_H

#include "h2o-analysis.h"
#include "so2-system-analysis.h"


namespace angle_analysis {

	using namespace md_analysis;


	template <typename T>
		class AngleHelper {
			public:
				typedef Analyzer<T> system_t;

				AngleHelper (system_t * t,
						double min1, double max1, double res1,
						double min2, double max2, double res2,
						std::string alphafile = std::string("alpha.dat"), 
						std::string betafile = std::string("beta.dat"))
					:
						_system(t), 
						_alpha(alphafile, min1, max1, res1, min2, max2, res2),
						_beta(betafile, min1, max1, res1, min2, max2, res2) 
						{
							double r[9] = {1,0,0,0,0,1,0,-1,0};
							_rotation = MatR(r);
						}


				virtual ~AngleHelper () { }

				void Alpha (const double val1, const double val2) { _alpha (val1, val2); }
				void Beta (const double val1, const double val2) { _beta (val1, val2); }

				double Alpha_TotalCount() const { return _alpha.TotalCount(); }
				double Beta_TotalCount() const { return _beta.TotalCount(); }

				virtual void DataOutput() {
					_alpha.OutputData();
					_beta.OutputData();
				}

			protected:

				system_t *	_system;

				// histograms are indexed by [position, angle]
				Histogram2DAgent		_alpha;
				Histogram2DAgent		_beta;

				MatR _rotation;		// rotates a vector from the system x-y-z axes to the analysis frame where
				MatR _dcm;

		};	 // angle helper



	template <typename T>
		class DistanceAngleHelper : public AngleHelper<T> {
			public:
				typedef Analyzer<T> system_t;

				DistanceAngleHelper (system_t * t)
					:	AngleHelper<T>(t, 
							this->_system->posmin, this->_system->posmax, this->_system->posres, 
							system_t::angmin, system_t::angmax, system_t::angres) { }

				virtual ~DistanceAngleHelper () { }
		};

	template <typename T>
		class AngleAngleHelper : public AngleHelper<T> {
			public:
				typedef Analyzer<T> system_t;

				AngleAngleHelper (system_t * t)
					:	AngleHelper<T>(t, 
							system_t::angmin, system_t::angmax, system_t::angres,
							system_t::angmin, system_t::angmax, system_t::angres) { }

				virtual ~AngleAngleHelper () { }
		};







	/************** H2O Angle Analysis **********************/
	/* An analysis to determine the spatial and angular distribution of waters in a slab system.
	 * Each water has an internal coordinate frame such that the molecular bisector points in the positive Z direction, and the molecular plane is in the molecular x-z plane.
	 * Additionally, a water's orientation is described by two angles - alpha & beta.
	 * Theta:
	 *		The angle formed between the water's molecular bisector and the system's positive Z-axis.
	 *		Theta ranges from 0 to 180 degrees, but is reported as cos(alpha) values ranging from +1 to -1
	 * Phi:
	 *		The angle of twist around the system Z-axis (referenced with the system X-axis as the 0 degree)
	 *		The angle is that one formed by the projection of the water's bisector onto the system's x-y plane and the positive system x axis.
	 *		Phi ranges from -180 to +180 degrees, but is reported in radians.
	 *
	 */

	template <typename T>
		class H2OAngleAnalysis : public AnalysisSet<T> {
			public:
				typedef Analyzer<T> system_t;

				H2OAngleAnalysis (system_t * t) :
					AnalysisSet<T>(t,
							std::string ("H2O Angle Analysis"),
							std::string ("")),
					h2os(t),
					angles(t) { 
						h2os.ReferencePoint(WaterSystem<T>::SystemParameterLookup("analysis.reference-location"));
					}

				virtual void BinAngles (MolPtr mol);
				virtual void Analysis ();
				void DataOutput () { angles.DataOutput(); }

			protected:
				h2o_analysis::H2OSystemManipulator<T>	h2os;
				DistanceAngleHelper<T>											angles;
		};



	template <typename T>
		void H2OAngleAnalysis<T>::BinAngles (MolPtr mol) {

			Water * wat = new Water(mol);
			wat->SetOrderAxes();
			double distance = system_t::Position(wat->ReferencePoint()) - h2os.SurfaceLocation();

			angles.Alpha(distance, wat->Bisector() < VecR::UnitY());
			// get the projection of the molecule's x-axis onto the system's x-y plane
			angles.Beta(distance, fabs(wat->Y() < VecR::UnitY()));

			// switch the lab frame axes to make the Y-axis the reference one instead of Z
			//r = _rotation * wat->Z();
			//r = wat->Z();
			//_dcm = wat->DCMToLab().transpose();
			//coordinate_conversion::DCM2EulerAngles_ZXZ (&_dcm(0,0), angles);

			delete wat;

			return;
		} // bin water angles



	template <typename T>
		void H2OAngleAnalysis<T>::Analysis () {

			h2os.FindWaterSurfaceLocation();
			h2os.Reload();

			for (Wat_it wat = h2os.begin(); wat != h2os.end(); wat++) {
				BinAngles(*wat);
			}
		}


	// *********** water OH angle analysis ****************
	template <typename T>
		class WaterOHAngleAnalysis : public AnalysisSet<T> {
			public:
				typedef Analyzer<T> system_t;

				WaterOHAngleAnalysis (system_t * t) :
					AnalysisSet<T>(t,
							std::string ("Water OH Angle Analysis - via SO2 transit"),
							std::string ("")),
					h2os(t),
					so2s(t),
					_alpha("alpha.dat", system_t::angmin, system_t::angmax, system_t::angres),
					_beta("beta.dat", system_t::angmin, system_t::angmax, system_t::angres),
					oh_calculator(VecR::UnitY()) { 
						h2os.ReferencePoint(WaterSystem<T>::SystemParameterLookup("analysis.reference-location"));
					}

				~WaterOHAngleAnalysis () { }
				virtual void Analysis ();
				void DataOutput () { _alpha.OutputData(); _beta.OutputData(); }

			protected:
				h2o_analysis::H2OSystemManipulator<T>	h2os;
				so2_analysis::SO2SystemManipulator<T>	so2s;
				md_analysis::Histogram1DAgent										_alpha;
				md_analysis::Histogram1DAgent										_beta;
				h2o_analysis::OHAngleCalculator											oh_calculator;
		};	// water oh angle analysis

	template <typename T>
		void WaterOHAngleAnalysis<T>::Analysis () {

			h2os.FindWaterSurfaceLocation();
			std::pair<double,double> p;

			if (h2os.TopSurface()) {
				for (Wat_rit wat = h2os.rbegin(); wat != h2os.rbegin()+70; wat++) {
					// for each water, find the angle that the oh bonds make with the surface normal.
					p = oh_calculator(*wat);

					_alpha (p.first);
					_beta (p.second);
				}
			}
			else {
				for (Wat_it wat = h2os.begin(); wat != h2os.begin()+70; wat++) {
					// for each water, find the angle that the oh bonds make with the surface normal.
					p = oh_calculator(*wat);

					_alpha (p.first);
					_beta (p.second);
				}
			}
		}



	// ************* so2 angle analysis *****************

	template <typename T>
		class SO2AngleAnalysis : public AnalysisSet<T> {
			protected:
				h2o_analysis::H2OSystemManipulator<T>	h2os;
				so2_analysis::SO2SystemManipulator<T>	so2s;
				DistanceAngleHelper<T>											angles;

			public:

				typedef Analyzer<T> system_t;

				SO2AngleAnalysis (system_t * t, 
						std::string description = std::string("SO2 angle analysis"), 
						std::string fn = std::string("")) :
					AnalysisSet<T> (t, description, fn),
					h2os(t), so2s(t), angles(t) { 
						h2os.ReferencePoint(WaterSystem<T>::SystemParameterLookup("analysis.reference-location"));
					}

				virtual ~SO2AngleAnalysis () { } 

				virtual void Analysis ();
				virtual void DataOutput () { angles.DataOutput(); }
				virtual void BinAngles (SulfurDioxide * so2);

		};	// so2 angle analysis



	template <typename T>
		void SO2AngleAnalysis<T>::BinAngles (SulfurDioxide * so2) {
			so2->SetOrderAxes();

			// calculate the position of the so2 relative to the water surface
			double distance;
			if (h2os.TopSurface()) {
				distance = system_t::Position(so2->ReferencePoint()) - h2os.SurfaceLocation();
				// get the value of theta: molecular bisector angle with system reference axis
				angles.Alpha(distance, so2->Bisector() < VecR::UnitY());
			} else {
				distance = h2os.SurfaceLocation() - system_t::Position(so2->ReferencePoint());		// bottom surface
				angles.Alpha(distance, -(so2->Bisector() < VecR::UnitY()));	// bottom surface
			}

			// get the value of phi - the molecular normal angle with the system ref
			angles.Beta(distance, fabs(so2->Y() < VecR::UnitY()));

			return;
		}



	template <typename T>
		void SO2AngleAnalysis<T>::Analysis () {

			h2os.Reload();
			h2os.FindWaterSurfaceLocation();

			for (std::vector<SulfurDioxide *>::iterator so2 = so2s.begin(); so2 != so2s.end(); so2++) {
				BinAngles(*so2);
			}
		}





	/************ so2 transit first-binding water angle analysis ************/

	template <typename T>
		class SO2AdsorptionWaterAngleAnalysis : public AnalysisSet<T> {
			public:
				typedef Analyzer<T> system_t;

				SO2AdsorptionWaterAngleAnalysis(system_t * t) 
					: 
						AnalysisSet<T> (t,
								std::string("Analysis of waters near an adsorbing so2"),
								std::string("first-adsorption-water.dat")),
						h2os(t), so2s(t), nm(t), first_bound_water ((WaterPtr)NULL), second_pass(false)
																																													 //angle_histo ("angle.dat", 0.0, 10.0, 0.1, -1.0, 1.0, 0.05)	// distance from 0 to 10 angstroms
			{ }
				~SO2AdsorptionWaterAngleAnalysis () {
					if (!first_bound_water) delete first_bound_water;
				}

				virtual void Analysis ();
				//virtual void DataOutput () { angle_histo.OutputData(); }

				void FindInteractions ();

			protected:
				h2o_analysis::H2OSystemManipulator<T>	h2os;
				so2_analysis::SO2SystemManipulator<T>	so2s;
				neighbor_analysis::NeighborManipulator<T>	nm;

				bondgraph::BondGraph	graph;
				std::vector<double>		so2_distances;
				std::vector<double>		water_angles;
				WaterPtr							first_bound_water;
				Atom_ptr_vec					bonded_atoms;
				Atom_ptr_vec					analysis_atoms;
				bool									second_pass;

				// predicate for checking atom residue types
				static bool residue_eq (const AtomPtr a, const std::string& res);
		};

	template <typename T>
		bool SO2AdsorptionWaterAngleAnalysis<T>::residue_eq (const AtomPtr a, const std::string& res) {
			if (a->Residue() == res) 
				std::cout << "found one!" << std::endl;
			return a->Residue() == res;
		}

	template <typename T>
		void SO2AdsorptionWaterAngleAnalysis<T>::FindInteractions () {
			// first sort the waters to find those closest to the so2 S
			nm.OrderAtomsByDistance (so2s.S());
			// then graph the closest several waters for analysis
			analysis_atoms.clear();
			Atom_it close_it = nm.closest(Atom::O);
			for (int i = 0; i < 10; i++) {
				analysis_atoms.push_back(*close_it);
				nm.next_closest(close_it, Atom::O);
			}
			analysis_atoms.push_back(so2s.S());
			// build a graph out of those atoms to find the interactions (if any) to the S
			graph.UpdateGraph(this->analysis_atoms); 

			// copy all the atoms bound to the S and the two Os
			bonded_atoms.clear();
			Atom_ptr_vec bound_atoms;

			bound_atoms = graph.BondedAtoms (so2s.S(), bondgraph::interaction);
			std::copy (bound_atoms.begin(), bound_atoms.end(), std::back_inserter(bonded_atoms));
			bound_atoms = graph.BondedAtoms (so2s.O1(), bondgraph::hbond);
			std::copy (bound_atoms.begin(), bound_atoms.end(), std::back_inserter(bonded_atoms));
			bound_atoms = graph.BondedAtoms (so2s.O2(), bondgraph::hbond);
			std::copy (bound_atoms.begin(), bound_atoms.end(), std::back_inserter(bonded_atoms));
			// remove duplicates
			std::sort(bonded_atoms.begin(), bonded_atoms.end(), std::ptr_fun(&Atom::id_cmp));
			Atom_it uniq_it = std::unique (bonded_atoms.begin(), bonded_atoms.end(), std::ptr_fun(&Atom::id_eq));
			bonded_atoms.resize(uniq_it - bonded_atoms.begin());

			// only keep waters that are bound, not other types of molecules
			bound_atoms.clear();
			for (Atom_it it = bonded_atoms.begin(); it != bonded_atoms.end(); it++) {
				if ((*it)->Residue() == "h2o")
					bound_atoms.push_back(*it);
			}
			bonded_atoms.clear();
			std::copy(bound_atoms.begin(), bound_atoms.end(), std::back_inserter(bonded_atoms));

			/*
				 bonded_atoms.erase(
				 std::remove_if(
				 bonded_atoms.begin(), 
				 bonded_atoms.end(), 
			//std::not1(
			std::bind2nd(std::ptr_fun(&SO2AdsorptionWaterAngleAnalysis<T>::residue_eq), "h2o")), 
			bonded_atoms.end());
			 */
		}	



	template <typename T>
		void SO2AdsorptionWaterAngleAnalysis<T>::Analysis () {

			if (!second_pass) {
				FindInteractions ();			// load up bonded_atoms with any atoms that are bound to the ref-atom

				if (bonded_atoms.size() && !first_bound_water) {
					std::cout << std::endl << "found the first interacting water at timestep " << this->_system->Timestep() << std::endl;
					first_bound_water = new Water(bonded_atoms[0]->ParentMolecule());	// Use the first of the bound atoms as it may be the closest... ?
					first_bound_water->Print();
					first_bound_water->SetAtoms();
					so2s.S()->Print();

					std::cout << "now rewinding and rerunning the analysis" << std::endl;
					this->_system->Rewind();
					std::cout << "system is rewound - starting over..." << std::endl;
					second_pass = true;
				}
			}

			else {
				// find distance of so2 to the water surface
				h2os.Reload();
				h2os.FindWaterSurfaceLocation();

				// output the posisition of the so2 wrt the water surface
				double so2_location;
				if (h2os.TopSurface()) {
					so2_location = system_t::Position(so2s.S()) - h2os.SurfaceLocation();
				} else {
					so2_location = h2os.SurfaceLocation() - system_t::Position(so2s.S())  ;
				}
				fprintf (this->output, " % 8.3f ", so2_location);

				// also the distance from the S to the O
				fprintf (this->output, " % 8.3f", MDSystem::Distance(first_bound_water->O(), so2s.S()).norm());

				// the standard deviation of distance from the surface of waters used in calculating the surface
				fprintf (this->output, " % 8.3f", h2os.SurfaceWidth());

				// distance of the reference water to the surface
				fprintf (this->output, " % 8.3f", system_t::Position(first_bound_water) - h2os.SurfaceLocation());

				// calculate the angle of the water of interest with the surface normal
				fprintf (this->output, " % 8.3f", first_bound_water->Bisector() < VecR::UnitY());

				// the (cos) angle between the h2o and so2 bisectors
				fprintf (this->output, " % 8.3f", first_bound_water->Bisector() < so2s.SO2()->Bisector());


				fprintf (this->output, "\n");
			}


			//std::cout << std::endl;
			//std::transform (h2os.begin(), h2os.begin()+20, std::ostream_iterator<double>(std::cout, " "), system_t::molecule_distance_generator(so2s.S()));
			//std::cout << std::endl;

		}	// analysis

}	// namespace angle analysis

#endif
