#ifndef ANGLE_ANALYSIS_H
#define ANGLE_ANALYSIS_H

#include "h2o-analysis.h"
#include "so2-system-analysis.h"


namespace angle_analysis {

	using namespace md_analysis;


	template <typename T>
		class AngleAnalyzer {
			public:
				typedef Analyzer<T> system_t;

				AngleAnalyzer (system_t * t,
						std::string alphafile = std::string("alpha.dat"), 
						std::string betafile = std::string("beta.dat"))
					:
					_system(t), 
					_alpha(alphafile, 
							_system->posmin, _system->posmax, _system->posres, 
							system_t::angmin, system_t::angmax, system_t::angres),
					_beta(betafile, 
							_system->posmin, _system->posmax, _system->posres, 
							system_t::angmin, system_t::angmax, system_t::angres) 
					{
						double r[9] = {1,0,0,0,0,1,0,-1,0};
						_rotation = MatR(r);
					}


				virtual ~AngleAnalyzer () { }

				void Alpha (const double distance, const double angle) { _alpha (distance, angle); }
				void Beta (const double distance, const double angle) { _beta (distance, angle); }

				virtual void DataOutput() {
					_alpha.OutputData();
					_beta.OutputData();
				}

			protected:

				system_t *	_system;

				double angles[3];

				// histograms are indexed by [position, angle]
				Histogram2DAgent		_alpha;
				Histogram2DAgent		_beta;

				MatR _rotation;		// rotates a vector from the system x-y-z axes to the analysis frame where
				MatR _dcm;
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
					angles(t) 
			{ 
				h2os.ReferencePoint(WaterSystem<T>::SystemParameterLookup("analysis.reference-location"));
			}

				virtual void BinAngles (MolPtr mol);
				void Analysis ();
				void DataOutput () { angles.DataOutput(); }

			protected:
				h2o_analysis::H2OSystemManipulator<T>	h2os;
				AngleAnalyzer<T>											angles;
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







	template <typename T>
		class SO2AngleAnalysis : public AnalysisSet<T> {
			protected:
				h2o_analysis::H2OSystemManipulator<T>	h2os;
				so2_analysis::SO2SystemManipulator<T>	so2s;
				AngleAnalyzer<T>											angles;

			public:

				typedef Analyzer<T> system_t;

				SO2AngleAnalysis (system_t * t, std::string description = std::string("SO2 angle analysis"), std::string fn = std::string("")) :
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
			double distance = h2os.SurfaceLocation() - system_t::Position(so2->ReferencePoint());		// bottom surface
			//double distance = system_t::Position(so2->ReferencePoint()) - h2os.SurfaceLocation();

			// get the value of theta: molecular bisector angle with system reference axis
			angles.Alpha(distance, -(so2->Bisector() < VecR::UnitY()));	// bottom surface
			//angles.Alpha(distance, so2->Bisector() < VecR::UnitY());

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


	/************ so2 transit angle analysis ******************/


	template <typename T>
		class SO2TransitAngleAnalysis : public SO2AngleAnalysis<T> {

			public:

				typedef Analyzer<T> system_t;

				SO2TransitAngleAnalysis (system_t * t) :
					SO2AngleAnalysis<T> (t,
							std::string("SO2 transit angle analysis")) { 
						double pos = system_t::Position(this->so2s.SO2()->ReferencePoint()); 
						this->h2os.ReferencePoint(pos);
					}

				~SO2TransitAngleAnalysis () { } 

				virtual void Analysis ();
		};

	template <typename T>
		void SO2TransitAngleAnalysis<T>::Analysis () {
			this->h2os.Reload();
			this->h2os.FindWaterSurfaceLocation();
			this->BinAngles(this->so2s.SO2());
		}

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

				void FindInteractions (const AtomPtr);

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
		};

	template <typename T>
		void SO2AdsorptionWaterAngleAnalysis<T>::FindInteractions (const AtomPtr ref_atom) {
			// first sort the waters to find those closest to the so2 S
			nm.OrderAtomsByDistance (ref_atom);
			// then graph the closest several waters for analysis
			analysis_atoms.clear();
			Atom_it it = nm.closest(Atom::O);
			for (int i = 0; i < 10; i++) {
				analysis_atoms.push_back(*it);
				nm.next_closest(it, Atom::O);
			}
			analysis_atoms.push_back(ref_atom);
			// build a graph out of those atoms to find the interactions (if any) to the S
			graph.UpdateGraph(this->analysis_atoms); 
			bonded_atoms.clear();
			bonded_atoms = graph.BondedAtoms (ref_atom, bondgraph::interaction);
		}	



	template <typename T>
		void SO2AdsorptionWaterAngleAnalysis<T>::Analysis () {

			if (!second_pass) {
				FindInteractions (so2s.S());			// load up bonded_atoms with any atoms that are bound to the ref-atom

				if (bonded_atoms.size() && !first_bound_water) {
					std::cout << std::endl << "found the first interacting water at timestep " << this->_system->Timestep() << std::endl;
					first_bound_water = new Water(bonded_atoms[0]->ParentMolecule());
					first_bound_water->SetAtoms();
					first_bound_water->Print();
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
				fprintf (this->output, " % 8.3f ", system_t::Position(so2s.S()) - h2os.SurfaceLocation());

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
