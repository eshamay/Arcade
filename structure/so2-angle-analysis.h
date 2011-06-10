#ifndef SO2_ANGLE_ANALYSIS_H
#define SO2_ANGLE_ANALYSIS_H

#include "histogram-analysis.h"
#include "h2o-analysis.h"
#include "so2-system-analysis.h"

namespace so2_angle_analysis {

	using namespace md_analysis;

		// ************* so2 angle analysis *****************

		template <typename T>
			class SO2AngleAnalyzer : public AnalysisSet<T> {
				protected:
					h2o_analysis::H2OSystemManipulator<T>	h2os;
					so2_analysis::XYZSO2Manipulator<T>		so2s;
					VecR ref_ax;

				public:

					typedef Analyzer<T> system_t;

					SO2AngleAnalyzer (system_t * t, 
							std::string description,
							std::string fn = std::string("")) :
						AnalysisSet<T> (t, description, fn),
						h2os(t,10), so2s(t), ref_ax(VecR::UnitZ()) {
							so2s.Initialize();
						}

					virtual ~SO2AngleAnalyzer () { } 

					virtual void Analysis () = 0;
					virtual void DataOutput () = 0;
					virtual void BinAngles () = 0;

					void SO2SetupAndBin ();
					void SetupAllAndBin ();
					void SO2AngleBinner ();
					virtual void UpdateHistograms (double * values) { }
			};	// so2 angle analysis




		template <typename T>
			void SO2AngleAnalyzer<T>::SO2SetupAndBin () {

				this->so2s.UpdateSO2();
				this->so2s.SO2()->SetOrderAxes();

				// now get the vector pointing from the center of mass to the so2
				this->BinAngles();
			}

		template <typename T>
			void SO2AngleAnalyzer<T>::SetupAllAndBin () {
				this->h2os.Reload();
				this->SO2SetupAndBin();
			}

		template <typename T>
			void SO2AngleAnalyzer<T>::SO2AngleBinner () {

				SulfurDioxide * so2 = this->so2s.SO2();

				double values[2];

				// find the angle of the so2 theta and phi wrt the reference axis
				double angle = so2->Z() < ref_ax;
				angle = 180. * acos(angle) / M_PI;
				values[0] = angle;

				// same for the phi angle
				double angle2 = fabs(so2->Y() < ref_ax);
				angle2 = 180. * acos(angle2) / M_PI;
				if (angle2 > 90.0)
					angle2 = 180.0 - angle2;
				values[1] = angle2;

				this->UpdateHistograms(values);
				return;
			}





		template <typename T>
			class SO2Angles2D : public SO2AngleAnalyzer<T> {
				protected:
					Histogram2DAgent	histo;

				public:

					typedef Analyzer<T> system_t;

					SO2Angles2D (system_t * t) : 
						SO2AngleAnalyzer<T> (t, std::string ("2D so2 angle analysis")),
						histo (std::string ("so2-angles.theta+phi.2d.dat"), 0.0, 180.0, 1.0, 0.0, 90.0, 0.5) { }

					virtual ~SO2Angles2D () { } 

					virtual void Analysis () { this->SO2SetupAndBin(); }
					virtual void DataOutput () { histo.OutputData(); }
					virtual void BinAngles () { this->SO2AngleBinner(); }
					void UpdateHistograms (double * values) { histo(values[0], values[1]); }

			};	// so2 2D angle analysis



		template <typename T>
			class SO2Angles1D : public SO2AngleAnalyzer<T> {
				protected:
					Histogram1DAgent		theta, phi;

				public:

					typedef Analyzer<T> system_t;

					SO2Angles1D (system_t * t) : 
						SO2AngleAnalyzer<T> (t, std::string ("1D so2 angle analysis")),
						theta (std::string ("so2-theta.z.dat"), 0.0, 180.0, 1.0),
						phi (std::string ("so2-phi.z.dat"), 0.0, 90.0, 0.5) { }
						

					virtual ~SO2Angles1D () { } 

					virtual void Analysis () { this->SO2SetupAndBin(); }
					virtual void DataOutput () { theta.OutputData(); phi.OutputData(); }
					virtual void BinAngles () { this->SO2AngleBinner(); }

					virtual void UpdateHistograms (double * values) {
						theta(values[0]);
						phi(values[1]);
					}

			};	// so2 2D angle analysis



		// same ADF as Moin2011, but with the 2nd dimension set to the O---H bondlength
		template <typename T>
			class SO2_H2O_Angles2D : public SO2AngleAnalyzer<T> {
				protected:
					Histogram2DAgent	histo;

				public:

					typedef Analyzer<T> system_t;

					SO2_H2O_Angles2D (system_t * t) : 
						SO2AngleAnalyzer<T> (t, std::string ("2D so2-h2o hbond angle analysis")),
						histo (std::string ("so2-angles.Os-Hw-Ow+Os-Hw.2d.dat"), 0.0, 180.0, 1.0, 0.0, 3.0, 0.1) { }

					virtual ~SO2_H2O_Angles2D () { } 

					virtual void Analysis () { this->SetupAllAndBin(); }
					virtual void DataOutput () { histo.OutputData(); }
					virtual void BinAngles ();
					void UpdateHistograms (double * values) { return; }

			};	// so2 2D angle analysis

		template <typename T>
			void SO2_H2O_Angles2D<T>::BinAngles () {

				bondgraph::BondGraph& graph = this->_system->BondGraph();
				graph.UpdateGraph(this->begin(), this->end());

				// find all the waters that are H-bonded to the so2 through the oxygens
				VecR OsH, HOw;
				WaterPtr wat;
				double angle;
				double distance;

				AtomPtr o = this->so2s.O1();
				Atom_ptr_vec bonded_atoms = graph.BondedAtoms (o, bondgraph::hbond);
				for (Atom_it it = bonded_atoms.begin(); it != bonded_atoms.end(); it++) {
					wat = new Water((*it)->ParentMolecule());
					wat->SetAtoms();

					// for each of the H-bonded waters, find the Os-Hw vector, and the Hw-Ow vector, and the angle between them
					HOw = MDSystem::Distance (*it, wat->O());
					OsH = MDSystem::Distance (*it, o);

					angle = acos(OsH < HOw) * 180.0 / M_PI;
					distance = OsH.Magnitude();

					histo(angle, distance);

					delete wat;
				}

				o = this->so2s.O2();
				bonded_atoms = graph.BondedAtoms (o, bondgraph::hbond);
				for (Atom_it it = bonded_atoms.begin(); it != bonded_atoms.end(); it++) {
					wat = new Water((*it)->ParentMolecule());
					wat->SetAtoms();

					// for each of the H-bonded waters, find the Os-Hw vector, and the Hw-Ow vector, and the angle between them
					HOw = MDSystem::Distance (*it, wat->O());
					OsH = MDSystem::Distance (*it, o);

					angle = acos(OsH < HOw) * 180.0 / M_PI;
					distance = OsH.Magnitude();

					histo(angle, distance);

					delete wat;
				}
			}


} // namespace so2 angle analysis
#endif
