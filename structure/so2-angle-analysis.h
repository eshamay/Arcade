#ifndef SO2_ANGLE_ANALYSIS_H
#define SO2_ANGLE_ANALYSIS_H

#include "histogram-analysis.h"
#include "h2o-analysis.h"
#include "so2-system-analysis.h"
#include "molecule-analysis.h"

namespace so2_angle_analysis {

	using namespace md_analysis;

	// ************* so2 angle analysis *****************
	class SO2ThetaPhiAnalyzer : public molecule_analysis::SO2Analysis {
		protected:
			double theta, phi;
			Multi2DHistogramAgent	histos;
			VecR axis, v1, v2, v3;

		public:

			SO2ThetaPhiAnalyzer (system_t * t) :
				SO2Analysis (t,
						std::string("SO2 theta-phi 2d angle analysis - slices by depth"),
						std::string ("temp")),
				axis(VecR::UnitY()),
				histos (
						-10.0, 8.0, 2.0,
						//-5.0, 4.0, 1.0,
						5.0,175.0,2.5,
						0.0,90.0,1.0,
						std::string("./theta-phi/so2-theta-phi."),
						std::string(".dat")) { }

			void MoleculeCalculation ();

			void DataOutput () {
				DivideByLeftSineDegrees func;
				histos.DataOutput(func);
			}
	};


	class SO2ThetaAnalyzer : public molecule_analysis::SO2Analysis {
		protected:
			double theta;
			Histogram2DAgent	histo;
			VecR axis, v1;

		public:

			SO2ThetaAnalyzer (system_t * t) :
				SO2Analysis (t,
						std::string("SO2 Theta analysis"),
						std::string ("temp")),
				axis(VecR::UnitY()),
				histo (std::string ("so2-theta.bulk.dat"),
						-8.0, 5.0, 0.1,
						5.0,175.0,2.5) { }

			void MoleculeCalculation ();

			void DataOutput () {
				DivideByRightSineDegrees func;
				histo.OutputData(func);
			}
	};



	class SO2AngleAnalyzer : public AnalysisSet {
		protected:
			h2o_analysis::H2OSystemManipulator	h2os;
			so2_analysis::XYZSO2Manipulator		so2s;
			VecR ref_ax;

		public:

			typedef Analyzer system_t;

			SO2AngleAnalyzer (system_t * t, 
					std::string description,
					std::string fn = std::string("")) :
				AnalysisSet (t, description, fn),
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
			void SO2AngleCOMDistanceBinner ();
			virtual void UpdateHistograms (double * values) { }
	};	// so2 angle analysis







	class SO2Angles2D : public SO2AngleAnalyzer {
		protected:
			Histogram2DAgent	histo;

		public:

			typedef Analyzer system_t;

			SO2Angles2D (system_t * t) : 
				SO2AngleAnalyzer (t, std::string ("2D so2 angle analysis")),
				histo (std::string ("so2-angle-distance.phi.2d.dat"), 4.0, 16.0, 0.5, 0.0, 90.0, 1.0) { }

			virtual ~SO2Angles2D () { } 

			virtual void Analysis () { this->SetupAllAndBin(); }
			virtual void DataOutput () { histo.OutputData(); }
			virtual void BinAngles () { this->SO2AngleCOMDistanceBinner(); }
			void UpdateHistograms (double * values) { histo(values[0], values[1]); fflush(this->output); }

	};	// so2 2D angle analysis



	class SO2Angles1D : public SO2AngleAnalyzer {
		protected:
			Histogram1DAgent		theta, phi;

		public:

			typedef Analyzer system_t;

			SO2Angles1D (system_t * t) : 
				SO2AngleAnalyzer (t, std::string ("1D so2 angle analysis")),
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
	class SO2_H2O_Angles2D : public SO2AngleAnalyzer {
		protected:
			Histogram2DAgent	histo;

		public:

			typedef Analyzer system_t;

			SO2_H2O_Angles2D (system_t * t) : 
				SO2AngleAnalyzer (t, std::string ("2D so2-h2o hbond angle analysis")),
				histo (std::string ("so2-angles.Os-Hw-Ow+Os-Hw.2d.dat"), 0.0, 180.0, 1.0, 0.0, 3.0, 0.1) { }

			virtual ~SO2_H2O_Angles2D () { } 

			virtual void Analysis () { this->SetupAllAndBin(); }
			virtual void DataOutput () { histo.OutputData(); }
			virtual void BinAngles ();
			void UpdateHistograms (double * values) { return; }

	};	// so2 2D angle analysis



} // namespace so2 angle analysis
#endif
