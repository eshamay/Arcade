#ifndef H2O_ANALYSIS_H
#define H2O_ANALYSIS_H

#include "analysis.h"
#include "manipulators.h"
#include "histogram-analysis.h"
#include "molecule-analysis.h"

namespace h2o_analysis {

	using namespace md_analysis;

	class WaterThetaPhiAnalysis : public molecule_analysis::H2OAnalysis {
		protected:
			Multi2DHistogramAgent	histos;
			VecR axis, v1, v2, v3;
			double theta, phi;

		public:
			WaterThetaPhiAnalysis (system_t * t) :
				molecule_analysis::H2OAnalysis (t,
						std::string("Water theta-phi 2d angle analysis"),
						std::string ("temp")),
				axis(VecR::UnitY()),
				histos (
						-8.0, 10.0, 2.0,
						10.0,170.0,2.5,
						0.0,90.0,1.0,
						std::string("./data/h2o-theta-phi."),
						std::string(".dat")) { }

			void MoleculeCalculation ();

			void DataOutput () {
				DivideByLeftSineDegrees func;
				histos.DataOutput(func);
			}
	};


	class WaterDipoleZComponentAnalysis : public molecule_analysis::H2OAnalysis {
		public:
			WaterDipoleZComponentAnalysis (Analyzer * t) :
				molecule_analysis::H2OAnalysis (t,
						std::string("Water dipole z-component analysis"),
						std::string ("WaterOrientation-z.dat")),
				min (-10.0), max(7.0), res(0.05),
				histo (int((max-min)/res), 0.0),
				counts (int((max-min)/res), 0.0),
				axis(VecR::UnitY()) { }

			void MoleculeCalculation ();

			void DataOutput ();

		protected:
			double min, max, res;
			std::vector<double> histo, counts;
			VecR axis, v1, bisector;
			double cos_sq, dep, avg;
	};


	class DistanceAngleAnalysis : public molecule_analysis::H2OAnalysis {

		protected:
			Histogram2DAgent	histo;
			double angle;
			VecR axis, v1;

		public:

			DistanceAngleAnalysis (Analyzer * t) :
				molecule_analysis::H2OAnalysis (t, 
						std::string ("H2O - Distance-angle analysis"),
						std::string ("")),
				histo(std::string ("WaterOrientation.dat"),
						-12.0, 5.0, 0.1,
						5.0, 175.0, 2.5),
				axis(VecR::UnitY()) { }

			void MoleculeCalculation ();

			void DataOutput () {
				DivideByRightSineDegrees func;
				histo.OutputData(func);
			}
	};



}	// namespace h2o analysis



#endif
