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




	class DistanceAngleAnalysis : public AnalysisSet {

		protected:
			H2ODoubleSurfaceManipulator h2os;
			Histogram2DAgent	histo;
			h2o_analysis::surface_distance_t	distance;
			double angle;
			VecR axis;

		public:
			typedef Analyzer system_t;

			DistanceAngleAnalysis (
					system_t * t, std::string desc, std::string output_name,
					double min1, double max1, double res1,
					double min2, double max2, double res2) :
				AnalysisSet(t, desc, std::string("")),
				h2os(t),
				histo(output_name,
						min1,max1,res1,
						min2,max2,res2) { }

			virtual void Setup () { h2os.Reload(); }

			virtual void PreAnalysis () { }
			virtual void Analysis ();
			virtual void PostAnalysis () { }

			virtual void MainWaterCalculation (WaterPtr wat) = 0;

			virtual void DataOutput () = 0;
	};



	class BisectorAnalysis : public DistanceAngleAnalysis {

		public:
			typedef Analyzer system_t;

			BisectorAnalysis (system_t * t) :
				DistanceAngleAnalysis(t,
						std::string ("Water Bisector Angle Analysis"),
						std::string("water-bisector-angle.dat"),
						-12.0,4.0,1.0,
						5.0,175.0,2.5) { }

			void MainWaterCalculation (WaterPtr wat);

			void DataOutput () {
				DivideByRightSineDegrees func;
				histo.OutputData(func);
			}
	};

	class TwistAnalysis : public DistanceAngleAnalysis {

		public:
			typedef Analyzer system_t;

			TwistAnalysis (system_t * t) :
				DistanceAngleAnalysis(t,
						std::string ("Water Twist Angle Analysis"),
						std::string("water-twist-angle.dat"),
						-12.0,4.0,1.0,
						5.0,175.0,2.5) { }

			void MainWaterCalculation (WaterPtr wat);

			void DataOutput () {
				histo.OutputData();
			}
	};
}	// namespace h2o analysis



#endif
