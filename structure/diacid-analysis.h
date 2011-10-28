#ifndef DIACID_ANALYSIS_H_
#define DIACID_ANALYSIS_H_

#include "molecule-analysis.h"
#include "angle-analysis.h"

namespace diacid {

	using namespace md_system;

	class Test : public molecule_analysis::DiacidAnalysis {
		public:
			Test (system_t * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid test"),
						std::string ("temp")) { }

			void MoleculeCalculation ();
	};


	class CarbonylThetaPhiAnalysis : public molecule_analysis::DiacidAnalysis {
		protected:
			angle_analysis::ThetaPhiAgent angles;

		public:
			CarbonylThetaPhiAnalysis (system_t * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid carbonyl theta-phi angle depth-slice analysis"),
						std::string ("temp")),
				angles (
						std::string ("./theta-phi/theta-phi."), // filename prefix
						std::string (".dat"),	// suffix
						-14.0, 4.0, 2.0, // depths
						5.0, 175.0, 2.5, // theta range
						0.0, 180.0, 2.0) // phi range
		{ }

			void MoleculeCalculation ();
			void DataOutput () { angles.DataOutput(); }
	};

	class MethylThetaPhiAnalysis : public molecule_analysis::DiacidAnalysis {
		protected:
			angle_analysis::ThetaPhiAgent angles;

		public:
			MethylThetaPhiAnalysis (system_t * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid methyl theta-phi angle depth-slice analysis"),
						std::string ("temp")),
				angles (
						std::string ("./methyl-theta-phi/theta-phi."), // filename prefix
						std::string (".dat"),	// suffix
						-14.0, 4.0, 2.0, // depths
						5.0, 175.0, 2.5, // theta range
						0.0, 90.0, 1.0) // phi range
		{ }

			void MoleculeCalculation ();
			void DataOutput () { angles.DataOutput(); }
	};


} // namespace 

#endif
