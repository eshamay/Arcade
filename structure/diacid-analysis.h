#ifndef DIACID_ANALYSIS_H_
#define DIACID_ANALYSIS_H_

#include "molecule-analysis.h"
#include "angle-analysis.h"
#include "rdf-analysis.h"
#include "dimergraph.h"

namespace diacid {

	using namespace md_system;
	using namespace md_analysis;

	class Test : public molecule_analysis::DiacidAnalysis {
		public:
			Test (system_t * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid test"),
						std::string ("temp")) { }

			void MoleculeCalculation ();
	};


	class CarboxylicThetaPhiAnalysis : public molecule_analysis::DiacidAnalysis {
		protected:
			angle_analysis::ThetaPhiAgent angles;

		public:
			CarboxylicThetaPhiAnalysis (system_t * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid carbonyl theta-phi angle depth-slice analysis"),
						std::string ("temp")),
				angles (
						std::string ("./carbonyl-theta-phi/theta-phi."), // filename prefix
						std::string (".dat"),	// suffix
						-14.0, 4.0, 2.0, // depths
						5.0, 175.0, 2.5, // theta range
						0.0, 180.0, 2.0) // phi range
		{ }

			void MoleculeCalculation ();
			void DataOutput () { angles.DataOutput(); }
	};


	class Dimers : public molecule_analysis::DiacidAnalysis {
		public:
			Dimers (system_t * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Look at the bonding between the dicarboxylic acids"),
						std::string ("DiacidDimers.dat")),
				initialized(false) { }

			void PreCalculation ();
			void MoleculeCalculation ();
			void PostCalculation ();

		protected:
			DimerGraph	_graph;
			bool initialized;
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


	class RDF : public molecule_analysis::DiacidAnalysis {
		protected:
			RDFAgent	rdf;
			double distance;
		public:
			RDF (Analyzer * t);
			void MoleculeCalculation ();
			void DataOutput () { rdf.OutputData(); }
	};




	class CarbonylThetaThetaAnalysis : public molecule_analysis::DiacidAnalysis {
		protected:
			angle_analysis::ThetaThetaAgent angles;

		public:
			CarbonylThetaThetaAnalysis (Analyzer * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("carbonyl C=O theta1-theta2 angle depth-slice analysis"),
						std::string ("temp")),
				angles (
						std::string ("./carbonyl-theta-theta/theta-theta."), // filename prefix
						std::string (".dat"),	// suffix
						-14.0, 4.0, 2.0, // depths
						5.0, 175.0, 2.5, // theta range
						5.0, 175.0, 2.5) // phi range
		{ }

			void MoleculeCalculation ();
			void DataOutput () { angles.DataOutput(); }
	};

} // namespace 

#endif
