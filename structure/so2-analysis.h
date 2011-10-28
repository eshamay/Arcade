#ifndef SO2_ANALYSIS_H_
#define SO2_ANALYSIS_H_

#include "rdf-analysis.h"
#include "molecule-analysis.h"

namespace so2_analysis {

	using namespace md_analysis;

	class SO2Test : public molecule_analysis::SO2XYZAnalysis {
		public:
			SO2Test (system_t * t) :
				SO2XYZAnalysis (t,
						std::string("SO2 test run"),
						std::string ("temp")) { }

			void MoleculeCalculation () { }
	};

	class SO2RDFAnalysis : public molecule_analysis::SO2Analysis {
		protected:
			RDFAgent	rdf;
			double distance;

		public:

			SO2RDFAnalysis (system_t * t) :
				SO2Analysis (t,
						std::string("SO2 - H2O RDF"),
						std::string ("temp")),
				//rdf (std::string ("rdf.so2-O.h2o-H.dat"), 0.5, 8.0, 0.05) { }
				rdf (std::string ("rdf.so2-S.h2o-O.dat"), 0.5, 8.0, 0.05) { }

			void MoleculeCalculation ();

			void DataOutput () { rdf.OutputData (); }

	};

}	// namespace so2 analysis

#endif
