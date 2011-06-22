#ifndef __RDF_ANALYSIS_H
#define __RDF_ANALYSIS_H

#include "../analysis.h"

namespace md_analysis {

		class RDFAnalyzer : public AnalysisSet {
			public:
			typedef Analyzer system_t;

				RDFAnalyzer (system_t * t) :
					AnalysisSet (t,
						std::string("RDF Analysis"),
						std::string("rdf.wat-H.wat-H.dat")),
					histo(0.5, 15.0, 0.05) { }

				void Analysis ();
				void DataOutput ();

			protected:

				histogram_utilities::Histogram1D<double> histo;
				
		};


} // namespace md_analysis
#endif
