#ifndef MOLECULE_ANALYSIS_H_
#define MOLECULE_ANALYSIS_H_

#include "analysis.h"
#include "manipulators.h"
#include "histogram-analysis.h"

namespace molecule_analysis {

	using namespace md_analysis;

	class SuccinicAcidAnalysis : public AnalysisSet {

		protected:
			h2o_analysis::H2ODoubleSurfaceManipulator	h2os;
			double com;
			h2o_analysis::surface_distance_t position;

		public:
			typedef Analyzer system_t;
			SuccinicAcidAnalysis (system_t * t, std::string desc, std::string fn) : 
				AnalysisSet (t, desc, fn),
				h2os(t) { }

			virtual void Setup () { h2os.Reload(); }
			virtual void Analysis ();

			virtual void SuccinicAcidCalculation (alkane::SuccinicAcid *) = 0;
			virtual void PreCalculation () { return; }
			virtual void PostCalculation () { return; }

	}; // succinic analysis set

} // namespace molecule analysis

#endif
