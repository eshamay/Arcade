#ifndef MALONIC_ANALYSIS_H_
#define MALONIC_ANALYSIS_H_

#include "analysis.h"

namespace malonic_analysis {

	using namespace md_analysis;

	class MalonicTest : public AnalysisSet {

		public:
			MalonicTest (Analyzer * t) :
				AnalysisSet (t, std::string ("Malonic Tester"), std::string("")) { }

			virtual void Analysis ();

	}; // malonic test class


} // namespace malonic analysis



#endif

