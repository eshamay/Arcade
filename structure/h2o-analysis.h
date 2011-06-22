#ifndef H2O_ANALYSIS_H
#define H2O_ANALYSIS_H

#include "analysis.h"
#include "manipulators.h"

namespace h2o_analysis {

	using namespace md_analysis;


	// functor takes a water and returns the values of the cos(angles) formed between the two oh-vectors. The first value of the pair is always the greater (magnitude) of the two values.
	class OHAngleCalculator : public std::unary_function <WaterPtr,std::pair<double,double> > {
		private:
			VecR axis;	// the reference axis to which the angles will be formed

		public:
			OHAngleCalculator (const VecR ax) : axis(ax) { }
			std::pair<double,double> operator() (const WaterPtr& wat);
	};


}	// namespace h2o analysis



#endif
