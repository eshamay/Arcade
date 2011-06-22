#include "h2o-analysis.h"

namespace h2o_analysis {

	std::pair<double,double> OHAngleCalculator::operator() (const WaterPtr& wat) {
		double angle1 = wat->OH1() < axis;
		double angle2 = wat->OH2() < axis;
		std::pair<double,double> p = (fabs(angle1) > fabs(angle2)) 
			? std::make_pair(angle1,angle2) 
			: std::make_pair(angle2,angle1);
		return p;
	}

}
