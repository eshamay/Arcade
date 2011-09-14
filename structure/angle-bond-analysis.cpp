#include "angle-bond-analysis.h"

namespace md_analysis {


	H2OAngleBondAnalysis::H2OAngleBondAnalysis (system_t * t) :
		AnalysisSet (t,
				std::string ("H2O H-O-H angle and O-H bondlength histograms"),
				std::string ("")),
		h2os(t),
		// outputting to two separate files.
		bonds("bondlength.h2o.O-H.dat", 
				WaterSystem::posmin, WaterSystem::posmax, Analyzer::posres, 
				WaterSystem::SystemParameterLookup("analysis.angle-bond-histogram.bondlength-min"), 
				WaterSystem::SystemParameterLookup("analysis.angle-bond-histogram.bondlength-max"), 
				WaterSystem::SystemParameterLookup("analysis.angle-bond-histogram.bondlength-res")),

		angles("angle.h2o.H-O-H.dat", 
				WaterSystem::posmin, WaterSystem::posmax, Analyzer::posres, 
				WaterSystem::SystemParameterLookup("analysis.angle-bond-histogram.angle-min"), 
				WaterSystem::SystemParameterLookup("analysis.angle-bond-histogram.angle-max"), 
				WaterSystem::SystemParameterLookup("analysis.angle-bond-histogram.angle-res")) { }


		void H2OAngleBondAnalysis::Analysis () {
			h2os.Reload();
			h2os.FindWaterSurfaceLocation ();

			for (Wat_it it = h2os.begin(); it != h2os.end(); it++) {
				double distance = system_t::Position((*it)->ReferencePoint()) - h2os.SurfaceLocation();

				// calculate and bin the water OH bond lengths
				bonds (distance, (*it)->OH1().norm());
				bonds (distance, (*it)->OH2().norm());

				// calculate the angle of the water molecule and bin it
				double angle = (*it)->Angle();
				angles (distance, acos(angle)*180.0/M_PI);
			}

			return;
		}

}
