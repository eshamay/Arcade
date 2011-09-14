#ifndef SO2_ANALYSIS_H_
#define SO2_ANALYSIS_H_

#include "analysis.h"
#include "manipulators.h"

namespace so2_analysis {

	using namespace md_analysis;



	class SO2PositionRecorder : public AnalysisSet {

		public:
			typedef Analyzer system_t;
			SO2PositionRecorder (system_t * t) : 
				AnalysisSet (t, 
						std::string ("Record position of so2 relative to the surface, and the surface location"),
						std::string ("so2-position.dat")),
				so2s(t), h2os(t) { }

			void Analysis ();

		private:
			SO2SystemManipulator	so2s;
			h2o_analysis::H2OSystemManipulator	h2os;
	};


	class SO2BondLengthAnalyzer : public AnalysisSet {

		public:
			typedef Analyzer system_t;
			SO2BondLengthAnalyzer (system_t * t) :
				AnalysisSet (t,
						std::string ("so2 bondlength analysis"),
						std::string ("so2-bondlengths.dat")),
				so2s(t) { so2s.Initialize(); }

			void Analysis ();

		protected:
			XYZSO2Manipulator	so2s;
	};


	class SO2AngleAnalyzer : public AnalysisSet {

		public:
			typedef Analyzer system_t;
			SO2AngleAnalyzer (system_t * t) :
				AnalysisSet (t,
						std::string ("so2 angle analysis"),
						std::string ("so2-angles.dat")),
				so2s(t) { so2s.Initialize(); }

			void Analysis ();

		protected:
			XYZSO2Manipulator	so2s;
	};



	class ClosestWaterBondlengths : public AnalysisSet {
		public:
			typedef Analyzer system_t;
			ClosestWaterBondlengths (system_t * t) :
				AnalysisSet (t,
						std::string ("so2's bound water oh bondlengths"),
						std::string ("top-water-bondlengths.dat")),
				so2s(t) { so2s.Initialize(); }

			void Analysis ();

		protected:
			XYZSO2Manipulator	so2s;
	};	// closest water bondlengths 



	class WaterAngleAnalyzer : public AnalysisSet {
		public:
			typedef Analyzer system_t;
			WaterAngleAnalyzer (system_t * t) : 
				AnalysisSet (t, 
						std::string ("Water's HOH Angle"),
						std::string ("h2o-angle.dat")),
				so2s(t) { so2s.Initialize(); }

			void Analysis ();

		private:
			XYZSO2Manipulator	so2s;
	};


	/*
		 class so2dipoleanalyzer : public analysisset<t> {

		 public:
		 typedef analyzer<t> system_t;
		 so2dipoleanalyzer (system_t * t) :
		 analysisset<t> (t,
		 std::string ("so2 dipole analysis"),
		 std::string ("so2-dipole.dat")),
		 so2s(t) { }

		 void analysis ();
		 void UpdateDipole (AtomPtr atom, const int num, VecR& dipole, const VecR& ref);

		 protected:
		 XYZSO2Manipulator	so2s;
		 };

		 void SO2DipoleAnalyzer::Analysis () {

		 so2s.UpdateSO2();
		 so2s.SO2()->UpdateCenterOfMass();
		 VecR com = so2s.SO2()->CenterOfMass();

		 VecR dipole;
		 UpdateDipole (so2s.S(), 1, dipole, com);
		 UpdateDipole (so2s.O1(), 4, dipole, com);
		 UpdateDipole (so2s.O2(), 4, dipole, com);

		 fprintf (this->output, "% 12.8f % 12.8f % 12.8f\n", dipole[x], dipole[y], dipole[z]);
		 }

		 void SO2DipoleAnalyzer::UpdateDipole (AtomPtr atom, const int num, VecR& dipole, const VecR& ref) {
		 vector_map_vec wans = so2s.GetWanniers (atom,num);
		 for (vector_map_it it = wans.begin(); it != wans.end(); it++) {
		 VecR r (MDSystem::Distance(ref, (*it)));
		 r *= (-2.0);
		 dipole += r;
		 }
		 dipole += MDSystem::Distance(ref, atom->Position()) * atom->Charge();
		 }
		 */

}	// namespace so2_analysis


#endif
