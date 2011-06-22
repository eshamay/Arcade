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

	// write out the position of the so2, the position of the surface, and the difference between the two
	void SO2PositionRecorder::Analysis () {
		h2os.FindWaterSurfaceLocation();
		double so2_pos, surface, distance;
		so2_pos = system_t::Position(so2s.S());
		surface = h2os.SurfaceLocation();
		if (h2os.TopSurface())
			distance = so2_pos - surface;
		else
			distance = surface - so2_pos;

		fprintf (this->output, "% 12.7f % 12.7f % 12.7f\n", so2_pos, surface, distance);
	}


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

	void SO2BondLengthAnalyzer::Analysis () {
		so2s.UpdateSO2();

		double so1 = so2s.SO2()->SO1().Magnitude();
		double so2 = so2s.SO2()->SO2().Magnitude();

		fprintf (this->output, "% 9.7f % 9.7f\n", so1, so2);
	}

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

	void SO2AngleAnalyzer::Analysis () {
		so2s.UpdateSO2();

		double angle = so2s.SO2()->Angle();


		fprintf (this->output, "% 12.7f\n", angle);
	}


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

	void ClosestWaterBondlengths::Analysis () { 
		so2s.UpdateSO2();
		this->LoadAll();

		Water_ptr_vec wats;
		for (Mol_it it = this->begin_mols(); it != this->end_mols(); it++) {
			if ((*it)->MolType() == Molecule::H2O) {
				WaterPtr wat = static_cast<WaterPtr>(*it);
				wat->SetAtoms();
				wats.push_back(wat);
			}
		}
		//std::sort (wats.begin(), wats.end(), WaterToSO2Distance_cmp (so2s.SO2()));
		//std::sort (wats.begin(), wats.end(), WaterToSO2Distance_cmp (so2s.SO2()));
		// sort the waters by position along the reference axis - first waters are lowest, last are highest
		std::sort (wats.begin(), wats.end(), typename system_t::molecule_position_pred(Atom::O));

		VecR oh1, oh2;
		double tally = 0.0;
		for (int i = 0; i < 10; i++) {
			WaterPtr wat = wats[i];
			// calculate the component of the oh bonds in the Z-direction
			oh1 = wat->OH1();
			oh2 = wat->OH2();

			tally += oh1[z];
			tally += oh2[z];

			//fprintf (this->output, "% 12.7f % 12.7f", oh1, oh2);
		}
		fprintf (this->output, "% 12.7f\n", tally);
	}




	// functor takes a water and returns the values of the cos(angles) formed between the two oh-vectors. The first value of the pair is always the greater (magnitude) of the two values.
	class SOAngleCalculator : public std::unary_function <SulfurDioxide*,std::pair<double,double> > {
		private:
			VecR axis;	// the reference axis to which the angles will be formed
		public:
			SOAngleCalculator (const VecR ax) : axis(ax) { }
			std::pair<double,double> operator() (const SulfurDioxide* so2) {
				double angle1 = so2->SO1() < axis;
				double angle2 = so2->SO2() < axis;
				std::pair<double,double> p = (fabs(angle1) > fabs(angle2)) 
					? std::make_pair(angle1,angle2) 
					: std::make_pair(angle2,angle1);
				return p;
			}
	};

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

	void WaterAngleAnalyzer::Analysis () {
		so2s.UpdateSO2();
		this->LoadAll();

		Water_ptr_vec wats;
		for (Mol_it it = this->begin_mols(); it != this->end_mols(); it++) {
			if ((*it)->MolType() == Molecule::H2O) {
				WaterPtr wat = static_cast<WaterPtr>(*it);
				wat->SetAtoms();
				wats.push_back(wat);
			}
		}
		std::sort (wats.begin(), wats.end(), WaterToSO2Distance_cmp (so2s.SO2()));

		double angle;
		for (int i = 0; i < 3; i++) {
			WaterPtr wat = wats[i];
			angle = wat->Angle();
			fprintf (this->output, "% 12.7f", angle);
		}
		fprintf (this->output, "\n");
	}

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
