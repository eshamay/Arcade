#include "so2-system-analysis.h"

namespace so2_analysis {

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


	void SO2BondLengthAnalyzer::Analysis () {
		so2s.UpdateSO2();

		double so1 = so2s.SO2()->SO1().Magnitude();
		double so2 = so2s.SO2()->SO2().Magnitude();

		fprintf (this->output, "% 9.7f % 9.7f\n", so1, so2);
	}


	void SO2AngleAnalyzer::Analysis () {
		so2s.UpdateSO2();

		double angle = so2s.SO2()->Angle();


		fprintf (this->output, "% 12.7f\n", angle);
	}



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

	std::pair<double,double> SOAngleCalculator::operator() (const SulfurDioxide* so2) {
		double angle1 = so2->SO1() < axis;
		double angle2 = so2->SO2() < axis;
		std::pair<double,double> p = (fabs(angle1) > fabs(angle2)) 
			? std::make_pair(angle1,angle2) 
			: std::make_pair(angle2,angle1);
		return p;
	}


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
		std::sort (wats.begin(), wats.end(), MoleculeToReferenceDistance_cmp (so2s.SO2()));

		double angle;
		for (int i = 0; i < 3; i++) {
			WaterPtr wat = wats[i];
			angle = wat->Angle();
			fprintf (this->output, "% 12.7f", angle);
		}
		fprintf (this->output, "\n");
	}


}
