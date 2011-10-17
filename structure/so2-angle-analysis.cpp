#include "so2-angle-analysis.h"

namespace so2_angle_analysis {


	void SO2ThetaPhiAnalyzer::MoleculeCalculation () {

		// find the center of mass location of the succinic acid
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);
		this->mol->SetOrderAxes();

		v1 = axis;// the reference axis - perp to the surface
		if (!(this->position.first))
			v1 = -v1;

		v2 = this->mol->Z();
		v3 = this->mol->SO1();

		theta = acos(v2 < v1) * 180.0 / M_PI;
		//phi = fabs(acos(this->so2->Y() < axis)) * 180.0 / M_PI;
		phi = Dihedral::Angle(v1,v2,v3) * 180.0 / M_PI;
		phi = fabs(phi);
		if (phi > 90.0)
			phi = 180.0 - phi;

		histos (this->position.second, theta, fabs(phi));
	}


	void SO2ThetaAnalyzer::MoleculeCalculation () {

		// find the center of mass location of the succinic acid
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);
		this->mol->SetOrderAxes();

		v1 = axis;// the reference axis - perp to the surface
		if (!(this->position.first))
			v1 = -v1;

		theta = acos(this->mol->Z() < v1) * 180.0 / M_PI;

		histo (this->position.second, theta);
	}






	void SO2AngleAnalyzer::SO2SetupAndBin () {

		this->so2s.UpdateSO2();
		this->so2s.SO2()->SetOrderAxes();

		// now get the vector pointing from the center of mass to the so2
		this->BinAngles();
	}

	void SO2AngleAnalyzer::SetupAllAndBin () {
		this->h2os.Reload();
		this->SO2SetupAndBin();
	}

	void SO2AngleAnalyzer::SO2AngleBinner () {

		SulfurDioxide * so2 = this->so2s.SO2();

		double values[2];

		// find the angle of the so2 theta and phi wrt the reference axis
		double angle = so2->Z() < ref_ax;
		angle = 180. * acos(angle) / M_PI;
		values[0] = angle;

		// same for the phi angle
		double angle2 = fabs(so2->Y() < ref_ax);
		angle2 = 180. * acos(angle2) / M_PI;
		if (angle2 > 90.0)
			angle2 = 180.0 - angle2;
		values[1] = angle2;

		this->UpdateHistograms(values);
		return;
	}

	void SO2AngleAnalyzer::SO2AngleCOMDistanceBinner () {

		SulfurDioxide * so2 = this->so2s.SO2();
		VecR h2o_com = h2os.CenterOfMass ();
		double distance = so2->S()->Position()[z] - h2o_com[z];

		double values[2];
		values[0] = distance;

		// find the angle of the so2 theta and phi wrt the reference axis
		double angle = so2->Z() < ref_ax;
		angle = 180. * acos(angle) / M_PI;
		//	values[1] = angle;

		// same for the phi angle
		double angle2 = fabs(so2->Y() < ref_ax);
		angle2 = 180. * acos(angle2) / M_PI;
		if (angle2 > 90.0)
			angle2 = 180.0 - angle2;

		values[1] = angle2;

		this->UpdateHistograms(values);
	}


	void SO2_H2O_Angles2D::BinAngles () {

		bondgraph::BondGraph graph;
		graph.UpdateGraph(this->begin(), this->end());

		// find all the waters that are H-bonded to the so2 through the oxygens
		VecR OsH, HOw;
		WaterPtr wat;
		double angle;
		double distance;

		AtomPtr o = this->so2s.O1();
		Atom_ptr_vec bonded_atoms = graph.BondedAtoms (o, bondgraph::hbond);
		for (Atom_it it = bonded_atoms.begin(); it != bonded_atoms.end(); it++) {
			wat = new Water((*it)->ParentMolecule());
			wat->SetAtoms();

			// for each of the H-bonded waters, find the Os-Hw vector, and the Hw-Ow vector, and the angle between them
			HOw = MDSystem::Distance (*it, wat->O());
			OsH = MDSystem::Distance (*it, o);

			angle = acos(OsH < HOw) * 180.0 / M_PI;
			distance = OsH.Magnitude();

			histo(angle, distance);

			delete wat;
		}

		o = this->so2s.O2();
		bonded_atoms = graph.BondedAtoms (o, bondgraph::hbond);
		for (Atom_it it = bonded_atoms.begin(); it != bonded_atoms.end(); it++) {
			wat = new Water((*it)->ParentMolecule());
			wat->SetAtoms();

			// for each of the H-bonded waters, find the Os-Hw vector, and the Hw-Ow vector, and the angle between them
			HOw = MDSystem::Distance (*it, wat->O());
			OsH = MDSystem::Distance (*it, o);

			angle = acos(OsH < HOw) * 180.0 / M_PI;
			distance = OsH.Magnitude();

			histo(angle, distance);

			delete wat;
		}
	}
} // so2 angle analysis namespace
