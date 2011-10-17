#include "succinic-analysis.h"


namespace succinic {

	void DensityDistribution::MoleculeCalculation () {
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);
		histo(this->position.second);
	}


	void SuccinicAcidCarbonChainDihedralAngleAnalysis::MoleculeCalculation () {
		// calculate the dihedral angle of the carbon chain.
		// the method sets the atoms of the carbon chain for the dihedral calculation...
		this->angle = this->mol->CalculateDihedralAngle() * 180.0/M_PI;

		// find the center of mass location of the succinic acid
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		this->histo(this->position.second, fabs(this->angle));
	}


	void SuccinicAcidCarbonylDihedralAngleAnalysis::MoleculeCalculation () {
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		this->mol->SetDihedralAtoms();
		// Aliphatic carbons are C2 and C3
		// Carbonyl carbons are C1 and C4
		// alcohol oxygens are O1 and O3
		// carbonyl oxygens are O2 and O4
		DihedralCalculation (this->mol->DihedralAtom(0), this->mol->GetAtom("O2"), this->mol->GetAtom("O1"));
		DihedralCalculation (this->mol->DihedralAtom(3), this->mol->GetAtom("O4"), this->mol->GetAtom("O3"));
	}

	// This calculates the two dihedrals of the y-axis, the aliphatic carbon, carbonyl carbon, and carbonyl oxygen.
	// The twist/dihedral value will tell us where the carbonyl oxygen is pointing, and if the O-C-O of the carbonyl group
	// is flat or perpendicular to the interface
	void SuccinicAcidCarbonylDihedralAngleAnalysis::DihedralCalculation (AtomPtr carbon, AtomPtr carbonyl, AtomPtr alcohol) {

		// calculate the dihedral of the ref-axis - aiphatic - carbonyl - oxygen set of vectors
		v1 = axis;// the reference axis - perp to the surface
		if (!(this->position.first))
			v1 = -v1;

		// the bisector of the two C-O bonds
		v2 = ((carbonyl->Position() - carbon->Position()) + (alcohol->Position() - carbon->Position())).normalized();
		// v3 points from carbonyl carbon to the carbonyl oxygen
		v3 = carbonyl->Position() - carbon->Position();

		twist = Dihedral::Angle(v1,v2,v3) * 180.0 / M_PI;

		this->histo (this->position.second,fabs(twist));
	}





	//// //// TILT TWIST ///// ///// 

	void SuccinicAcidCarbonylTiltTwistAnglesAnalysis::MoleculeCalculation () {
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		this->mol->SetDihedralAtoms();
		// Aliphatic carbons are C2 and C3
		// Carbonyl carbons are C1 and C4
		// alcohol oxygens are O1 and O3
		// carbonyl oxygens are O2 and O4
		DihedralCalculation (this->mol->DihedralAtom(0), this->mol->GetAtom("O2"), this->mol->GetAtom("O1"));
		DihedralCalculation (this->mol->DihedralAtom(3), this->mol->GetAtom("O4"), this->mol->GetAtom("O3"));
	}


	// This calculates the two dihedrals of the y-axis, the aliphatic carbon, carbonyl carbon, and carbonyl oxygen.
	// The twist/dihedral value will tell us where the carbonyl oxygen is pointing, and if the O-C-O of the carbonyl group
	// is flat or perpendicular to the interface
	void SuccinicAcidCarbonylTiltTwistAnglesAnalysis::DihedralCalculation (AtomPtr carbon, AtomPtr carbonyl, AtomPtr alcohol) {

		// calculate the dihedral of the ref-axis - aiphatic - carbonyl - oxygen set of vectors
		v1 = axis;// the reference axis - perp to the surface
		if (!(this->position.first))
			v1 = -v1;

		// the bisector of the two C-O bonds
		v2 = ((carbonyl->Position() - carbon->Position()) + (alcohol->Position() - carbon->Position())).normalized();
		// v3 points from carbonyl carbon to the carbonyl oxygen
		v3 = carbonyl->Position() - carbon->Position();

		twist = Dihedral::Angle(v1,v2,v3) * 180.0 / M_PI;
		tilt = acos(v2 < v1) * 180.0/M_PI;

		histos(this->position.second, tilt, fabs(twist));

		return;
	}







	void SuccinicAcidBondAngleAnalysis::MoleculeCalculation () {
		// Set the atoms in the succinic acid molecule
		this->mol->SetDihedralAtoms();
		// get the molecule's position, and the surface it's closest to
		com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];

		// check the alphatic C-C bond orientation
		// Aliphatic carbons are C2 and C3
		// Carbonyl carbons are C1 and C4
		//AngleDistanceCalculation (succ->DihedralAtom(1), succ->DihedralAtom(0));
		//AngleDistanceCalculation (succ->DihedralAtom(2), succ->DihedralAtom(3));

		// grab the aliphatic carbon and then the alcohol oxygen
		AngleDistanceCalculation (this->mol->DihedralAtom(0), this->mol->GetAtom("O1"));
		AngleDistanceCalculation (this->mol->DihedralAtom(3), this->mol->GetAtom("O3"));
		return;
	}



	void SuccinicAcidBondAngleAnalysis::AngleDistanceCalculation (AtomPtr aliphatic, AtomPtr carbonyl) {
		this->position = this->h2os.TopOrBottom(com);

		// the vector points from C2 to C1 (from the aliphatic towards the carbonyl)
		bond = MDSystem::Distance (aliphatic, carbonyl);

		angle = bond < axis;
		// if closer to the top surface, then no problem. Otherwise, invert the angle value
		if (!position.first)
			angle = -angle;
		angle = acos(angle) * 180.0 / M_PI;

		histo(this->position.second, angle);

		return;
	}



	void NeighboringWaterOrientation::MoleculeCalculation () {
		// first find the depth of the acid
		//this->com = succ->UpdateCenterOfMass()[WaterSystem::axis];

		this->mol->SetDihedralAtoms();

		// then we cycle through each water in the system and find out it's distance to the acid
		for (Wat_it wat = this->h2os.begin(); wat != this->h2os.end(); wat++) {

			AngleDepthCalculation(this->mol->GetAtom("O1"), *wat); // alcohol oxygens
			AngleDepthCalculation(this->mol->GetAtom("O3"), *wat); 

			//AngleDepthCalculation(succ->GetAtom("O2"), *wat); // carbonyl oxygens
			//AngleDepthCalculation(succ->GetAtom("O4"), *wat); 
		}
	}

	void NeighboringWaterOrientation::AngleDepthCalculation (AtomPtr oxygen, WaterPtr wat) {
			// for the given carboxylic acid oxygen, calculate the distance to the water
			oo_distance = MDSystem::Distance (oxygen, wat->O()).Magnitude();
			
			// if the distance is close enough
			VecR ax = axis;
			if (oo_distance < 10.0) {

				// we find the depth of the acid oxygen
				this->position = this->h2os.TopOrBottom(oxygen->Position()[WaterSystem::axis]);
				depth = this->position.second;

				// and also calculate its tilt wrt the water surface
				if (!this->position.first)
					ax = -ax;

				tilt = wat->Bisector() < axis;
				tilt = acos(tilt) * 180.0/M_PI;

				histos(depth, oo_distance, tilt);
			}
	}

	void CarbonylGroupDistance::MoleculeCalculation () {
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		this->mol->SetDihedralAtoms();

		//distance = MDSystem::Distance(succ->DihedralAtom(0), succ->DihedralAtom(3)).Magnitude();
		distance = this->mol->DihedralAtom(0)->Position()[WaterSystem::axis] - this->mol->DihedralAtom(3)->Position()[WaterSystem::axis];
		distance = fabs(distance);
			
		histo (this->position.second, distance);
	}


	void MethyleneBisectorTilt::MoleculeCalculation () {
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		ax = VecR::UnitY();
		if (!this->position.first)
			ax = -ax;

		this->mol->SetMethyleneBisectors();

		//angle = acos(succ->CH2_1() < succ->CH2_2())*180.0/M_PI;
		//histo (this->position.second, angle);
		histo (this->position.second, acos(this->mol->CH2_1() < ax)*180.0/M_PI);
		histo (this->position.second, acos(this->mol->CH2_2() < ax)*180.0/M_PI);
		return;
	}



} // namespace succinic
