#include "succinic-analysis.h"


namespace succinic {



	void DensityDistribution::SuccinicAcidCalculation (alkane::SuccinicAcid * succ) {
		this->com = succ->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);
		histo(this->position.second);
	}


	void SuccinicAcidCarbonChainDihedralAngleAnalysis::SuccinicAcidCalculation (alkane::SuccinicAcid * succ) {
		// calculate the dihedral angle of the carbon chain.
		// the method sets the atoms of the carbon chain for the dihedral calculation...
		this->angle = succ->CalculateDihedralAngle() * 180.0/M_PI;

		// find the center of mass location of the succinic acid
		this->com = succ->UpdateCenterOfMass() [WaterSystem::axis];
		this->distance = this->h2os.TopOrBottom(com);

		this->histo(this->distance.second, fabs(this->angle));
	}


	void SuccinicAcidCarbonylDihedralAngleAnalysis::SuccinicAcidCalculation (alkane::SuccinicAcid * succ) {
		this->com = succ->UpdateCenterOfMass() [WaterSystem::axis];
		this->distance = this->h2os.TopOrBottom(com);

		succ->SetDihedralAtoms();
		// Aliphatic carbons are C2 and C3
		// Carbonyl carbons are C1 and C4
		// alcohol oxygens are O1 and O3
		// carbonyl oxygens are O2 and O4
		DihedralCalculation (succ->DihedralAtom(0), succ->GetAtom("O2"), succ->GetAtom("O1"));
		DihedralCalculation (succ->DihedralAtom(3), succ->GetAtom("O4"), succ->GetAtom("O3"));
	}

	// This calculates the two dihedrals of the y-axis, the aliphatic carbon, carbonyl carbon, and carbonyl oxygen.
	// The twist/dihedral value will tell us where the carbonyl oxygen is pointing, and if the O-C-O of the carbonyl group
	// is flat or perpendicular to the interface
	void SuccinicAcidCarbonylDihedralAngleAnalysis::DihedralCalculation (AtomPtr carbon, AtomPtr carbonyl, AtomPtr alcohol) {

		// calculate the dihedral of the ref-axis - aiphatic - carbonyl - oxygen set of vectors
		v1 = axis;// the reference axis - perp to the surface
		if (!(this->distance.first))
			v1 = -v1;

		// the bisector of the two C-O bonds
		v2 = ((carbonyl->Position() - carbon->Position()) + (alcohol->Position() - carbon->Position())).normalized();
		// v3 points from carbonyl carbon to the carbonyl oxygen
		v3 = carbonyl->Position() - carbon->Position();

		twist = Dihedral::Angle(v1,v2,v3) * 180.0 / M_PI;

		this->histo (distance.second,fabs(twist));
	}





	//// //// TILT TWIST ///// ///// 

	void SuccinicAcidCarbonylTiltTwistAnglesAnalysis::SuccinicAcidCalculation (alkane::SuccinicAcid * succ) {
		this->com = succ->UpdateCenterOfMass() [WaterSystem::axis];
		this->distance = this->h2os.TopOrBottom(com);

		succ->SetDihedralAtoms();
		// Aliphatic carbons are C2 and C3
		// Carbonyl carbons are C1 and C4
		// alcohol oxygens are O1 and O3
		// carbonyl oxygens are O2 and O4
		DihedralCalculation (succ->DihedralAtom(0), succ->GetAtom("O2"), succ->GetAtom("O1"));
		DihedralCalculation (succ->DihedralAtom(3), succ->GetAtom("O4"), succ->GetAtom("O3"));
	}


	// This calculates the two dihedrals of the y-axis, the aliphatic carbon, carbonyl carbon, and carbonyl oxygen.
	// The twist/dihedral value will tell us where the carbonyl oxygen is pointing, and if the O-C-O of the carbonyl group
	// is flat or perpendicular to the interface
	void SuccinicAcidCarbonylTiltTwistAnglesAnalysis::DihedralCalculation (AtomPtr carbon, AtomPtr carbonyl, AtomPtr alcohol) {

		// calculate the dihedral of the ref-axis - aiphatic - carbonyl - oxygen set of vectors
		v1 = axis;// the reference axis - perp to the surface
		if (!(this->distance.first))
			v1 = -v1;

		// the bisector of the two C-O bonds
		v2 = ((carbonyl->Position() - carbon->Position()) + (alcohol->Position() - carbon->Position())).normalized();
		// v3 points from carbonyl carbon to the carbonyl oxygen
		v3 = carbonyl->Position() - carbon->Position();

		twist = Dihedral::Angle(v1,v2,v3) * 180.0 / M_PI;
		tilt = acos(v2 < v1) * 180.0/M_PI;

		histos(distance.second, tilt, fabs(twist));

		return;
	}







	void SuccinicAcidBondAngleAnalysis::SuccinicAcidCalculation (alkane::SuccinicAcid * succ) {
		// Set the atoms in the succinic acid molecule
		succ->SetDihedralAtoms();
		// get the molecule's position, and the surface it's closest to
		com = succ->UpdateCenterOfMass() [WaterSystem::axis];

		// check the alphatic C-C bond orientation
		// Aliphatic carbons are C2 and C3
		// Carbonyl carbons are C1 and C4
		//AngleDistanceCalculation (succ->DihedralAtom(1), succ->DihedralAtom(0));
		//AngleDistanceCalculation (succ->DihedralAtom(2), succ->DihedralAtom(3));

		// grab the aliphatic carbon and then the alcohol oxygen
		AngleDistanceCalculation (succ->DihedralAtom(0), succ->GetAtom("O1"));
		AngleDistanceCalculation (succ->DihedralAtom(3), succ->GetAtom("O3"));
		return;
	}



	void SuccinicAcidBondAngleAnalysis::AngleDistanceCalculation (AtomPtr aliphatic, AtomPtr carbonyl) {
		distance = this->h2os.TopOrBottom(com);

		// the vector points from C2 to C1 (from the aliphatic towards the carbonyl)
		bond = MDSystem::Distance (aliphatic, carbonyl);

		angle = bond < axis;
		// if closer to the top surface, then no problem. Otherwise, invert the angle value
		if (!distance.first)
			angle = -angle;
		angle = acos(angle) * 180.0 / M_PI;

		histo(distance.second, angle);

		return;
	}



	void NeighboringWaterOrientation::SuccinicAcidCalculation (alkane::SuccinicAcid * succ) {
		// first find the depth of the acid
		//this->com = succ->UpdateCenterOfMass()[WaterSystem::axis];

		succ->SetDihedralAtoms();

		// then we cycle through each water in the system and find out it's distance to the acid
		for (Wat_it wat = this->h2os.begin(); wat != this->h2os.end(); wat++) {

			AngleDepthCalculation(succ->GetAtom("O1"), *wat); // alcohol oxygens
			AngleDepthCalculation(succ->GetAtom("O3"), *wat); 

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
				this->distance = this->h2os.TopOrBottom(oxygen->Position()[WaterSystem::axis]);
				depth = distance.second;

				// and also calculate its tilt wrt the water surface
				if (!this->distance.first)
					ax = -ax;

				tilt = wat->Bisector() < axis;
				tilt = acos(tilt) * 180.0/M_PI;

				histos(depth, oo_distance, tilt);
			}
	}


} // namespace succinic
