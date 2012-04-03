#include "diacid-analysis.h"

namespace diacid {

	BondLengths::BondLengths (Analyzer * t) :
		DiacidAnalysis (t,
				std::string ("Intramolecular Bondlengths"),
				std::string("IntraMolecularBondLengths.dat")),
			min(0.0), max(7.0), res(0.005) {
			lengths[c1o1] = new histogram_utilities::Histogram1D<double> (min, max, res);
			lengths[c1oh1] = new histogram_utilities::Histogram1D<double> (min, max, res);
			lengths[c2o2] = new histogram_utilities::Histogram1D<double> (min, max, res);
			lengths[c2oh2] = new histogram_utilities::Histogram1D<double> (min, max, res);
			lengths[h1o2] = new histogram_utilities::Histogram1D<double> (min, max, res);
			lengths[h2o1] = new histogram_utilities::Histogram1D<double> (min, max, res);
			lengths[h1oh1] = new histogram_utilities::Histogram1D<double> (min, max, res);
			lengths[h2oh2] = new histogram_utilities::Histogram1D<double> (min, max, res);
			lengths[h1oh2] = new histogram_utilities::Histogram1D<double> (min, max, res);
			lengths[h2oh1] = new histogram_utilities::Histogram1D<double> (min, max, res);
			//lengths[o1waterh] = new histogram_utilities::Histogram1D<double> (min, max, res);
			//lengths[o2waterh] = new histogram_utilities::Histogram1D<double> (min, max, res);
		}

	void BondLengths::CalcDistance (AtomPtr atom1, AtomPtr atom2, bond_t bond) {
		distance = MDSystem::Distance(atom1, atom2).norm();
		lengths[bond]->operator()(distance);
	}

	void BondLengths::MoleculeCalculation () {
		this->mol->SetAtoms();

		CalcDistance(this->mol->C1(), this->mol->O1(), c1o1);
		CalcDistance(this->mol->C2(), this->mol->O2(), c2o2);
		CalcDistance(this->mol->C1(), this->mol->OH1(), c1oh1);
		CalcDistance(this->mol->C2(), this->mol->OH2(), c2oh2);

		CalcDistance(this->mol->H1(), this->mol->OH1(), h1oh1);
		CalcDistance(this->mol->H1(), this->mol->OH2(), h1oh2);
		CalcDistance(this->mol->H1(), this->mol->O2(), h1o2);

		CalcDistance(this->mol->H2(), this->mol->OH2(), h2oh2);
		CalcDistance(this->mol->H2(), this->mol->OH1(), h2oh1);
		CalcDistance(this->mol->H2(), this->mol->O1(), h2o1);
	}

	void BondLengths::OutputDataPoint (bond_t bond, double position) {
		fprintf (this->output, " %6.4f", lengths[bond]->Population(position));
	}

	void BondLengths::DataOutput () {
		rewind (this->output);

		fprintf (this->output, "position c1o1 c2o2 c1oh1 c2oh2 h1oh1 h1oh2 h1o2 h2oh2 h2oh1 h2o1\n");
		for (double pos = min; pos < max; pos += res) {
			fprintf (this->output, "% 8.5f ", pos);
			OutputDataPoint (c1o1, pos);
			OutputDataPoint (c2o2, pos);
			OutputDataPoint (c1oh1, pos);
			OutputDataPoint (c2oh2, pos);
			OutputDataPoint (h1oh1, pos);
			OutputDataPoint (h1oh2, pos);
			OutputDataPoint (h1o2, pos);
			OutputDataPoint (h2oh2, pos);
			OutputDataPoint (h2oh1, pos);
			OutputDataPoint (h2o1, pos);
			fprintf (this->output, "\n");
		}
	}


	void Dimers::PreCalculation () {
		if (!initialized) {
			molecule_analysis::SingleMoleculeAnalysis<alkane::Diacid>::PreCalculation ();

			std::vector<alkane::Diacid *> mols;
			for (Mol_it it = this->analysis_mols.begin(); it != this->analysis_mols.end(); it++)
				mols.push_back(static_cast<alkane::Diacid *>(*it));

			this->_graph.Initialize(mols.begin(), mols.end());

			std::for_each(mols.begin(), mols.end(), std::mem_fun(&alkane::Diacid::SetAtoms));
			initialized = true;
		} 

		else {
			this->_graph.RecalculateBonds();

		}
	}

	void Dimers::MoleculeCalculation () {

		/*
			 this->mol->LoadAtomGroups();
			 this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
			 this->position = this->h2os.TopOrBottom(com);

			 alkane::Diacid * dia;
			 int n=0;
		// run the calculation between the current molecule and all the other diacids
		for (Mol_it it = this->analysis_mols.begin(); it != this->analysis_mols.end(); it++) {
		if (*it == this->mol) continue;	// don't check against itself

		dia = static_cast<alkane::Diacid *>(*it);
		//n += NumberOfHBonds(this->mol, dia);

		// calculate the number of H-bonds formed between the two diacids
		}

		fprintf (this->output, "%d\n", n);
		fflush(this->output);
		*/
	}

	void Dimers::PostCalculation () {

		DimerGraph::hbond_data_t nums = this->_graph.NumHBonds();
		fprintf (this->output, "%5d %5d %5d %5d\n", 
				nums.acids, nums.methyls, nums.carbonyls, this->_graph.InterAcidConnections());
		fflush(this->output);
	}


	RDF::RDF (Analyzer * t) :
		molecule_analysis::DiacidAnalysis (t,
				std::string ("Malonic RDFs"),
				std::string("")),
		rdf_alc (std::string("MalonicRDF.alcO-H.surface.dat"), 
				0.5, 10.0, 0.05),
		rdf_carb (std::string("MalonicRDF.carbO-H.surface.dat"), 
				0.5, 10.0, 0.05) { }

	void RDF::MoleculeCalculation () {
		this->mol->SetAtoms();
		//this->mol->LoadAtomGroups();
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		WaterPtr wat;
		AtomPtr o, oh;
		for (alkane::Diacid::atom_group_list::const_iterator coo = this->mol->carbonyls_begin(); coo != this->mol->carbonyls_end(); coo++) {
			o = coo->Right();	// carbonyl
			oh = coo->Left();	// alcohol

			for (Wat_it it = this->h2os.begin(); it != this->h2os.end(); it++) 
			{
				wat = *it;
				distance = MDSystem::Distance (oh, wat->H1()).norm();
				rdf_alc(distance);
				distance = MDSystem::Distance (oh, wat->H2()).norm();
				rdf_alc(distance);

				distance = MDSystem::Distance (o, wat->H1()).norm();
				rdf_carb(distance);
				distance = MDSystem::Distance (o, wat->H2()).norm();
				rdf_carb(distance);
			}
		}
	}



	void Test::MoleculeCalculation () {
		this->mol->SetAtoms();
		this->mol->Print();

		std::cout << "methyl Hs" << std::endl;
		Atom_ptr_vec Hs = this->mol->methyl_hydrogens();
		std::for_each(Hs.begin(), Hs.end(), std::mem_fun(&Atom::Print));

		std::cout << "carbonyl Os" << std::endl;
		Atom_ptr_vec Os = this->mol->carbonyl_oxygens();
		std::for_each(Os.begin(), Os.end(), std::mem_fun(&Atom::Print));

		std::cout << "carbonyl Hs" << std::endl;
		Hs = this->mol->carbonyl_hydrogens();
		std::for_each(Hs.begin(), Hs.end(), std::mem_fun(&Atom::Print));
	}




	void CarboxylicThetaPhiAnalysis::MoleculeCalculation () {

		// find the center of mass location of the succinic acid
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);
		this->mol->SetAtoms();

		// get the dihedral angle
		angles (this->mol->CarbonylBisector1(), this->mol->CO1(), this->position);
		angles (this->mol->CarbonylBisector2(), this->mol->CO2(), this->position);

	}

	void MethylThetaPhiAnalysis::MoleculeCalculation () {

		this->mol->SetAtoms();

		// find the center of mass location of the succinic acid
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		// run through each methyl group and grab the angles needed
		for (alkane::Diacid::atom_group_list::const_iterator it = this->mol->methyls_begin(); 
				it != this->mol->methyls_end(); it++) {
			angles (it->Bisector(), it->Bond1(), this->position);
		}
	}



	void CarbonBackboneThetaCarboxylicDihedral::MoleculeCalculation () {
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);
		this->mol->SetAtoms();

		v1 = axis;// the reference axis - perp to the surface
		if (!(position.first))
			v1 = -v1;

		ccc.SetAtoms(this->mol->GetAtom("C1"), this->mol->GetAtom("C2"), this->mol->GetAtom("C3"));

		theta = acos(ccc.Bisector() < v1) * 180.0 / M_PI;

		psi = alkane::Diacid::MalonicDihedralAngle (this->mol);

		//printf ("\ntheta = %f\nphi = %f\nposition = %f\n", theta, phi, position.second);
		// bin the dihedral angles of both the carboxylic groups against the backbone tilt
		angles.Override(this->position.second, theta, psi.first);
		angles.Override(this->position.second, theta, psi.second);
	}



	void CarbonBackboneThetaPhi::MoleculeCalculation () {
		// find the center of mass location of the succinic acid
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);
		this->mol->SetAtoms();

		// find the bisector axis
		// it's splits the two vectors C2->C1 and C2->C3
		// This is specific to malonic acid
		ccc.SetAtoms(this->mol->GetAtom("C1"), this->mol->GetAtom("C2"), this->mol->GetAtom("C3"));

		// get the angles binned
		angles (ccc.Bisector(), ccc.Bond1(), this->position);
	}

	void COTheta::MoleculeCalculation () {
		this->mol->SetAtoms();
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		v1 = axis;// the reference axis - perp to the surface
		if (!(this->position.first))
			v1 = -v1;

		double theta1 = acos(this->mol->CO1() < v1) * 180.0/M_PI;
		double theta2 = acos(this->mol->CO2() < v1) * 180.0/M_PI;

		angles (this->position.second, theta1);
		angles (this->position.second, theta2);
	}

	void CHTheta::MoleculeCalculation () {
		this->mol->SetAtoms();
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		v1 = axis;// the reference axis - perp to the surface
		if (!(this->position.first))
			v1 = -v1;

		VecR ch1, ch2;
		// run through each methyl group and grab the angles needed
		for (alkane::Diacid::atom_group_list::const_iterator it = this->mol->methyls_begin(); 
				it != this->mol->methyls_end(); it++) {

			ch1 = it->Bond1();
			ch2 = it->Bond2();

			double theta1 = acos(ch1 < v1) * 180.0/M_PI;
			double theta2 = acos(ch2 < v1) * 180.0/M_PI;

			angles (this->position.second, theta1);
			angles (this->position.second, theta2);
		}

	}

	void CarboxylicDihedralPsiPsi::MoleculeCalculation () {
		this->mol->SetAtoms();
		//this->mol->LoadAtomGroups();
		this->com = this->mol->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		std::pair<double,double> psi = alkane::Diacid::MalonicDihedralAngle (this->mol);

		angles.Override(this->position.second, fabs(psi.first), fabs(psi.second));

		return;
	}


} // namespace diacid
