#include "angle-analysis.h"

namespace angle_analysis {

	std::pair<double,double> OHAngleCalculator::operator() (const WaterPtr& wat) {
		double angle1 = wat->OH1() < axis;
		double angle2 = wat->OH2() < axis;
		std::pair<double,double> p = (fabs(angle1) > fabs(angle2)) 
			? std::make_pair(angle1,angle2) 
			: std::make_pair(angle2,angle1);
		return p;
	}

	void H2OAngleAnalysis::BinAngles (MolPtr mol) {

		Water * wat = new Water(mol);
		wat->SetOrderAxes();
		double distance;
		if (h2os.TopSurface()) {
			distance = system_t::Position(wat->ReferencePoint()) - h2os.SurfaceLocation();
			angles.Alpha(distance, wat->Bisector() < VecR::UnitY());
		}
		else {
			distance = h2os.SurfaceLocation() - system_t::Position(wat->ReferencePoint()); 
			angles.Alpha(distance, -(wat->Bisector() < VecR::UnitY()));
			// get the projection of the molecule's x-axis onto the system's x-y plane
		}

		// get the projection of the molecule's x-axis onto the system's x-y plane
		angles.Beta(distance, fabs(wat->Y() < VecR::UnitY()));

		// switch the lab frame axes to make the Y-axis the reference one instead of Z
		//r = _rotation * wat->Z();
		//r = wat->Z();
		//_dcm = wat->DCMToLab().transpose();
		//coordinate_conversion::DCM2EulerAngles_ZXZ (&_dcm(0,0), angles);

		delete wat;

		return;
	} // bin water angles



	void H2OAngleAnalysis::Analysis () {

		h2os.Reload();
		h2os.FindWaterSurfaceLocation();

		for (Wat_it wat = h2os.begin(); wat != h2os.end(); wat++) {
			BinAngles(*wat);
		}
	}



	void OHAngleAnalysis::Analysis () {

		h2os.FindWaterSurfaceLocation();
		h2os.UpdateAnalysisWaters();
		OHAngleCalculator::angle_pair_t p;

		double oh1, oh2;

		for (Wat_it wat = h2os.begin(); wat != h2os.end(); wat++) {

			if (h2os.TopSurface()) {
				distance = system_t::Position(*wat) - h2os.SurfaceLocation();
			} else {
				distance = h2os.SurfaceLocation() - system_t::Position(*wat);
			}

			(*wat)->SetAtoms();
			oh1 = (*wat)->OH1() < VecR::UnitY();
			oh2 = (*wat)->OH2() < VecR::UnitY();


			if (h2os.TopSurface()) {
				_alpha (distance, oh1);
				_alpha (distance, oh2);
			}
			else {
				_alpha (distance, -oh1);
				_alpha (distance, -oh2);
			}

		}
	}




	void WaterOHAngleAnalysis::Analysis () {

		h2os.FindWaterSurfaceLocation();
		std::pair<double,double> p;

		for (Wat_it wat = h2os.begin(); wat != h2os.end(); wat++) {
			distance = system_t::Position(*wat) - h2os.SurfaceLocation();

			// for each water, find the angle that the oh bonds make with the surface normal.
			p = oh_calculator(*wat);

			_alpha (distance, p.first);
			//_alpha (distance, p.second);
		}
	}


	std::pair<double,double> SOAngleCalculator::operator() (const SulfurDioxide* so2) {
		double angle1 = so2->SO1() < axis;
		double angle2 = so2->SO2() < axis;
		std::pair<double,double> p = (fabs(angle1) > fabs(angle2)) 
			? std::make_pair(angle1,angle2) 
			: std::make_pair(angle2,angle1);
		return p;
	}

	void ReferenceSO2AngleAnalysis::BinAngles (SulfurDioxide * so2) {
		so2->SetOrderAxes();

		// calculate the position of the so2 relative to the water surface
		double distance;
		if (h2os.TopSurface()) {
			distance = system_t::Position(so2->ReferencePoint()) - h2os.SurfaceLocation();
			// get the value of theta: molecular bisector angle with system reference axis
			angles.Alpha(distance, so2->Bisector() < VecR::UnitY());
		} else {
			distance = h2os.SurfaceLocation() - system_t::Position(so2->ReferencePoint());		// bottom surface
			angles.Alpha(distance, -(so2->Bisector() < VecR::UnitY()));	// bottom surface
		}

		// get the value of phi - the molecular normal angle with the system ref
		angles.Beta(distance, fabs(so2->Y() < VecR::UnitY()));

		return;
	}

	void ReferenceSO2AngleAnalysis::Analysis () {
		h2os.FindWaterSurfaceLocation();
		BinAngles(so2s.SO2());
	}


	void SOAngleAnalysis::Analysis () {

		h2os.FindWaterSurfaceLocation();
		std::pair<double,double> p;

		if (h2os.TopSurface()) {
			for (so2_analysis::so2_it so2 = so2s.begin(); so2 != so2s.end(); so2++) {
				distance = system_t::Position(*so2) - h2os.SurfaceLocation();

				(*so2)->SetAtoms();
				// for each water, find the angle that the oh bonds make with the surface normal.
				p = so_calculator(*so2);

				_alpha (distance, p.first);
				_alpha (distance, p.second);
			}
		}
		else {
			for (so2_analysis::so2_it so2 = so2s.begin(); so2 != so2s.end(); so2++) {
				distance = h2os.SurfaceLocation() - system_t::Position(*so2);

				// for each water, find the angle that the oh bonds make with the surface normal.
				(*so2)->SetAtoms();
				p = so_calculator(*so2);

				_alpha (distance, -p.first);
				_alpha (distance, -p.second);
			}
		}
	} // analysis


	void WaterOrientationNearSO2::Analysis () {

		this->LoadAll();

		Water_ptr_vec wats;
		wats.clear();

		for (Mol_it it = this->begin_mols(); it != this->end_mols(); it++) {
			if ((*it)->MolType() != Molecule::H2O) continue;
			WaterPtr wat = static_cast<Water *>(*it);
			wat->SetAtoms();
			wats.push_back(wat);
		}

		// grab the so2 of the system
		MolPtr mol = Molecule::FindByType(this->begin_mols(), this->end_mols(), Molecule::SO2);
		SulfurDioxide * so2 = static_cast<SulfurDioxide *>(mol);
		so2->SetAtoms();

		// a) for each water, find the axis between the water-oxygen and the so2-sulfur.
		// b) calculate the distance between the two molecules
		// c) calculate the angle between the two bisectors (are they aligned, anti-aligned, etc)
		// d) histogram the angle by distance

		VecR axis;	// pointing from the water-O to the so2-S
		double distance, angle;
		for (Wat_it wat = wats.begin(); wat != wats.end(); wat++) {

			axis = MDSystem::Distance ((*wat)->O(), so2->S());
			distance = axis.norm();
			angle = acos((*wat)->Bisector() < so2->Bisector()) * 180.0/M_PI;

			angles.Alpha(distance, angle);
			angles.Beta(distance, 180.0/M_PI * acos((*wat)->Bisector() < axis));
		}

	}	// analysis
	bool SO2AdsorptionWaterAngleAnalysis::residue_eq (const AtomPtr a, const std::string& res) {
		if (a->Residue() == res) 
			std::cout << "found one!" << std::endl;
		return a->Residue() == res;
	}

	void SO2AdsorptionWaterAngleAnalysis::FindInteractions () {
		// first sort the waters to find those closest to the so2 S
		nm.OrderAtomsByDistance (so2s.S());
		// then graph the closest several waters for analysis
		analysis_atoms.clear();
		Atom_it close_it = nm.closest(Atom::O);
		for (int i = 0; i < 10; i++) {
			analysis_atoms.push_back(*close_it);
			nm.next_closest(close_it, Atom::O);
		}
		analysis_atoms.push_back(so2s.S());
		// build a graph out of those atoms to find the interactions (if any) to the S
		graph.UpdateGraph(this->analysis_atoms); 

		// copy all the atoms bound to the S and the two Os
		bonded_atoms.clear();
		Atom_ptr_vec bound_atoms;

		bound_atoms = graph.BondedAtoms (so2s.S(), bondgraph::interaction);
		std::copy (bound_atoms.begin(), bound_atoms.end(), std::back_inserter(bonded_atoms));
		bound_atoms = graph.BondedAtoms (so2s.O1(), bondgraph::hbond);
		std::copy (bound_atoms.begin(), bound_atoms.end(), std::back_inserter(bonded_atoms));
		bound_atoms = graph.BondedAtoms (so2s.O2(), bondgraph::hbond);
		std::copy (bound_atoms.begin(), bound_atoms.end(), std::back_inserter(bonded_atoms));
		// remove duplicates
		std::sort(bonded_atoms.begin(), bonded_atoms.end(), std::ptr_fun(&Atom::id_cmp));
		Atom_it uniq_it = std::unique (bonded_atoms.begin(), bonded_atoms.end(), std::ptr_fun(&Atom::id_eq));
		bonded_atoms.resize(uniq_it - bonded_atoms.begin());

		// only keep waters that are bound, not other types of molecules
		bound_atoms.clear();
		for (Atom_it it = bonded_atoms.begin(); it != bonded_atoms.end(); it++) {
			if ((*it)->Residue() == "h2o")
				bound_atoms.push_back(*it);
		}
		bonded_atoms.clear();
		std::copy(bound_atoms.begin(), bound_atoms.end(), std::back_inserter(bonded_atoms));

		/*
			 bonded_atoms.erase(
			 std::remove_if(
			 bonded_atoms.begin(), 
			 bonded_atoms.end(), 
		//std::not1(
		std::bind2nd(std::ptr_fun(&SO2AdsorptionWaterAngleAnalysis::residue_eq), "h2o")), 
		bonded_atoms.end());
		*/
	}	



	void SO2AdsorptionWaterAngleAnalysis::Analysis () {

		if (!second_pass) {
			FindInteractions ();			// load up bonded_atoms with any atoms that are bound to the ref-atom

			if (bonded_atoms.size() && !first_bound_water) {
				std::cout << std::endl << "found the first interacting water at timestep " << this->_system->Timestep() << std::endl;
				first_bound_water = new Water(bonded_atoms[0]->ParentMolecule());	// Use the first of the bound atoms as it may be the closest... ?
				first_bound_water->Print();
				first_bound_water->SetAtoms();
				so2s.S()->Print();

				std::cout << "now rewinding and rerunning the analysis" << std::endl;
				this->_system->Rewind();
				std::cout << "system is rewound - starting over..." << std::endl;
				second_pass = true;
			}
		}

		else {
			// find distance of so2 to the water surface
			h2os.Reload();
			h2os.FindWaterSurfaceLocation();

			// output the posisition of the so2 wrt the water surface
			double so2_location;
			if (h2os.TopSurface()) {
				so2_location = system_t::Position(so2s.S()) - h2os.SurfaceLocation();
			} else {
				so2_location = h2os.SurfaceLocation() - system_t::Position(so2s.S())  ;
			}
			fprintf (this->output, " % 8.3f ", so2_location);

			// also the distance from the S to the O
			fprintf (this->output, " % 8.3f", MDSystem::Distance(first_bound_water->O(), so2s.S()).norm());

			// the standard deviation of distance from the surface of waters used in calculating the surface
			fprintf (this->output, " % 8.3f", h2os.SurfaceWidth());

			// distance of the reference water to the surface
			fprintf (this->output, " % 8.3f", system_t::Position(first_bound_water) - h2os.SurfaceLocation());

			// calculate the angle of the water of interest with the surface normal
			fprintf (this->output, " % 8.3f", first_bound_water->Bisector() < VecR::UnitY());

			// the (cos) angle between the h2o and so2 bisectors
			fprintf (this->output, " % 8.3f", first_bound_water->Bisector() < so2s.SO2()->Bisector());


			fprintf (this->output, "\n");
		}


	}	// analysis


	void ThetaThetaAgent::operator() 
		(VecR vec1, VecR vec2, h2o_analysis::surface_distance_t position) {
			v1 = axis;// the reference axis - perp to the surface
			if (!(position.first))
				v1 = -v1;

			theta1 = acos(vec1 < v1) * 180.0 / M_PI;
			theta2 = acos(vec2 < v1) * 180.0 / M_PI;

			//if (theta1 > theta2)
				histos (position.second, theta1, theta2);
			//else
				//histos (position.second, theta2, theta1);
		}

	void ThetaPhiAgent::operator() 
		( VecR bisector, VecR ref_bond, 
			h2o_analysis::surface_distance_t position) {

			v1 = axis;// the reference axis - perp to the surface
			if (!(position.first))
				v1 = -v1;

			//bisector.Print();
			//ref_bond.Print();

			theta = acos(bisector < v1) * 180.0 / M_PI;
			phi = Dihedral::Angle(v1,bisector,ref_bond) * 180.0 / M_PI;
			phi = fabs(phi);
			//if (phi > 90.0)
				//phi = 180.0 - phi;

			//printf ("\ntheta = %f\nphi = %f\nposition = %f\n", theta, phi, position.second);
			histos (position.second, theta, fabs(phi));
		}


} // namespace angle_analysis
