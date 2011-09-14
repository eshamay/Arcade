#include "bond-analysis.h"

namespace bond_analysis {

	void BondLengthAnalyzer::Analysis () {
		this->LoadAll();

		bondgraph::BondGraph graph;
		bondgraph::WaterCoordination_pred p (bondgraph::OH, graph);

		Water_ptr_vec wats;
		VecR_vec ohs;
		VecR z_ax = VecR::UnitZ();
		for (Mol_it it = this->begin_mols(); it != this->end_mols(); it++) {
			if ((*it)->MolType() != Molecule::H2O) continue;
			WaterPtr wat = static_cast<Water *>(*it);

			//if ((wat->OH1() < z_ax) > 0.9659258)
			if (fabs(wat->OH1() < z_ax) < 0.25882)
				ohs.push_back(wat->OH1());
			//if ((wat->OH2() < z_ax) > 0.9659258)
			if (fabs(wat->OH2() < z_ax) < 0.25882)
				ohs.push_back(wat->OH2());
		}

		/*
		// grab the molecule

		//std::cout << std::endl << wats.size() << " --> ";
		wats.erase(std::remove_if(wats.begin(), wats.end(), std::not1(p)), wats.end());
		if (wats.empty())
		std::cerr << std::endl << "didn't find any" << std::endl;
		//std::cout << wats.size() << std::endl;
		*/

		//MolPtr mol = Molecule::FindByType(this->begin_mols(), this->end_mols(), Molecule::SO2);
		//SulfurDioxide * so2 = static_cast<SulfurDioxide *>(mol);
		//so2->SetAtoms();
		/*
			 int N = 5;
		//std::sort (wats.begin(), wats.end(), WaterToSO2Distance_cmp (so2));
		std::sort (wats.begin(), wats.end(), typename system_t::molecule_position_pred(Atom::O));
		// take the top 10 waters and sort by distance to the so2
		Water_ptr_vec close_wats;
		std::copy (wats.rbegin(), wats.rbegin()+(2*N), std::back_inserter(close_wats));
		std::sort (close_wats.begin(), close_wats.end(), WaterToSO2Distance_cmp(so2)); 
		*/

		if (ohs.size() > 0) {
			double sum = 0.0;
			for (VecR_it it = ohs.begin(); it != ohs.end(); it++) {
				sum += it->Magnitude();
			}
			fprintf (this->output, "% 9.4f\n", sum/(double)ohs.size());
		}
		/*
			 if (wats.size() > 0) {
			 double sym = 0.0;
			 double antisym = 0.0;
			 for (Wat_it it = wats.begin(); it != wats.end(); it++) {
			 double oh1 = ((*it)->OH1()).Magnitude();
			 double oh2 = ((*it)->OH2()).Magnitude();

			 sym += oh1 + oh2;
			 antisym += oh1 - oh2;
			 }

			 fprintf (this->output, "% 9.4f % 9.4f\n", sym/(double)wats.size(), antisym/(double)wats.size());
			 }
			 */

		//fprintf (this->output, "%7.4f %7.4f %7.4f\n", so2->SO1().Magnitude(), so2->SO2().Magnitude(), acos(so2->Angle())*180.0/M_PI);
		//
		//
		//
		//alkane::Alkane * mal = static_cast<alkane::Alkane *>(mol);

		/*
		// print out some bond info
		double co, ch1, ch2;
		fprintf (this->output, " %12.5f %12.5f %12.5f\n",
		form->CH1().Magnitude(),
		form->CH2().Magnitude(),
		form->CO().Magnitude()
		);

		delete form;
		*/
	}

	void SO2BondingAnalyzer::FindCoordination () {
		this->LoadAll();
		this->so2s.UpdateSO2();

		graph.UpdateGraph(this->begin(), this->end());

		this->so2 = so2s.SO2();

		// for each of the three so2 atoms, find the bonds made to them
		AtomPtr atom = this->so2s.O1();	// set the atom we're interested in
		o1_bonds.clear();
		o1_bonds = graph.BondedAtoms(atom, bondgraph::hbond);
		o1 = o1_bonds.size();

		atom = this->so2s.O2();	
		o2_bonds.clear();
		o2_bonds = graph.BondedAtoms(atom, bondgraph::hbond);
		o2 = o2_bonds.size();

		atom = this->so2s.S();	
		s_bonds.clear();
		s_bonds = graph.InteractingAtoms (atom);
		s = s_bonds.size();

		coordination = coordination_t(s*10 + o1 + o2);
	}




	void SO2CoordinationAnalyzer::Analysis () {
		this->FindCoordination ();
		fprintf (this->output, "%d\n", this->s*100 + this->o1 * 10 + this->o2);
	} // analysis






	void SO2CoordinationAngleAnalyzer::Analysis () {
		this->FindCoordination ();
		printf ("%d, %d\n", this->coordination, this->coordination % 10);
		if (this->coordination / 10 == 2) {
			SulfurDioxide * so2 = this->so2s.SO2();
			so2->SetOrderAxes();

			theta (acos(so2->Bisector() < VecR::UnitZ()) * 180. / M_PI);

			double phi_v = acos(fabs(so2->Y() < VecR::UnitZ())) * 180. / M_PI;
			phi (phi_v);
		}
	}





}
