#ifndef BOND_ANALYSIS_H
#define BOND_ANALYSIS_H

#include "analysis.h"
#include "dipole-analysis.h"

namespace bond_analysis {

	using namespace md_system;
	using namespace md_analysis;

		class BondLengthAnalyzer : public AnalysisSet {

			public:
				typedef Analyzer system_t;
				BondLengthAnalyzer (system_t * t) :
					AnalysisSet (t,
							std::string ("bondlength analysis"),
							std::string ("so2-bond+angles.normal_modes.dat")) { }

				void Analysis ();

		};	 // bond length analyzer class

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



	// find the coordination of the so2 molecule for the 3 atoms separately
		class SO2BondingAnalyzer : public AnalysisSet {
			protected:
				so2_analysis::XYZSO2Manipulator		so2s;

			public:
				typedef Analyzer system_t;

				typedef struct {
					AtomPtr atom;
					double bondlength;
				} bond;
					
				SO2BondingAnalyzer (system_t * t) :
					AnalysisSet (t,
							std::string ("so2 Coordination analyzer"),
							std::string ("so2-coordination.dat")),
					so2s(t) {
						so2s.Initialize();
					}

				void Analysis ();

				// comparator for bondlengths
				class Bond_cmp : public std::binary_function<bond,bond,bool> {
					public:
						bool operator() (const bond& left, const bond& right) const {
							return left.bondlength < right.bondlength;
						}
				};
		};	 // bond length analyzer class

		void SO2BondingAnalyzer::Analysis () {
			this->LoadAll();
			this->so2s.UpdateSO2();

			bondgraph::BondGraph graph;
			//bondgraph::WaterCoordination_pred p (bondgraph::OH, graph);

			graph.UpdateGraph(this->begin(), this->end());

			int o1, o2, s;
			// for each of the three so2 atoms, find the bonds made to them
			AtomPtr atom = this->so2s.O1();	// set the atom we're interested in
			o1 = graph.NumHBonds(atom);

			atom = this->so2s.O2();	
			o2 = graph.NumHBonds(atom);

			atom = this->so2s.S();	
			s = graph.NumInteractions(atom);

			fprintf (this->output, "%d\n", s*100 + o1*10+ o2);
		} // analysis



}	// namespace bond analysis


#endif
