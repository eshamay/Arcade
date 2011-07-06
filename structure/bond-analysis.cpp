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

			coordination = s*10 + o1 + o2;
		}


		void SO2CoordinationAnalyzer::Analysis () {
			this->FindCoordination ();
			fprintf (this->output, "%d\n", this->s*100 + this->o1*10+ this->o2);
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





		void SO2CycleCoordinationAnalyzer::Analysis () {
			this->FindCoordination ();


			cm.SetReferenceAtom(so2->S());
			cm.SetCycleSize (30);
			cm.BuildGraph();
			cm.ParseCycles();

			if (this->coordination / 10 >= 1 && this->coordination % 10 >= 1) {
				CheckCycles ();
				++total;
			}
		}



		void SO2CycleCoordinationAnalyzer::CheckCycles () {

			CycleManipulator::cycle_type_it cycle_type = cm.cycle_type_begin();
			if (*cycle_type == CycleManipulator::HALFBRIDGE) {

				for (CycleManipulator::cycle_list_it cycle = cm.cycle_begin(); cycle != cm.cycle_end(); cycle++) {
					CycleManipulator::cycle_it _first = cycle->begin(); _first++;	
					AtomPtr atom1 = *_first;
					AtomPtr atom2 = cycle->back();

					Molecule::Molecule_t mol1_t, mol2_t;
					mol1_t = atom1->ParentMolecule()->MolType(); 
					mol2_t = atom2->ParentMolecule()->MolType(); 
					if ((mol1_t == Molecule::SO2 && mol2_t != Molecule::SO2) 
							|| (mol1_t != Molecule::SO2 && mol2_t == Molecule::SO2)) {

						cm.FindUniqueMembers(*cycle);

						int num_mol = cm.NumUniqueMoleculesInCycle(); 
						int num_atoms = cm.NumUniqueAtomsInCycle(); 

						if (num_mol == 2 && num_atoms == 4)
							++single_cycles;
						else if (num_mol == 3)
							++double_cycles;

						else if (num_mol == 4 && num_atoms == 8) {
							std::pair<int,int> minmax = CountMoleculeAtoms(cycle);
							if (minmax.first == 1 && minmax.second == 3)
								++type_1_triple_cycles;
							else if (minmax.first == 2 && minmax.second == 2)
								++type_2_triple_cycles;
							else {
								std::cerr << "check it -- " << minmax.first << "," << minmax.second << "  !" << std::endl;
								std::for_each (cycle->begin(), cycle->end(), std::mem_fun(&Atom::Print));
								exit(1);
							}
						}
						else if (num_mol == 4)
							++other_triple_cycles;

						else if (num_mol <= 4) {
							std::cerr << "hey boo!" << std::endl;
							std::cerr << num_mol << " mols and " << num_atoms << " atoms" << std::endl;
							std::for_each (cycle->begin(), cycle->end(), std::mem_fun(&Atom::Print));
							exit(1);
						}
					}

					cycle_type++;
				}

			}

			return;
		}

		std::pair<int,int> SO2CycleCoordinationAnalyzer::CountMoleculeAtoms (CycleManipulator::cycle_list_it& cycle) {

			typedef std::vector<int> int_vec;
			typedef int_vec::iterator int_it;

			int_vec mol_ids;
			std::transform (cycle->begin(), cycle->end(), std::back_inserter(mol_ids), std::mem_fun<int>(&Atom::MolID));
			std::sort(mol_ids.begin(), mol_ids.end());

			int_vec unique_ids;
			std::copy(mol_ids.begin(), mol_ids.end(), std::back_inserter(unique_ids));
			int_it it = std::unique(unique_ids.begin(), unique_ids.end());
			unique_ids.resize(it - unique_ids.begin());

			int_vec count;

			for (int_it jt = unique_ids.begin(); jt != unique_ids.end(); jt++) {
				count.push_back((int)std::count(mol_ids.begin(), mol_ids.end(), *jt));
			}

			int_it min = std::min_element(count.begin(), count.end());
			int_it max = std::max_element(count.begin(), count.end());
			return std::make_pair (*min,*max);
		}


		void SO2CycleCoordinationAnalyzer::DataOutput () {
			rewind (this->output);
			fprintf (this->output, "%d %d %d %d %d %d\n", 
					total, single_cycles, double_cycles, type_1_triple_cycles, type_2_triple_cycles, other_triple_cycles);
			fflush (this->output);
		}
}
