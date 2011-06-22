#ifndef DIPOLE_ANALYSIS_H_
#define DIPOLE_ANALYSIS_H_

#include "analysis.h"

namespace md_analysis {

	using namespace so2_analysis;

	template <typename T>
		class SystemDipoleAnalyzer : public AnalysisSet {

			public:
				typedef Analyzer system_t;
				SystemDipoleAnalyzer (system_t * t) :
					AnalysisSet (t,
							std::string ("System Dipole Analysis"),
							std::string ("system-dipole.dat")) { }

				void Analysis ();

				VecR_vec dipoles;

				typedef Atom_ptr_vec::iterator	Atom_ncit;	// non-const iterators



				class FixAtomCharges : public std::unary_function <AtomPtr, void> {
					public:
						void operator() (AtomPtr a) const {
							switch (a->Element()) {
								case Atom::H : a->Charge(0.365); break;
								case Atom::O : a->Charge(-0.73); break;
								case Atom::S : a->Charge(0.47); break;
								default: break;
							}
						}
				};

		}; // system dipole analyzer

	template <>
		void SystemDipoleAnalyzer<AmberSystem>::Analysis () {
			dipoles.clear();
			this->LoadAll();
			this->_system->LoadWaters();

			//for (wannier_it wi = begin_wanniers(); wi != end_wanniers(); wi++) {
				//wi->Print();
			//}

			std::for_each (this->begin(), this->end(), FixAtomCharges());

			std::transform (this->begin_wats(), this->end_wats(), std::back_inserter(dipoles), std::ptr_fun(&MDSystem::CalcClassicDipole));
			VecR dipole = std::accumulate (dipoles.begin(), dipoles.end(), VecR(0.0,0.0,0.0), vecr_add());

			//dipole.Print();
			fprintf (this->output, "% 12.8f % 12.8f % 12.8f\n", dipole[x], dipole[y], dipole[z]);

		}	// analysis for amber systems



	template <>
		void SystemDipoleAnalyzer<XYZSystem>::Analysis () {
			dipoles.clear();
			this->LoadAll();

			//std::for_each (this->begin_mols(), this->end_mols(), std::mem_fun(&Molecule::Print));

			// set the system waters
			Water_ptr_vec wats;
			for (Mol_it it = this->begin_mols(); it != this->end_mols(); it++) {
				if ((*it)->MolType() == Molecule::H2O) {
					WaterPtr wat = static_cast<WaterPtr>(*it);
					wat->SetAtoms();
					wats.push_back(wat);
				}
			}

			// grab the so2
			SulfurDioxide * so2;
			for (Mol_it it = this->begin_mols(); it != this->end_mols(); it++) {
				if ((*it)->MolType() == Molecule::SO2) {
					so2 = static_cast<SulfurDioxide *>(*it);
					so2->SetAtoms();
				}
			}
			
			// calculate the dipole of each molecule
			std::transform (this->begin_mols(), this->end_mols(), std::back_inserter(dipoles), std::ptr_fun(&MDSystem::CalcWannierDipole));

			// then grab the 5 waters nearest the so2
			// by sorting them according to the distance to the so2
			std::sort (wats.begin(), wats.end(), MoleculeToReferenceDistance_cmp (so2));

			//VecR dipole = std::accumulate (dipoles.begin(), dipoles.end(), VecR(0.0,0.0,0.0), vecr_add());
			VecR dipole (0.,0.,0.);
			dipole += so2->Dipole();
			for (int i = 0; i < 5; i++) {
					dipole += wats[i]->Dipole();
			}
			fprintf (this->output, "% 12.8f % 12.8f % 12.8f\n", dipole[x], dipole[y], dipole[z]);
		}	// analysis	for XYZ systems

}	// namespace md_analysis

#endif
