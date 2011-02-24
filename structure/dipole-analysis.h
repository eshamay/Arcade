#ifndef DIPOLE_ANALYSIS_H_
#define DIPOLE_ANALYSIS_H_

#include "analysis.h"
#include "histogram-analysis.h"

namespace md_analysis {

	/*
	template <typename T>
	class h2o_dipole_magnitude_histogram_analyzer : public histogram_analyzer<T> {
		public:
			typedef typename histogram_analyzer<T>::system_t system_t;

			h2o_dipole_magnitude_histogram_analyzer () :
				histogram_analyzer<T> (
						std::string ("Generate a histogram of H2O dipole moment magnitudes"),
						std::string ("h2o-dipole-magnitude-histogram.dat")) { }

			void Analysis (system_t& t);
	};	// H2O dipole magnitude histogram


	template <typename T>
	void h2o_dipole_magnitude_histogram_analyzer<T>::Analysis (system_t& t) {
		t.LoadWaters();

		std::for_each (t.int_wats.begin(), t.int_wats.end(), MDSystem::CalcWannierDipole);

		for (Mol_it it = t.int_wats.begin(); it != t.int_wats.end(); it++) {
			this->values.push_back((*it)->Dipole().Magnitude());
		}

		return;
	}

	*/
	template <typename T>
		class SystemDipoleAnalyzer : public AnalysisSet<T> {

			public:
				typedef Analyzer<T> system_t;
				SystemDipoleAnalyzer (system_t * t) :
					AnalysisSet<T> (t,
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

			std::for_each (WaterSystem<AmberSystem>::sys_atoms.begin(), WaterSystem<AmberSystem>::sys_atoms.end(), FixAtomCharges());

			std::transform (WaterSystem<AmberSystem>::int_wats.begin(), WaterSystem<AmberSystem>::int_wats.end(), std::back_inserter(dipoles), std::ptr_fun(&MDSystem::CalcClassicDipole));
			VecR dipole = std::accumulate (dipoles.begin(), dipoles.end(), VecR(0.0,0.0,0.0), vecr_add());

			//dipole.Print();
			fprintf (this->output, "% 12.8f % 12.8f % 12.8f\n", dipole[x], dipole[y], dipole[z]);

		}	// analysis for amber systems



	template <>
		void SystemDipoleAnalyzer<XYZSystem>::Analysis () {
			dipoles.clear();
			this->LoadAll();

			std::transform (this->begin_mols(), this->end_mols(), std::back_inserter(dipoles), std::ptr_fun(&MDSystem::CalcWannierDipole));
			VecR dipole = std::accumulate (dipoles.begin(), dipoles.end(), VecR(0.0,0.0,0.0), vecr_add());

			//dipole.Print();
			fprintf (this->output, "% 12.8f % 12.8f % 12.8f\n", dipole[x], dipole[y], dipole[z]);
		}	// analysis	for XYZ systems

}	// namespace md_analysis

#endif
