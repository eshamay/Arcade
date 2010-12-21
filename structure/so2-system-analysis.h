#ifndef SO2_ANALYSIS_H_
#define SO2_ANALYSIS_H_

#include "h2o-analysis.h"

namespace so2_analysis {

	using namespace md_analysis;
	//using namespace md_system;
	//using namespace md_files;


	// a convenience class for working with systems comprised of at least 1 SO2 molecule and a whole bunch of waters
	template <typename T>
		class SO2SystemAnalyzer : public h2o_analysis::H2OSystemAnalyzer<T> {
			public:
				typedef Analyzer<T> system_t;

				SO2SystemAnalyzer (std::string desc, std::string filename) :
					h2o_analysis::H2OSystemAnalyzer<T> (desc, filename) { }

				virtual ~SO2SystemAnalyzer () { 
					delete this->so2;
				}

				// The method by which the SO2 of interest is found in the system
				virtual void FindSO2 (system_t& t) {
					// find the so2 in the system and set some pointers up
					MolPtr mol = Molecule::FindByType(t.sys_mols, Molecule::SO2);
					this->so2 = new SulfurDioxide(mol);
				}

				// things to do after the standard setup
				virtual void PostSetup (system_t& t) { 
					this->FindSO2 (t);

					this->so2->SetAtoms();
					this->s = so2->S();
					this->o1 = so2->O1();
					this->o2 = so2->O2();

					// grab the first location of the so2 as a reference for later analyses
					this->reference_point = system_t::Position(this->s);
				}

			protected:
				SulfurDioxide * so2;	// the sulfur dioxide of interest
				AtomPtr s, o1, o2;	// sulfur dioxide's atoms
		};

}	// namespace so2_analysis
#endif
