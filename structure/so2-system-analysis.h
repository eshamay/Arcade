#ifndef SO2_ANALYSIS_H_
#define SO2_ANALYSIS_H_

#include "h2o-analysis.h"

namespace so2_analysis {

	using namespace md_analysis;


	// a convenience class for working with systems comprised of at least 1 SO2 molecule and a whole bunch of waters
	template <typename T>
		class SO2SystemManipulator : public SystemManipulator<T> {
			public:
				typedef Analyzer<T> system_t;

				SO2SystemManipulator (system_t * t) : SystemManipulator<T>(t) { 

					this->_system->LoadAll();

					this->FindSO2 ();

					this->so2->SetAtoms();
					this->s = so2->S();
					this->o1 = so2->O1();
					this->o2 = so2->O2();

					so2s.clear();
					// load all the system so2s into the container - in case
					for (Mol_it it = this->_system->sys_mols.begin(); it != this->_system->sys_mols.end(); it++) {
						if ((*it)->Name() == "so2") {
							SulfurDioxide * mol = new SulfurDioxide(*it);
							mol->SetAtoms();
							so2s.push_back(mol);
						}
					}
					printf ("\nFound %zu Sulfur Dioxides in the system\n", so2s.size());
				}

				virtual ~SO2SystemManipulator () { 
					delete this->so2;
					for (std::vector<SulfurDioxide *>::iterator it = so2s.begin(); it != so2s.end(); it++) {
						delete *it;
					}
				}

				// The method by which the SO2 of interest is found in the system
				virtual void FindSO2 () {
					// find the so2 in the system and set some pointers up
					int id = WaterSystem<T>::SystemParameterLookup("analysis.reference-molecule-id");
					MolPtr mol = this->_system->sys_mols[id];
					//MolPtr mol = Molecule::FindByType(this->_system->sys_mols, Molecule::SO2);
					this->so2 = new SulfurDioxide(mol);
				}

				SulfurDioxide * SO2 () { return so2; }
				AtomPtr S () { return s; }
				AtomPtr O1 () { return o1; }
				AtomPtr O2 () { return o2; }


				std::vector<SulfurDioxide *>::iterator	begin () { return so2s.begin(); }
				std::vector<SulfurDioxide *>::iterator	end () { return so2s.end(); }

			protected:
				SulfurDioxide * so2;	// the sulfur dioxide of interest
				AtomPtr s, o1, o2;	// sulfur dioxide's atoms
				std::vector<SulfurDioxide *> so2s;
		};


}	// namespace so2_analysis


#endif
