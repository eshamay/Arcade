#ifndef MANIPULATORS_H_
#define MANIPULATORS_H_

namespace md_analysis {

		// manipulator superclass
		class SystemManipulator {
			public:

				SystemManipulator (Analyzer * sys) : _system(sys) { 
					_system->LoadAll();
					std::copy(WaterSystem::sys_atoms.begin(), WaterSystem::sys_atoms.end(), std::back_inserter(all_atoms));
					std::copy(WaterSystem::sys_mols.begin(), WaterSystem::sys_mols.end(), std::back_inserter(all_mols));
					Reload();
				}

				virtual ~SystemManipulator () { }

				// reload all the analysis atoms and molecules
				virtual void Reload () {
					analysis_atoms.clear();
					analysis_mols.clear();
					std::copy(all_atoms.begin(), all_atoms.end(), std::back_inserter(analysis_atoms));
					std::copy(all_mols.begin(), all_mols.end(), std::back_inserter(analysis_mols));
				}

				Atom_it atoms_begin() { return analysis_atoms.begin(); }
				Atom_it atoms_end() { return analysis_atoms.end(); }

				Mol_it mols_begin() { return analysis_mols.begin(); }
				Mol_it mols_end() { return analysis_mols.end(); }

			protected:
				Analyzer * _system;

				Atom_ptr_vec		all_atoms;
				Mol_ptr_vec			all_mols;

				Atom_ptr_vec		analysis_atoms;
				Mol_ptr_vec			analysis_mols;
		};	// system manipulator

}


#endif
