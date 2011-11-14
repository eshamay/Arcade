#ifndef MOLECULE_ANALYSIS_H_
#define MOLECULE_ANALYSIS_H_

#include "analysis.h"
#include "manipulators.h"
#include "histogram-analysis.h"

namespace molecule_analysis {

	using namespace md_analysis;


	//////////////// SINGLE MOLECULE ANALYSIS ///////////////////
	template <typename T>
	class SingleMoleculeAnalysis : public AnalysisSet {
		protected:
			Mol_ptr_vec analysis_mols;
			Molecule::Molecule_t moltype;
			T * mol;

		public:

			typedef Analyzer system_t;
			SingleMoleculeAnalysis (system_t * t, std::string desc, std::string fn) : 
				AnalysisSet (t, desc, fn) { }

			virtual void PreCalculation () {
				analysis_mols.clear();
				algorithm_extra::copy_if (this->begin_mols(), this->end_mols(), std::back_inserter(analysis_mols), member_functional::mem_fun_eq (&Molecule::MolType, moltype));
			}

			virtual void MoleculeCalculation () = 0;
			virtual void PostCalculation () { return; }

			void Analysis ();
	};


	template <typename T>
		void SingleMoleculeAnalysis<T>::Analysis () {
			// load up the molecules to be analyzed
			// or do what needs to be done first
			PreCalculation ();

			// run through the analysis
			for (Mol_it it = analysis_mols.begin(); it != analysis_mols.end(); it++) {
				mol = static_cast<T *>(*it);
				MoleculeCalculation ();
			}
			PostCalculation ();
			return;
		}





	//////////////// XYZ MOLECULE ANALYSIS ///////////////////
	//
	//
	template <typename T>
		class SingleMoleculeXYZAnalysis : public SingleMoleculeAnalysis<T> {
			public:
				typedef Analyzer system_t;
				SingleMoleculeXYZAnalysis (system_t * t, std::string desc, std::string fn) : 
					SingleMoleculeAnalysis<T> (t, desc, fn) { }

				virtual void PreCalculation () {
					this->LoadAll();
				}
		};



	//////////////// AMBER MOLECULE ANALYSIS ///////////////////
	//
	//
	template <typename T>
		class SingleMoleculeAmberAnalysis : public SingleMoleculeAnalysis<T> {

			protected:
				h2o_analysis::H2ODoubleSurfaceManipulator	h2os;
				double com;
				h2o_analysis::surface_distance_t position;

			public:
				typedef Analyzer system_t;
				SingleMoleculeAmberAnalysis (system_t * t, std::string desc, std::string fn) : 
					SingleMoleculeAnalysis<T> (t, desc, fn),
					h2os(t) { }

				virtual void Setup () { 
					h2os.Reload(); 
					AnalysisSet::Setup(); 
				}

				virtual void PreCalculation () {
					h2os.FindWaterSurfaceLocation();
					SingleMoleculeAnalysis<T>::PreCalculation ();
					//this->analysis_mols.clear();
					//algorithm_extra::copy_if (this->begin_mols(), this->end_mols(), std::back_inserter(this->analysis_mols), member_functional::mem_fun_eq (&Molecule::MolType, moltype));
				}
		}; 



	class SO2Analysis : public SingleMoleculeAmberAnalysis<SulfurDioxide> {
		public:
			SO2Analysis (Analyzer * t, std::string desc, std::string fn) : 
				SingleMoleculeAmberAnalysis<SulfurDioxide> (t, desc, fn) { this->moltype = Molecule::SO2; }
	}; // so2 analysis

	class SO2XYZAnalysis : public SingleMoleculeXYZAnalysis<SulfurDioxide> {
		public:
			SO2XYZAnalysis (Analyzer * t, std::string desc, std::string fn) : 
				SingleMoleculeXYZAnalysis<SulfurDioxide> (t, desc, fn) { this->moltype = Molecule::SO2; }
	}; // so2 xyz analysis



	class DiacidAnalysis : public SingleMoleculeAmberAnalysis<alkane::Diacid> {
		public:
			DiacidAnalysis (Analyzer * t, std::string desc, std::string fn) : 
				SingleMoleculeAmberAnalysis<alkane::Diacid> (t, desc, fn) { this->moltype = Molecule::DIACID; }
	}; // diacid analysis



	class SuccinicAcidAnalysis : public SingleMoleculeAmberAnalysis<alkane::SuccinicAcid> {
		public:
			SuccinicAcidAnalysis (Analyzer * t, std::string desc, std::string fn) : 
				SingleMoleculeAmberAnalysis<alkane::SuccinicAcid> (t, desc, fn) { this->moltype = Molecule::SUCCINIC; }
	}; // succinic analysis

	class H2OAnalysis : public SingleMoleculeAmberAnalysis<Water> {
		public:
			H2OAnalysis (Analyzer * t, std::string desc, std::string fn) : 
				SingleMoleculeAmberAnalysis<Water> (t, desc, fn) { this->moltype = Molecule::H2O; }
	}; // succinic analysis




	class MalonicAnalysis : public SingleMoleculeXYZAnalysis<alkane::MalonicAcid> {
		public:
			MalonicAnalysis (Analyzer * t, std::string desc, std::string fn) : 
				SingleMoleculeXYZAnalysis<alkane::MalonicAcid> (t, desc, fn) { }

			virtual void PreCalculation () {
				SingleMoleculeXYZAnalysis<alkane::MalonicAcid>::PreCalculation();
				// copy in all the different types of malonic acid
				analysis_mols.clear();
				algorithm_extra::copy_if (this->begin_mols(), this->end_mols(), std::back_inserter(analysis_mols), member_functional::mem_fun_eq (&Molecule::MolType, Molecule::MALONIC));
				algorithm_extra::copy_if (this->begin_mols(), this->end_mols(), std::back_inserter(analysis_mols), member_functional::mem_fun_eq (&Molecule::MolType, Molecule::MALONATE));
				algorithm_extra::copy_if (this->begin_mols(), this->end_mols(), std::back_inserter(analysis_mols), member_functional::mem_fun_eq (&Molecule::MolType, Molecule::DIMALONATE));
			}

	}; // malonic xyz analysis



} // namespace molecule analysis



#endif
