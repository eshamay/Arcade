#ifndef MOLECULE_ANALYSIS_H_
#define MOLECULE_ANALYSIS_H_

#include "analysis.h"
#include "manipulators.h"
#include "histogram-analysis.h"

namespace molecule_analysis {

	using namespace md_analysis;

	template <typename T>
	class SingleMoleculeAmberAnalysis : public AnalysisSet {

		protected:
			h2o_analysis::H2ODoubleSurfaceManipulator	h2os;
			double com;
			h2o_analysis::surface_distance_t position;
			Mol_ptr_vec analysis_mols;
			Molecule::Molecule_t moltype;
			T * mol;

		public:
			typedef Analyzer system_t;
			SingleMoleculeAmberAnalysis (system_t * t, std::string desc, std::string fn) : 
				AnalysisSet (t, desc, fn),
				h2os(t) { }

			virtual void Setup () { h2os.Reload(); }
			virtual void Analysis ();

			virtual void MoleculeCalculation () = 0;

			virtual void PreCalculation () {
				algorithm_extra::copy_if (this->begin_mols(), this->end_mols(), std::back_inserter(analysis_mols), member_functional::mem_fun_eq (&Molecule::MolType, moltype));
			}

			virtual void PostCalculation () { return; }
	}; 


	template <typename T>
	void SingleMoleculeAmberAnalysis<T>::Analysis () {
		h2os.FindWaterSurfaceLocation();

		// load up the molecules to be analyzed
		analysis_mols.clear();
		PreCalculation ();

		for (Mol_it it = analysis_mols.begin(); it != analysis_mols.end(); it++) {
			mol = static_cast<T *>(*it);
			MoleculeCalculation ();
		}

		PostCalculation ();

		return;
	}




	class SO2Analysis : public SingleMoleculeAmberAnalysis<SulfurDioxide> {
		public:
			SO2Analysis (Analyzer * t, std::string desc, std::string fn) : 
				SingleMoleculeAmberAnalysis<SulfurDioxide> (t, desc, fn) { this->moltype = Molecule::SO2; }
	}; // so2 analysis



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

} // namespace molecule analysis

#endif
