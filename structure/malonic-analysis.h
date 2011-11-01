#ifndef MALONIC_ANALYSIS_H_
#define MALONIC_ANALYSIS_H_

#include "manipulators.h"
#include "molecule-analysis.h"
#include "histogram-analysis.h"

namespace malonic {

	using namespace md_analysis;

	class MalonicTest : public molecule_analysis::MalonicAnalysis {

		public:
			MalonicTest (Analyzer * t) :
				MalonicAnalysis (t, 
						std::string ("Malonic Tester"), 
						std::string("")) { }

			virtual void MoleculeCalculation ();

	}; // malonic test class


	class MalonicBondLengthAnalysis : public AnalysisSet {
		public: 
			MalonicBondLengthAnalysis (Analyzer * t) :
				AnalysisSet (t, 
						std::string ("Malonic bond-length analyzer"), 
						std::string("malonic.bond-lengths.dat")) 
		{ 
			// create the pairs of atom ids of the various atom groups
			int num = 6;
			int first[] = {0,4,1,5,6,2};
			int second[] = {2,6,8,9,8,9};

			atom_ids.clear();
			for (int i = 0; i < num; i++) {
				atom_ids.push_back (std::make_pair(first[i],second[i]));
			}
		}

			void Analysis ();
			double BondLength (const int id1, const int id2);

		private:
			double distance;
			AtomPtr atom1, atom2;
			std::vector< std::pair<int,int> >	atom_ids;
	};




	class MethyleneTilt : public AnalysisSet {

		protected:
			Histogram1DAgent	histo;
			VecR bisector;
			double angle;
			alkane::MalonicAcid * mal;

		public:
			MethyleneTilt (Analyzer * t) :
				AnalysisSet (t, 
						std::string ("Malonic CH2 tilt"),
						std::string("")),
				histo (std::string ("ch2-tilt.dat"),
						5.0, 175.0, 2.5) { }

			void Analysis ();
			void DataOutput () { histo.OutputData(); }
	};




} // namespace malonic analysis



#endif

