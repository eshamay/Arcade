#ifndef MALONIC_ANALYSIS_H_
#define MALONIC_ANALYSIS_H_

#include "analysis.h"

namespace malonic_analysis {

	using namespace md_analysis;

	class MalonicTest : public AnalysisSet {

		public:
			MalonicTest (Analyzer * t) :
				AnalysisSet (t, std::string ("Malonic Tester"), std::string("")) { }

			virtual void Analysis ();

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

} // namespace malonic analysis



#endif

