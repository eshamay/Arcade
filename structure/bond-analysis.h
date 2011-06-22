#ifndef BOND_ANALYSIS_H
#define BOND_ANALYSIS_H

#include "analysis.h"
#include "manipulators.h"

namespace bond_analysis {

	using namespace md_system;
	using namespace md_analysis;

		class BondLengthAnalyzer : public AnalysisSet {

			public:
				typedef Analyzer system_t;
				BondLengthAnalyzer (system_t * t) :
					AnalysisSet (t,
							std::string ("bondlength analysis"),
							std::string ("so2-bond+angles.normal_modes.dat")) { }

				void Analysis ();

		};	 // bond length analyzer class




	// find the coordination of the so2 molecule for the 3 atoms separately
		class SO2BondingAnalyzer : public AnalysisSet {
			protected:
				so2_analysis::XYZSO2Manipulator		so2s;

			public:
				typedef Analyzer system_t;

				typedef struct {
					AtomPtr atom;
					double bondlength;
				} bond;
					
				SO2BondingAnalyzer (system_t * t) :
					AnalysisSet (t,
							std::string ("so2 Coordination analyzer"),
							std::string ("so2-coordination.dat")),
					so2s(t) {
						so2s.Initialize();
					}

				void Analysis ();

				// comparator for bondlengths
				class Bond_cmp : public std::binary_function<bond,bond,bool> {
					public:
						bool operator() (const bond& left, const bond& right) const {
							return left.bondlength < right.bondlength;
						}
				};
		};	 // bond length analyzer class




}	// namespace bond analysis


#endif
