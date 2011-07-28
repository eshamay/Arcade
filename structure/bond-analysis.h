#ifndef BOND_ANALYSIS_H
#define BOND_ANALYSIS_H

#include "analysis.h"
#include "manipulators.h"
#include "histogram-analysis.h"

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


	typedef struct {
		AtomPtr atom;
		double bondlength;
	} bond_t;

	typedef int atom_coordination_t;

	typedef enum {
		unbound = 0,
		O = 1, OO = 2, OOO = 3, OOOO=4,
		S=10, SO = 11, SOO=12, SOOO=13,
		SS=20, SSO=21, SSOO=22, SSOOO=23,
		SSS=30, SSSO=31, SSSOO=32, SSSOOO=33
	} coordination_t;

	// find the coordination of the so2 molecule for the 3 atoms separately
	class SO2BondingAnalyzer : public AnalysisSet {

		public:
			typedef Analyzer system_t;

			SO2BondingAnalyzer (system_t * t, const std::string& desc, const std::string& fn) :
				AnalysisSet (t, desc, fn),
				so2s(t) {
					so2s.Initialize();
				}
			virtual ~SO2BondingAnalyzer () { }

			virtual void Analysis () = 0;
			void FindCoordination ();

		protected:
			so2_analysis::XYZSO2Manipulator		so2s;
			SulfurDioxide * so2;
			atom_coordination_t	coordination, s, o1, o2;
			Atom_ptr_vec				s_bonds, o1_bonds, o2_bonds;
			bondgraph::BondGraph graph;

	};	 // bond length analyzer class


	// comparator for bondlengths
	class Bond_cmp : public std::binary_function<bond_t,bond_t,bool> {
		public:
			bool operator() (const bond_t& left, const bond_t& right) const {
				return left.bondlength < right.bondlength;
			}
	};

	class SO2CoordinationAnalyzer : public SO2BondingAnalyzer {
		public:

			typedef enum {
				UNBOUND = 0, O = 1, OO = 2, OOO = 3, OOOO = 4, OOOOO = 5,
				S = 10, SO = 11, SOO = 12, SOOO = 13, SOOOO = 14,
				SS = 20, SSO = 21, SSOO = 22, SSOOO = 23, SSOOOO = 24,
				SSS = 30, SSSO = 31, SSSOO = 32
			} coordination_t;

			SO2CoordinationAnalyzer (Analyzer * t) :
				SO2BondingAnalyzer (t, 
						std::string ("so2 Coordination analyzer"),
						std::string ("so2-coordination.dat")) { }

			void Analysis ();
	};


	class SO2CoordinationAngleAnalyzer : public SO2BondingAnalyzer {
		protected:
			Histogram1DAgent theta;
			Histogram1DAgent phi;
			atom_coordination_t ref_coord;

		public:
			SO2CoordinationAngleAnalyzer (Analyzer * t) :
				SO2BondingAnalyzer (t, 
						std::string ("so2 Angle by Coordination"),
						std::string ("")),
				theta ("so2.theta.OO.dat",
						0.0, 180.0, 2.0),
				phi ("so2.phi.OO.dat",
						-180.0, 180.0, 1.0),
				ref_coord(1)
	 	{ }

			void Analysis ();
			void DataOutput () { theta.OutputData(); phi.OutputData(); return; }
	};

}	// namespace bond analysis


#endif
