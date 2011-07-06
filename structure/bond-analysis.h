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

	typedef int so2_coordination_t;

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
			so2_coordination_t	coordination, s, o1, o2;
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
			so2_coordination_t ref_coord;

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

	// Find the percentage of SO-coordination structures that are also cyclic from the S to the O through a single water, or through 2 waters.
	class SO2CycleCoordinationAnalyzer : public SO2BondingAnalyzer {
		private:
			int single_cycles, double_cycles, type_1_triple_cycles, type_2_triple_cycles, other_triple_cycles;
			int total;

			// h* and o* are the atoms on the connected waters
			// o_ref is the so2-oxygen through which the coordination happens
			AtomPtr ha, oa, hb, ob, o_ref;

			WaterPtr h2o_a, h2o_b;

			Atom_ptr_vec o_bonds, // the list of atoms bound to the o_ref
									 ob_bonds;	// all the Hs h-bonded to the h2o_b oxygen

			std::pair<int,int> CountMoleculeAtoms (CycleManipulator::cycle_list_it& cycle);
			void CheckCycles ();

			CycleManipulator cm;

		public:
			SO2CycleCoordinationAnalyzer (Analyzer * t) :
				SO2BondingAnalyzer (t, 
						std::string ("so2 cyclic-SO coordination analyzer"),
						std::string ("so2.cyclic-SO.dat")),
				single_cycles(0), double_cycles(0), 
				type_1_triple_cycles(0), type_2_triple_cycles(0), other_triple_cycles(0), total(0),
				cm(t)
	 	{ }

			void Analysis ();
			void DataOutput ();
	};
}	// namespace bond analysis


#endif
