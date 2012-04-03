#ifndef DIACID_ANALYSIS_H_
#define DIACID_ANALYSIS_H_

#include "molecule-analysis.h"
#include "angle-analysis.h"
#include "rdf-analysis.h"
#include "dimergraph.h"

namespace diacid {

	using namespace md_system;
	using namespace md_analysis;

	class Test : public molecule_analysis::DiacidAnalysis {
		public:
			Test (system_t * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid test"),
						std::string ("")) { }

			void MoleculeCalculation ();
	};

	typedef enum {
		c1o1, c1oh1,
		c2o2, c2oh2,
		h1o2, h2o1,
		h1oh1, h2oh2,
		h1oh2, h2oh1,
		o1waterh, o2waterh
	} bond_t;

	class BondLengths : public molecule_analysis::DiacidAnalysis {
		private:
			typedef std::map< bond_t, histogram_utilities::Histogram1D<double>* > bondlength_map;
			bondlength_map lengths;
			double min,max,res;
			double distance;

			void CalcDistance (AtomPtr atom1, AtomPtr atom2, bond_t bond);
			void OutputDataPoint (bond_t bond, double position);

		public:
			BondLengths (Analyzer * t);

			void MoleculeCalculation ();
			void DataOutput ();
	};



	class CarbonBackboneThetaPhi : public molecule_analysis::DiacidAnalysis {
		protected:
			angle_analysis::ThetaPhiAgent angles;
			ThreeAtomGroup ccc;

		public:
			CarbonBackboneThetaPhi (system_t * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid carbon C-C-C backbone theta-phi angle depth-slice analysis"),
						std::string ("")),
				angles (
						std::string ("./carbonbackbone-theta-phi/theta-phi."), // filename prefix
						std::string (".dat"),	// suffix
						-8.0, 4.0, 2.0, // depths
						1.0, 179.0, 2.5, // theta range
						0.0, 90.0, 1.0) // phi range
		{ }

			void MoleculeCalculation ();
			void DataOutput () { angles.DataOutput(); }
	};


	class CarboxylicThetaPhiAnalysis : public molecule_analysis::DiacidAnalysis {
		protected:
			angle_analysis::ThetaPhiAgent angles;

		public:
			CarboxylicThetaPhiAnalysis (system_t * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid carbonyl theta-phi angle depth-slice analysis"),
						std::string ("")),
				angles (
						std::string ("./carbonyl-theta-phi/theta-phi."), // filename prefix
						std::string (".dat"),	// suffix
						-14.0, 4.0, 2.0, // depths
						5.0, 175.0, 2.5, // theta range
						0.0, 180.0, 2.0) // phi range
		{ }

			void MoleculeCalculation ();
			void DataOutput () { angles.DataOutput(); }
	};


	class Dimers : public molecule_analysis::DiacidAnalysis {
		public:
			Dimers (system_t * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Look at the bonding between the dicarboxylic acids"),
						std::string ("DiacidDimers.dat")),
				initialized(false) { }

			void PreCalculation ();
			void MoleculeCalculation ();
			void PostCalculation ();

		protected:
			DimerGraph	_graph;
			bool initialized;
	};




	class MethylThetaPhiAnalysis : public molecule_analysis::DiacidAnalysis {
		protected:
			angle_analysis::ThetaPhiAgent angles;

		public:
			MethylThetaPhiAnalysis (system_t * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid methyl theta-phi angle depth-slice analysis"),
						std::string ("")),
				angles (
						std::string ("./methyl-theta-phi/theta-phi."), // filename prefix
						std::string (".dat"),	// suffix
						-14.0, 4.0, 2.0, // depths
						5.0, 175.0, 2.5, // theta range
						0.0, 90.0, 1.0) // phi range
		{ }

			void MoleculeCalculation ();
			void DataOutput () { angles.DataOutput(); }
	};


	class RDF : public molecule_analysis::DiacidAnalysis {
		protected:
			RDFAgent	rdf_alc;
			RDFAgent	rdf_carb;
			double distance;
		public:
			RDF (Analyzer * t);
			void MoleculeCalculation ();
			void DataOutput () { 
				rdf_alc.OutputData(); 
				rdf_carb.OutputData();
			}
	};


	class CarbonBackboneThetaCarboxylicDihedral : public molecule_analysis::DiacidAnalysis {

		protected:
			angle_analysis::ThetaPhiAgent angles;
			double theta;
			ThreeAtomGroup ccc;
			VecR v1, axis;
			std::pair<double,double> psi;


		public:
			CarbonBackboneThetaCarboxylicDihedral (Analyzer * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid O=C-C-C dihedral v C-C-C theta angle depth-slice analysis"),
						std::string ("")),
				angles (
						std::string ("./carbonbackbone-theta-carbonyl-psi/theta-psi."), // filename prefix
						std::string (".dat"),	// suffix
						-14.0, 4.0, 2.0, // depths
						5.0, 175.0, 2.5, // theta range
						-180.0, 180.0, 4.0), // dihedral range
				axis(VecR::UnitY())
		{ }

			void MoleculeCalculation ();
			void DataOutput () { angles.DataOutput(); }
	};


	class CHTheta : public molecule_analysis::DiacidAnalysis {
		private:
			angle_analysis::PositionThetaAgent	angles;
			VecR v1;
			const VecR axis;
			double theta1, theta2;

		public:
			CHTheta (Analyzer * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid Methyl C-H theta vs distance in water"),
						std::string ("")),
				angles (
						std::string ("MethylCHDistanceTheta.dat"),	
						-10.0, 4.0, 0.2, // depths
						5.0, 175.0, 2.0), // theta
				axis(VecR::UnitY()) { }

			void MoleculeCalculation ();

			void DataOutput () {
				DivideByRightSineDegrees func;
				//DoNothing2D func;
				angles.OutputData(func);
			}
	};


	class COTheta : public molecule_analysis::DiacidAnalysis {
		private:
			angle_analysis::PositionThetaAgent	angles;
			VecR v1;
			const VecR axis;
			double theta1, theta2;

		public:
			COTheta (Analyzer * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid Carbonyl C=O theta vs distance in water"),
						std::string ("")),
				angles (
						std::string ("CarbonylDistanceTheta.dat"),	
						-10.0, 4.0, 0.2, // depths
						//5.0, 175.0, 2.5), // theta
						5.0, 175.0, 2.0), // theta
				axis(VecR::UnitY()) { }

			void MoleculeCalculation ();

			void DataOutput () {
				DivideByRightSineDegrees func;
				//DoNothing2D func;
				angles.OutputData(func);
			}
	};

	// looking at the two dihedrals of malonic acid
	class CarboxylicDihedralPsiPsi : public molecule_analysis::DiacidAnalysis {
		protected:
			angle_analysis::PsiPsiAgent angles;
			VecR v1,v2,v3;
			double theta1, theta2;
			unsigned int molcount;

		public:
			CarboxylicDihedralPsiPsi (Analyzer * t) :
				molecule_analysis::DiacidAnalysis (t,
						std::string("Diacid O=C-C-C psi1-psi2 dihedral angle depth-slice analysis"),
						std::string ("")),
				angles (
						std::string ("./carboxylic-dihedral-psi-psi/psi-psi."), // filename prefix
						std::string (".dat"),	// suffix
						-8.0, 4.0, 2.0, // depths
						0.0, 180.0, 2.5, // psi1
						0.0, 180.0, 2.5), // psi2
				molcount(0)
		{ }

			void MoleculeCalculation ();
			void DataOutput () { angles.DataOutput(); }
	};

} // namespace 

#endif
