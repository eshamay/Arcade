#ifndef MALONIC_ANALYSIS_H_
#define MALONIC_ANALYSIS_H_

#include "manipulators.h"
#include "molecule-analysis.h"
#include "histogram-analysis.h"
#include "rdf-analysis.h"

namespace malonic {

	using namespace md_analysis;

	typedef enum {
		c1o1, c1oh1,
		c2o2, c2oh2,
		h1o2, h2o1,
		h1oh1, h2oh2,
		h1oh2, h2oh1,
		o1waterh, o2waterh
	} bond_t;


	class CarbonBackboneThetaPhi : public molecule_analysis::MalonicAnalysis {
		protected:
			Histogram2DAgent	angles;
			ThreeAtomGroup ccc;
			VecR axis, v1, bisector;
			double theta, phi;

		public:
			CarbonBackboneThetaPhi (system_t * t) :
				molecule_analysis::MalonicAnalysis (t,
						std::string("Diacid carbon C-C-C backbone theta-phi angle analysis"),
						std::string ("")),
				angles (
						std::string ("CarbonBackbone-Theta-Phi.dat"),	// filename
						1.0, 179.0, 2.5, // theta range
						0.0, 90.0, 1.0), // phi range
				axis(VecR::UnitZ())
		{ }

			void MoleculeCalculation ();
			void DataOutput () { 
				DoNothing2D func;
				//DivideByLeftSineDegrees func;
				angles.OutputData(func); 
			}
	};

	class COTheta : public molecule_analysis::MalonicAnalysis {
		private:
			Histogram1DAgent	angles;
			VecR axis, co_bond;
			double theta;

		public:
			COTheta (Analyzer * t) :
				molecule_analysis::MalonicAnalysis (t,
						std::string("Diacid Carbonyl C=O theta vs distance in water"),
						std::string ("")),
				angles (
						std::string ("Carbonyl-DistanceTheta.dat"),	
						5.0, 175.0, 2.0), // theta
				axis(VecR::UnitZ()) { }

			void MoleculeCalculation ();

			void DataOutput () {
				angles.OutputData();
			}
	};


	// looking at the two dihedrals of malonic acid
	class CarboxylicDihedralPsiPsi : public molecule_analysis::MalonicAnalysis {
		protected:
			Histogram2DAgent	angles;
			ThreeAtomGroup ccc;
			VecR axis, v1, bisector;
			double psi1, psi2;

		public:
			CarboxylicDihedralPsiPsi (Analyzer * t) :
				molecule_analysis::MalonicAnalysis (t,
						std::string("Malonic O=C-C-C psi1-psi2 dihedral angle analysis"),
						std::string ("")),
				angles (
						std::string ("Carboxylic-Psi-Psi.dat"), // filename 
						0.0, 180.0, 2.5, // psi1
						0.0, 180.0, 2.5), // psi2
				axis(VecR::UnitZ())
		{ }

			void MoleculeCalculation ();
			void DataOutput () {
				DoNothing2D func;
				//DivideByLeftSineDegrees func;
				angles.OutputData(func); 
			}
	};

	class RDF : public molecule_analysis::MalonicAnalysis {
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

	class MolecularDipole : public molecule_analysis::MalonicAnalysis {
		public:
			MolecularDipole (Analyzer * t);

			void MoleculeCalculation ();
			void DataOutput () { }

	};


	class BondLengths : public molecule_analysis::MalonicAnalysis {

		private:
			typedef std::map< bond_t, std::vector<double> > bondlength_map;
			bondlength_map lengths;

			void CalcDistance (AtomPtr atom1, AtomPtr atom2, bond_t bond);
			void OutputDataPoint (bond_t bond, int timestep);

		public:
			BondLengths (Analyzer * t);

			void MoleculeCalculation ();
			void DataOutput ();

	}; // malonic test class

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

