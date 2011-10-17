#ifndef SUCCINIC_ANALYSIS_H
#define SUCCINIC_ANALYSIS_H

#include "analysis.h"
#include "manipulators.h"
#include "histogram-analysis.h"
#include "neighbor-analysis.h"
#include "molecule-analysis.h"

namespace succinic {

	using namespace md_analysis;
	// 2D analysis of the dihedral angle of the carbon chain in succinic acid as a function
	// of position relative to the water surface.
	class SuccinicAcidDihedralAngleAnalysis : public molecule_analysis::SuccinicAcidAnalysis {

		protected:
			Histogram2DAgent		histo;
			double angle;

		public:
			typedef Analyzer system_t;
			SuccinicAcidDihedralAngleAnalysis (system_t * t, std::string desc, std::string fn) : 
				molecule_analysis::SuccinicAcidAnalysis (t, desc, std::string("")),
				histo (fn,
						//-20.0, 10.0, 0.2,
						0.0, 180.0, 1.0,
						0.0,180.0,1.0)	// dihedral angle values are folded to be between 0 and 180.
				// negative angle values don't count because of molecular symmetry
		{ }

			virtual void DataOutput () { histo.OutputData(); }

	};	// succinic dihedral analysis



	// 2D analysis of the dihedral angle of the carbon chain in succinic acid as a function
	// of position relative to the water surface.
	class SuccinicAcidCarbonChainDihedralAngleAnalysis : public SuccinicAcidDihedralAngleAnalysis {

		public:
			typedef Analyzer system_t;
			SuccinicAcidCarbonChainDihedralAngleAnalysis (system_t * t) :
				SuccinicAcidDihedralAngleAnalysis (t,
						std::string("succinic acid carbon-chain dihedral vs distance to surface"),
						std::string ("dihedrals.v.distance.dat")) { }

			void MoleculeCalculation ();

	};	// succinic dihedral analysis



	// 2D analysis of the dihedral angle of the carbonyl group in succinic acid as a function
	// of position relative to the water surface.
	class SuccinicAcidCarbonylDihedralAngleAnalysis : public SuccinicAcidDihedralAngleAnalysis {
		protected:
			VecR axis, v1, v2, v3;
			double twist;

		public:
			typedef Analyzer system_t;
			SuccinicAcidCarbonylDihedralAngleAnalysis (system_t * t) :
				SuccinicAcidDihedralAngleAnalysis (t,
						std::string("succinic acid carbonyl dihedral vs distance to surface"),
						std::string ("carbonyl-bisector-dihedral.v.distance.dat")),
				axis(VecR::UnitY()) { }

			void DihedralCalculation (AtomPtr aliphatic, AtomPtr carbonyl, AtomPtr oxygen);
			void MoleculeCalculation ();

	};	// succinic dihedral analysis


	class DensityDistribution : public molecule_analysis::SuccinicAcidAnalysis {
		protected:
			Histogram1DAgent histo;

		public:
			typedef Analyzer system_t;
			DensityDistribution (system_t * t) :
				SuccinicAcidAnalysis (t,
						std::string("density distribution of succinic acid"),
						std::string ("")),
				histo (
						std::string ("density-distribution.dat"),
						-15.0, 7.0, 0.1) { }


			void MoleculeCalculation ();
			void DataOutput () { histo.OutputData(); }

	};


	// cuts the surface into slices and gets the tilt-twist histogram for each depth
	// tilt of the carboxylic acid O-C-O bisector vector with the surface normal
	// twist of the O-C-O about the bisector
	class SuccinicAcidCarbonylTiltTwistAnglesAnalysis : public SuccinicAcidDihedralAngleAnalysis {
		protected:
			VecR axis, v1, v2, v3;
			double tilt, twist;
			Multi2DHistogramAgent	histos;

		public:
			typedef Analyzer system_t;
			SuccinicAcidCarbonylTiltTwistAnglesAnalysis (system_t * t) :
				SuccinicAcidDihedralAngleAnalysis (t,
						std::string("succinic acid carbonyl bisector twist vs dihedral twist"),
						std::string ("temp")),
				axis(VecR::UnitY()),
				histos (
						-14.0, 4.0, 1.0,
						5.0,175.0,2.5,
						5.0,175.0,2.5,
						std::string("./carbonyl-tilt-twist-histos/carbonyl-tilt-twist."),
						std::string(".dat")) { }

			void DataOutput () {
				DivideByLeftSineDegrees func;
				histos.DataOutput(func);
			}

			void DihedralCalculation (AtomPtr aliphatic, AtomPtr carbonyl, AtomPtr oxygen);
			void MoleculeCalculation ();

	};	// succinic dihedral analysis

	// cuts the surface into slices and gets the angle of neighboring waters around the carbonyl groups
	// dimensions are: depth under surface, distance of water to the carbonyl group, tilt of water molecule wrt surface
	class NeighboringWaterOrientation : public SuccinicAcidDihedralAngleAnalysis {
		protected:
			double oo_distance, tilt, depth;
			Multi2DHistogramAgent	histos;
			AtomPtr oxygen;
			VecR axis;

		public:
			typedef Analyzer system_t;
			NeighboringWaterOrientation (system_t * t) :
				SuccinicAcidDihedralAngleAnalysis (t,
						std::string("orientation of waters around succinic acid carboxylic acid headgroups"),
						std::string ("temp")),
				histos (
						-7.0, 3.0, 1.0,	// depths
						1.0, 8.0, 0.1,	// distance from carboxy to water
						5.0,175.0,2.5,	// tilt of water
						std::string("./nearby-water-orientation/alcohol-oxygen."),
						std::string(".dat")),
				axis(VecR::UnitY()) { }

			void DataOutput () {
				DivideByRightSineDegrees func;
				histos.DataOutput(func);
			}

			void MoleculeCalculation ();
			void AngleDepthCalculation (AtomPtr oxygen, WaterPtr wat);

	};	// succinic dihedral analysis

	// angle of a bond vector relative to the surface normal
	class SuccinicAcidBondAngleAnalysis : public molecule_analysis::SuccinicAcidAnalysis {

		private:
			Histogram2DAgent		histo;
			double angle;
			VecR bond, axis;

			void AngleDistanceCalculation (AtomPtr, AtomPtr);	// does the calculation for each bond

		public:
			typedef Analyzer system_t;
			SuccinicAcidBondAngleAnalysis (system_t * t) : 
				molecule_analysis::SuccinicAcidAnalysis (t,
						std::string("succinic acid bond-angle analysis"),
						std::string("")),
				histo (std::string ("CO-alcohol-bond-angle.v.distance.dat"), 
						-15.0, 4.0, 0.2,	// distance to interface
						5.0,175.0,2.0),	// angle
				axis(VecR::UnitY())
		{ }

			void DataOutput () { 
				DivideByRightSineDegrees func;
				histo.OutputData(func); 
			}
			void MoleculeCalculation ();

	};	// succinic bond angle

	class CarbonylGroupDistance : public molecule_analysis::SuccinicAcidAnalysis {
		protected:
			Histogram2DAgent histo;
			double distance;

		public:
			typedef Analyzer system_t;
			CarbonylGroupDistance (system_t * t) :
				SuccinicAcidAnalysis (t,
						std::string("Distribution of distances between carboxylic head groups in succinic acid"),
						std::string ("")),
				histo (
						std::string ("head-group-distance.dat"),
						-12.0, 5.0, 0.1,
						0.0, 5.0, 0.02) { }


			void MoleculeCalculation ();
			void DataOutput () { histo.OutputData(); }
	};

	class MethyleneBisectorTilt : public molecule_analysis::SuccinicAcidAnalysis {
		protected:
			Histogram2DAgent histo;
			double angle;
			VecR ax;

		public:
			typedef Analyzer system_t;
			MethyleneBisectorTilt (system_t * t) :
				SuccinicAcidAnalysis (t,
						std::string("Succinic acid CH2 bisector tilt distribution"),
						std::string ("")),
				histo (
						std::string ("succinic-CH2.dat"),
						-12.0, 3.0, 0.5,
						5.0, 175.0, 2.5) { }


			void MoleculeCalculation ();
			void DataOutput () { 
				DivideByRightSineDegrees func;
				histo.OutputData(func); 
				//histo.OutputData(); 
			}
	};
}	 // namespace succinic

#endif
