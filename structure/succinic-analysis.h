#ifndef SUCCINIC_ANALYSIS_H
#define SUCCINIC_ANALYSIS_H

#include "analysis.h"
#include "manipulators.h"
#include "histogram-analysis.h"
#include "neighbor-analysis.h"
#include "molecule-analysis.h"

namespace succinic {

	// 2D analysis of the dihedral angle of the carbon chain in succinic acid as a function
	// of position relative to the water surface.
	class SuccinicAcidDihedralAngleAnalysis : public molecule_analysis::SuccinicAcidAnalysis {

		protected:
			Histogram2DAgent		histo;
			double angle;
			double com;
			h2o_analysis::surface_distance_t	distance;

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
			virtual void SuccinicAcidCalculation (alkane::SuccinicAcid *) = 0;

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

			void SuccinicAcidCalculation (alkane::SuccinicAcid *);

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
			void SuccinicAcidCalculation (alkane::SuccinicAcid *);

	};	// succinic dihedral analysis




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
			void SuccinicAcidCalculation (alkane::SuccinicAcid *);

	};	// succinic dihedral analysis


	// angle of a bond vector relative to the surface normal
	class SuccinicAcidBondAngleAnalysis : public molecule_analysis::SuccinicAcidAnalysis {

		private:
			Histogram2DAgent		histo;
			double angle;
			double com;
			h2o_analysis::surface_distance_t	distance;
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
			void SuccinicAcidCalculation (alkane::SuccinicAcid *);

	};	// succinic bond angle


}	 // namespace succinic
