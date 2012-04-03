#ifndef TEST_H_
#define TEST_H_

#include "analysis.h"
#include "utility.h"
#include "angle-bond-analysis.h"
#include "dipole-analysis.h"
#include "neighbor-analysis.h"
#include "density-analysis.h"
#include "rdf-analysis.h"
#include "angle-analysis.h"
#include "so2-analysis.h"
#include "so2-angle-analysis.h"
#include "so2-system-analysis.h"
#include "bond-analysis.h"
#include "cycle-analysis.h"
#include "h2o-analysis.h"
#include "diacid-analysis.h"
#include "malonic-analysis.h"


typedef std::vector<double> double_vec;
typedef double_vec::const_iterator double_it;

namespace md_analysis {

	using namespace md_system;
	using namespace md_files;

	typedef enum {AMBER=0, XYZ, TRR, XTC} system_type;


	template <typename T>
		class StructureAnalyzer {
			public:
				typedef std::vector<AnalysisSet *>	analysis_vec;

				StructureAnalyzer (const int choice = -1) : _analysis_choice(choice) {
					LoadSystemAnalyses ();
					PromptForAnalysisFunction(); 
				}

				virtual ~StructureAnalyzer () {
					for (typename analysis_vec::iterator it = analyses.begin(); it != analyses.end(); it++) {
						delete *it;
					}
					delete analyzer;
					delete sys;
				}

				void SystemAnalysis (AnalysisSet& an);

			protected:
				WaterSystem * sys;
				Analyzer * analyzer;

				int _analysis_choice;
				// output a bit of text to stdout and ask the user for a choice as to which type of analysis to perform - then do it.
				void PromptForAnalysisFunction ();
				//! Loads all the analyses that are compatible with the given system-type
				void LoadSystemAnalyses ();
				// The set of possible analyses to perform on a given system
				analysis_vec analyses;
		};


	template <typename T>
	void StructureAnalyzer<T>::SystemAnalysis (AnalysisSet& an) 
	{
		// do some initial setup
		an.Setup();

		// start the analysis - run through each timestep
		for (Analyzer::timestep = 0; Analyzer::timestep < Analyzer::timesteps; Analyzer::timestep++) 
		{
			// Perform the main loop analysis that works on every timestep of the simulation
			an.Analysis ();

			// output the status of the analysis (to the screen or somewhere useful)
			analyzer->OutputStatus ();
			// Output the actual data being collected to a file or something for processing later
			if (analyzer->ReadyToOutputData())
				an.DataOutput();

			// load the next timestep
			analyzer->LoadNext();
		}
		// do one final data output to push out the finalized data set
		an.DataOutput();

		// do a little work after the main analysis loop (normalization of a histogram? etc.)
		an.PostAnalysis ();
		return;
	} // System Analysis w/ analysis set


	//! Loads all the system analyses that can be performed on Amber systems
	template <>
		void StructureAnalyzer<AmberSystem>::LoadSystemAnalyses () {
			sys = new AmberWaterSystem ();
			analyzer = new Analyzer (sys);

			analyses.push_back (new H2OAngleBondAnalysis(analyzer));									
			analyses.push_back (new density::MolecularDensityDistribution(analyzer));							
			//analyses.push_back (new SystemDensitiesAnalysis(analyzer));					
			analyses.push_back (new angle_analysis::H2OAngleAnalysis(analyzer));
			analyses.push_back (new angle_analysis::ReferenceSO2AngleAnalysis(analyzer));			
			analyses.push_back (new neighbor_analysis::SO2BondingCycleAnalysis(analyzer));	
			analyses.push_back (new neighbor_analysis::SO2HBondingAnalysis(analyzer));		
			//analyses.push_back (new md_analysis::H2OSurfaceStatisticsAnalysis(analyzer));		
			analyses.push_back (new angle_analysis::SO2AdsorptionWaterAngleAnalysis(analyzer));	
			analyses.push_back (new angle_analysis::WaterOHAngleAnalysis(analyzer));				
			analyses.push_back (new angle_analysis::OHAngleAnalysis(analyzer));				
			//analyses.push_back (new so2_analysis::SO2PositionRecorder(analyzer));				
			analyses.push_back (new angle_analysis::SOAngleAnalysis(analyzer));				
			analyses.push_back (new neighbor_analysis::SO2NearestNeighborAnalysis(analyzer));
			//analyses.push_back (new RDFByDistanceAnalyzer(analyzer));
			//analyses.push_back (new succinic::DensityDistribution(analyzer));
			//analyses.push_back (new succinic::SuccinicAcidBondAngleAnalysis(analyzer));
			//analyses.push_back (new succinic::SuccinicAcidCarbonChainDihedralAngleAnalysis(analyzer));
			//analyses.push_back (new succinic::SuccinicAcidCarbonylDihedralAngleAnalysis(analyzer));
			//analyses.push_back (new succinic::SuccinicAcidCarbonylTiltTwistAnglesAnalysis(analyzer));
			//analyses.push_back (new succinic::NeighboringWaterOrientation(analyzer));
			//analyses.push_back (new succinic::CarbonylGroupDistance(analyzer));
			//analyses.push_back (new succinic::MethyleneBisectorTilt(analyzer));
			analyses.push_back (new h2o_analysis::WaterDipoleZComponentAnalysis(analyzer));
			analyses.push_back (new h2o_analysis::DistanceAngleAnalysis(analyzer));
			analyses.push_back (new so2_angle_analysis::SO2ThetaPhiAnalyzer(analyzer));
			analyses.push_back (new so2_angle_analysis::SO2ThetaAnalyzer(analyzer));
			analyses.push_back (new h2o_analysis::WaterThetaPhiAnalysis(analyzer));
			analyses.push_back (new so2_analysis::SO2RDFAnalysis(analyzer));
			analyses.push_back (new diacid::CarboxylicDihedralPsiPsi(analyzer));
			analyses.push_back (new diacid::CarbonBackboneThetaCarboxylicDihedral(analyzer));
			analyses.push_back (new diacid::CarbonBackboneThetaPhi(analyzer));
			analyses.push_back (new diacid::MethylThetaPhiAnalysis(analyzer));
			analyses.push_back (new diacid::RDF(analyzer));
			analyses.push_back (new diacid::Dimers(analyzer));
			analyses.push_back (new diacid::Test(analyzer));
			analyses.push_back (new diacid::COTheta(analyzer));
			analyses.push_back (new diacid::CHTheta(analyzer));
			analyses.push_back (new diacid::BondLengths(analyzer));
		}


	template <>
		void StructureAnalyzer<XYZSystem>::LoadSystemAnalyses () {
			sys = new XYZWaterSystem ();
			analyzer = new Analyzer (sys);
			AnalysisSet * a;

			//analyses.push_back(new md_analysis::SystemDipoleAnalyzer<XYZSystem>(analyzer));
			//analyses.push_back(new bond_analysis::BondLengthAnalyzer (analyzer));
			analyses.push_back(new bond_analysis::SO2CoordinationAnalyzer(analyzer));
			analyses.push_back(new cycle_analysis::SO2CycleCoordinationAnalyzer(analyzer));
			analyses.push_back(new cycle_analysis::SO2CycleLifespanAnalyzer(analyzer));
			analyses.push_back(new md_analysis::RDFAnalyzer(analyzer));
			/*
			analyses.push_back(new bond_analysis::SO2CoordinationAngleAnalyzer(analyzer));
			analyses.push_back(new cycle_analysis::SO2CycleCoordinationAnalyzer(analyzer));
			analyses.push_back(new cycle_analysis::SO2CycleLifespanAnalyzer(analyzer));
			analyses.push_back(new cycle_analysis::SO2CycleCoordinateWriter(analyzer));
			analyses.push_back(new malonic::MalonicTest(analyzer));
			analyses.push_back(new malonic::BondLengths(analyzer));
			analyses.push_back(new malonic::MolecularDipole(analyzer));
			analyses.push_back(new malonic::RDF(analyzer));
			analyses.push_back(new malonic::CarbonBackboneThetaPhi(analyzer));
			analyses.push_back(new malonic::CarboxylicDihedralPsiPsi(analyzer));
			analyses.push_back(new malonic::COTheta(analyzer));
			*/

		}


	template <typename T>
		void StructureAnalyzer<T>::PromptForAnalysisFunction () {

			if (_analysis_choice < 0) {
				printf ("Choose the system analysis to perform from the list below\n\n");

				int choice = 0;
				for (typename analysis_vec::iterator it = analyses.begin(); it != analyses.end(); it++) {
					printf ("\t%d) %s\n", choice, (*it)->Description().c_str());
					++choice;
				}
				//printf ("\n\nperforming analysis (%d) using output filename \"%s\"\n", choice, analyses[choice-1]->Filename().c_str());
				exit(1);
			}

			else if (_analysis_choice < 0 || _analysis_choice >= analyses.size()) {
				std::cerr << "Analysis choice is out of range." << std::endl;
				exit(1);
			}
			else 
				SystemAnalysis(*analyses[_analysis_choice]);

			return;
		}

}	// namespace md_analysis

#endif

