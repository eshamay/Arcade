#include "analysis.h"

namespace md_analysis {
	double	Analyzer::posres;
	int			Analyzer::posbins;

	double	Analyzer::angmin;
	double	Analyzer::angmax;
	double	Analyzer::angres;
	int			Analyzer::angbins;

	int			Analyzer::timestep;
	int 		Analyzer::timesteps;
	int			Analyzer::restart;

	Analyzer::Analyzer (WaterSystem * water_sys) :
		sys(water_sys),
		output_freq(WaterSystem::SystemParameterLookup("analysis.output-frequency")) 
	
	{ 

		try {
			sys->Initialize();
		}
		catch (const libconfig::SettingTypeException &stex) {
			std::cerr << "WaterSystem::_InitializeSystem() -- Something wrong with initializing the system. Try checking filenames in the system.cfg" << std::endl;
			exit(EXIT_FAILURE);
		}
		Analyzer::posres = WaterSystem::SystemParameterLookup("analysis.position-range")[2];
		Analyzer::posbins = int((WaterSystem::posmax - WaterSystem::posmin)/posres);

		Analyzer::angmin = WaterSystem::SystemParameterLookup("analysis.angle-range")[0];
		Analyzer::angmax = WaterSystem::SystemParameterLookup("analysis.angle-range")[1];
		Analyzer::angres = WaterSystem::SystemParameterLookup("analysis.angle-range")[2];
		Analyzer::angbins = int((angmax - angmin)/angres);

		Analyzer::timestep = 0;
		Analyzer::timesteps = WaterSystem::SystemParameterLookup("system.timesteps");
		Analyzer::restart = WaterSystem::SystemParameterLookup("analysis.restart-time");

		status_updater.Set (output_freq, timesteps, 0);
		this->registerObserver(&status_updater);

		this->_OutputHeader();
	} // Analyzer ctor



	void Analyzer::_OutputHeader () const {

		printf ("Analysis Parameters:\n\tScreen output frequency = 1/%d\n\n\tPosition extents for analysis:\n\t\tMin = % 8.3f\n\t\tMax = % 8.3f\n\t\tPosition Resolution = % 8.3f\n\n\tPrimary Axis = %d\nNumber of timesteps to be analyzed = %d\n",
				output_freq, WaterSystem::posmin, WaterSystem::posmax, Analyzer::posres, int(WaterSystem::axis), Analyzer::timesteps);

#ifdef ANALYZER_SURFACE_AVG
		printf ("\n\nThe analysis is averaging about the two interfaces located as:\n\tLow  = % 8.3f\n\tHigh = % 8.3f\n\n", int_low, int_high);
#endif
		return;
	}

	Analyzer::~Analyzer () { return; }

	void Analyzer::OutputStatus () {
		this->notifyObservers ();
		return;
	}


	void Analyzer::LoadNext () {
		this->sys->LoadNext();
		return;
	}





	/* Find the periodic-boundary-satistfying location of a molecule, atom, vector, or raw coordinate along the reference axis */
	double Analyzer::Position (const MolPtr mol) {
		return Analyzer::Position(mol->ReferencePoint());
	}

	double Analyzer::Position (const AtomPtr patom) {
		//return Analyzer::Position(patom->Position());
		return WaterSystem::AxisPosition(patom);
	}

	double Analyzer::Position (const VecR& v) {
		double position = v[WaterSystem::axis];
		return Analyzer::Position(position);
	}

	double Analyzer::Position (const double d) {
		double pos = d;
		if (pos < WaterSystem::pbcflip) pos += MDSystem::Dimensions()[WaterSystem::axis];
		return pos;
	}

	void AnalysisSet::OpenDataOutputFile () {

		output = (FILE *)NULL;
		if (filename == "") {
			printf ("\nAnalysis:: No filename specified for dataoutput.\n");
		}
		else {
			output = fopen(filename.c_str(), "w");

			if (output == (FILE *)NULL) {
				printf ("AnalysisSet::_OpenDataOutputFile() - couldn't open the data output file, \"%s\", given in the analysis set!\n", filename.c_str());
				exit(1);
			}

			printf ("\nOutputting data to \"%s\"\n", filename.c_str());
		}
		return;
	}

} // namespace md_analysis
