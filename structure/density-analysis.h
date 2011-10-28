#ifndef DENSITY_HISTOGRAM_H
#define DENSITY_HISTOGRAM_H

#include "analysis.h"
#include "manipulators.h"
#include "utility.h"
#include <map>
#include <set>

namespace density {

	using namespace md_analysis;

	class MolecularDensityDistribution : public AnalysisSet {
		protected:
			typedef histogram_utilities::Histogram1D<double> histo_t;
			typedef std::map <Molecule::Molecule_t,histo_t> histo_map;
			histo_map histos;

			double min,max,res;

			h2o_analysis::H2ODoubleSurfaceManipulator	h2os;

		public:
			MolecularDensityDistribution (Analyzer * t) :
				AnalysisSet (t, 
						std::string("Molecular Density Analysis"),
						std::string("MolecularDensities.dat")),
				min(-15.0), max(10.0), res(0.1),
				h2os(t) { }

			void Setup ();
			void Analysis ();
			void DataOutput ();
	};





/*
 * This analysis package calculates atomic densities throughout a system along the 3 different axes of the system (X,Y,Z). A system.cfg file is required in the working directory, and the input files from the system must have specific names (or links to the files with specific names). For Amber Systems there must be a files.prmtop file (or link) and a files.mdcrd file, using those names exactly as lsited in system.cfg.
 *
 * In the system.cfg file, aside from all the other usual (mandatory) sections, another section is required for this analysis:
 *
 * analysis.density.atom-names
 *		This section is a list of the atom names (as defined in the xyz or prmtop files) that will be analyzed and output to the data file.
 *		example: 
 *				density:
					{
						atom-names = ( "O", "H1", "H2" );
					};

 * Output will be written to the file specific in system.cfg's analysis.density.filename and will overwrite previously written data if not backed up.
 *
 * Data in the file will be a series of columns, 3 column-pairs for x,y,z directions for each atom listed. Each column pair contains the position along the direction, and the number density of the atoms. Thus, for the example listed above, the columns would be:
 *		X-pos, O_x_population, Y-pos, O_y_population, Z-pos, O_z_population, X-pos, H1_x_population, etc.
 *
 *		Analysis will use the position extents listed in system.cfg's analysis.position-range +/- 10-angstroms, with a slice resolution of analysis.resolution.position.
 */


	/*
		class AtomicDensityAnalysis : public AnalysisSet {
			public:
				typedef Analyzer system_t;

				AtomicDensityAnalysis (system_t * t) : 
					AnalysisSet (t,
							std::string("An analysis of the density of atoms in a system relative to surface position"),
							std::string("atomic-density.dat")),
					h2os(t) { } 

				virtual ~AtomicDensityAnalysis () { }

				void Setup ();
				void Analysis ();
				// For each atom type (name) in the system, the histograms in each direction will be output
				void DataOutput ();

			protected:
				std::vector<std::string> atom_name_list;
				// Every atom-name will have its own histogram of positions in the system. Each position is held as a vector to the atom site.
				typedef histogram_utilities::Histogram1D<double>	histogram_t;
				typedef std::pair<std::string, histogram_t>				histogram_map_elmt;
				typedef std::map<std::string, histogram_t>				histogram_map;
				histogram_map																			histograms;
				h2o_analysis::H2OSystemManipulator							h2os;

				class atomic_position_binner : public std::unary_function<AtomPtr,void> {
					private:
						double surface_location;
						histogram_map * histos;
						bool top_surface;
					public:
						atomic_position_binner (const double surface, histogram_map * hs, bool top) 
							: surface_location (surface), histos(hs), top_surface(top) { }
						void operator() (AtomPtr atom) const {
							histogram_map::iterator it = histos->find(atom->Name());
							if (it == histos->end()) std::cout << "couldn't find the atom named: " << atom->Name() << std::endl;
							//histogram_t* hs = &(histos->operator[] (atom->Name()));
							else {
								double distance;
								if (top_surface)
									distance = system_t::Position(atom) - surface_location;	// distance above the water surface
								else 
									distance = surface_location - system_t::Position(atom);	
								(it->second.operator()) (distance);
							}
						}
				};	// atomic density binner

		};	// atomic density class


		void AtomicDensityAnalysis::Setup () {

			AnalysisSet::Setup();

			this->_system->LoadAll();

			// grab the list of atomic names/types that will be used for the analysis and create the vector-histograms
			libconfig::Setting &atom_names = WaterSystem::SystemParameterLookup("analysis.density.atom-names");
			for (int i = 0; i < atom_names.getLength(); i++)
			{
				std::string atom_name = atom_names[i];
				atom_name_list.push_back(atom_name);

				histogram_t hs (histogram_t(WaterSystem::posmin, WaterSystem::posmax, Analyzer::posres));
				histograms.insert(histogram_map_elmt(atom_name, hs));
			}

			// narrow down the system atoms to just those with names we're looking for
			md_name_utilities::KeepByNames (this->Atoms(), atom_name_list);

		}	// Setup


		void AtomicDensityAnalysis::Analysis () { 

			h2os.FindWaterSurfaceLocation();

			atomic_position_binner binner (h2os.SurfaceLocation(), &histograms, h2os.TopSurface());
			std::for_each (this->begin(), this->end(), binner);
		}

		void AtomicDensityAnalysis::DataOutput () {

			rewind(this->output);

			// first output the header of all the atom-names
			fprintf (this->output, "position ");
			for (std::vector<std::string>::const_iterator it = atom_name_list.begin(); it != atom_name_list.end(); it++) {
				fprintf (this->output, " %s ", it->c_str());
			}
			fprintf (this->output, "\n");

			// output the data from the histograms
			double dr = system_t::posres;
			double min = WaterSystem::posmin;
			double max = WaterSystem::posmax;

			for (double r = min; r < max; r+=dr) {
				fprintf (this->output, "% 8.4f ", r);	// print the position

				for (std::vector<std::string>::const_iterator name = atom_name_list.begin(); name != atom_name_list.end(); name++) {

					histogram_t * hs = &histograms.find(*name)->second;
					fprintf (this->output, "% 8.3f ", hs->Population(r)/this->_system->timestep);
				}
				fprintf (this->output, "\n");
			}

			fflush(this->output);

		}	// Data Output
		*/

	//********************* Atomic System Density - not correlated to surface location ********************/
	/*
		class SystemDensitiesAnalysis : public AnalysisSet {
			public:
				typedef Analyzer system_t;

				SystemDensitiesAnalysis (system_t * t) : 
					AnalysisSet (t,
							std::string("An analysis of the density of atoms in a system based on atomic position"),
							std::string("system-densities.dat")) { }

				virtual ~SystemDensitiesAnalysis () { }

				void Setup ();
				void Analysis ();
				// For each atom type (name) in the system, the histograms in each direction will be output
				void DataOutput ();

			protected:
				std::vector<std::string> atom_name_list;
				// Every atom-name will have its own histogram of positions in the system. Each position is held as a vector to the atom site.
				typedef histogram_utilities::Histogram1D<double>	histogram_t;
				typedef std::pair<std::string, histogram_t>				histogram_map_elmt;
				typedef std::map<std::string, histogram_t>				histogram_map;
				histogram_map																			histograms;

				class atomic_position_binner : public std::unary_function<AtomPtr,void> {
					private:
						histogram_map * histos;
					public:
						atomic_position_binner (histogram_map * hs) : histos(hs) { }

						void operator() (AtomPtr atom) const {
							histogram_map::iterator it = histos->find(atom->Name());
							if (it == histos->end()) std::cout << "couldn't find the atom named: " << atom->Name() << std::endl;
							else { 
								double position = system_t::Position (atom);
								it->second.operator()(position);
							}
						}
				};	// atomic density binner



		};	// atomic density class


		void SystemDensitiesAnalysis::Setup () {

			AnalysisSet::Setup();

			this->_system->LoadAll();

			// grab the list of atomic names/types that will be used for the analysis and create the vector-histograms
			libconfig::Setting &atom_names = WaterSystem::SystemParameterLookup("analysis.density.atom-names");
			for (int i = 0; i < atom_names.getLength(); i++)
			{
				std::string atom_name = atom_names[i];
				atom_name_list.push_back(atom_name);

				histogram_t hs (histogram_t(WaterSystem::posmin, WaterSystem::posmax, system_t::posres));
				histograms.insert(histogram_map_elmt(atom_name, hs));
			}

			// narrow down the system atoms to just those with names we're looking for
			md_name_utilities::KeepByNames (this->Atoms(), atom_name_list);

		}	// Setup


		void SystemDensitiesAnalysis::Analysis () { 

			this->_system->LoadAll();
			std::for_each (this->begin(), this->end(), atomic_position_binner (&histograms));
		}

		void SystemDensitiesAnalysis::DataOutput () {

			rewind(this->output);

			// first output the header of all the atom-names
			fprintf (this->output, "position ");
			for (std::vector<std::string>::const_iterator it = atom_name_list.begin(); it != atom_name_list.end(); it++) {
				fprintf (this->output, " %s ", it->c_str());
			}
			fprintf (this->output, "\n");

			// output the data from the histograms
			double dr = system_t::posres;
			double min = WaterSystem::posmin;
			double max = WaterSystem::posmax;

			for (double r = min; r < max; r+=dr) {
				fprintf (this->output, "% 8.4f ", r);	// print the position

				for (std::vector<std::string>::const_iterator name = atom_name_list.begin(); name != atom_name_list.end(); name++) {

					histogram_t * hs = &histograms.find(*name)->second;
					fprintf (this->output, "% 8.3f ", hs->Population(r)/this->_system->timestep);
				}
				fprintf (this->output, "\n");
			}

			fflush(this->output);

		}	// Data Output

		*/



	//********************** H2O Surface Statistics ***************/
	/*

	// output some statistics about the water surface - i.e. standard deviation of the position of the top several waters
		class H2OSurfaceStatisticsAnalysis : public AnalysisSet {
			public:
				typedef Analyzer system_t;

				H2OSurfaceStatisticsAnalysis (system_t * t) : 
					AnalysisSet (t,
													std::string("Surface water statistical analysis"),
													std::string("h2o-surface-statistics.dat")),
					h2os(t) { 
						//h2os.ReferencePoint(WaterSystem::SystemParameterLookup("analysis.reference-location"));
						h2os.ReferencePoint(45.0);
					}

				void Analysis ();

			private:
				double	std;	// standard deviation
				h2o_analysis::H2OSystemManipulator	h2os;
				std::vector<double>		surface_locations;


		};	// h2o surface statistics analysis

		void H2OSurfaceStatisticsAnalysis::Analysis () {
			h2os.Reload();
			h2os.FindWaterSurfaceLocation();

			// find the standard deviation of the overall surface locations
			surface_locations.push_back(h2os.SurfaceLocation());
			// record the mean location and the standard deviation of the water positions, and the running SD of the surface widths
		  fprintf (this->output, "% 8.3f % 8.3f % 8.3f\n", surface_locations.back(), h2os.SurfaceWidth(), gsl_stats_sd(&surface_locations[0], 1, surface_locations.size()));
		}	// analysis

		*/





	/*
	 * This analysis is for an SO2/H2O system where the SO2 starts outside of the water and then adsorbs to the surface and is absorbed into the water phase. To watch the adsorption rate, the water phase is analyzed and the topmost and lowest water (along the main axis) is determined and then all the SO2 molecules located within the water phase (and a fudge-factor if desired) are counted up and printed out.
	 * The output will be 1 column in a datafile:
	 *		# of SO2s in the water phase
	 */
	/*
		 class so2_uptake_analysis : public AnalysisSet< Analyzer > {
		 public:
		 typedef Analyzer system_t;

		 so2_uptake_analysis () : 
		 AnalysisSet<system_t> (
		 std::string("An analysis of the adsorption of SO2 into a water phase"),
	//std::string("3d-atomic-density.dat")) { }
	std::string("so2-adsorption.dat")) { }

	virtual ~so2_uptake_analysis () { }

	void Analysis (system_t& t);
	protected:
	MolPtr high_water, low_water;
	double high_position, low_position;
	int numAdsorbed;

	void FindHighAndLowWaters (system_t& t);
	void CountSO2InWater (system_t& t);
	};

	void so2_uptake_analysis::FindHighAndLowWaters (system_t& t) {

	t.LoadWaters();

	// sort the waters according to position in the slab
	std::sort(t.int_wats.begin(), t.int_wats.end(), Analyzer<AmberSystem>::molecule_position_pred(Atom::O));
	// find the highest and the lowest at the extents of the water slab
	high_water = t.int_wats.back();
	high_position = Analyzer::Position(high_water->GetAtom(Atom::O));
	low_water = t.int_wats.front();
	low_position = Analyzer::Position(low_water->GetAtom(Atom::O));

	}

	void so2_uptake_analysis::CountSO2InWater (system_t& t) {

	// get the listing of the SO2s in the system
	t.LoadAll();
	t.int_mols.clear();
	for (Mol_it so2 = t.sys_mols.begin(); so2 != t.sys_mols.end(); so2++) {
	if ((*so2)->MolType() == Molecule::SO2)
	t.int_mols.push_back(*so2);
	}

	Atom * s;
	numAdsorbed = 0;
	double pos;
	for (Mol_it so2 = t.int_mols.begin(); so2 != t.int_mols.end(); so2++) {
	pos = (*so2)->GetAtom(Atom::O)->Position()[WaterSystem::axis];
	if (pos < high_position && pos > low_position) {
	++numAdsorbed;
	}
	}
	}

	void so2_uptake_analysis::Analysis (system_t& t) {
	this->FindHighAndLowWaters(t);
	this->CountSO2InWater(t);

	fprintf (t.Output(), "%12d %8d\n", t.Timestep(), this->numAdsorbed);
	}
	 */

}	// namespace density




#endif
