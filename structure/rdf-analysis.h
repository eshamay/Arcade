#ifndef __RDF_ANALYSIS_H
#define __RDF_ANALYSIS_H

#include "analysis.h"
#include "threading.h"


namespace md_analysis {

	template <class T>
		class rdf_analysis : public AnalysisSet< Analyzer<T> > {
			public:
				typedef Analyzer<T> system_t;
				rdf_analysis () :
					AnalysisSet<system_t> (
							std::string("RDF Analysis"),
							std::string("rdf.dat")) { }

				virtual ~rdf_analysis () { }

				void Setup (system_t& t);
				void Analysis (system_t& t);
				void DataOutput (system_t& t);

			protected:
				std::vector<boost::thread>	_threads;
				std::vector<int>					_histogram;

				double min, max, res;

				Atom::Element_t		elmt1, elmt2;
		};

	template <typename T>
		void rdf_analysis<T>::Setup(system_t& t) {
			// get the rdf parameters from the configuration file
			min = t.SystemParameterLookup("analysis.rdf.minimum");
			max = t.SystemParameterLookup("analysis.rdf.maximum");
			res = t.SystemParameterLookup("analysis.rdf.resolution");

			_histogram.resize ((int)((max-min)/res), 0);

			// form all the histograms and the mapping between element pairs with the histograms to be used for rdf binning
			libconfig::Setting &atompairs = t.SystemParameterLookup("analysis.rdf.atom-pairs");
			elmt1 = Atom::String2Element (atompairs[0][0]);
			elmt2 = Atom::String2Element (atompairs[0][1]);


			// load all the atoms to work with for the duration of the analysis
			t.LoadAll();

		} // rdf analysis setup

	template <typename T>
		void rdf_analysis<T>::Analysis(system_t& t) {
			/*
			double atomic_distance;
			// go through each atom pair and 
			for (Atom_it it = t.sys_atoms.begin(); it != t.sys_atoms.end() - 1; it++) {
				for (Atom_it jt = it+1; jt != t.sys_atoms.end(); jt++) {
					// first check if the pair of atoms is one of the pairs to analyze - is the pair in the list of pairs?
					element_pair_t ep = std::make_pair((*it)->Element(), (*jt)->Element());
					element_pair_list::iterator epl_it = pair_utility::PairListMember (ep, element_pairs.begin(), element_pairs.end());

					if (epl_it == element_pairs.end()) continue;	// if the pair isn't in the list, then don't bin it!

					// for each pair in the list, bin the distance between the atoms into the proper histogram
					atomic_distance = MDSystem::Distance(*it, *jt).norm();
					histograms.find(ep)->second(atomic_distance);
				}
			}
			*/

			/*
			// specific analysis to bag only the O1 and O2 rdfs for so2
			t.LoadWaters();
			Atom::KeepByElement (t.int_atoms, Atom::H);
			AtomPtr o1 = so2->O1();
			AtomPtr o2 = so2->O2();
			for (Atom_it it = t.int_atoms.begin(); it != t.int_atoms.end(); it++) {
			atomic_distance = MDSystem::Distance(*it, o1).norm();
			histograms.begin()->second(atomic_distance);
			atomic_distance = MDSystem::Distance(*it, o2).norm();
			histograms.begin()->second(atomic_distance);
			}
			*/

		}	// rdf analysis analysis



	template <typename T>
		void rdf_analysis<T>::DataOutput(system_t& t) {
			/*
			rewind (t.Output());

			double minimum = t.SystemParameterLookup("analysis.rdf.minimum");
			double maximum = t.SystemParameterLookup("analysis.rdf.maximum");
			double resolution = t.SystemParameterLookup("analysis.rdf.resolution");

			double dV, n, N;
			double total_volume = 4.0/3.0 * M_PI * pow(maximum, 3);

			// go through each value in the RDF position range
			for (double r = minimum; r < maximum; r += resolution) {
				// print the position to the first column
				fprintf (t.Output(), "%12.3f ", r);
				dV = 4.0 * M_PI * pow(r, 2) * resolution;	// volume of the shell being considered

				// then for each element pair analyzed, print a column with the value of g(r) for the given position
				for (element_histogram_map_t::iterator it = histograms.begin(); it != histograms.end(); it++) {
					n = it->second.Population(r);	// the population of the shell
					N = it->second.Count();				// total population of the histogram

					fprintf (t.Output(), "%12.3f ", n * total_volume / dV / N);
				}
				fprintf(t.Output(), "\n");
			}

			fflush(t.Output());
			return;
			*/
		}	// rdf analysis data output

} // namespace md_analysis
#endif
