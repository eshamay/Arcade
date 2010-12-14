#ifndef __RDF_ANALYSIS_H
#define __RDF_ANALYSIS_H

#include "analysis.h"
#include "threading.h"

#define NUMTHREADS	16

namespace md_analysis {

	void * build_histo (void * thread_data);

	struct thread_data_t {
		Atom **atoms;									// system atom array
		double min, max, res;					// rdf parameters
		Atom::Element_t	elmt1, elmt2;	// elements to process
		int * big_histo;
		int thread_id;
		int blocklow, blockhigh;

		thread_data_t (
				Atom ** _atoms,
				double _min, double _max, double _res, 
				Atom::Element_t _elmt1, Atom::Element_t _elmt2,
				int * _big_histo, int id, int low, int high)
			: 
				atoms(_atoms),
				min(_min), max(_max), res(_res),
				elmt1(_elmt1), elmt2(_elmt2),
				big_histo(_big_histo), thread_id(id),
				blocklow(low), blockhigh(high)
	 	{ }

	};

	pthread_mutex_t histogram_mutex;

	template <class T>
		class rdf_analysis : public AnalysisSet< Analyzer<T> > {
			public:
				typedef Analyzer<T> system_t;
				rdf_analysis () :
					AnalysisSet<system_t> (
							std::string("RDF Analysis"),
							std::string("rdf.dat")) { }

				virtual ~rdf_analysis () { 
					for (int i = 0; i < NUMTHREADS; i++)
						delete _thread_data[i];
					pthread_mutex_destroy(&histogram_mutex);
					pthread_exit(NULL);
				}

				void Setup (system_t& t);
				void Analysis (system_t& t);
				void DataOutput (system_t& t);

				typedef std::vector<int>  _histogram_t;
				_histogram_t							_histogram;

				thread_data_t			*_thread_data[NUMTHREADS];
				pthread_t					_thread[NUMTHREADS];

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

			//printf ("total size = %zu\n", t.sys_atoms.size());
			Atom ** atoms = &t.sys_atoms[0];
			for (int i = 0; i < NUMTHREADS; i++) {
				int low = threads::block_low(i, NUMTHREADS, t.sys_atoms.size());
				int high = threads::block_high(i, NUMTHREADS, t.sys_atoms.size()); 
				//printf ("%d) %d -- %d\n", i, low, high);
				_thread_data[i] = new thread_data_t (
						atoms,
						min, max, res, 
						elmt1, elmt2,
						&_histogram[0], i, 
						low, high);
			}

		} // rdf analysis setup



	void * build_histo (void * thread_data) {
		thread_data_t * data = (thread_data_t *) thread_data;
		//printf ("%d starting\n", data->thread_id);

		double atomic_distance;
		std::vector<int>	histo;
		histo.clear();
		histo.resize((int)((data->max-data->min)/data->res), 0);
	 	

		// go through each atom pair and process the rdf histogram
		Atom *ai, *aj;
		for (int i = data->blocklow; i < data->blockhigh; i++) {
			ai = data->atoms[i];

			for (int j = i+1; j <= data->blockhigh; j++) {
				aj = data->atoms[j];

				// first check if the pair of atoms is one of the pairs to analyze - is the pair in the list of pairs?
				if ( ((ai->Element() == data->elmt1) && (aj->Element() == data->elmt2)) ||
						((ai->Element() == data->elmt2) && (aj->Element() == data->elmt1)) ) {

					// for each pair in the list, bin the distance between the atoms into the proper histogram
					atomic_distance = MDSystem::Distance(ai, aj).norm();
					if (atomic_distance > data->max || atomic_distance < data->min) continue;

					++histo[(int)((atomic_distance-data->min)/data->res)];
				}
			}
		}

		pthread_mutex_lock (&histogram_mutex);
		for (int i = 0; i < histo.size(); i++) {
			data->big_histo[i] += histo[i];
		}
		pthread_mutex_unlock (&histogram_mutex);

		pthread_exit(NULL);
		return;
	}




	template <typename T>
		void rdf_analysis<T>::Analysis(system_t& t) {
			pthread_attr_t attr;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

			int rc;
			for (int i = 0; i < NUMTHREADS; i++) {
				rc = pthread_create(&_thread[i], &attr, build_histo, (void *)_thread_data[i]);
				if (rc) {
					printf("ERROR; return code from pthread_create() is %d\n", rc);
					exit(-1);
				}
			}

			pthread_attr_destroy(&attr);
			void * status;
			for (int i = 0; i < NUMTHREADS; i++) {
				rc = pthread_join(_thread[i], &status);
				if (rc) {
					printf("ERROR; return code from pthread_join() is %d\n", rc);
					exit(-1);
				}
			}
		}	// rdf analysis analysis



	template <typename T>
		void rdf_analysis<T>::DataOutput(system_t& t) {
			rewind (t.Output());

			double dV, n, N;
			double total_volume = 4.0/3.0 * M_PI * pow(max, 3);
			N = std::accumulate (_histogram.begin(), _histogram.end(), 0);

			// go through each value in the RDF position range
			for (double r = min; r < max; r += res) {
				dV = 4.0 * M_PI * pow(r, 2) * res;	// volume of the shell being considered
				n = _histogram[(int)((r-min)/res)];	// the population of the shell

				fprintf (t.Output(), "%12.3f %12.3f\n", r, n * total_volume / dV / N);
			}

			fflush(t.Output());
			return;
		}	// rdf analysis data output

} // namespace md_analysis
#endif
