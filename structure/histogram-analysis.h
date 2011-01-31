#ifndef HISTOGRAM_ANALYSIS_H_
#define HISTOGRAM_ANALYSIS_H_

#include "utility.h"

namespace md_analysis {


	// used for analyses that have histograms  - should help with organization and output
	class Histogram1DAgent {
		protected:
			histogram_utilities::Histogram1D<double>	histogram;
			std::string				filename;
			FILE *						output; 

		public:

			Histogram1DAgent (const std::string file, const double min, const double max, const double res)
				: 
					histogram (min, max, res),
					filename (file),
					output ((FILE *)NULL) { }

			~Histogram1DAgent () { }

			virtual void OutputData ();

			double max () const { return histogram.Max(); }
			double min () const { return histogram.Min(); }
			double res () const { return histogram.Resolution(); }
			int size () const { return histogram.Size(); }

			void operator() (const double i) { histogram(i); }

			double Count () const { return histogram.Count(); }
			double Population (const double i) const { return histogram.Population(i); }

	}; // histogram 1d agent

	void Histogram1DAgent::OutputData () {

		output = fopen(filename.c_str(), "w");
		if (!output) {
			std::cerr << "couldn't open the file for output --> \" " << filename << " \"" << std::endl;
			exit(1);
		}
		rewind(output);

		for (double alpha = histogram.Min(); alpha < histogram.Max(); alpha += histogram.Resolution()) {
			fprintf (output, "% 12.3f % 12.3f\n", alpha, histogram.Population(alpha));
		}
		fflush(output);
		fclose(output);
		return;
	}



	// used for analyses that have histograms to help with output
	class Histogram2DAgent {
		protected:
			histogram_utilities::Histogram2D<double>	histogram;
			std::string				filename;
			FILE *						output; 

		public:

			typedef histogram_utilities::Histogram2D<double>::pair_t	pair_t;

			Histogram2DAgent (const std::string file,
					const double min_1, const double max_1, const double res_1,
					const double min_2, const double max_2, const double res_2)
				: 
					histogram (
							std::make_pair(min_1,min_2),
							std::make_pair(max_1,max_2),
							std::make_pair(res_1,res_2)),
					filename (file),
					output ((FILE *)NULL) { }

			~Histogram2DAgent () { }

			virtual void OutputData ();
			virtual void OutputDataMatrix ();

			pair_t max () const { return histogram.max; }
			pair_t min () const { return histogram.min; }
			pair_t res () const { return histogram.resolution; }
			pair_t size () const { return histogram.size; }

			void operator() (const double i, const double j) { histogram(i,j); }

			double Count (const double i) const { return histogram.Count(i); }
			double TotalCount() const { return histogram.TotalCount(); }
			double Population (const double i, const double j) const { return histogram.Population(i,j); }

	}; // histogram 2d agent

	void Histogram2DAgent::OutputData () {

		output = fopen(filename.c_str(), "w");
		if (!output) {
			std::cerr << "couldn't open the file for output --> \" " << filename << " \"" << std::endl;
			exit(1);
		}
		rewind(output);

		for (double alpha = histogram.min.first; alpha < histogram.max.first; alpha += histogram.resolution.first) {
			for (double beta = histogram.min.second; beta < histogram.max.second; beta += histogram.resolution.second) {
				fprintf (output, "% 12.3f % 12.3f % 12f\n", alpha, beta, histogram.Population(alpha, beta));
			}
		}
		fflush(output);
		fclose(output);
		return;
	}


	// rows are first index, columns are second
	void Histogram2DAgent::OutputDataMatrix () {

		output = fopen(filename.c_str(), "w");
		if (!output) {
			std::cerr << "couldn't open the file for output --> \" " << filename << " \"" << std::endl;
			exit(1);
		}
		rewind(output);

		for (double alpha = histogram.min.first; alpha < histogram.max.first; alpha += histogram.resolution.first) {
			for (double beta = histogram.min.second; beta < histogram.max.second; beta += histogram.resolution.second) {
				fprintf (output, " % 12.3f ", histogram.Population(alpha, beta));
			}
			printf ("\n");
		}
		fflush(output);
		fclose(output);
		return;
	}



	/*
		 template <typename T, typename U=double>
		 class histogram_analysis : public AnalysisSet< Analyzer<T> > {
		 public:
		 typedef Analyzer<T> system_t;

		 histogram_analysis (
		 std::string description, std::string filename,
		 const U min, const U max, const U resolution)
		 : 
		 AnalysisSet<system_t> (description, filename),
		 histogram (min, max, resolution) { }

		 virtual ~histogram_analysis() { }

		 protected:
		 typedef histogram_utilities::Histogram1D<U>	histogram_t;
		 histogram_t histogram;
		 };



		 template <typename T>
		 class double_histogram_analysis : public histogram_analysis<T> {
		 public:
		 typedef typename histogram_analysis<T>::system_t system_t;

		 double_histogram_analysis (
		 analysisstd::string desc, std::string fn,
		 const double min, const double max, 
		 const double min_2, const double max_2, const int number_of_bins)
		 : 
		 histogram_analysis<T> (desc, fn, min, max, (max-min)/double(number_of_bins)),
		 histogram_2 (min_2, max_2, (max_2-min_2)/double(number_of_bins)) { }

		 virtual ~double_histogram_analysis() { }
		 virtual void DataOutput (system_t&);

		 protected:
		 typedef typename histogram_analysis<T>::histogram_t histogram_t;
		 histogram_t histogram_2;

		 };



		 template <typename T>
		 void double_histogram_analysis<T>::DataOutput(system_t& t) {

		 rewind(t.Output());

		 for (double val_1 = this->histogram.Min(), val_2 = this->histogram_2.Min(); val_1 < this->histogram.Max(); 
		 val_1 += this->histogram.Resolution(), val_2 += this->histogram_2.Resolution()) {
		 fprintf (t.Output(), "% 8.3f % 8.3f % 8.3f % 8.3f\n", val_1, double(this->histogram.Population(val_1))/double(t.Timestep()), val_2, double(this->histogram_2.Population(val_2))/double(t.Timestep()));
		 }

		 fflush(t.Output());
		 }



		 template <typename T>
	// an analyzer for generating a single histogram (of doubles) of a system property
	class histogram_analyzer : public AnalysisSet< Analyzer<T> > {
	public:
	typedef Analyzer<T> system_t;

	histogram_analyzer (std::string desc, std::string fn) : AnalysisSet<system_t> (desc, fn) { }
	virtual ~histogram_analyzer() { }

	virtual void DataOutput (system_t& t);
protected:
typedef std::vector<double> double_vec;
double_vec values;
};	// single histogram analyzer


template <typename T>
// an analyzer for generating histograms (of doubles) of two system properties
class double_histogram_analyzer : public histogram_analyzer<T> {
	public:
		typedef typename histogram_analyzer<T>::system_t system_t;

		double_histogram_analyzer (std::string desc, std::string fn) : histogram_analyzer<T> (desc,fn) { }
		virtual ~double_histogram_analyzer() { }
		void DataOutput (system_t& t);
	protected:
		typedef typename histogram_analyzer<T>::double_vec double_vec;
		double_vec second_values;
};	// double histogram analyzer




template <typename T>
void histogram_analyzer<T>::DataOutput (system_t& t) {

	rewind(t.Output());

	histogram_utilities::histogram_t first_histo = histogram_utilities::Histogram (values.begin(), values.end(), 100);
	int first_histo_max = histogram_utilities::MaxPopulation (first_histo.begin(), first_histo.end());

	for (unsigned i = 0; i < first_histo.size(); i++) {
		fprintf (t.Output(), "% 8.4f % 8.4f\n",
				first_histo[i].first, (double)first_histo[i].second/first_histo_max);
	}
	return;
}


template <typename T>
void double_histogram_analyzer<T>::DataOutput (system_t& t) {

	rewind(t.Output());

	histogram_utilities::histogram_t first_histo = histogram_utilities::Histogram (this->values.begin(), this->values.end(), 200);
	int first_histo_max = histogram_utilities::MaxPopulation (first_histo.begin(), first_histo.end());

	histogram_utilities::histogram_t second_histo = histogram_utilities::Histogram (this->second_values.begin(), this->second_values.end(), 200);
	int second_histo_max = histogram_utilities::MaxPopulation (second_histo.begin(), second_histo.end());

	for (unsigned i = 0; i < first_histo.size(); i++) {
		fprintf (t.Output(), "% 8.4f % 8.4f % 8.4f % 8.4f\n", 
				first_histo[i].first, (double)first_histo[i].second/first_histo_max, 
				second_histo[i].first, (double)second_histo[i].second/second_histo_max);
	}
	return;
}
*/

}	// namespace md_analysis

#endif
