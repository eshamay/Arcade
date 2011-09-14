#include "histogram-analysis.h"

namespace md_analysis { 

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

	void Histogram2DAgent::OutputData (const DataOutput2DFunction& func) {

		output = fopen(filename.c_str(), "w");
		if (!output) {
			std::cerr << "couldn't open the file for output --> \" " << filename << " \"" << std::endl;
			exit(1);
		}
		rewind(output);

		for (double alpha = histogram.min.first; alpha < histogram.max.first; alpha += histogram.resolution.first) {
			for (double beta = histogram.min.second; beta < histogram.max.second; beta += histogram.resolution.second) {
				fprintf (output, "% 12.3f % 12.3f % 12f\n", alpha, beta, func(alpha, beta,histogram.Population(alpha, beta)));
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




	///////////// Multi 2D Histogram Agent //////////////
	Multi2DHistogram::Multi2DHistogram (
			const double minimum1, const double maximum1, const double resolution1
			const double minimum2, const double maximum2, const double resolution2
			const double minimum3, const double maximum3, const double resolution3,
			std::string prefix, std::string suffix) :

		min (minimum1), max(maximum1), res(resolution1) // extents of the analysis and thickness of the slices
		{
			// set up all the histograms
			histos.clear();
			histos.resize (
					int((max-min)/res),
					Histogram2DAgent (std::string (""), 
						minimum2, maximum2, resolution2,
						minimum3, maximum3, resolution3));

			// set the name for each of the histograms
			double pos;
			for (int i = 0; i < histos.size(); i++) {
				pos = res * i + min;
				std::stringstream sstr;
				sstr.clear();
				std::string filenum;
				filenum.clear();
				sstr << pos;
				filenum = sstr.str();
				//std::string filepath (std::string("./alcohol-oxygen-water-hydrogen.distance-rdfs/rdf.") + filenum + ".dat");
				std::string filepath (prefix + filenum + suffix);
				histos[i].SetOutputFilename (filepath);
			}
		}


	Histogram2DAgent* Multi2DHistogramAgent::FindHistogram (const double val) {
		// given the max/min/res parameters for chopping up the system/surface into slices
		// this will return the corrent histogram to use
		Histogram2DAgent* histo;

		if (val <= min)
			histo = &histos.front();

		else if (val >= max)
			histo = &histos.back();

		else {
			int r = int((val-min)/res);
			histo = &histos[r];
		}

		return histo;
	}

	void Multi2DHistogramAgent::operator() (const double val1, const double val2, const double val3) {
		Histogram2DAgent* histo = FindHistogram(val1);
		histo->operator() (val2,val3);
		return;
	}


	void Multi2DHistogramAgent::DataOutput (DataOutput2DFunction func) {
		for (std::vector<Histogram2DAgent>::iterator hist = histos.begin(); hist != histos.end(); hist++) {
			hist->OutputData(func);
		}
	}

}
