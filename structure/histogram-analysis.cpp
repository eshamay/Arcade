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


}
