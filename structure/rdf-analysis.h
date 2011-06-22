#ifndef __RDF_ANALYSIS_H
#define __RDF_ANALYSIS_H

#include "../analysis.h"

namespace md_analysis {


		class RDFAnalyzer : public AnalysisSet {
			public:
			typedef Analyzer system_t;

				RDFAnalyzer (system_t * t) :
					AnalysisSet (t,
						std::string("RDF Analysis"),
						std::string("rdf.wat-H.wat-H.dat")),
					histo(0.5, 15.0, 0.05) { }

				void Analysis ();
				void DataOutput ();

			protected:

				histogram_utilities::Histogram1D<double> histo;
				
		};

	void RDFAnalyzer::Analysis () {

		this->LoadAll();

		// find the so2
		MolPtr mol = Molecule::FindByType(this->begin_mols(), this->end_mols(), Molecule::SO2);
		SulfurDioxide * so2 = static_cast<SulfurDioxide *>(mol);
		so2->SetAtoms();

		// find the waters
		this->LoadWaters();

		WaterPtr wat, wat2;
		double distance;
		for (Mol_it mol = this->begin_wats(); mol != this->end_wats(); mol++) {
		for (Mol_it mol2 = mol+1; mol2 != this->end_wats(); mol2++) {

			wat = static_cast<WaterPtr>(*mol);
			wat2 = static_cast<WaterPtr>(*mol2);

			// grab the distances from the so2 sulfur to all the water oxygens
			distance = MDSystem::Distance (wat->H1(), wat2->H1()).Magnitude();
			histo(distance);
			distance = MDSystem::Distance (wat->H2(), wat2->H1()).Magnitude();
			histo(distance);
			distance = MDSystem::Distance (wat->H1(), wat2->H2()).Magnitude();
			histo(distance);
			distance = MDSystem::Distance (wat->H2(), wat2->H2()).Magnitude();
			histo(distance);
		}}

	}

		void RDFAnalyzer::DataOutput () {
			rewind (this->output);

			//double minimum = WaterSystem::SystemParameterLookup("analysis.rdf.minimum");
			//double maximum = WaterSystem::SystemParameterLookup("analysis.rdf.maximum");
			//double resolution = WaterSystem::SystemParameterLookup("analysis.rdf.resolution");

			double min = histo.Min();
			double max = histo.Max();
			double res = histo.Resolution();

			double dV, n, N;
			double total_volume = 4.0/3.0 * M_PI * pow(max, 3);

			// go through each value in the RDF position range
			for (double r = min; r < max; r += res) {
				// print the position to the first column
				fprintf (this->output, "%12.3f ", r);
				dV = 4.0 * M_PI * pow(r, 2) * res;	// volume of the shell being considered

				n = histo.Population(r);
				N = histo.Count();
				fprintf (this->output, "%12.3f\n", n * total_volume / dV / N);
			}

			fflush(this->output);
			return;
		}

} // namespace md_analysis
#endif
