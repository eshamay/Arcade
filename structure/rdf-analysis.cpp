#include "rdf-analysis.h"

namespace md_analysis {

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
			//for (Mol_it mol2 = mol+1; mol2 != this->end_wats(); mol2++) {

				wat = static_cast<WaterPtr>(*mol);
				//wat2 = static_cast<WaterPtr>(*mol2);

				// grab the distances from the so2 sulfur to all the water oxygens
				distance = MDSystem::Distance (so2->O1(), wat->H1()).Magnitude();
				histo(distance);
				distance = MDSystem::Distance (so2->O1(), wat->H2()).Magnitude();
				histo(distance);
				distance = MDSystem::Distance (so2->O2(), wat->H1()).Magnitude();
				histo(distance);
				distance = MDSystem::Distance (so2->O2(), wat->H2()).Magnitude();
				histo(distance);
			}//}

	}

	// see http://www.physics.emory.edu/~weeks/idl/gofr2b.html for how this is done here
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
			//fprintf (this->output, "%12.3f\n", n / dV / N);
		}

		fflush(this->output);
		return;
	}

	void RDFAgent::OutputData () {
		FILE * fout = fopen (filename.c_str(), "w");

		std::vector<double> output = CalcRDF();
		double pos;
		for (int i = 0; i < output.size(); i++) {
			pos = histo.min() + i * histo.res();
			fprintf (fout, "%.3f %.4f\n", pos, output[i]);
		}

		fclose(fout);
	}

	std::vector<double> RDFAgent::CalcRDF () {
		double dV, n, N, g_r;
		double total_volume = 4.0/3.0 * M_PI * pow(histo.max(), 3);
		std::vector<double> output;

		// go through each value in the RDF position range
		for (double r = histo.min(); r < histo.max(); r += histo.res()) {
			// print the position to the first column
			//fprintf (this->output, "%12.3f ", r);
			dV = 4.0 * M_PI * pow(r, 2) * histo.res();	// volume of the shell being considered

			n = histo.Population(r);
			N = histo.Count();

			g_r = n * total_volume / dV / N;
			output.push_back(g_r);
		}

		return output;
	}


	void RDFByDistanceAnalyzer::SuccinicAcidCalculation(alkane::SuccinicAcid * succ) {
		this->com = succ->UpdateCenterOfMass() [WaterSystem::axis];
		this->position = this->h2os.TopOrBottom(com);

		RDFAgent * rdf = FindRDFAgent (position.second);
		//succ->SetDihedralAtoms();

		double distance;
		// the carbonyl oxygens
		AtomPtr o1 = succ->GetAtom ("O2");
		AtomPtr o2 = succ->GetAtom ("O3");

		for (Wat_it wat = h2os.begin(); wat != h2os.end(); wat++) {
			distance = MDSystem::Distance (o1, (*wat)->O()).norm();
			rdf->operator()(distance);
			distance = MDSystem::Distance (o2, (*wat)->O()).norm();
			rdf->operator()(distance);
		}

		return;
	}

	RDFAgent* RDFByDistanceAnalyzer::FindRDFAgent (const double pos) {
		// given the max/min/res parameters for chopping up the system/surface into slices
		// this will return the corrent histogram to use
		RDFAgent* rdf;

		if (pos <= posmin)
			rdf = &rdfs.front();

		else if (pos >= posmax)
			rdf = &rdfs.back();

		else {
			int r = int((pos-posmin)/posres);
			rdf = &rdfs[r];
		}

		return rdf;
	}

} // namespace
