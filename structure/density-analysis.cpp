#include "density-analysis.h"

namespace density {

	void MolecularDensityDistribution::Setup () {

		AnalysisSet::Setup();

		h2os.Reload(); 

		// gather all the different moltypes in the system
		std::set<Molecule::Molecule_t> moltypes;
		std::transform (this->begin_mols(), this->end_mols(), std::inserter(moltypes,moltypes.end()), std::mem_fun<Molecule::Molecule_t>(&Molecule::MolType));

		// create a histogram for each molecule type
		histos.clear();
		for (std::set<Molecule::Molecule_t>::iterator it = moltypes.begin(); it != moltypes.end(); it++) {
			histos.insert (std::pair<Molecule::Molecule_t,histo_t> (
						*it, 
						histo_t (min, max, res)));
		}

		// set up all the waters
		for (Mol_it mol = this->begin_mols(); mol != this->end_mols(); mol++) {
			if ((*mol)->MolType() == Molecule::H2O)
				(*mol)->SetAtoms();
		}
	}


	void MolecularDensityDistribution::Analysis () {
		h2os.FindWaterSurfaceLocation();
		h2o_analysis::surface_distance_t surface_distance;
		double com;

		for (Mol_it mol = this->begin_mols(); mol != this->end_mols(); mol++) {
			com = (*mol)->UpdateCenterOfMass()[WaterSystem::axis];
			surface_distance = h2os.TopOrBottom(com);
			//if ((*mol)->MolType() == Molecule::DIACID) {
				//printf ("%f %f %f %f\n", com, h2os.TopSurfaceLocation(), h2os.BottomSurfaceLocation(), surface_distance.second);
			//}
			histo_map::iterator it = histos.find((*mol)->MolType());
			it->second(surface_distance.second);
		}
	}

	void MolecularDensityDistribution::DataOutput () {
		rewind (this->output);

		// output a header line with all the molecule types (comma delimited)
		for (histo_map::iterator it = histos.begin(); it != histos.end(); it++) {
			std::string type = Molecule::Moltype2String(it->first);
			fprintf (this->output, "%s ", type.c_str());
		}
		fprintf (this->output, "\n");


		// output comma delimited rows containing the position 
		// followed by histogram populations for each molecule type
		for (double pos = min; pos < max; pos += res) {
			fprintf (this->output, "%f ", pos);	// the position is the first entry in each row
			// each row then is followed by the populations in the histograms
			for (histo_map::iterator it = histos.begin(); it != histos.end(); it++) {
				fprintf (this->output, "%f ", it->second.Population(pos));
			}
			fprintf (this->output, "\n");
		}

	}

} // namespace density
