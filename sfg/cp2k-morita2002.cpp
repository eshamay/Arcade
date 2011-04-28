#include "cp2k-morita2002.h"

namespace morita {

	// for the small cp2k systems, just use all the waters
	void CP2KMorita2008Analysis::SelectAnalysisWaters () { 
		this->analysis_wats.clear();
		//std::sort (this->all_wats.begin(), this->all_wats.end(), typename system_t::molecule_position_pred(Atom::O));
		// testing a 2-water system
		//std::copy(this->all_wats.begin(), this->all_wats.begin()+2, back_inserter(this->analysis_wats));
		this->analysis_wats.push_back(this->all_wats[0]);
		return;
	}

	/*
		 void CP2KMorita2008Analysis::SetAnalysisWaterDipoleMoments () {
		 this->CalcWannierDipoles ();
		 return;
		 }

		 void CP2KMorita2008Analysis::SetAnalysisWaterPolarizability () {
		 this->MoritaH2OPolarizabilities();
		 return;
		 }
		 */

} // namespace morita




// used to calculate the SFG spectrum based on the morita/hynes 2008 method
int main (int argc, char **argv) {

	md_analysis::Analyzer<md_system::XYZSystem> analyzer;
	//morita::CP2KMorita2008Analysis analysis ("cp2k-morita2008.dat");
	morita::CP2KMorita2008Analysis analysis (&analyzer);
	analyzer.SystemAnalysis(analysis);

	return 0;
}
