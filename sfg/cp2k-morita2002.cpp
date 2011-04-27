#include "cp2k-morita2002.h"

namespace morita {

	// for the small cp2k systems, just use all the waters
	void CP2KMorita2008Analysis::SelectAnalysisWaters () { return; }

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
