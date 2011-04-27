#pragma once
#ifndef CP2K_MORITA2008_H_
#define CP2K_MORITA2008_H_

#include "morita2002.h"

/*!
 * A derived Morita2008 analysis routine. This performs the SFG spectral analysis (J. Phys. Chem. B 2008, 106, 673-685) on an MD data set produced by the CP2K (Quickstep) ab initio (DFT) molecular dynamics suite (http://cp2k.berlios.de/). This analysis assumes that the input MD data has been calculated along with the corresponding wannier localization centers in order to reproduce molecular and system dipole moments.
 */
namespace morita {
	using namespace md_system;
	using namespace md_analysis;

  class CP2KMorita2008Analysis : public Morita2008LookupAnalysis<XYZSystem> {
	public:

		typedef Analyzer<XYZSystem> system_t;
		CP2KMorita2008Analysis (system_t * t) :
			Morita2008LookupAnalysis<XYZSystem>(t) { }

	protected:
		//! Selects the waters to be analyzed by loading the entire set - doesn't discriminate on location because there are so few water molecules
		void SelectAnalysisWaters ();
		//! sets the dipole moments of all the system waters by means of the wannier localization centers calculated by CP2K
		//void SetAnalysisWaterDipoleMoments ();
		//! Calculate the polarizability of the given water molecule
		//void SetAnalysisWaterPolarizability ();
	};	// class CP2K sfg analyzer

} // namespace morita

#endif

