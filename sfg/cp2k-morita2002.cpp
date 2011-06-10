#include "cp2k-morita2002.h"

namespace morita {

	// for the small cp2k systems, just use all the waters
	void CP2KMorita2008Analysis::SelectAnalysisWaters () { 

		//std::cout << "pre-water selection" << std::endl;
		this->analysis_wats.clear();
		std::sort(this->all_wats.begin(), this->all_wats.end(), system_t::molecule_position_pred(Atom::O));
		std::copy(this->all_wats.rbegin(), this->all_wats.rbegin()+10, std::back_inserter(this->analysis_wats));

		/*
		VecR z_ax = VecR::UnitZ();
		for (Morita_it it = this->all_wats.begin(); it != this->all_wats.end(); it++) {
				if (((*it)->OH1() < z_ax) > cos(35.0*M_PI/180.0))
					this->analysis_wats.push_back (*it);
				if (((*it)->OH2() < z_ax) > cos(25.0*M_PI/180.0))
					this->analysis_wats.push_back (*it);
		}
		*/
		//std::cout << std::endl << this->analysis_wats.size() << std::endl;
		//std::cout << "post-water selection" << std::endl;


		/*
		bondgraph::BondGraph graph;
		Atom_ptr_vec atoms;
		for (Morita_it it = this->analysis_wats.begin(); it != this->analysis_wats.end(); it++) {
			std::copy((*it)->begin(), (*it)->end(), std::back_inserter(atoms));
		}
		graph.UpdateGraph(atoms.begin(), atoms.end());
		bondgraph::coordination coord = bondgraph::OH;
		bondgraph::WaterCoordination_pred p (coord, graph);

		this->analysis_wats.clear();
		for (Morita_it it = this->all_wats.begin(); it != this->all_wats.end(); it++) {
			//bondgraph::coordination crd = graph.WaterCoordination(*it);
			if (p(*it))
				this->analysis_wats.push_back(*it);
		}
		*/

		/*
		printf ("\n%zu --> ", this->analysis_wats.size());
		Morita_ptr_vec::iterator pend = std::remove_if(this->analysis_wats.begin(), this->analysis_wats.end(), std::not1(p));
		this->analysis_wats.erase(pend, this->analysis_wats.end());
		*/

		/*
		Morita_ptr_vec new_wats;
		std::copy (this->analysis_wats.begin(), this->analysis_wats.end(), std::back_inserter(new_wats));
		this->analysis_wats.clear();
		std::copy (new_wats.begin(), new_wats.end(), std::back_inserter(this->analysis_wats));
		*/
		//printf ("%zu\n", this->analysis_wats.size()); fflush(stdout);


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
