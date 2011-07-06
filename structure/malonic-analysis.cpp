#include "malonic-analysis.h"

namespace malonic_analysis {

	void MalonicTest::Analysis () {
		this->LoadAll();

		int numalk = alkane::Alkane::numAlkanes; 
		int nummal = alkane::MalonicAcid::numMalonicAcid;
		int numate = alkane::MalonicAcid::numMalonate;
		int numdiate = alkane::MalonicAcid::numDimalonate;

		//std::cout << "num alkanes: " << numalk << std::endl;
		//std::cout << "num malonic: " << nummal << std::endl;
		//std::cout << "num malonate: " << numate << std::endl;
		//std::cout << "num dimalonate: " << numdiate << std::endl;

		//if (numate) {
			//alkane::MalonicAcid * mal = static_cast<alkane::MalonicAcid *>( Molecule::FindByType (this->begin_mols(), this->end_mols(), Molecule::MALONATE));

			//mal->Print();

		//}
	}

}
