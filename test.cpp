#include "test.h"

using namespace md_system;

int main () {

	analysis::Analyzer<AmberSystem> sys;
	
	sys.LoadAll();
	for (int i = 0; i < 10; i++) {
		MolPtr a = sys.int_mols[0];
		a->Print();
		sys.LoadNext();
	}
	
	return 0;
}
