#include "test.h"

using namespace md_system;

int main () {

	md_system::WaterSystem<md_files::AmberSystem> sys ("system.cfg");
	//md_files::AmberSystem sys ("prmtop", "mdcrd");

	sys.LoadAll();
	for (int i = 0; i < 10; i++) {
		AtomPtr a = sys.int_atoms[5];
		a->Print();
		sys.LoadNext();
	}
	
	return 0;
}
