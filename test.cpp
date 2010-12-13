#include "ambersystem.h"

int main () {

	md_files::AmberSystem sys ("prmtop", "mdcrd");

	for (int i = 0; i < 10; i++) {
		sys.Molecules(5)->Print();
		sys.LoadNext();
	}
	
	return 0;
}
