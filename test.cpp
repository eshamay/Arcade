#include "crdfile.h"

int main () {

	md_files::CRDFile crd ("mdcrd", 2703);

	for (int i = 0; i < 10; i++) {
		crd(2702).Print();
		crd.LoadNext();
	}
	
	return 0;
}
