#ifndef MOLECULE_FACTORY_H_
#define MOLECULE_FACTORY_H_

#include "molecule.h"
#include "h2o.h"
#include "oh.h"
#include "h.h"
#include "h3o.h"
#include "hno3.h"
#include "so2.h"
//#include "ctc.h"
#include "alkane.h"

namespace md_system {

	MolPtr MoleculeFactory (const std::string name);

}	// namespace molecule

#endif
