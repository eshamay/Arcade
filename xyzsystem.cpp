#include "xyzsystem.h"

namespace md_files {

	using namespace md_system;

	XYZSystem::XYZSystem (const std::string& filepath, const VecR& size, const std::string& wannierpath) :
		_xyzfile(filepath),
		_wanniers(wannierpath),
		_reparse_limit(1),	// initially set to parse everything everytime
		_reparse_step(0)
	{
		MDSystem::Dimensions (size);
		//this->LoadNext();
	}


	XYZSystem::~XYZSystem () {
		for (Mol_it it = _mols.begin(); it != _mols.end(); it++)
			delete *it;
	}



	void XYZSystem::_ParseMolecules () {

		/***********************************************************************************
		 * This is the top-level parsing routine to give the overall idea of what's going on
		 * *********************************************************************************/
		// first things first - we need the interatomic distances and bonding information - atomic bonding graph
		try { graph.UpdateGraph (_xyzfile.Atoms()); }

		catch (bondgraph::BondGraph::graphex& ex) {
			std::cout << "Caught an exception while updating the bond graph" << std::endl;
		}

		this->_InitializeSystemAtoms();

		this->_FindMolecules();

		/*
			 this->_ParseSimpleMolecule<Hydronium> (Atom::O, Atom::H, 3);
			 this->_ParseSimpleMolecule<Water> (Atom::O, Atom::H, 2);
		// parse NO3- ions
		this->_ParseSimpleMolecule<Nitrate> (Atom::N, Atom::O, 3);
		// turn any NO3- ions that have a covalently bound H into HNO3 - nitric acid
		this->_ParseNitricAcids ();
		this->_ParseSimpleMolecule<SulfurDioxide> (Atom::S, Atom::O, 2);
		this->_ParseSimpleMolecule<Hydroxide> (Atom::O, Atom::H, 1);
		this->_ParseProtons ();
		this->_ParseAlkanes ();
		*/

		/*
			 try {
			 this->_CheckForUnparsedAtoms ();
			 }
			 catch (xyzsysex& ex) {
			 std::cerr << "Exception thrown while checking for unparsed atoms in the system - some atoms were not parsed into molecules" << std::endl;
			 throw;
			 }
			 */

		return;
	}	// parse molecules



	void XYZSystem::FindMoleculesByMoleculeGraph () {
		// track which atoms have already been parsed in the system
		_unparsed.clear();
		std::copy (_xyzfile.begin(), _xyzfile.end(), std::back_inserter(_unparsed));

		// now do the work of parsing out the molecules. The idea is to grab sub-graphs of the atomic bondgraph that was built that represent covalently bound atom groups.
		// For each subgraph we tease out the type of molecule it is, then we set the molid and add it to the list of molecules.
		std::vector<AtomPtr> parsed;
		while (!_unparsed.empty()) {
			// build the molecule graph of whatever is at the front of the unparsed list
			molgraph::MoleculeGraph molgraph;
			molgraph.Initialize(_unparsed.front(), graph);

			// then generate the molecule from the molecule graph
			MolPtr newmol = molgraph::MoleculeGraph2Molecule (molgraph);

			parsed = molgraph.Atoms();
			this->_UpdateUnparsedList (parsed);

			newmol->MolID((int)_mols.size());
			newmol->FixAtoms();
			_mols.push_back(newmol);
			//newmol->SetAtoms();
		}
	}	// find molecules by molecule graph



	void TopologyXYZSystem::_FindMolecules () {

		if (!_parsed) {

			int id = 0;
			for (MolecularTopologyFile::mol_topology_it top = _topology.begin(); top != _topology.end(); top++) {
				MolPtr mol = MoleculeFactory (top->name);
				mol->MolID(id++);

				for (std::vector<int>::iterator atom_id = top->atoms.begin(); atom_id != top->atoms.end(); atom_id++) {
					mol->AddAtom(_xyzfile[*atom_id]);
				}

				this->_mols.push_back(mol);
			}	// for each molecular topology
			_parsed = true;
		}	// if not parsed

	}	// topology xyz find molecules

	/*
		 void XYZSystem::_ParseNitricAcids () { }

		 for (Mol_it mol = _mols.begin(); mol != _mols.end(); mol++) {
		 if ((*mol)->MolType() != Molecule::NO3) continue;

		 AtomPtr N = (*mol)->GetAtom(Atom::N);
// find any Hs bound to oxygens of the NO3
Atom_ptr_vec Os = graph.BondedAtoms(N, bondgraph::covalent, Atom::O);
Atom_ptr_vec Hs; 
for (Atom_it O = Os.begin(); O != Os.end(); O++) {
Atom_ptr_vec H_check = graph.BondedAtoms(*O, bondgraph::covalent, Atom::H);
if (H_check.size() == 1)
Hs.push_back(H_check[0]);
if (H_check.size() > 1) {
printf ("A nitrate Oxygen has %zu covalently bound hydrogens! what's going on here?\n", H_check.size());
(*O)->Print();
exit(1);
}
}

if (Hs.size() == 1)
(*mol)->AddAtom (Hs[0]);
else if (Hs.size() > 1) {
printf ("This nitric acid has %zu covalently bound hydrogens!\n", Hs.size());
(*mol)->Print();
for (Atom_it ai = Hs.begin(); ai != Hs.end(); ai++) {
(*ai)->Print();
}
}
}
return;
}	// Parse Nitric acids
*/


void XYZSystem::_UpdateUnparsedList (Atom_ptr_vec& parsed) {
	// fix up the list for tracking atoms that have already been added into molecules.
	// _unparsed should contain only atoms that are not in the parsed list
	// these atom vectors get sorted according the the atom's ID

	_unparsed.erase(
			std::remove_if (_unparsed.begin(), _unparsed.end(), std::bind2nd(AtomPtr_In_List<Atom_ptr_vec>(), parsed)), _unparsed.end());
	return;
} 

void XYZSystem::_CheckForUnparsedAtoms () const {

	if (!_unparsed.empty()) {

		std::cout << "The following atoms were found unparsed into molecules after all molecules had been formed" << std::endl;

		// print out every atom that hasn't been parsed
		for (Atom_it it = _unparsed.begin(); it != _unparsed.end(); it++) {
			std::cout << std::endl;
			(*it)->Print();	// show all remaining atoms

			// and all the atoms to which it is bound
			Atom_ptr_vec bound (graph.BondedAtoms(*it));
			for (Atom_it jt = bound.begin(); jt != bound.end(); jt++) { 
				// and the distance between them
				std::cout << "^--~ (" << graph.Distance(*jt, *it) << ")  ";
				(*jt)->Print();
			}
		}
		throw (unaccountedex());
	}

	return;
}	// check for unparsed atoms



/******************************
 * Add in the wannier centers
 *
 * What we'll do here is run through each of the wannier centers and find which molecule they belong to. For each center we'll check its distance against all the atoms and when we find the one it's bound to we'll shove it into the parent molecule.
 * ****************************/
void XYZSystem::_ParseWanniers () {
	// we've got to clear out all the wanniers already loaded into the molecules
	std::for_each (_mols.begin(), _mols.end(), std::mem_fun(&Molecule::ClearWanniers));
	//std::for_each (_mols.begin(), _mols.end(), std::mem_fun(&Molecule::SetAtoms));

	int num;
	std::map<Molecule::Molecule_t, int>::iterator mapend = WannierFile::numWanniers.end();
	std::map<Molecule::Molecule_t, int>::iterator it;
	for (Mol_it mol = _mols.begin(); mol != _mols.end(); mol++) {

		it = WannierFile::numWanniers.find((*mol)->MolType());

		if (it == mapend) {
			printf ("Couldn't find the number of wanniers for the molecule: %s. Add it into wanniers.cpp.\n", (*mol)->Name().c_str());
			(*mol)->Print();
			fflush (stdout);
			exit(1);
		}

		num = it->second;

		AddWanniers (*mol, num);
	}

	/*
	// keep track of the wanniers that have been parsed
	for (Atom_it at = _xyzfile.begin(); at != _xyzfile.end(); at++) {
// find every oxygen and sulfur atom (those are the ones that contain wanniers
if ((*at)->Element() != Atom::O && (*at)->Element() != Atom::S) continue;
//(*at)->Print();

MolPtr mol = (*at)->ParentMolecule();

// sort the wanniers by distance to the given atom
//std::sort (_wanniers.begin(), _wanniers.end(), vecr_distance_cmp((*at)->Position()));

// add the appropriate number into the molecules
int num_wans = 0;
if ((*at)->Element() == Atom::O) num_wans = 4;
if ((*at)->Element() == Atom::S) num_wans = 1;

double distance;
int i = 0;
vector_map_vec parsed_ws;
while (num_wans) {

distance = MDSystem::Distance ((*at)->Position(), _wanniers[i]).Magnitude();
if (distance < WANNIER_BOND) {
mol->AddWannier(_wanniers[i]);
parsed_ws.push_back (_wanniers[i]);
--num_wans;
}
++i;
}
*/

/*
// and then go through all the wanniers in order to find the ones attached to the oxygen
for (WannierFile::Wannier_it wi = _wanniers.begin(); wi != _wanniers.end(); wi++) {

// we find the distance to the oxygen from the wannier center
double distance = MDSystem::Distance ((*at)->Position(), *wi).Magnitude();

// if it's close enough, then that oxygen gets it
if (distance < WANNIER_BOND) {
mol->AddWannier(*wi);
parsed_ws.push_back (*wi);
}
}
}	// for-each atom
*/

//printf ("%d parsed wanniers\n", (int)parsed_ws.size());
// now find which wanniers have not been parsed
//std::cout << " found " << _wanniers.size() - parsed_ws.size() << " wanniers unparsed" << std::endl;

return;
}	// Parse Wanniers

void XYZSystem::AddWanniers (MolPtr mol, const int num) {
	mol->SetAtoms();
	mol->UpdateCenterOfMass();
	// sort all the wannier centers by their distance to the given atom location
	std::sort (_wanniers.begin(), _wanniers.end(), vecr_distance_cmp(mol->ReferencePoint()));

	// add in the given num of wannier centers to the molecule
	for (vector_map_it it = _wanniers.begin(); it != _wanniers.begin() + num; it++)
		mol->AddWannier(*it);
}

void XYZSystem::Rewind () {
	_xyzfile.Rewind();
	this->LoadNext();
}	// rewind


void XYZSystem::LoadNext () {
	_xyzfile.LoadNext();

	if (_wanniers.Loaded()) {
		_wanniers.LoadNext();
	}
	//try {
	//if (++_reparse_step == _reparse_limit) {
	this->_ParseMolecules();
	if (_wanniers.Loaded())
		this->_ParseWanniers();
	//_reparse_step = 0;
	//}
	//} catch (xyzsysex& ex) {
	//std::cout << "Exception caught while parsing the molecules of the XYZ system" << std::endl;
	//throw;
	//}
} // Load Next




// This will calculate the total dipole moment of the system based on atom locations and wannier center positions
// The origin is shifted to the center of the system in order to get closest images (wrapped into the box) of all the atoms/wanniers
VecR XYZSystem::SystemDipole () {

	VecR dipole;
	dipole.Set(0.0,0.0,0.0);

	for (Atom_it it = _xyzfile.begin(); it != _xyzfile.end(); it++) {
		VecR ri = (*it)->Position();
		dipole += ri * (*it)->Charge();
	}

	for (WannierFile::Wannier_it it = _wanniers.begin(); it != _wanniers.end(); it++) {
		dipole -= (*it) * 2.0;
	}

	return (dipole);
}

/*
	 void XYZSystem::_ParseAlkanes () {

	 for (Atom_it it = _xyzfile.begin(); it != _xyzfile.end(); it++) {
// form all the remaining atoms, find a carbon atom that it likely part of an alkane
if ((*it)->Element() != Atom::C) continue;
// check that it's still not parsed
if (!_Unparsed(*it)) continue;

int molIndex = (int)_mols.size();	// set the molecule index
alkane::Alkane * alk = new alkane::Alkane ();
alk->MolID (molIndex);
alk->Initialize (*it, graph);	// build the alkane molecule
_mols.push_back(alk);

// remove all the atoms in the alkane from the unparsed list
Atom_ptr_vec parsed;
std::copy (alk->begin(), alk->end(), std::back_inserter(parsed));
this->_UpdateUnparsedList (parsed);
}
}	// parse alkane
*/

bool XYZSystem::_Unparsed (const AtomPtr atom) const {
	bool ret;
	ret = (std::find (_unparsed.begin(), _unparsed.end(), atom) != _unparsed.end());
	return ret;
}

/*
	 void XYZSystem::_ParseProtons () {

	 Atom_ptr_vec parsed;

	 for (Atom_it H = _unparsed.begin(); H != _unparsed.end(); H++) {
// form all the remaining un-bound hydrogens into their own "proton" molecules
if ((*H)->Element() != Atom::H) continue;

Atom_ptr_vec cov = graph.BondedAtoms(*H, bondgraph::covalent);
if (cov.empty()) {

int molIndex = (int)_mols.size();	// set the molecule index
_mols.push_back (new Proton ());

MolPtr newmol = _mols[molIndex];
newmol->MolID (molIndex);
newmol->AddAtom (*H);

parsed.push_back (*H);
}
// otherwise, if there are covalent bonds, AND the hydrogen is not parsed into a molecule... something is wrong
else if (!(*H)->ParentMolecule()) {
printf ("Found a hydrogen that is not parsed into a molecule, but has covalent bonds!\n");
(*H)->Print();
printf ("It is covalently bound to:\n");
for (Atom_it it = cov.begin(); it != cov.end(); it++) {
(*it)->Print();
}
printf ("And it is h-bonded to:\n");
Atom_ptr_vec hbonds = graph.BondedAtoms(*H, bondgraph::hbond);
for (Atom_it it = hbonds.begin(); it != hbonds.end(); it++) {
(*it)->Print();
}

exit(1);
}
}

_UpdateUnparsedList (parsed);
return;
}
*/

} // namespace md_files

