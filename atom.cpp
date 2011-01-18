#include "atom.h"

namespace md_system {

	Atom::Atom (const std::string& name, const double * position, const double * force) :
		_name(name), _residue(""), 
		_ID(-1), _molid(-1),
		_pmolecule((Molecule *)NULL),
		_position(Eigen::Map<VecR> (position)), _force(Eigen::Map<VecR> (force))
	{ this->SetAtomProperties (); }

	// copy constructor will pull over all the member values, and copy over the vector values
	Atom::Atom (const Atom& oldAtom) :
		_name(oldAtom._name), _residue(oldAtom._residue), _ID(oldAtom._ID), _molid(oldAtom._molid),
		_pmolecule(oldAtom._pmolecule),
		_mass(oldAtom._mass), _charge(oldAtom._charge),
		_position(oldAtom._position), _force(oldAtom._force),
		_element(oldAtom._element)	
	{ this->SetAtomProperties (); }


	Atom::~Atom () { }

	void Atom::Position (const vector_base& position) { 
		_position = position;
	}

	void Atom::Position (const double X, const double Y, const double Z) { 
		_position.Set(X, Y, Z); 
	}

	void Atom::Position (coord const axis, double const value) { 
		_position.Set (axis, value); 
	}

	double Atom::operator[] (const coord index) const {
		double pos = -1.0;

		switch (index) {
			case x: pos = _position.x();
			case y: pos = _position.y();
			case z: pos = _position.z();
		}
		return pos;
	}

	double Atom::operator- (const Atom& input) const {
		return (_position - input._position).norm();
	}

	bool Atom::operator< (const AtomPtr& rhs) const {
		return this->_ID < rhs->ID();
	}
	bool Atom::operator< (const Atom& rhs) const {
		return this->_ID < rhs.ID();
	}

	void Atom::Print () const {
		printf ("%5s (ID:%6d nuc:%4d)\t%4s (molID:%5d)\t% f\t% f\t% f\n", _name.c_str(), _ID, _element, _residue.c_str(), _molid, _position.x(), _position.y(), _position.z());
	}

	void Atom::SetAtomProperties () {

		_mass = 0.0;
		_charge = 0.0;
		_element = NO_ELEMENT;

		//const double AMU2KG	= 1.6762158e-27;		// conversion for amu -> kg
		if (_name.find("DW") != std::string::npos) {
			_mass = 0.0;
			_charge = 0.0;
			_element = DW;
		}
		else if (_name.find("SW") != std::string::npos) {
			_mass = 0.0;
			_charge = 0.0;
			_element = SW;
		}
		else if (_name.find("Cl") != std::string::npos) {
			_mass = 35.453;
			_charge = 7.0;
			_element = Cl;
		}
		else if (_name.find("H") != std::string::npos) {
			_mass = 1.0078250321;
			_charge = 1.0;
			_element = H;
		}
		else if (_name.find("O") != std::string::npos) {
			_mass = 15.9949146221;
			_charge = 6.0;
			_element = O;
		}
		else if (_name.find("N") != std::string::npos) {
			_mass = 14.0030740052;
			_charge = 5.0;
			_element = N;
		}
		else if (_name.find("D") != std::string::npos) {
			_mass = 2.0156500641;
			_charge = 1.0;
			_element = D;
		}
		else if (_name.find("C") != std::string::npos) {
			_mass = 12.0000000;
			_charge = 4.0;
			_element = C;
		}
		else if (_name.find("Si") != std::string::npos) {
			_mass = 28.0855;
			_charge = 4.0;
			_element = Si;
		}
		else if (_name.find("S") != std::string::npos) {
			_mass = 32.065;
			_charge = 6.0;
			_element = S;
		}
		else if (_name.find("F") != std::string::npos) {
			_mass = 18.9984;
			_charge = 7.0;
			_element = F;
		}
		return;
	}

	Atom::Element_t Atom::String2Element (const std::string& name) {
		const double * pos=(double *)NULL, * frc=(double *)NULL;
		AtomPtr atom = new Atom(name,pos,frc);
		Atom::Element_t elmt = atom->Element();
		delete atom;
		return elmt;
	}

	std::string Atom::Element2String (Atom::Element_t elmt) {

		std::string name;
		switch (elmt) {
			case(0) :	name = "NO_ELEMENT"; break;
			case(1) :	name = "H"; break;
			case(2) :	name = "He"; break;
			case(5) :	name = "B"; break;
			case(6) :	name = "C"; break;
			case(7) :	name = "N"; break;
			case(8) :	name = "O"; break;
			case(9) :	name = "F"; break;
			case(10) :	name = "Ne"; break;
			case(11) :	name = "Na"; break;
			case(12) :	name = "Mg"; break;
			case(13) :	name = "Al"; break;
			case(14) :	name = "Si"; break;
			case(15) :	name = "P"; break;
			case(16) :	name = "S"; break;
			case(17) :	name = "Cl"; break;
			case(18) :	name = "Ar"; break;
			case(19) :	name = "L"; break;
			case(20) :	name = "Ca"; break;
			case(53) :	name = "I"; break;
			case(55) :	name = "Cs"; break;
			case(56) :	name = "Ba"; break;
			default :	name = "NO_ELEMENT_DEFAULT"; break;
		}

		return name;
	}


	void Atom::KeepByElement (Atom_ptr_vec& u, const Element_t& elmt) {
		u.erase(
				remove_if(u.begin(), u.end(), std::not1(member_functional::ptr_mem_fun_eq<Atom,Atom::Element_t>(&Atom::Element, elmt))), u.end());
		return;
	}

	bool Atom::ElementCombo (const AtomPtr& ai, const AtomPtr& aj, const Element_t element_a, const Element_t element_b) {
		return 
			((ai->Element() == element_a && aj->Element() == element_b)
			 ||
			 (ai->Element() == element_b && aj->Element() == element_a));
	}


	bool Atom::element_eq (const AtomPtr& first, const AtomPtr& second) {
		return first->Element() == second->Element();
	}


	bool Atom::id_cmp (const AtomPtr& first, const AtomPtr& second) {
		return first->ID() < second->ID();
	}

	bool Atom::id_eq (const AtomPtr& first, const AtomPtr& second) {
		return first->ID() == second->ID();
	}

	AtomPtr Atom::FindByElement (const Atom_ptr_vec& apv, Element_t elmt) {
		Atom_it a = std::find_if (apv.begin(), apv.end(), member_functional::mem_fun_eq(&Atom::Element,elmt));
		return *a;
	}

	AtomPtr Atom::FindByID (const Atom_ptr_vec& apv, int id) {
		Atom_it a = std::find_if (apv.begin(), apv.end(), member_functional::mem_fun_eq(&Atom::ID,id));
		return *a;
	}

}	// namespace md system
