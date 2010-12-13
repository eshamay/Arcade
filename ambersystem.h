#pragma once
#ifndef AMBERSYSTEM_H_
#define AMBERSYSTEM_H_

#include "mdsystem.h"
#include "crdfile.h"
//#include "forcefile.h"
#include "topfile.h"
//#include "graph.h"

namespace md_files {

	class AmberSystem : public md_system::MDSystem {

		protected:
			//AmberSystem () { }
			//Atom_ptr_vec	_atoms;
			TOPFile		_topfile;
			CRDFile		_coords;
			//ForceFile	_forces;
			bool _periodic;


			void _ParseAtomInformation ();
			void _ParseMolecules ();

		public:
			// constructors
			AmberSystem (const std::string& prmtop, const std::string& mdcrd, const bool periodic=true);//, const std::string& mdvel = "");
			virtual ~AmberSystem ();

			// Controller & Calculation methods
			void LoadNext ();	 					// Update the system to the next timestep
			void LoadFirst ();

			bool eof () const { return _coords.eof(); }

			// Output
			const VecR&	Dimensions () 		const 	{ return _coords.Dimensions(); }		// returns the system size.

			//int 	Current ()		const 	{ return _coords.Current(); }

			//void PrintCRDFile () const;						// to output a frame of the system in .crd format

			//bondgraph::BondGraph graph;

	};

}
#endif
