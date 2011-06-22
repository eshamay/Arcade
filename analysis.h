#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "watersystem.h"
#include "utility.h"
#include "patterns.h"
#include "dataoutput.h"


namespace md_analysis {

	using namespace md_system;
	using namespace md_files;


	/*********** ANALYZER **************/


	class Analyzer : public patterns::observer::observable {

		protected:

			WaterSystem * sys;

			int	output_freq;

			void _OutputHeader () const;
			//md_analysis::StarStatusBarUpdater	status_updater;
			md_analysis::PercentProgressBar	status_updater;

		public:
			Analyzer (WaterSystem * water_sys);
			virtual ~Analyzer ();

			// position boundaries and bin widths for gathering histogram data
			static double	posres;
			static int		posbins;
			static double angmin, angmax, angres;
			static int		angbins;
			static int		timestep;
			int						Timestep () const { return timestep; }
			static int 		timesteps;
			static int		restart;

			static double Position (const MolPtr);
			static double Position (const AtomPtr);
			static double Position (const VecR&);
			static double Position (const double);

			void LoadAll () { sys->LoadAll(); }
			void LoadNext ();
			void Rewind() { 
				this->sys->Rewind();
				timestep = 1;
			}

			void LoadWaters () { sys->LoadWaters(); }

			void OutputStatus ();
			bool ReadyToOutputData () const { 
				return Analyzer::timestep % (Analyzer::output_freq * 10) || Analyzer::timestep;
			}	

			//Atom_ptr_vec& Atoms () { return WaterSystem::int_atoms; } 
			//Mol_ptr_vec& Molecules () { return WaterSystem::int_mols; }
			//Mol_ptr_vec& Waters () { return WaterSystem::int_wats; }

			// calculate the system's center of mass
			template <typename Iter> static VecR CenterOfMass (Iter first, Iter last);

			//! Predicate for sorting a container of molecules based on position along the main axis of the system, and using a specific element type to determine molecular position. i.e. sort a container of waters based on the O position, or sort a container of NO3s based on the N position, etc.
			class molecule_position_pred; 		
			class atomic_distance_cmp;
			class molecule_distance_cmp;
			class molecule_distance_generator;

			class MoleculeAbovePosition;
			class MoleculeBelowPosition;

	};	// Analyzer



	template <class Iter>	// Has to be iterators to a container of molecules
		VecR Analyzer::CenterOfMass (Iter first, Iter last)
		{
			double mass = 0.0;
			VecR com;
			com.setZero();

			typedef typename std::iterator_traits<Iter>::value_type val_t;

			for (Iter it = first; it != last; it++) {
				for (Atom_it jt = (*it)->begin(); jt != (*it)->end(); jt++) {
					mass += (*jt)->Mass();
					com += (*jt)->Position() * (*jt)->Mass();
				}
			}
			com /= mass;
			return com;
		}

	//! Predicate for sorting molecules based on their positions along the system reference axis. The position of the element supplied (elmt) is used. e.g. if elmt = Atom::O, then the first oxygen of the molecule will be used
	class Analyzer::molecule_position_pred : public std::binary_function <Molecule*,Molecule*,bool> {
		private:
			Atom::Element_t _elmt;	// determines the element in a molecule to use for position comparison
		public:
			//! upon instantiation, the element to be used for specifying molecular position is provided
			molecule_position_pred (const Atom::Element_t elmt) : _elmt(elmt) { }

			bool operator()(const Molecule* left, const Molecule* right) const {
				AtomPtr left_o = left->GetAtom(_elmt);
				AtomPtr right_o = right->GetAtom(_elmt);
				double left_pos = Analyzer::Position(left_o);
				double right_pos = Analyzer::Position(right_o);

				return left_pos < right_pos;
			}
	};



	// this predicate is used for distance calculations/sorting between atoms given a reference atom or position
	class Analyzer::atomic_distance_cmp : public std::binary_function <AtomPtr,AtomPtr,bool> {
		private:
			VecR _v;	// the molecule that will act as the reference point for the comparison
		public:
			atomic_distance_cmp (const AtomPtr refatom) : _v (refatom->Position()) { }
			atomic_distance_cmp (const VecR v) : _v (v) { }
			// return the distance between the two molecules and the reference mol
			bool operator()(const AtomPtr left, const AtomPtr right) const {
				double left_dist = MDSystem::Distance(left->Position(), _v).norm();
				double right_dist = MDSystem::Distance(right->Position(), _v).norm();
				return left_dist < right_dist;
			}
	};

	// given a reference point, this returns a molecule's distance to that point
	class Analyzer::molecule_distance_generator : public std::unary_function <MolPtr,double> {
		private:
			VecR _v;	// the molecule that will act as the reference point for the comparison

		public:
			molecule_distance_generator (const MolPtr refmol) : _v(refmol->ReferencePoint()) { }
			molecule_distance_generator (const AtomPtr refatom) : _v(refatom->Position()) { }
			molecule_distance_generator (const VecR v) : _v(v) { }

			double operator()(const MolPtr mol) const {
				return MDSystem::Distance(mol->ReferencePoint(), _v).norm();
			}
	}; // molecule distance generator

	class Analyzer::molecule_distance_cmp : public std::binary_function <MolPtr,MolPtr,bool> {
		private:
			VecR _v;	// the molecule that will act as the reference point for the comparison
		public:
			molecule_distance_cmp (const MolPtr refmol) : _v(refmol->ReferencePoint()) { }
			molecule_distance_cmp (const AtomPtr refatom) : _v(refatom->Position()) { }
			molecule_distance_cmp (const VecR v) : _v(v) { }

			// return the distance between the two molecules and the reference
			bool operator()(const MolPtr left, const MolPtr right) const {
				/*
					 left->Print();
					 left->SetAtoms();
					 left->ReferencePoint().Print();
					 _v.Print();
					 */
				double left_dist = MDSystem::Distance(left->ReferencePoint(), _v).norm();
				double right_dist = MDSystem::Distance(right->ReferencePoint(), _v).norm();
				return left_dist < right_dist;
			}
	};

	// predicate tells if a molecule's reference point along a given axis is above a given value
	class Analyzer::MoleculeAbovePosition : public std::unary_function <MolPtr,bool> {
		private:
			double position;
			coord axis;
		public:
			MoleculeAbovePosition (const double pos, const coord ax) : position(pos), axis(ax) { }
			bool operator() (const MolPtr mol) {
				return Analyzer::Position(mol->ReferencePoint()) > position;
			}
	};

	// predicate tells if a molecule's reference point along a given axis is above a given value
	class Analyzer::MoleculeBelowPosition : public std::unary_function <MolPtr,bool> {
		private:
			double position;
			coord axis;
		public:
			MoleculeBelowPosition (const double pos, const coord ax) : position(pos), axis(ax) { }
			bool operator() (const MolPtr mol) {
				return Analyzer::Position(mol->ReferencePoint()) < position;
			}
	};




	// An analysis that will be performed on a system by an analyzer
	class AnalysisSet {
		public:

			virtual ~AnalysisSet () {
				if (output != (FILE *)NULL)
					fclose(output);
			}

			AnalysisSet (Analyzer * sys, std::string desc, std::string fn) 
				: 
					_system(sys),
					description (desc), filename(fn),
					output((FILE *)NULL) { }

			// default setup
			virtual void Setup () {
				OpenDataOutputFile ();
				_system->LoadAll();
				return;
			}

			// each analyzer has to have an analysis function to do some number crunching
			virtual void Analysis () = 0;
			// normally this can be done in the analysis section, but just for style we can have something different defined here
			virtual void DataOutput () { }
			virtual void PostAnalysis () { }

			void OpenDataOutputFile ();

			std::string& Description () { return description; }
			std::string& Filename () { return filename; }

			void LoadAll () const { this->_system->LoadAll(); }
			void LoadWaters () const { this->_system->LoadWaters(); }

			Atom_it_non_const begin () const { return WaterSystem::int_atoms.begin(); }
			Atom_it_non_const end () const { return WaterSystem::int_atoms.end(); }
			Atom_ptr_vec& Atoms () { return WaterSystem::int_atoms; }

			Mol_it begin_mols () const { return WaterSystem::int_mols.begin(); }
			Mol_it end_mols () const { return WaterSystem::int_mols.end(); }
			Mol_ptr_vec& Mols () { return WaterSystem::int_mols; }

			Mol_it begin_wats () const { return WaterSystem::int_wats.begin(); }
			Mol_it end_wats () const { return WaterSystem::int_wats.end(); }
			Mol_ptr_vec& IntWats () { return WaterSystem::int_wats; }

			//wannier_it begin_wanniers () const { return this->_system->begin_wanniers(); }
			//wannier_it end_wanniers () const { return this->_system->end_wanniers(); }

		protected:

			Analyzer * _system;
			std::string description;	// describes the analysis that is performed
			std::string filename;		// filename to use for data output
			FILE * output;

	};  // class AnalysisSet



}
#endif
