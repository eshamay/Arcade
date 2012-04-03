#ifndef	WATERSYSTEM_H_
#define	WATERSYSTEM_H_

#include "mdsystem.h"
#include "ambersystem.h"
#include "xyzsystem.h"
//#include "gmxsystem.h"

#include "utility.h"

#include <libconfig.h++>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <numeric>

namespace md_system { 

	using namespace md_files;

	class WaterSystem {

		public:

			WaterSystem (const std::string configuration_filename);
			virtual ~WaterSystem ();

			static libconfig::Config * config_file;	/* Configuration file */

			static libconfig::Setting& SystemParameterLookup (std::string param) {
				try {
					return WaterSystem::config_file->lookup(param);
				}
				catch(const libconfig::SettingTypeException &stex) {
					std::cerr << "Something is wrong with the " << param << " setting in the system.cfg configuration file." << std::endl;
					exit(EXIT_FAILURE);
				}
				catch(const libconfig::SettingNotFoundException &snfex) {
					std::cerr << "Couldn't find the " << param << " setting in the system.cfg configuration file." << std::endl;
					exit(EXIT_FAILURE);
				}
			}


			static double	posmin, posmax;
			static double	pbcflip;			// location to flip about periodic boundaries
			static coord axis;				// axis normal to the interface
			static VecR ref_axis;			// vector representation of the reference axis
			static double int_low, int_high, middle;	// the positions of analysis cutoffs

			static Atom_ptr_vec	sys_atoms;		// all atoms/mols in the system
			static Mol_ptr_vec	sys_mols;

			static Atom_ptr_vec	int_atoms;		// interfacial water atoms (or as above)
			static Mol_ptr_vec 	int_mols;

			static Mol_ptr_vec	int_wats;		// interfacial waters, or just all the waters in the system depending on the function call

			wannier_it begin_wanniers() const;
			wannier_it end_wanniers() const;

			void OpenFile ();


			static double AxisPosition (const AtomPtr a) {
				double pos = a->Position()[axis];
				pos = (pos > pbcflip) ? pos : pos + MDSystem::Dimensions()[axis];
				return pos;
			}

			typedef std::pair<double,double> Double_pair;
			// quick way to make a pair for the oft-used extents std::pair defaulting to the posmin/posmax in the config file
			Double_pair ExtentPair (
					const double low = WaterSystem::posmin,
					const double high = WaterSystem::posmax) const {
				return std::make_pair<double,double> (low, high);
			}

			void LoadAll ();							// Loads all the molecules and atoms in the system into the containers
			//void SliceWaterCoordination (const bondgraph::coordination c);

			// predicate determines if an atom sits within a particular slice of the system
			class AtomPositionInSlice : public std::binary_function<AtomPtr, Double_pair&, bool> {
				public:
					bool operator() (const AtomPtr atom, const Double_pair& extents) const
					{
						double pos = AxisPosition (atom);
						return pos > extents.first && pos < extents.second;
					}
			};

			// remove all the atoms in the defined slice (extents=(min,max)) of the slab
			void RemoveAtomsInSlice (Atom_ptr_vec& atoms, Double_pair& extents) {
				atoms.erase(
						remove_if(atoms.begin(), atoms.end(), std::bind2nd(AtomPositionInSlice(), extents)), atoms.end());
				return;
			}

			void KeepAtomsInSlice (Atom_ptr_vec& atoms, Double_pair& extents) {
				atoms.erase(
						remove_if(atoms.begin(), atoms.end(), std::not1(std::bind2nd(AtomPositionInSlice(), extents))), atoms.end());
				return;
			}

			// predicate to find if a water molecule is within a region (based on the position of the oxygen atom)
			template <typename U>
				class WaterInSlice : public std::binary_function<U, Double_pair, bool> {
					private:
						AtomPositionInSlice apis;
					public:
						bool operator() (const U wat, const Double_pair& extents) const
						{
							return apis (wat->GetAtom("O"), extents);
						}
				};

			// This slices the water molecules and leaves only those within a given location of the slab. However, it does not do any loading or unloading of the atoms from int_atoms or anything else...
			template <typename U>
				static void SliceWaters (Mol_ptr_vec& mols, Double_pair& extents) {
					mols.erase(
							remove_if(mols.begin(), mols.end(), std::not1(std::bind2nd(WaterInSlice<U>(), extents))), mols.end());

					//this->UpdateAtoms(int_wats, int_atoms);
					return;
				}


			/* loads the int_wats and int_atoms with only waters and water atoms */
			void LoadWaters () {
				LoadAll();
				int_wats.clear();

				// copy over all the water molecules into the int_wats container
				algorithm_extra::copy_if (sys_mols.begin(), sys_mols.end(), std::back_inserter(int_wats), member_functional::mem_fun_eq(&Molecule::MolType, Molecule::H2O));

				// load in the water atoms to int_atoms
				this->UpdateAtoms (int_wats, int_atoms);

				return;
			}


			/*
			// Predicate to test if a water molecule has a given coordination (H-bonding pattern)
			class WaterCoordination_p : public std::binary_function<Water *, bondgraph::coordination, bool> {
			public:
			bool operator() (const Water * wat, const bondgraph::coordination c) const {
			return sys->graph.WaterCoordination(wat) == c;
			}
			};


			void KeepWaterByCoordination (Mol_ptr_vec& mols, const bondgraph::coordination c) {
			mols.erase(
			remove_if(mols.begin(), mols.end(), std::not1(std::bind2nd(WaterCoordination_p(), c))), mols.end());
			return;
			}
			*/



			// Sets the atom container to hold only the atoms of the set of molecules given
			void UpdateAtoms (const Mol_ptr_vec& mols, Atom_ptr_vec& atoms) {
				atoms.clear();
				for (Mol_it it = mols.begin(); it != mols.end(); it++) {
					std::copy ((*it)->begin(), (*it)->end(), std::back_inserter(atoms));
				}
				return;
			}

			//void UpdateGraph () { sys->graph.UpdateGraph (int_atoms); }
			//void UpdateGraph (const Atom_ptr_vec& atoms) { sys->graph.UpdateGraph (atoms); }
			//bondgraph::BondGraph& Graph () const { return sys->graph; }

			virtual void Initialize () = 0;
			void LoadNext() const { sys->LoadNext(); }
			virtual void Rewind() const { sys->Rewind(); }

		protected:
			MDSystem * sys;	/* System coordinate & files */

	};	// class watersystem



	class AmberWaterSystem : public WaterSystem {
		public:
			AmberWaterSystem (const std::string configuration_filename = std::string ("system.cfg")) : 
				WaterSystem (configuration_filename) { }
			~AmberWaterSystem () { delete this->sys; }

			virtual void Initialize () {
				try {
					std::string prmtop = this->SystemParameterLookup("system.files.prmtop");
					std::string mdcrd = this->SystemParameterLookup("system.files.mdcrd");
					bool periodic = this->SystemParameterLookup("system.periodic");
					this->sys = new AmberSystem(prmtop, mdcrd, periodic);
					printf ("\n\tSystem Files::\n\t\tprmtop = %s\n\t\tmdcrd = %s\n", prmtop.c_str(), mdcrd.c_str());
				}
				catch (const libconfig::SettingNotFoundException &snfex) {
					std::cerr << "Couldn't find the Amber system filenames listed in the configuration file" << std::endl;
					exit(EXIT_FAILURE);
				}
				return;
			}
	};	// amber water system


	class XYZWaterSystem : public WaterSystem {
		public:
			XYZWaterSystem (const std::string configuration_filename = std::string ("system.cfg")) : 
				WaterSystem (configuration_filename) { }
			~XYZWaterSystem () { delete this->sys; }

			//wannier_it begin_wanniers() const { return this->sys->begin_wanniers(); }
			//wannier_it end_wanniers() const { return this->sys->end_wanniers(); }

			virtual void Initialize () {

				try {
					std::string filepath = this->SystemParameterLookup("system.files.xyzfile");

					double a,b,c;
					a = this->SystemParameterLookup("system.dimensions")[0];
					b = this->SystemParameterLookup("system.dimensions")[1];
					c = this->SystemParameterLookup("system.dimensions")[2];

					std::string wanniers = SystemParameterLookup("system.files.wanniers");
					VecR dims(a,b,c);
					printf ("system dimensions are: ");
					dims.Print();

					std::cout << "new xyz system" << std::endl;
					this->sys = new XYZSystem(filepath, dims, wanniers);
				}
				catch (const libconfig::SettingNotFoundException &snfex) {
					std::cerr << "Couldn't find the xyz system parameters in the configuration file" << std::endl;
					exit(EXIT_FAILURE);
				}
				//catch (std::exception &e) {
					//std::cerr << "something else went wrong while starting the xyzsystem" << std::endl;
				//}


				return;
			}

	};	// xyz water system


}	// namespace

#endif
