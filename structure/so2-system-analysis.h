#ifndef SO2_ANALYSIS_H_
#define SO2_ANALYSIS_H_

#include "h2o-analysis.h"

namespace so2_analysis {

	using namespace md_analysis;
	typedef SulfurDioxide*								so2_ptr;
	typedef std::vector<SulfurDioxide *>	so2_vec;
	typedef so2_vec::iterator							so2_it;
	typedef so2_vec::reverse_iterator			so2_rit;


	// a convenience class for working with systems comprised of at least 1 SO2 molecule and a whole bunch of waters
	template <typename T>
		class SO2SystemManipulator : public SystemManipulator<T> {
			public:
				typedef Analyzer<T> system_t;

				SO2SystemManipulator (system_t * t) : SystemManipulator<T>(t) { }

				virtual ~SO2SystemManipulator () { 
					delete this->so2;
					for (std::vector<SulfurDioxide *>::iterator it = so2s.begin(); it != so2s.end(); it++) {
						delete *it;
					}
				}

				virtual void Initialize () {
					this->_system->LoadAll();

					this->FindSO2 ();

					this->so2->SetAtoms();

					//printf ("\nFound %zu Sulfur Dioxides in the system\n", so2s.size());
					this->FindAllSO2s();
				}

				// The method by which the SO2 of interest is found in the system
				virtual void FindSO2 () {
					// find the so2 in the system and set some pointers up
					int id = WaterSystem<T>::SystemParameterLookup("analysis.reference-molecule-id");
					MolPtr mol = this->_system->sys_mols[id];
					//MolPtr mol = Molecule::FindByType(this->_system->sys_mols, Molecule::SO2);
					this->so2 = new SulfurDioxide(mol);
				}

				virtual void UpdateSO2 () {
					delete this->so2;
					for (std::vector<SulfurDioxide *>::iterator it = so2s.begin(); it != so2s.end(); it++) {
						delete *it;
					}
					this->Initialize();
				}
					
				virtual void FindAllSO2s () {
					so2s.clear();
					// load all the system so2s into the container - in case
					for (Mol_it it = this->_system->sys_mols.begin(); it != this->_system->sys_mols.end(); it++) {
						if ((*it)->Name() == "so2") {
							SulfurDioxide * mol = new SulfurDioxide(*it);
							mol->SetAtoms();
							so2s.push_back(mol);
						}
					}
				}

				SulfurDioxide * SO2 () { return so2; }
				AtomPtr S () { return so2->S(); }
				AtomPtr O1 () { return so2->O1(); }
				AtomPtr O2 () { return so2->O2(); }


				so2_it begin () { return so2s.begin(); }
				so2_it end ()		{ return so2s.end(); }
				so2_rit rbegin () { return so2s.rbegin(); }
				so2_rit rend () { return so2s.rend(); }

			protected:
				SulfurDioxide * so2;	// the sulfur dioxide of interest
				std::vector<SulfurDioxide *> so2s;
		}; // so2 system manipulator


	template <typename T>
	class XYZSO2Manipulator : public SO2SystemManipulator<T> {

		public:
			typedef Analyzer<T> system_t;

			XYZSO2Manipulator (system_t * t) : SO2SystemManipulator<T>(t) { }

				virtual void Initialize () {
					this->_system->LoadAll();

					this->FindSO2 ();
				}

			virtual void FindSO2 () {
				for (Mol_it it = this->_system->sys_mols.begin(); it != this->_system->sys_mols.end(); it++) {
					if ((*it)->MolType() == Molecule::SO2) {
						this->so2 = new SulfurDioxide(*it);
						this->so2->SetAtoms();
					}
				}
			}

			void UpdateSO2 () {
				delete this->so2;
				this->Initialize();
			}

			// get the wannier centers associated with a particular atom
			vector_map_vec GetWanniers (const AtomPtr& atom, const int num) const {
				std::sort (this->so2->wanniers_begin(), this->so2->wanniers_end(), vecr_distance_cmp(atom->Position()));
				vector_map_vec ret;
				std::copy (this->so2->wanniers_begin(), this->so2->wanniers_begin()+num, std::back_inserter(ret));
				return ret;
			}

	};	// xyz so2 manipulator



	template <typename T>
		class SO2PositionRecorder : public AnalysisSet<T> {

			public:
				typedef Analyzer<T> system_t;
				SO2PositionRecorder (system_t * t) : 
					AnalysisSet<T> (t, 
							std::string ("Record position of so2 relative to the surface, and the surface location"),
							std::string ("so2-position.dat")),
					so2s(t), h2os(t) { }

				void Analysis ();

			private:
				SO2SystemManipulator<T>	so2s;
				h2o_analysis::H2OSystemManipulator<T>	h2os;
		};

	// write out the position of the so2, the position of the surface, and the difference between the two
	template <typename T>
		void SO2PositionRecorder<T>::Analysis () {
			h2os.FindWaterSurfaceLocation();
			double so2_pos, surface, distance;
			so2_pos = system_t::Position(so2s.S());
			surface = h2os.SurfaceLocation();
			if (h2os.TopSurface())
				distance = so2_pos - surface;
			else
				distance = surface - so2_pos;

			fprintf (this->output, "% 12.7f % 12.7f % 12.7f\n", so2_pos, surface, distance);
		}



	// functor takes a water and returns the values of the cos(angles) formed between the two oh-vectors. The first value of the pair is always the greater (magnitude) of the two values.
	class SOAngleCalculator : public std::unary_function <SulfurDioxide*,std::pair<double,double> > {
		private:
			VecR axis;	// the reference axis to which the angles will be formed
		public:
			SOAngleCalculator (const VecR ax) : axis(ax) { }
			std::pair<double,double> operator() (const SulfurDioxide* so2) {
				double angle1 = so2->SO1() < axis;
				double angle2 = so2->SO2() < axis;
				std::pair<double,double> p = (fabs(angle1) > fabs(angle2)) 
					? std::make_pair(angle1,angle2) 
					: std::make_pair(angle2,angle1);
				return p;
			}
	};


	template <typename T>
		class SO2DipoleAnalyzer : public AnalysisSet<T> {

			public:
				typedef Analyzer<T> system_t;
				SO2DipoleAnalyzer (system_t * t) :
					AnalysisSet<T> (t,
							std::string ("SO2 Dipole Analysis"),
							std::string ("so2-dipole.dat")),
					so2s(t) { }

				void Analysis ();
				void UpdateDipole (AtomPtr atom, const int num, VecR& dipole, const VecR& ref);

			protected:
				XYZSO2Manipulator<T>	so2s;
		};

	template <typename T>
		void SO2DipoleAnalyzer<T>::Analysis () {

			so2s.UpdateSO2();
			so2s.SO2()->UpdateCenterOfMass();
			VecR com = so2s.SO2()->CenterOfMass();

			VecR dipole;
			UpdateDipole (so2s.S(), 1, dipole, com);
			UpdateDipole (so2s.O1(), 4, dipole, com);
			UpdateDipole (so2s.O2(), 4, dipole, com);

			fprintf (this->output, "% 12.8f % 12.8f % 12.8f\n", dipole[x], dipole[y], dipole[z]);
		}

	template <typename T>
		void SO2DipoleAnalyzer<T>::UpdateDipole (AtomPtr atom, const int num, VecR& dipole, const VecR& ref) {
			vector_map_vec wans = so2s.GetWanniers (atom,num);
			for (vector_map_it it = wans.begin(); it != wans.end(); it++) {
				VecR r (MDSystem::Distance(ref, (*it)));
				r *= (-2.0);
				dipole += r;
			}
			dipole += MDSystem::Distance(ref, atom->Position()) * atom->Charge();
		}

}	// namespace so2_analysis


#endif
