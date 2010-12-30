#ifndef __ANGLE_BOND_ANALYSIS_H_
#define __ANGLE_BOND_ANALYSIS_H_

#include "histogram-analysis.h"
#include "h2o-analysis.h"

namespace md_analysis {


	template <typename T>
	class H2OAngleBondAnalysis : public AnalysisSet<T> {
		public:
			typedef Analyzer<T>	system_t;

			H2OAngleBondAnalysis (system_t * t) :
				AnalysisSet<T> (t,
						std::string ("H2O H-O-H angle and O-H bondlength histograms"),
						std::string ("")),
				h2os(t),
				// outputting to two separate files.
				bonds("bondlength.h2o.O-H.dat", 
						t->posmin, t->posmax, t->posres, 
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.bondlength-min"), 
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.bondlength-max"), 
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.bondlength-res")),

				angles("angle.h2o.H-O-H.dat", 
						t->posmin, t->posmax, t->posres, 
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.angle-min"), 
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.angle-max"), 
						WaterSystem<T>::SystemParameterLookup("analysis.angle-bond-histogram.angle-res")) { }

			void Analysis ();

			virtual void DataOutput () { 
				bonds.OutputData();
				angles.OutputData();
			}

		protected:
			h2o_analysis::H2OSystemManipulator<T>	h2os;
			Histogram2DAgent	bonds;
			Histogram2DAgent	angles;

	};	 // h2o angle & bond histogram analyzer

	template <typename T>
		void H2OAngleBondAnalysis<T>::Analysis () {
			h2os.Reload();
			h2os.FindWaterSurfaceLocation ();

			for (Wat_it it = h2os.begin(); it != h2os.end(); it++) {
				double distance = system_t::Position((*it)->ReferencePoint()) - h2os.SurfaceLocation();

				// calculate and bin the water OH bond lengths
				bonds (distance, (*it)->OH1().norm());
				bonds (distance, (*it)->OH2().norm());

				// calculate the angle of the water molecule and bin it
				double angle = (*it)->Angle();
				angles (distance, acos(angle)*180.0/M_PI);
			}

			return;
		}


	/*

		 template <typename T>
		 class so2_angle_bond_analyzer : public so2_analysis::SO2SystemAnalyzer<T> {
		 public:
		 typedef typename so2_analysis::SO2SystemAnalyzer<T>::system_t system_t;

		 so2_angle_bond_analyzer () :
		 so2_analysis::SO2SystemAnalyzer<T> (
		 std::string("SO2 molecular angle and S-O bondlengths"),
		 std::string ("so2-angle+bonds.dat")) { }

		 virtual ~so2_angle_bond_analyzer () { }

		 void PostSetup (system_t& t) {
//AnalysisSet<system_t>::Setup(t);
// Ox-H-n == the distance from a given so2-O to the nth closest h2o-H
// S-O-n == distance from so2-S to the nth closest h2o-O
fprintf (t.Output(), "step location distance SO-1 SO-2 OSO-theta Bisector-theta S-O-1 S-O-2 S-O-3 O1-H-1 O1-H-2 O1-H-3 O2-H-1 O2-H-2 O2-H-3\n");

}

void Analysis (system_t& t);

protected:
double angle;
};




template <typename T>
void so2_angle_bond_analyzer<T>::Analysis (system_t& t) {

// calculate the OSO angle
double oso_angle = acos(this->so2->Angle())*180.0/M_PI;
VecR bisector = this->so2->Bisector();
VecR Y = Vector3d::UnitY();
// and the angle of the bisector to the system normal
double system_angle = acos(bisector < Y)*180.0/M_PI;

this->Reload ();
this->FindWaterSurfaceLocation();

// distance from the top of the water surface to the sulfur of the so2
double distance_to_location = system_t::Position(this->s) - this->surface_location;

// sort the waters by distance: h2o-O to the so2-S - the first waters in the vector will be closest to the SO2
Atom::KeepByElement(t.int_atoms, Atom::O);

std::sort(t.int_atoms.begin(), t.int_atoms.end(), 
Analyzer<T>::atomic_reference_distance_pred(this->so2->S()));

// grab the distances from so2-S to closest h2o-Os
double so_1 = MDSystem::Distance(this->so2->S(), t.int_atoms[0]).norm();
double so_2 = MDSystem::Distance(this->so2->S(), t.int_atoms[1]).norm();
double so_3 = MDSystem::Distance(this->so2->S(), t.int_atoms[2]).norm();

// sort the water atoms by distance h2o-H to so2-O1
t.LoadWaters();
Atom::KeepByElement(t.int_atoms, Atom::H);
std::sort(t.int_atoms.begin(), t.int_atoms.end(), Analyzer<T>::atomic_reference_distance_pred(this->so2->GetAtom("O1")));

// grab the distances from so2-S to closest h2o-Os
double oh1_1 = MDSystem::Distance(this->so2->O1(), t.int_atoms[0]).norm();
double oh1_2 = MDSystem::Distance(this->so2->O1(), t.int_atoms[1]).norm();
double oh1_3 = MDSystem::Distance(this->so2->O1(), t.int_atoms[2]).norm();

// sort the water atoms by distance h2o-H to so2-O1
t.LoadWaters();
Atom::KeepByElement(t.int_atoms, Atom::H);
std::sort(t.int_atoms.begin(), t.int_atoms.end(), Analyzer<T>::atomic_reference_distance_pred(this->so2->GetAtom("O2")));

// grab the distances from so2-S to closest h2o-Os
double oh2_1 = MDSystem::Distance(this->so2->O2(), t.int_atoms[0]).norm();
double oh2_2 = MDSystem::Distance(this->so2->O2(), t.int_atoms[1]).norm();
double oh2_3 = MDSystem::Distance(this->so2->O2(), t.int_atoms[2]).norm();

// find the fixed oxygen that is tethered for the steered MD
//AtomPtr fixed_o = Atom::FindByID(t.int_atoms, 1161);

// output the distance and the two S-O bondlengths and the SO2 oso_angle for each timestep
fprintf (t.Output(), "%d % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f % 12.4f\n", 
		t.Timestep(),
		this->surface_location,
		distance_to_location,
		this->so2->SO1().norm(), this->so2->SO2().norm(), 
		oso_angle, 
		system_angle,
		so_1, so_2, so_3,
		oh1_1, oh1_2, oh1_3,
		oh2_1, oh2_2, oh2_3);

return;
}


template <typename T>
class so2_angle_bond_histogram_analyzer : public double_histogram_analyzer<T> {
	public:
		typedef typename double_histogram_analyzer<T>::system_t system_t;

		virtual ~so2_angle_bond_histogram_analyzer() { }
		so2_angle_bond_histogram_analyzer() :
			double_histogram_analyzer<T> (
					std::string("SO2 O-S-O angle and S-O bondlength histograms"),
					std::string("so2-bond-angle-histograms.dat")) { }

		void Analysis (system_t& t);
	protected:
		SulfurDioxide * so2;
};






template <typename T>
void so2_angle_bond_histogram_analyzer<T>::Analysis (system_t& t) {
	t.LoadAll();

	MolPtr mol = Molecule::FindByType(t.sys_mols, Molecule::SO2);
	so2 = new SulfurDioxide(mol);
	so2->SetAtoms();

	this->values.push_back(so2->SO1().norm());
	this->values.push_back(so2->SO2().norm());
	double angle = so2->Angle();
	this->second_values.push_back(acos(angle)*180.0/M_PI);
	delete so2;

	return;
}
*/


}	// namespace md_analysis

#endif
