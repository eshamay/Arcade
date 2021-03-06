#ifndef MORITA2008_H_
#define MORITA2008_H_

#include "analysis.h"
#include "sfgunits.h"
#include "moritah2o.h"
#include "data-file-parser.h"

//#include <Eigen/LU>
#include <iostream>
#include <mkl_blas.h>
#include <mkl_lapack.h>
#include <functional>


namespace morita {

	using namespace md_system;
	using namespace md_analysis;

	typedef std::pair<double,double> sfgdata_pair_t;

		/*!
		 * A pure virtual base class for performing an SFG analysis based on the method of Morita/Hynes (J. Phys. Chem. B 2008, 106, 673-685). Depending on the specific MD system used (Amber, CP2K, etc) the pure virtual methods define the actions to be taken to alter how various components are defined, and customize the analysis.
		 */
		template <class T>
		class Morita2008Analysis : public AnalysisSet<T>{
			public:
				typedef Analyzer<T> system_t;

				Morita2008Analysis (system_t * t, std::string desc, std::string filename) :
					AnalysisSet<T> (t,
							std::string (desc),
							std::string (filename)) { }
						//_p(0), _alpha(0,0), _IDENT(0,0), _g(0,0), _g_trans(0,0) { }

				virtual ~Morita2008Analysis ();

				virtual void Analysis ();

			protected:
				typedef std::vector<MatR, Eigen::aligned_allocator<MatR> > mat_vec;
				typedef std::vector< mat_vec >	tensor;
				typedef std::vector<VecR, Eigen::aligned_allocator<VecR> > vec_vec; 
				//! all the waters in the MD system
				Morita_ptr_vec		all_wats;		
				//! Water molecules that will be analyzed by the morita 2008 routines
				Morita_ptr_vec		analysis_wats;	
				//! Total system dipole moment (3Nx1 tensor)
				vec_vec	_p;
				//! The sum-total of all molecular polarizabilities of the water molecules in the system
				tensor						_alpha;
				//! System dipole field tensor
				tensor						_T;	
				//! Local field correction tensor
				tensor						_g;	
				//! Local field correction tensor
				mat_vec	_f;	

				//! total system polarizability
				MatR			_A;		
				//! total system dipole moment (at time zero)
				VecR			_M;		

				/*!
				 * Copies all the waters in the system into a separate container in order to differentiate between those that will be analyzed and the rest. Also, the waters and converted into a special Morita water class that makes it cleaner to calculate the dipole and polarizability tensors
				 */
				void SetupSystemWaters ();

				/*!
				 * Determines which of the molecules (waters) in the system are to be used for the ensuing analysis. Only waters that have been copied into the analysis_wats container will be used.
				 */
				virtual void SelectAnalysisWaters () = 0;	// routine that determines which waters will be used in the analysis

				/*!
				 * Method of calculating the various tensors (molecular dipole moments, polarizabilities, dipole field tensor, etc) used in the ensuing analysis
				 */
				void CalculateTensors();

				//! Sets the dipole moment of each water used in the analysis.
				virtual void SetAnalysisWaterDipoleMoments () = 0;

				//! Sets the polarizability matrix of each water used in the analysis (the lab-frame projection of the matrix)
				virtual void SetAnalysisWaterPolarizability () = 0;

				// dipole field tensor 'T' as used in the morita&hynes paper
				// it's a square 3Nx3N matrix, where N = number of particles
				MatR DipoleFieldTensor (const MoritaH2O_ptr wat1, const MoritaH2O_ptr wat2);

				void CalculateTotalDipole ();
				void CalculateTotalPolarizability ();
				void CalculateLocalFieldCorrection ();

				void CalcWannierDipoles ();
				void MoritaH2OPolarizabilities ();
				void MoritaH2ODipoles ();
				void Dipole_PolarizabilityLookupTables ();
				// sfgdata_pair_t SSP_SPS_Result (const int s1, const int s2, const int p, const MatR& a) const;
		};	// sfg-analyzer


	template <class U> Morita2008Analysis<U>::~Morita2008Analysis () {
		for (Morita_it it = all_wats.begin(); it != all_wats.end(); it++) {
			delete *it;
		}
	}


	template <class U>
		void Morita2008Analysis<U>::SetupSystemWaters () {

			// load all the waters into the int_wats container
			this->_system->LoadWaters();
			//std::for_each(t.int_wats.begin(), t.int_wats.end(), std::mem_fun(&Molecule::Print));

			for (Morita_it it = all_wats.begin(); it != all_wats.end(); it++) {
				delete *it;
			}
			all_wats.clear();

			// load up the all_wats with new derived Morita waters that have some extra functionality
			for (Mol_it it = this->_system->int_wats.begin(); it != this->_system->int_wats.end(); it++) {
				MoritaH2O_ptr ptr (new MoritaH2O (*it));
				ptr->SetAtoms();
				all_wats.push_back(ptr);
			}

			// analysis_wats will hold only those waters that will by analyzed
			analysis_wats.clear();
			std::copy(all_wats.begin(), all_wats.end(), back_inserter(analysis_wats));
			// here filter out the waters to use for analysis
			this->SelectAnalysisWaters ();

			return;
		}

	template <class U>
		void Morita2008Analysis<U>::Analysis () {

			this->SetupSystemWaters ();

			// zero out all the matrices/vectors
			int N = analysis_wats.size();
			_p.clear();
			_p.resize(N, VecR::Zero());
			_alpha.clear();
			_alpha.resize(N, mat_vec(N, MatR::Zero()));
			_T.clear();
			_T.resize(N, mat_vec(N, MatR::Zero()));
			_g.clear();
			_g.resize(N, mat_vec(N, MatR::Identity())); // setting up the identity matrix early for ease later
			_f.clear();
			_f.resize(N, MatR::Zero());	// setting up the identity matrix early for ease later

			//std::cout << "tensor calcs" << std::endl;
			// Sets up the p, alpha, and T tensors by calculating through each water's dipole moment, polarizability, and dipole field contributions
			this->CalculateTensors();

			//std::cout << "local field corr" << std::endl;
			// determine the local field correction to each of the molecular polarizabilities
			this->CalculateLocalFieldCorrection ();

			//std::cout << "total dip" << std::endl;
			// sum all the molecular dipoles to get the total system value
			this->CalculateTotalDipole();

			//std::cout << "calc polariz" << std::endl;
			// sum all the molecular polarizabilities to get the total system value
			this->CalculateTotalPolarizability ();

			//std::cout << "output" << std::endl;
			for (unsigned int i = 0; i < 3; i++) {
				fprintf (this->output, "% 12.7f ", _M(i));
			} 
			// output in row-major order
			for (unsigned int i = 0; i < 3; i++) {
				for (unsigned int j = 0; j < 3; j++) {
					fprintf (this->output, "% 12.7f ", _A(i,j));
				}
			}
			fprintf (this->output,"\n");
			fflush(this->output);
			//std::cout << "post output" << std::endl;

		} // Analysis


	// calculates the ssp and sps value of the autocorrelation function for a given set of pqr space-frame coordinates
	/*
		 template <class U>
		 sfgdata_pair_t Morita2008Analysis<U>::SSP_SPS_Result (const int s1, const int s2, const int p, const MatR& a) const {
	// for now, only output the SSP and SPS components of the correlation function
	double a_sps, m_sps;
	double a_ssp, m_ssp;

	a_sps = a(s1,p);
	m_sps = _M(s1);

	a_ssp = a(s1,s1);
	m_ssp = _M(p);
	//a_sps = (a(s1,p) + a(s2,p))/2.0;
	//a_ssp = (a(s1,s1) + a(s2,s2))/2.0;

	//m_sps = (_M(s1) + _M(s2))/2.0;
	//m_ssp = _M(p);

	return std::make_pair (a_ssp*m_ssp, a_sps*m_sps);
	}
	*/


	//template <class U>
	//void Morita2008Analysis<U>::DataOutput () {


	// try this - for each timestep, output the vector dipole, and the matrix polarizability of the system

	// *********** time-domain output work *********** //
	// Output will be 6 columns. 2 each for the X,Y, and Z axes as the "P" axis, and each pair of columns will have ssp, sps data
	/*
		 sfgdata_pair_t datapair;
		 for (MatR_it it = _vA.begin(); it != _vA.end(); it++) {
	//datapair = SSP_SPS_Result(1,2,0,*it);	// for the X-axis
	//fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
	//datapair = SSP_SPS_Result(0,2,1,*it);	// for the Y-axis
	//fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
	datapair = SSP_SPS_Result(1,2,0,*it);	// for the Z-axis
	fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
	datapair = SSP_SPS_Result(2,1,0,*it);	// for the Z-axis
	fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
	datapair = SSP_SPS_Result(0,2,1,*it);	// for the Z-axis
	fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
	datapair = SSP_SPS_Result(2,0,1,*it);	// for the Z-axis
	fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
	datapair = SSP_SPS_Result(0,1,2,*it);	// for the Z-axis
	fprintf (this->output, "% 12.6f % 12.6f", datapair.first, datapair.second);
	datapair = SSP_SPS_Result(1,0,2,*it);	// for the Z-axis
	fprintf (this->output, "% 12.6f % 12.6f\n", datapair.first, datapair.second);
	}
	*/

	//fflush(t.Output());
	//}


	// ********  fftw work ******** //
	/*
		 double *in;
		 fftw_complex *sps_fft, *ssp_fft;
		 fftw_plan p;

		 in = (double*) fftw_malloc(sizeof(double)*sps.size());

	// sps transform
	std::copy(sps.begin(), sps.end(), in);
	sps_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*sps.size());
	p = fftw_plan_dft_r2c_1d(sps.size(),in,sps_fft,FFTW_DESTROY_INPUT);
	fftw_execute(p);

	// ssp transform
	std::copy(ssp.begin(), ssp.end(), in);
	ssp_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*ssp.size());
	p = fftw_plan_dft_r2c_1d(ssp.size(),in,ssp_fft,FFTW_DESTROY_INPUT);
	fftw_execute(p);

	// output one column each for
	// time-domain data, real, imaginary, magnitude
	rewind(output);
	for (int i = 0; i < sps.size(); i++) {
	fprintf (output, "% 23.8f % 23.8f % 23.8f % 23.8f % 23.8f % 23.8f % 23.8f % 23.8f\n",
	sps[i], sps_fft[i][0], sps_fft[i][1], sqrt(sps_fft[i][0]*sps_fft[i][0] + sps_fft[i][1]*sps_fft[i][1]),
	ssp[i], ssp_fft[i][0], ssp_fft[i][1], sqrt(ssp_fft[i][0]*ssp_fft[i][0] + ssp_fft[i][1]*ssp_fft[i][1])
	);
	}

	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(sps_fft);
	fftw_free(ssp_fft);


	if (!(timestep == timesteps))
	this->Setup();
	*/

	//} // Data Output

template <class U>
void Morita2008Analysis<U>::CalculateTensors() {

	// Calculate the dipole moment of each water
	this->SetAnalysisWaterDipoleMoments ();

	// Calculate the polarizability of each water molecule
	this->SetAnalysisWaterPolarizability ();

	int N = analysis_wats.size();

	for (int i = 0; i < N; i++) {

		// set the value for p_not 
		_p[i] = analysis_wats[i]->Dipole() * sfg_units::ANG2BOHR;	// these are all in (a.u. charge) * bohr after conversion

		_alpha[i][i] = analysis_wats[i]->Polarizability();	// these are in atomic units (bohr)

		int j = i+1;
		while (j < N) {
			// Calculate the tensor 'T' which is formed of 3x3 matrix elements
			MatR dft (DipoleFieldTensor(analysis_wats[i], analysis_wats[j]));

			_T[i][j] = dft;
			_T[j][i] = dft;

			j++;
		}

	}

	return;
}	// Calculate Tensors


template <class U>
MatR Morita2008Analysis<U>::DipoleFieldTensor (const MoritaH2O_ptr wat1, const MoritaH2O_ptr wat2) {

	VecR r = MDSystem::Distance(wat1->GetAtom(Atom::O), wat2->GetAtom(Atom::O));	// this is in angstroms
	r = r * sfg_units::ANG2BOHR; // work in atomic units (bohr)

	double distance = r.Magnitude();	// units of bohr (a.u.)
	double ir3 = 1.0/pow(distance,3.0);	// in units of 1/a_0^3 - a.u. (inverse bohr^3)
	double ir5 = 1.0/pow(distance,5.0);	// in units of 1/a_0^3 - a.u. (inverse bohr^3)

	// calculate T as in eq. 10 of the Morita/Hynes 2008 paper
	MatR dft = Matrix3d::Identity();

	double tmp;
	for (unsigned i = 0; i < 3; i++) {
		for (unsigned j = 0; j < 3; j++) {
			tmp = (dft(i,j)*ir3) - (3.0*r[i]*r[j]*ir5);
			dft(i,j) = tmp;
		}
	}

	return dft;
} //dipole field tensor


// take care of the polarization calculations for the first timestep
template <class U>
void Morita2008Analysis<U>::CalculateTotalDipole () {

	// first calculate the corrected (for local field) _p tensor
	// This is done as with the polarizability, alpha, by means of the dgemm blas routine
	// p_corrected = trans(g) * p
	// where p_corrected is p corrected for the local field effect of the neighboring waters

	/*

		 int n = 1;
		 char transa = 'T';
		 char transb = 'N';
		 double scalea = 1.0;
		 double scaleb = 0.0;

		 dgemm (&transa, &transb, &N, &n, &N, &scalea, &_g(0,0), &N, &_p(0,0), &N, &scaleb, &_p(0,0), &N);
	// _p now holds the local-field corrected dipole moments
	*/

	int N = analysis_wats.size();
	_M.setZero();
	//for (vec_vec::iterator it = _p.begin(); it != _p.end(); it++) {
		for (int p = 0; p < N; ++p) {
		//for (int q = 0; q < N; ++q) {
		//_g_t = _g[p][q].transpose(); 
		//_M += _g_t * _p[q];
		_M += _f[p].transpose() * _p[p];
		//_M += _p[p];
		//_M += *it;
	}//}

	// perform the summation of the total system dipole moment
	/*
		 double tally;
		 VecR mu;
		 MatR g;
		 for (int p = 0; p < N; ++p) {
		 for (int q = 0; q < N; ++q) {
		 mu = _p.block(3*q,0,3,1);
		 g = _g_trans.block(3*p,3*q,3,3);
		 for (int i = 0; i < 3; i++) {
		 tally = 0.0;
		 for (int j = 0; j < 3; j++) {
		 tally += mu[j] * g[i][j];
		 } 
		 _M[i] += tally;
		 }
		 }
		 }
		 */

	//for (int i = 0; i < N; i++) {
	//VecR mu = _p.block(3*i,0,3,1);
	//_M += mu;
	//_M += _p.block(3*i,0,3,1);
	//}

	/*
		 VecR_vec dipoles;
		 std::transform (analysis_wats.begin(), analysis_wats.end(), std::back_inserter(dipoles), std::mem_fun<VecR,Molecule>(&Molecule::Dipole));

	// sum all the molecular dipoles together to get the total system dipole
	_M = std::accumulate (dipoles.begin(), dipoles.end(), VecR());

	// in case it's needed - keep the dipole for every timestep in a running list
	_vM.push_back(_M);
	*/

	}	// calculate total dipole


	template <class U>
		void Morita2008Analysis<U>::CalculateLocalFieldCorrection () {
			// following equation 23 - (1 + T*alpha) f = h.
			// first, set up _g = (1 + T*alpha)
			// then solve the equation for f with a lapack routine
			//
			// note: g = inv(1+T*alpha)

			int N = analysis_wats.size();

			//for (int p = 0; p < N; ++p) {
			//_g[p][p] = MatR::Identity();
			//}
			MatR tt, aa;
			for (int p = 0; p < N; ++p) {
				for (int q = 0; q < N; ++q) {
					aa = _alpha[q][q];
					tt = _T[p][q];
					_g[p][q] += tt * aa;
				}}//}
			/*
				 for (int p = 0; p < N; ++p) {
				 for (int q = 0; q < N; ++q) {
				 _g[p][q].Print();
				 }}
				 */

			// now perform the inverse on each piece of g
			MatR g_inv;
			for (int p = 0; p < N; ++p) {
				for (int q = 0; q < N; ++q) {
					g_inv = _g[p][q].inverse();	
					_g[p][q] = g_inv; /** problem here **/
				}}

			for (int p = 0; p < N; ++p) {
				for (int q = 0; q < N; ++q) {
					_f[p] += _g[p][q];
				}
			}

			/*
				 for (int p = 0; p < N; ++p) {
				 for (int q = 0; q < N; ++q) {
				 _g[p][q].Print();
				 }}
				 */


			// first step: g = I

			//char trans = 'N';
			//double scale = 1.0;

			//boost::timer t;
			//t.restart();
			// set g = T*alpha + I,
			// though the blas routine is really doing a g = T*alpha + g, since g already = I
			//dgemm (&trans, &trans, &N, &N, &N, &scale, &_T(0,0), &N, &_alpha(0,0), &N, &scale, &_g(0,0), &N);

			//_IDENT.setIdentity(N,N);
			//_g = _T;
			//_g *= _alpha;
			//_g += _IDENT;

			//_g = _IDENT + _T*_alpha;

			//std::cout << "blas matrix mult. " << t.elapsed() << std::endl;
			//
			// now g = T*alpha + I
			// At this point _g is still the inverse of what we're really looking for. It's painful to calculate the inverse of a tensor directly... however, the method is to use the LU decomposition to get the inverse

			// for now, f is 'h', the 3Nx3 block identity tensor
			//_h.setZero(N,3);
			//tensor::BlockIdentity(_h,3);

			/********** now solve for f in the equation g*f = h using the lapack dsgesv *********/
			/*
				 int nrhs = 3;
				 int ipiv[N];
				 for (int i = 0; i < N; i++) ipiv[i] = 0;
				 int info = 0;

				 int lwork = N*nrhs;
				 double work[lwork];
				 float swork[N*(N+nrhs)];
				 int iter;
				 _f.setZero(N,3);

			//t.restart();
			dsgesv (&N, &nrhs, &_g(0,0), &N, ipiv, &_h(0,0), &N, &_f(0,0), &N, work, swork, &iter, &info);
			if (info != 0) {
			std::cout << "DSGESV.info parameter had a value of " << info << " meaning that something went wrong!" << std::endl;
			exit(1);
			}
			//std::cout << "DSGESV (iterative) system solve:  " << t.elapsed() << std::endl;
			*/

			//_f = _g;
			/******** Solve for the inverse of _g through LU decomposition and direct application of the lapack _getri routine *******/

			/*
				 int ipiv[N+1];
				 int info = 0;
			// perform LU decomposition
			dgetrf (&N, &N, &_g(0,0), &N, ipiv, &info);


			// before continuing, find out the best working memory chunk size
			double lwork_query[10];
			int lwork = -1;	// tell dgetri to perform a query of the optimal workspace size, instead of performing the inversion
			// perform the work-size query
			dgetri (&N, &_g(0,0), &N, ipiv, lwork_query, &lwork, &info);
			lwork = (int)lwork_query[0];
			double work[lwork];



			// and do the inverse calculation to get the real g in Eq. 19 of the 2008 Morita/Ishiyama paper
			dgetri (&N, &_g(0,0), &N, ipiv, work, &lwork, &info);
			if (info != 0) {
			printf ("\nsomething was wrong with the calculation of the inverse: info = %d\n", info);
			}
			*/
			// now _g is as expected - the inverse of (1+T*alpha) - the value of Eq. 19 - 2008 Morita/Ishiyama


			//_f *= _g;
			/*
				 double scaleb = 0.0;
				 dgemm (&trans, &trans, &N, &N, &N, &scale, &_g(0,0), &N, &_f(0,0), &N, &scaleb, &_f(0,0), &N);
				 std::cout << "------" << std::endl;
				 std::cout << _g << std::endl;
				 std::cout << "------" << std::endl;
				 std::cout << "------" << std::endl;
				 std::cout << "------" << std::endl;
				 std::cout << "------" << std::endl;
				 std::cout << _f << std::endl;
				 std::cout << "------" << std::endl;
				 */
			//t.restart();
			//dgesv (&N, &nrhs, &_g(0,0), &N, ipiv, &_f(0,0), &N, &info);
			//std::cout << "DGESV system solve:  " << t.elapsed() << std::endl;

			//MatrixXd _x;
			//t.restart();
			//_g.lu().solve(_f, &_x);
			//std::cout << "Eigen LU solve:  " << t.elapsed() << std::endl;

}	// Local field correction


// calculates A for the current timestep
template <class U>
void Morita2008Analysis<U>::CalculateTotalPolarizability () {

	int N = analysis_wats.size();
	/*

	// The alpha tensor is calculated as: alpha_corrected = trans(g) * alpha * g
	// thus, matrix multiplication twice 

	char transa = 'T';
	char transb = 'N';
	double scale = 1.0;
	double scaleb = 0.0;
	*/

	// determine 'A', the summed (total) system polarizability.
	_A.setZero();
	//MatR _g_t;
	//for (int p = 0; p < N; ++p) {
	for (int q = 0; q < N; ++q) {
		//for (int r = 0; r < N; ++r) {
		// grab the piece of the local field tensor and find the inner product with alpha
		//_g_t = _g[p][q].transpose(); 
		//_A += _g_t * _alpha[q][q] * _g[q][r];
		//_A += (_f[q].transpose()) * _alpha[q][q] * _f[q];
		_A += _f[q].transpose() * _alpha[q][q] * _f[q];
		//_A += _alpha[q][q];
	}//}//}
}	// Calculate total polarizability


/*
// Using the LAPACK solver for some simple systems
extern "C" {

void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);

void pdgesv_ (int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv, double *B, int *ib, int *jb, int *descb, int *info); 

} // extern
*/

/************* helper class for morita analysis based on lookup tables *************/
template <class T>
class Morita2008LookupAnalysis : public Morita2008Analysis<T> {

	public:
		typedef Analyzer<T> system_t;
		Morita2008LookupAnalysis (system_t * t) : 
			Morita2008Analysis<T>(t,
					std::string ("SFG Analysis via M/I 2008 w/ polarizability lookup"),
					std::string ("sfg.dat")),
			pdf ("h2o_polarizability.dat") { }
		//ddf ("h2o_dipole.dat") { }

	protected:
		//! Selects the waters to be analyzed by loading the entire set of molecules above a particular cutoff location
		virtual void SelectAnalysisWaters () = 0;
		//! sets the dipole moments of all the system waters by means of the parameterized model of Morita/Hynes 2008
		virtual void SetAnalysisWaterDipoleMoments ();

		virtual void SetAnalysisWaterPolarizability ();

		//! Data file containing pre-calculated polarizabilities of different geometries of water molecule
		datafile_parsers::PolarizabilityDataFile pdf;	
		//datafile_parsers::DipoleDataFile ddf;	

};	// class lookup analysis



template <>
void Morita2008LookupAnalysis<XYZSystem>::SetAnalysisWaterDipoleMoments () {
	std::for_each (this->analysis_wats.begin(), this->analysis_wats.end(), MDSystem::CalcWannierDipole);
}

template <class T>
void Morita2008LookupAnalysis<T>::SetAnalysisWaterPolarizability () {

	MatR alpha;
	VecR mu;
	MatR dcm;
	double oh1, oh2, theta;

	// for every water to be analyzed:
	//		set the molecule up
	//		calculate the direction cosine matrix from the local to lab frame rotation
	//		find the polarizability from the lookup-table
	//		rotate the polarizability to the lab frame
	//		set the molecule's polarizability to the value
	for (Morita_it it = this->analysis_wats.begin(); it != this->analysis_wats.end(); it++) {

		//(*it)->SetAtoms();	// first get the atoms and bonds set in the water

		dcm = MatR::Zero();
		dcm = (*it)->DCMToLab();	// set up the direction cosine matrix for rotating the polarizability to the lab-frame

		// bondlengths are in angstroms
		oh1 = (*it)->OH1().Magnitude();
		oh2 = (*it)->OH2().Magnitude();
		// angle in degrees
		theta = acos((*it)->Angle()) * 180.0/M_PI;

		// lookup the molecular (frame) polarizability from the data file
		alpha = MatR::Zero();
		// The polarizability lookup value is given in atomic units (bohr^3)
		alpha = pdf.Value(oh1, oh2, theta);	// using the polarizability data file for lookup (pdf)
		// rotate the polarizability tensor into the lab-frame
		MatR dcm_t = dcm.transpose();
		alpha = dcm * alpha * dcm_t;
		(*it)->SetPolarizability (alpha);

		/*
			 mu = VecR::Zero();
		// dipole moment is given in Debye units
		mu = ddf.Value(oh1, oh2, theta);	// using the dipole data file (ddf)
		// rotate the polarizability tensor into the lab-frame
		mu = dcm * mu;
		(*it)->SetDipoleMoment (mu);
		//std::cout << mu.Magnitude() << std::endl;
		*/

		/*
			 (*it)->Print();
			 std::cout << std::endl << "oh's" << std::endl;
			 (*it)->OH1()->normalized().Print();
			 (*it)->OH2()->normalized().Print();
			 (*it)->Bisector().Print();
			 std::cout << "mu.unit" << std::endl;
			 mu.normalized().Print();
		//std::cout << mu.Magnitude() << std::endl; //<< ")  "; mu.Print(); std::cout << std::endl;
		//dcm.Print();
		//(*it)->Print();
		*/
	}

}	// lookup analysis - polarizability



/************* Methods for calculating dipoles and polarizabilities ****************/

template <class U>
void Morita2008Analysis<U>::CalcWannierDipoles () {
	std::for_each (analysis_wats.begin(), analysis_wats.end(), MDSystem::CalcWannierDipole);
}

template <class U>
void Morita2008Analysis<U>::MoritaH2OPolarizabilities () {
	std::for_each (analysis_wats.begin(), analysis_wats.end(), std::mem_fun<void,MoritaH2O>(&MoritaH2O::CalculateGeometricalPolarizability));
}

template <class U>
void Morita2008Analysis<U>::MoritaH2ODipoles () {
	std::for_each (analysis_wats.begin(), analysis_wats.end(), std::mem_fun<void,MoritaH2O>(&MoritaH2O::CalculateGeometricalDipoleMoment));
}



}	// namespace morita

#endif
