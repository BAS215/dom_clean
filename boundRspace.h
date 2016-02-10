#ifndef _BOUNDRSPACE_H_
#define _BOUNDRSPACE_H_
#include "numerovR.h"
#include <iostream>
#include <cmath>
#include "pot.h"
#include "whit.h"
#include "coul.h"
#include "eigen.h"
#include "types.h"
#include "meshes.h"
#include "density.h"
#include <list>
#include <utility>

// type for holding an eigenvalue and its corresponding eigenvector
typedef std::pair< double, std::vector<double > > eigen_t;

// type for holding the eigenvalues and corresponding eigenvector
// for a given lj combination
typedef std::vector< eigen_t > lj_eigen_t;

// type for holding energy mesh 
typedef std::vector< std::pair< double, double > > mesh_t;

// type for holding propagator as a function of energy and r and r'
// (or k and k'). 
typedef std::vector< cmatrix_t > prop_t;

// a vector of vectors
typedef std::vector< std::vector< double > > vector_of_vectors;

class boundRspace : public numerovR
{
 public:

  static double const kconstant;
  static double const e2;
  static double const pi;
  double rmax;
  int hamiltonian_pts;
  int Lmax;
  double Efermi;
  double ph_gap;
  double Emin;
  double Z;
  double A;
  double Zp;
  double mu;
  void init( double rmax, int ham_pts, double Efermi, double ph_gap, 
             int Lmax, double Z, double Zp, double A );
  pot * Pot;
  double Kwave2;
  double Kwave;
  double gamma;
  double rStart;
  double rStop;
  double rdelt; // grid spacing used in constructing the Hamiltonian
  double Ecm;
  double muhbar;
  double* WaveArray;
  int Nham; // number of grid points used for constructing the Hamiltonian
  int Nnumerov;
  int nWave;
  int mWave;
  int max_n;

  boundRspace();
  ~boundRspace();

  boundRspace( double rmax, int ham_pts, double Efermi, double ph_gap, 
               int Lmax, double Z, double Zp, double A, pot*Pot ); 

  int nodes();
  int searchNonLoc(double & Elower, double Eupper, double j, int l, int Ifine);
  int searchLoc(double & Elower, double Eupper, double j, int l, int Ifine);
  int searchFine(double& Elower, double Eupper, double j, int l);
  void newEnergy(double Ecm);
  void integrateWaveFunction(double j, int l);

  std::vector< double > make_rmesh(); 

 std::vector< double > make_rmesh_point( );
 
  matrix_t 
  re_hamiltonian( const std::vector< double > &rmesh, double Ecm, 
                int L, double J );
  cmatrix_t 
  c_hamiltonian( const std::vector< double > &rmesh, double Ecm, 
                int L, double J );

  cmatrix_t 
  propagator( const std::vector<double> &rmesh, double E, int l, double j ); 
  std::vector< double > real_eigvals( const matrix_t &Ham ); 
  std::vector< eigen_t > real_eigvecs( const matrix_t &Ham ); 

  double
  find_level( const std::vector<double> &rmesh, double Estart, 
              int N, int l, double j, double tol ); 

  std::vector< double >
  normalize( const std::vector< double > &rmesh, 
             const std::vector<double> &eigenvector );
  eigen_t
  find_boundstate( const std::vector<double> &rmesh, double Estart, 
                   int N, int l, double j, double tol ); 

  double 
  sfactor( const std::vector<double> &rmesh, double QPE, int l, 
           double xj, const std::vector<double> &QPF ); 

  double 
  rms_radius( const std::vector<double> &rmesh, 
              const std::vector<double> &QPF ); 

  double 
  occupation( const std::vector<double> &rmesh, const matrix_t &d_mtx, 
              const std::vector<double> &QPF ); 

  double
  spectral_function_k_space( const std::vector< double > &rmesh, 
                             double E, double k ); 

  std::pair< double, double >
  find_max( const std::vector< double > &x_values, 
            const std::vector< double > &y_values ); 

  double 
  approx_width( const std::vector<double> &rmesh, double QPE, 
                int L, double J, const std::vector< double > &QPF ); 

  double
  find_width( const std::vector< double > &e_values,
              const std::vector< double > &s_of_E_values ); 

  std::vector< lj_eigen_t > 
  get_bound_levels( const std::vector< double > &rmesh, double tol ); 

  std::vector< mesh_t >
  get_emeshes( const std::vector< double > &rmesh, double Emin, double Emax,
               const std::vector< lj_eigen_t > &bound_levels ); 

  int
  L_from_index( unsigned int index ); 

  double
  J_from_index( unsigned int index ); 

  int
  index_from_LJ( int L, double J ); 


  mesh_t
  get_lj_emesh( int L, double J, const std::vector< mesh_t > &emesh_vec ); 

  std::vector< double > 
  get_lj_energy_vector( int L, double J, 
                        const std::vector< mesh_t > &emesh_vec ); 

  std::vector< prop_t >
  get_propagators( const std::vector< double > &rmesh, 
                   const std::vector< mesh_t > &emesh_vec ); 

  std::vector< double >
  spectral_strength( int L, double J,
                     const std::vector< double > &rmesh, 
                     const std::vector< mesh_t > &emesh_vec, 
                     const std::vector< prop_t > &prop_vec ); 


 std::vector< std::vector< double > >
 spectral_strength_der( int L, double J,
                                const std::vector< double > &rmesh, 
                                const std::vector< mesh_t > &emesh_vec, 
                                const std::vector< prop_t > &prop_vec );

  std::vector< double >
  spectral_strength_QH( int L, double J, 
                        const std::vector< double > &rmesh, 
                        const std::vector< mesh_t > &emesh_vec, 
                        const std::vector< prop_t > &prop_vec,
                        const std::vector< double > &QHF ); 

  std::vector< double >
  charge_density( const std::vector< double > &rmesh, double Emax,
                  const std::vector< mesh_t > &emesh_vec,
                  const std::vector< prop_t > &prop_vec,
                  const std::vector< lj_eigen_t > &bound_levels ); 

  std::vector< double >
  point_distribution( const std::vector< double > &rmesh, double Emax,
                      const std::vector< mesh_t > &emesh_vec,
                      const std::vector< prop_t > &prop_vec,
                      const std::vector< lj_eigen_t > &bound_levels ); 
  double
  chd_rms_radius( const std::vector< double > &rmesh,
                  const std::vector< double > &chd ); 

  double 
  particle_number( const std::vector< double > &rmesh, double Emax,
                   const std::vector< mesh_t > &emesh_vec,
                   const std::vector< prop_t > &prop_vec,
                   const std::vector< lj_eigen_t > &bound_levels );

  double y_end;
  double dydr_end;

  int phase; //!< which phase of the nonlocal iteration
  double Elower0; //!< saves the bound energy from the first iteration
  bool generateVloc; //!<indicates the equiv. local potential needs to generated

  void generateLocEquiv(int l);
  void nonLocalArray1(bool generateVloc, int l);
  void nonLocalArray0(int l);
  double nonlocalFact(double r1, double r2, double beta, int l);
  void integrateWaveFunction0(double y, double dydr, int l);

  double mstar; // !< effective mass relative to nucleon mass 
  double ** LogDerArray;
  int Nl; 
  int LogDerMin;
  int LogDerMax;
  int NLogDer;
  double LogDerDifference(double j, int l);
  double exteriorLogDer(int l);
  void normalizeWF();
  void exteriorWaveFunct(int l);

};

bool 
compare_N( eigen_t pair1, eigen_t pair2 );



#endif // _BOUNDRSPACE_H_
