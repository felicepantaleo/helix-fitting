/*!
 * \file Track.h
 * \brief The track class is the list of RealVectors ...
 * plus their covariance matrix.
 * 
 * A simulated track can be generated by calling .Simulate()
 * A reference track can be generated by calling .ReferenceTrack()
 * Riemann fit can be applied to a track by calling .RiemannFit()
 * 
 * \author Wolfgang Waltenberger
 * \date Mon 09 Jul 2001 14:27:21 CEST
 */

#ifndef Track_H
#define Track_H

#include <string>
#include "DetElement.h"
#include "Detector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "Hit.h"
#include "Plane.h"
#include "FullVector.h"
#include "StateVector.h"
#include "Plane.h"
#include "defines.h"
#include <vector>

#ifdef GL
#include <Inventor/nodes/SoSelection.h>
#endif

// maybe good for gdb?
HepDouble ret_el ( const HepMatrix & m, int row, int col );
HepDouble ret_els ( const HepSymMatrix & m, int row, int col );

/// The track class.
class Track {
public:
	vector < Hit > hits;
private:
	int Num; ///< number of track vectors.
	HepSymMatrix _cov;
	/// we need to preserve the MS part as well 
	/// (RiemannFit on Discs)
	HepSymMatrix _cov_ms;
	vector < HepDouble > allR() const;
	HepVector allz() const;
	Detector *detector;
	#ifdef GL
	void GLinitMaterial ( );
	#endif
	// we keep copies of the propagated reference track.
	// else we would have to propagate all through for every fit.
	bool rt_computed;
	vector < StateVector > rt;
	// the arclengths of the StateVectors.
	HepVector rts;

	/// We store the 'updated' R coordinates separately.
	/// This is specific to the Riemann Fit on Discs -
	/// the other algos dont use this.
	bool R_pred_computed;
	HepVector R_pred;
	HepVector VarRPhi_pred;

	HepMatrix _mB; // B ( C * A' * G)
	HepSymMatrix _fv_cov; // the 5x5 covariance matrix of the detel-0_fv.
	
	/// compute the covariance matrix of a fv.
	/// B ( C * A' * G) is stored in _mB.
	HepSymMatrix computeFVCov( );
	int computed_cov, computed_fv_cov;
	bool last_fit_was_special;
	void initCov (); ///< Fills the diagonal of the 24x24-covariance matrix
	void initCov_special(); ///< Fills the diagonal of the 24x24-covariance matrix
	RealPlane RiemannFit_Plane ( bool MS , bool compute_cov );
	/// Computes the Covariance matrix.
	/// init_special=true: we use VarRPhi_pred for
	/// initialization.
	/// fullvar = true: we correct the NO intercept
	/// with the full covariance matrix.
	void force_computeCov ( int level );
	void computeCov ( int level , bool init_special, 
			bool fullvar, HepDouble kappa );
	void computeCov ( int level, bool init_special=false, bool fullvar=false );
	
	// void computeCov ( int level , bool init_special );
	// void computeCov ( int level );
	void computeCovMS ( Plane *, int level );
	FullVector KalmanDis ( int level );
	FullVector KalmanCyl ( int level );
	FullVector GlobalDis ( int level );
	FullVector GlobalCyl ( int level );
	#ifdef DEBUG
	/// for debugging purposes.
	void ReferenceTrackFromTrueTrack ( bool );
	/// we keep the 'true' StateVector;
	#endif
	StateVector truevec;
	inline void initR_pred();
	inline void force_initR_pred();
	/// Correct NonOrthogonal intercept in the covariance matrix.
	/// uses only the MS part of the covariance matrix for
	/// correction.
	inline void correct_NO ();
	/// Correct NonOrthogonal intercept in the covariance matrix.
	/// fullvar = true: use the full covariance matrix for
	/// correction.
	inline void correct_NO ( bool fullVar );
	inline void correct_NO ( bool fullVar, HepDouble );
public:
	Track () {};
	Track ( Detector * );
	~Track () {};
	
	Track & operator = ( Track &t );
	/// Simulates a 'regular' Track.
	/// Simulates a 'regular' Track.
	/// smear=0: Dont smear.
	/// smear=1: Do smear.
	/// smear=2: Do smear deterministically (only available in DEBUG mode)
	/// scatter=0: Dont scatter.
	/// scatter=1: Scatter randomly
	/// scatter=2: Scatter deterministically
	StateVector Simulate (char scatter=1, char smear=0, string="" );
	
	HepSymMatrix CovRPhi() const {return _cov.sub(1,Num); };
	HepSymMatrix CovPhi() const;
	// The phi quadrant of the MS part of the covariance matrix.
	HepSymMatrix CovMSPhi() const;
	// The MS part of the covariance matrix, given in Phi/R
	HepSymMatrix CovMSPhiR() const;
	HepSymMatrix Covz() const {return _cov.sub(Num+1,2*Num); }; // CYL
	HepSymMatrix CovR() const {return _cov.sub(Num+1,2*Num); }; // DIS
	HepSymMatrix Cov() const {return _cov; };
	void write() ;
	void write_raw( char * ) const;
	
	/// manually set the Reference Track.
	void setReferenceTrack ( StateVector, bool comp_der );

	/// Compute a reference track. Currently this 
	/// equals a riemannfit w/o multiscattering.
	void ReferenceTrack ( bool comp_der );
	/// Reference Track, but we keep a copy of the Plane object.
	void ReferenceTrack ( Plane  * , bool comp_der );
	
	// The algorithms.
	FullVector RiemannFit ( int level, int iterations = 1 );
	FullVector GlobalFit ( int level );
	FullVector KalmanFilter ( int level );
	#ifdef GL
	string Info( string );
	string MoreInfo( string );
	void GLdraw ( );
	#ifdef DEBUG
	void GLdrawTrueTrack();
	void GLdrawRefTrack();
	#endif
	#endif
	
//                                        --- more extensive documentation ---
/*! \fn void computeCov ( int level )
	\brief Computes the Covarianve matrix
    \param level

	possible levels are:
   - 0: no multiscattering.
   - 1: approximate multiscattering.
   - 2: exact multiscattering.
   - 3: exact maultscattering; z and RPhi components 'mix'.
 */

/*! \fn RealPlane RiemannFit ( int level , bool MS )
   \brief The Riemann fitting algorithm.

   See: 'Particle Tracks fitted on the Riemann sphere'
   R. Fr�hwirth, A. Strandlie, J. Wroldsen and B. Lillekjendlie.
   Computer Physics Communications 131 (2000) 95. 
   This method returns a RealPlane, instead of a FullVector.
 */

};
#endif /* Track_H */