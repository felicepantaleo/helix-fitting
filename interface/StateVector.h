/*!
 * \file StateVector.h
 * \brief StateVector: the class for 5vectors.
 *
 * \author Wolfgang Waltenberger
 * \date Tue 30 Oct 2001 11:00:53 CET
 */

#ifndef StateVector_H
#define StateVector_H

#include "DetElement.h"
#include "defines.h"
#include "Plane.h"
#include <CLHEP/config/iostream.h>
#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"

#ifdef GL
#include <Inventor/nodes/SoSelection.h>
#endif

class StateVector {
protected:
	#ifdef GL
	void GLinitMaterial(SoSelection *,HepDouble,HepDouble,HepDouble ) const;
	void drawonCyl( SoSelection *, string, HepDouble );
	void drawonDis( SoSelection *, string, HepDouble );
	SoSelection *sel;
	#endif
	
	DetElement _el;
	HepDouble _r;
	HepDouble _x;
	HepDouble _beta;
	HepDouble _y;
	HepDouble _z, _theta, _phi, _Phi, _kappa; // Delphi coordinates.
	HepDouble _a0, _psi, _kappa_cms;          // CMS coordinates.
	HepMatrix _Deriv, _Der;
	HepDouble mys; // arclength
	void CMS2Delphi ();
	void CMS2Delphi ( HepDouble );
	void Delphi2CMS ();
	bool useCMS;
	bool useDelphi;
	StateVector PropagateCyl ( const DetElement & nxt , bool Drvtives );
	StateVector PropagateDis ( const DetElement & nxt , bool Drvtives );

public:
	/*
	HepDouble s() { return mys; };
	void put_s( HepDouble s) { mys=s;};
	*/
	StateVector ( const DetElement & );
	StateVector ( const DetElement & , const HepDouble & , const HepDouble & , const HepDouble & , const HepDouble & , const HepDouble & );
	StateVector ();
	/// 'unwrap' Phi, given the previous StateVector.
	void unwrapPhi ( StateVector * );
	/// randomize.
	/// give the intervals either by hand, or specify a file.
	void StateVector::Randomize( const string & filename );
	void StateVector::Randomize( );
	void StateVector::Randomize( interval, interval, interval, 
			interval, interval );
	
	/// The multiscattering effects are represented by this method.
	/// det = true: scatter deterministically.
	void StateVector::Scatter ( );
	void StateVector::Scatter ( bool det );

	/// setDelphi ( Phi, beta, kappa ) for CYL
	/// setDelphi ( Phi , r , phi ) for DIS
	/// [ used to be setDelphi ( x , y , phi ) for DIS ]
	void setDelphi ( const HepDouble & , const HepDouble & , const HepDouble& );

	/// setDelphi ( Phi , z , theta , beta , kappa) for CYL
	/// setDelphi ( Phi , r , theta , phi , kappa)  for DIS
	/// [ used to be setDelphi ( x , y , theta , beta , kappa) for DIS ]
	void setDelphi ( const HepDouble &, const HepDouble &, const HepDouble & ,  const HepDouble &, const HepDouble & );
	void setCMS ( const HepDouble & a0, const HepDouble & Psi, const HepDouble & kappa );
	void setCMS ( const HepDouble & a0, const HepDouble & z, const HepDouble & theta, const HepDouble & Psi, const HepDouble & kappa );
	void clear ();

	// Coordinates, general
	HepDouble z()     const { return _z; };
	HepDouble x()     const { return _x; };
	HepDouble y()     const { return _y; };
	HepDouble theta() const { return _theta; };
	HepDouble R()     const;
	HepDouble Radius() const;
	//	HepDouble s() ;
	HepMatrix Deriv() const { return _Deriv; }; ///< 5x5 Derivative Matrix
	/// 5x5 Derivative Matrix, 'chained ruled' to the beginning
	HepMatrix Der() const { return _Der; }; 
	HepDouble zVth_fak() const;
	/// Delphi coordinates
	HepVector Delphi();
	HepDouble Phi();
	HepDouble phi();
	HepDouble beta();
	HepDouble kappa();
	HepDouble eff_thick();
	DetElement Element() const { return _el; };

	// CMS coordinates
	HepDouble a0();
	HepDouble psi(); // for kappa, R, and theta: see above
	HepDouble kappa_cms();
	// 'Circle' coordinates
	HepVector XY0();

	// momentum of particle
	HepDouble pr () ;

	void write_raw ();
	void write_raw ( char * );
	void write ();
	void writeParam();
	void writeCMS ();
	void writeDelphi ();
	string DelphiValues( StateVector );
	string CMSValues( StateVector );
	string DelphiValues( StateVector * );
	string CMSValues( StateVector * );
	string DelphiValues();
	string CMSValues();
	string DelphiValues( char precision );
	string CMSValues( char precision );

	// inserts into StateVector the parameter of the plane.
	void getPlane ( Plane & , const DetElement & );
	void getPlane ( Plane & , const DetElement &, const HepDouble &);
	
	StateVector Propagate ( const DetElement & nxt , bool Drvtives );
	
	#ifdef GL
	string CMSString();
	string DelphiShortString();
	string CMSShortString();
	string DelphiString();
	void GLinitMaterial ( SoSelection * ) const;
	void GLdraw( SoSelection * , HepDouble lngth );
	void GLdraw( SoSelection * , HepDouble red, HepDouble green, \
			HepDouble blue, HepDouble lngth );
	void GLdraw( SoSelection * , const string & , HepDouble lngth );
	void GLdraw( SoSelection * , const string & , HepDouble red, HepDouble green, \
			HepDouble blue , HepDouble lngth );
	string getInfo ( int ) ;
	#endif
	// identifier:
	string DelphiFixedString();
	string identifierDelphi();
	string identifierCMS();
	string identifier; ///< the 'base' identifier
	string short_id;
	/// 5 x 5 derivative matrix, referring to the previous DetEl.
	string identifierDeriv(); 
	/// 5 x 5 derivative matrix, referring to the first DetEl (chainruled)
	string identifierDer(); 
};

#endif
