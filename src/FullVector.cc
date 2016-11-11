// $Id: FullVector.cc,v 1.1 2001/11/16 10:29:55 wwalten Exp $

#include <CLHEP/config/iostream.h>
#include <CLHEP/config/CLHEP.h>
#include <stdlib.h>
#include <stdio.h>
#include "defines.h"
#include "FullVector.h"
#include "Plane.h"
#include "Track.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "any2str.h"

#ifdef GL
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoSelection.h>
#include <Inventor/nodes/SoMaterial.h>
#include "hepvis/SoHelicalTrack.hh"
#endif /* ifdef GL */

#ifdef GL
extern void addChild ( SoSelection *, SoNode *, string (*info)(int), int );
extern void addChild ( SoSelection *, SoNode *, string (*info)(int), int,bool);
extern void addChild ( SoSelection *, SoNode * );
extern int hlx_shows;
extern int matrix_prec;
vector <FullVector> fv_list;
#endif

// a nicer way to print large matrices.
void mycout ( HepSymMatrix x )
{
	cout << any2str((char *) "",x);
};

FullVector::FullVector()
{
	StateVector::StateVector();
	_CMSCov = HepSymMatrix ( 5, 0);
	_DelphiCov = HepSymMatrix ( 5,0);
	useCMSCov=false; useDelphiCov=false;
	#ifdef GL
	identifier="FullVector";
	#endif
};

FullVector::FullVector ( const DetElement & el, HepDouble Phi , HepDouble z , \
		HepDouble theta, HepDouble beta, HepDouble kappa )
{
	_CMSCov = HepSymMatrix ( 5, 0);
	_DelphiCov = HepSymMatrix ( 5,0);
	useCMSCov=false; useDelphiCov=false;
	_el=el;
	_a0=0;
	_psi=0;
	if ( _el.Type() == CYL ) {
		_r=_el.R();
		_Phi=Phi; 
		_theta=theta; 
		_beta=beta; 
		_kappa=kappa; 
		_z=z;
	} else {
		ERROR ("not implemented");
	};
	_kappa_cms=0;
	useCMS=false;
	useDelphi=true;
	_Deriv=HepSymMatrix(5,0);
	_Der=HepSymMatrix(5,1);
	#ifdef GL
	identifier="FullVector";
	#endif
};

FullVector::FullVector( const DetElement & elm )
{
	_CMSCov = HepSymMatrix ( 5,0 );
	_DelphiCov = HepSymMatrix ( 5,0 );
	useCMSCov=false; useDelphiCov=false;
	_el=elm;
	if (_el.Type() == DIS) {
		_r=0; _z=_el.Z();
	} else {
		_r=_el.R(); _z=0; 
	};
	_Phi=0; _theta=0; _beta=0; _kappa=0;
	_kappa_cms=0; _phi=0;
	_a0=0; _psi=0; useCMS=false; useDelphi=true;
	_Deriv=HepSymMatrix(5,0);
	_Der=HepSymMatrix(5,1);
	
	#ifdef GL
	identifier="FullVector";
	#endif
};

/*
void FullVector::getPlane ( RealPlane pl , DetElement el )
{
	StateVector::getPlane (pl , el);
};

void FullVector::getPlane ( RealPlane pl , DetElement el,  )
{
	StateVector::getPlane (pl , el, r);
};
*/
void FullVector::writeCMS()
{
	StateVector::writeCMS();
	// cout << "5Vector/CovCMS:" << endl;
};

void FullVector::writeDelphi()
{
	StateVector::writeDelphi();
	// cout << "5Vector/CovDelphi:" << endl;
} 

void FullVector::setDelphiCov ( const HepSymMatrix & x)
{
	useDelphiCov=true;
	_DelphiCov=x;
	useCMSCov=false;
};

void FullVector::setCMSCov ( const HepSymMatrix & x)
{
	useCMSCov=true;
	_CMSCov=x;
	useDelphiCov=false;
};

HepSymMatrix FullVector::DelphiCov ()
{
	if (!useDelphiCov) CMS2DelphiCov();
	return _DelphiCov;
};

HepSymMatrix FullVector::CMSCov ()
{
	if (!useCMSCov) Delphi2CMSCov();
	return _CMSCov;
};

void FullVector::CMS2DelphiCov ()
{
	#ifdef DEBUG
	if (!useCMSCov) {
		ERROR ("cannot convert covariance matrix");
//		exit (-1);
	};
	#endif
	ERROR("Cov Matrix conversion not implemented.");
	useDelphiCov=true;
	_DelphiCov=_CMSCov; // oh FIXME 
};

void FullVector::Delphi2CMSCov ()
{
	#ifdef DEBUG
	if (!useDelphiCov) {
		ERROR ("cannot convert covariance matrix");
//		exit (-1);
	};
	#endif
	ERROR("Cov Matrix conversion not implemented.");
	useCMSCov=true;
	_CMSCov=_DelphiCov; // oh FIXME 
};
//                   ------------- graphix ------------
#ifdef GL

string FullVector::getInfo ( int what )
{
	int num=track_cov.num_row(), fail;
//	char name[355];
	switch (what) {
	case 0: // CMS
//		snprintf(name,355,"%s\n    Cov=", identifierCMS());
		return any2str ( identifierCMS() + "\n      Cov=", 
				CMSCov(), matrix_prec+3 );
	case 3: // delphi
		return any2str ( identifierDelphi()+"\n     Cov=" , 
				DelphiCov(), matrix_prec+3 );
	case 1: // covrphi
		return any2str ( identifier + "\n     CovRPhi(hits)=" , \
				track_cov.sub(1,num / 2 ), matrix_prec+3 );
	case 2: // covz
		if ( Element().Type()==CYL ) {
		return any2str ( identifier + "\n     Covz(hits)=" , \
				track_cov.sub( 1+num / 2, num ), matrix_prec+3 );
		} else {
		return any2str ( identifier + "\n     CovR(hits)=" , \
				track_cov.sub( 1+num / 2, num ), matrix_prec+3 );
		};
	case 4: // covrphi^-1
		return any2str ( identifier + "\n     CovRPhi^-1=" , \
				track_cov.sub(1,num / 2 ).inverse( fail ), matrix_prec+3 );
	case 99:
		return identifierCMS();
	default:
		return "error";
	};
};

string getfvInfo (int a)
{
	return fv_list[a].getInfo( hlx_shows );
};

void FullVector::drawonDis ( SoSelection *root , const string & idnt, HepDouble lngth ) {
	SoHelicalTrack *newHelix = new SoHelicalTrack;
	if (!(psi()<M_2PI) || !(theta() < M_2PI) ) {
		// FIXME I actually wanted an (SoText3) question mark here,
		// but this did not work on the Red (Dread) Hat 6.1/6.2
		newHelix->inverseRadius= 0;
		newHelix->phi0=M_PI;
		newHelix->z0=0;
		newHelix->d0=0;
		newHelix->cotTheta= 0;
		newHelix->s1=lngth;
		fv_list.push_back(*this);
		identifier=idnt+" !illegal!!";
		addChild ( root, newHelix , (getfvInfo) , fv_list.size()-1 , true );
	} else {
		newHelix->inverseRadius= 1 / Radius();
		newHelix->phi0= psi(); // -M_PI;
		HepDouble diff=phi()-psi() ;
		if (diff > M_PI ) diff -= M_2PI;
		if (diff < - M_PI ) diff += M_2PI; 
		newHelix->z0= _el.Z() + diff * tan ( M_PI_2 - theta() ) / kappa_cms ();
		newHelix->d0= - a0();
		newHelix->cotTheta= Z_SCALE / tan(theta()); // yes, scaling is that easy
		newHelix->s1=lngth;
		fv_list.push_back(*this);
		identifier=idnt;
		addChild ( root, newHelix , (getfvInfo) , fv_list.size()-1 , true );
	};
};
void FullVector::drawonCyl ( SoSelection *root , const string & idnt, HepDouble lngth ) {
	SoHelicalTrack *newHelix = new SoHelicalTrack;
	if (!(psi()<M_2PI) || !(theta() < M_2PI) ) {
		// FIXME I actually wanted an (SoText3) question mark here,
		// but this did not work on the Red (Dread) Hat 6.1/6.2
		newHelix->inverseRadius= 0;
		newHelix->phi0=0;
		newHelix->z0=0;
		newHelix->d0=0;
		newHelix->cotTheta= 0;
		newHelix->s1=lngth;
		fv_list.push_back(*this);
		identifier=idnt+" !illegal!!";
		addChild ( root, newHelix , (getfvInfo) , fv_list.size()-1 , true );
	} else {
		newHelix->inverseRadius= 1 / Radius();
		newHelix->phi0=psi() - M_PI;
		// the tan theta is the correction for teh difference z(R=0) 
		// and z(R=R(DetEl1))
		// FIXME sometime we should implement a 
		// FullVector sitting at R=0
		newHelix->z0= (z() - _el.R() / tan(theta())) * Z_SCALE;
		newHelix->d0=a0();
		newHelix->cotTheta= Z_SCALE / tan(theta()); // yes, scaling is that easy
		newHelix->s1=lngth;
		fv_list.push_back(*this);
		identifier=idnt;
		addChild ( root, newHelix , (getfvInfo) , fv_list.size()-1 , true );
	};
};

void FullVector::GLdraw ( SoSelection *root ,HepDouble r, HepDouble g,\
		HepDouble b, HepDouble lngth ) {
	GLinitMaterial ( root, r , g , b );
	if ( _el.Type() == CYL ) {
		drawonCyl ( root , "Track.", lngth );
	} else {
		drawonDis ( root , "Track.", lngth );
	};
};

void FullVector::GLdraw ( SoSelection *root , HepDouble lngth ) {
//	GLinitMaterial ( root );
	if ( _el.Type() == CYL ) {
		drawonCyl ( root , "Track.", lngth );
	} else {
		drawonDis ( root , "Track.", lngth );
	};
};

	
void FullVector::GLdraw ( SoSelection *root ,const string & idn , HepDouble r, \
		HepDouble g, HepDouble b , HepDouble lngth ) {
	GLinitMaterial ( root, r , g , b );
	if ( _el.Type() == CYL ) {
		drawonCyl ( root , idn, lngth );
	} else {
		drawonDis ( root , idn, lngth );
	};
};

void FullVector::GLdraw ( SoSelection *root , const string & idn , HepDouble lngth ) {
//	GLinitMaterial ( root );
	if ( _el.Type() == CYL ) {
		drawonCyl ( root , idn, lngth );
	} else {
		drawonDis ( root , idn, lngth );
	};
};
#endif /* ifdef GL */
