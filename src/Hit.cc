//
// Detector Hits:
// Hits and Slices
// \file Hit.cc
//

#include <CLHEP/config/iostream.h>
#include <fstream>
#include <stdlib.h>
#include "Hit.h"
#include "StateVector.h"
#include "defines.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Random/RandGauss.h"

#ifdef GL
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoSelection.h>
#endif

// drawing the hits.
#define SIG_RPHI_FAKT_CYL 300
#define SIG_RPHI_FAKT_DIS 200
#define SIG_Z_FAKT_DIS 1
#define SIG_Z_FAKT_CYL 1
#define SIG_RPHI_MIN_CYL 0.000025
#define SIG_Z_MIN_CYL 0.005
#define SIG_RPHI_MIN_DIS 0.000025
#define SIG_Z_MIN_DIS 0.005

#ifdef GL
extern void addChild ( SoSelection *root, SoNode *node , string (*info)(int) , int );
extern void addChild ( SoSelection *root, SoNode *node );
// extern void addChild ( SoSelection *root, SoNode *node , Element *el );
#endif
extern int hit_shows;
extern int matrix_prec; // precision of the matrices.

vector <Hit> hit_list;

// constructor
Hit::Hit ()
{
	_RPhi=0; _z=0; _sig=HepSymMatrix ( 2, 0);
};

Hit::Hit (const DetElement & el, const HepDouble & RPhi, const HepDouble & z )
{
	_sig=HepSymMatrix ( 2, 0);
	_sig(1,1)=el.sigRPhi(); 
	_el = el;
	if ( _el.Type()==CYL ) {
		_RPhi = RPhi; _z = z;
		_sig(2,2)=el.sigZ();
	} else {
		_x = RPhi;
		_y = z;
		_sig(2,2)=el.sigR();
	};
};

Hit::Hit ( StateVector & sv )
{
	_sig=HepSymMatrix ( 2, 0);
	_el= sv.Element();
	_sig(1,1)=_el.sigRPhi();
	_sv=sv; 
	if ( _el.Type() == CYL ) {
		_RPhi=sv.Phi() * sv.R(); _z=sv.z();
		_sig(2,2)=_el.sigZ();
	} else {
		_x=sv.x(); _y = sv.y();
		_sig(2,2)=_el.sigR();
	};
};


HepDouble Hit::R()  const { 
	if ( _el.Type()==CYL ) {
		return _el.R(); 
	} else {
		return norm(x(),y());
	};
};

HepDouble Hit::RPhi()   const { 
	if ( _el.Type()==CYL ) {
		return _RPhi; 
	} else { // this is not needed, but anyways ...
		return R()*Phi();
	};
};

HepDouble Hit::z()   const { 
	if ( _el.Type()==CYL ) {
		return _z; 
	} else {
		return _el.Z();
	};
};

HepDouble Hit::x()   const { 
	if ( _el.Type()==CYL ) {
		return R() * cos ( Phi() );
	} else {
		return _x;
	};
};

HepDouble Hit::y()   const { 
	if ( _el.Type()==CYL ) {
		return R() * sin ( Phi() );
	} else {
		return _y;
	};
};

HepDouble Hit::Sigz() const {
	#ifdef DEBUG
	if ( _el.Type()==DIS ) {
		WARNING ( "sigz with Disc coordinates.");
		return 0;
	};
	#endif
	return _sig(2,2);
};

HepDouble Hit::Varx() const {
	#ifdef DEBUG
	if ( _el.Type() == CYL )
		ERROR ("used Varx with cylinder coordinates.");
	#endif
	return sqr(cos(Phi())) * VarR()+ sqr(R())*sqr(sin(Phi()))*VarPhi();
};

HepDouble Hit::Vary() const {
	#ifdef DEBUG
	if ( _el.Type() == CYL )
		ERROR ("used Vary with cylinder coordinates.");
	#endif
	return sqr(sin(Phi())) * VarR()+ sqr(R())*sqr(cos(Phi()))*VarPhi();
};

HepDouble Hit::Covxy() const {
	#ifdef DEBUG
	if ( _el.Type() == CYL )
		ERROR ("used Vovxy with cylinder coordinates.");
	#endif
	return sin(Phi()) * cos(Phi()) * VarR() - \
		sqr(R())*sin(Phi())*cos(Phi())*VarPhi();
};

HepDouble Hit::VarR() const {
	#ifdef DEBUG
	if ( _el.Type() == CYL )
		ERROR ("used VarR with cylinder coordinates.");
	#endif
	return sqr( SigR() );
};

HepDouble Hit::Varz() const {
	#ifdef DEBUG
	if ( _el.Type() == DIS )
		ERROR ("used VarR with disc coordinates.");
	#endif
	return sqr( Sigz() );
};

HepDouble Hit::VarPhi() const {
	#ifdef DEBUG
	if ( _el.Type() == CYL )
		ERROR ("used VarPhi with cylinder coordinates.");
	#endif
	return sqr ( SigRPhi() / /* correct: R() */ _sv.R() );
};


HepDouble Hit::VarRPhi() const {
	return sqr ( SigRPhi() );
};

HepDouble Hit::Sigx() const {
	#ifdef DEBUG
	if ( _el.Type()==CYL ) {
		WARNING ( "sigx with cyl coordinates.");
		return 0;
	};
	#endif
	return sqrt( Varx());
};

HepDouble Hit::Sigy() const {
	#ifdef DEBUG
	if ( _el.Type()==CYL ) {
		WARNING ( "sigy with cyl coordinates." );
		return 0;
	};
	#endif
	return sqrt( Vary() );
};


HepDouble Hit::SigR() const {
	#ifdef DEBUG
	if ( _el.Type()==CYL ) {
		WARNING ( "sigr with cyl coordinates.");
		return 0;
	};
	#endif
	return _sig(2,2);
};

HepDouble Hit::SigRPhi() const {
	/* #ifdef DEBUG
	if ( _el.Type()==DIS ) {
		WARNING ( "sigrphi with disc coordinates.");
		return 0;
	};
	#endif */
	return _sig(1,1);
};

HepDouble Hit::Phi() const {
	if ( _el.Type()==CYL ) {
		return _RPhi / _el.R();
	} else {
		return myatan ( y() , x() );
	};
};


void Hit::write_raw () const {
	if ( _el.Type() == CYL ) {
		cout << "Hit Pt" << _el.TypeName() << " " << _RPhi << " " << _z << endl;
	} else {
		cout << "Hit Pt" << _el.TypeName() << " " << Phi() << " " << R()<< endl;
	};
};

void Hit::write_raw (const char *file) const {
	if (!strcmp(file,"")) {
		Hit::write_raw(); // hm. How can we merge the two?
	} else {
		static ofstream fout(file,ios::app);
		if ( _el.Type() == CYL ) {
			fout << "Hit Pt" << _el.TypeName() << " " << _RPhi << " " << _z << endl;
		} else {
			fout << "Hit Pt" << _el.TypeName() << " " << Phi() << " " << R()<< endl;
		};
	};
};
void Hit::write() const {
	cout.setf(ios::floatfield);
	cout.setf(ios::internal);
	cout.precision(3);
	cout << "Hit " << _el.TypeName();
	if (_el.Type()==CYL) {
		cout << " r=" << _el.R();
		cout.precision(5);
		cout << " Phi=" << _RPhi/_el.R() << " z=" << _z;
	} else {
		cout << " z=" << _el.Z();
		cout.precision(5);
		cout << " x=" << x() << " y=" << y();
	};
};

/// Smears out the measurement. deterministically, if det==true
void Hit::Smear ( const bool det )
{
	if (!det) {
		Smear();
		return;
	};
	if ( _el.Type() == CYL ) {
		HepDouble ran_dz    = pow(-1,_el.num) * _el.sigZ();
		HepDouble ran_dRPhi = pow(-1,_el.num) * _el.sigRPhi();
		_RPhi += ran_dRPhi;
		_z += ran_dz;
	} else {
		HepDouble myr, myphi;
		myr= R()+pow(-1,_el.num) * _el.sigR();
		myphi= Phi()+pow(-1,_el.num) * _el.sigRPhi() / /* correct: myr */ R();
		_x = myr * cos (myphi);
		_y = myr * sin (myphi);
	};
}

/// Smears out the measurement
void Hit::Smear ( )
{
	if ( _el.Type() == CYL ) {
		HepDouble ran_dz   = RandGauss::shoot(0, _el.sigZ() );
		HepDouble ran_dRPhi = RandGauss::shoot(0, _el.sigRPhi() );
		_RPhi += ran_dRPhi;
		_z += ran_dz;
	} else {
		HepDouble ran_dr   = RandGauss::shoot(0, _el.sigR() );
		HepDouble ran_dphi = RandGauss::shoot(0, _el.sigRPhi() / R() );
		HepDouble myr= R()+ran_dr;
		HepDouble myphi= Phi()+ran_dphi;
		_x = myr * cos (myphi);
		_y = myr * sin (myphi);
	};
}

#ifdef GL
void Hit::GLinitMaterial ( SoSelection *root ) const {
	SoMaterial *hitMaterial = new SoMaterial;
	hitMaterial->emissiveColor.setValue(0.6, 0, 0);
	hitMaterial->diffuseColor.setValue(0.8, 0, 0);
	addChild(root,hitMaterial);
};

string Hit_getInfo( const int a) 
{
	return hit_list[a].getInfo();
};

string Hit::Coords() const
{
	char myname[139];
	if ( _el.Type() == CYL ) {
		snprintf(myname,139,
			"r=%.*f Phi=%.*f z=%.*f\n     SigRPhi=%.*f, Sigz=%.*f",
		matrix_prec+3, _el.R(),
		matrix_prec+3, _RPhi / _el.R() ,
		matrix_prec+3, _z,
		matrix_prec+3, _sig[0][0] ,
		matrix_prec+3, _sig[1][1] );
	} else {
		snprintf(myname,139,
			"x=%.*f y=%.*f R=%.*f Phi=%.*f z=%.*f\n     SigRPhi=%.*f, SigR=%.*f", 
		matrix_prec+3, x() ,
		matrix_prec+3, y() ,
		matrix_prec+3, R() ,
		matrix_prec+3, Phi() ,
		matrix_prec+3, _el.Z(),
		matrix_prec+3, _sig[0][0] ,
		matrix_prec+3, _sig[1][1] );
	};
	return (string) myname;
};

string Hit::getInfo () 
{
	switch (hit_shows) {
	case 1:
		return _sv.identifierDelphi();
	case 2:
		return _sv.identifierCMS();
	case 3:
		return _sv.identifierDeriv();
	case 4:
		return _sv.identifierDer();
	default:
		return "Hit\n     "+Coords();
	};
};

void Hit::GLdraw ( SoSelection *root ) {


	SoEllipsoid *newSphere = new SoEllipsoid();
	register HepDouble myx,myy,myz;
	if ( _el.Type()==CYL ) {
		myy= - _el.R()*sin( _RPhi / _el.R() );
		myx= -_el.R()*cos( _RPhi / _el.R() );
		myz= _z * Z_SCALE;
	} else {
		myx= x();
		myy= y();
		myz= _el.Z();
	};
	HepDouble myc_z, myc_rphi;
	HepDouble rot;
	if ( _el.Type() == CYL ) {
		myc_z= _sig[1][1] < SIG_Z_MIN_CYL ? SIG_Z_MIN_CYL : _sig[1][1];
		myc_rphi= _sig[0][0] < SIG_RPHI_MIN_CYL ? SIG_RPHI_MIN_CYL : _sig[0][0];
		rot= fmod ( _sv.Phi()+_sv.beta() , M_PI ) ;
		newSphere->eigenvalues.setValue(  .01 , SIG_RPHI_FAKT_CYL * myc_rphi , \
			SIG_Z_FAKT_CYL * myc_z );
	} else {
		myc_z= _sig[1][1] < SIG_Z_MIN_DIS ? SIG_Z_MIN_DIS : _sig[1][1];
		myc_rphi= _sig[0][0] < SIG_RPHI_MIN_DIS ? SIG_RPHI_MIN_DIS : _sig[0][0];
		rot= fmod ( _sv.Phi() , M_PI ) ;
		newSphere->eigenvalues.setValue(  .01 , SIG_RPHI_FAKT_DIS * myc_rphi , \
			SIG_Z_FAKT_DIS * myc_z );
	};
	SbVec3f vec(0,0,1);
	
	newSphere->rotation.setValue ( vec , rot );
	newSphere->center.setValue ( myx, myy , myz );
	
	hit_list.push_back(*this);
	addChild(root, newSphere, (Hit_getInfo) , hit_list.size()-1 );
};

void Hit::getHit( StateVector & sv )
{
	_sig=HepSymMatrix(2,0);
	_el = sv.Element();
	if ( _el.Type() == CYL ) {
		_RPhi=sv.Phi() * sv.R(); _z = sv.z();
		_sig(1,1)=_el.sigRPhi();
		_sig(2,2)=_el.sigZ();
	} else {
		_x=sv.x(); _y = sv.y();
		_sig(1,1)=_el.sigRPhi();
		_sig(2,2)=_el.sigR();
	};
	_sv=sv;
};

#endif
