// $Id: StateVector.cc,v 1.1 2001/11/16 10:29:55 wwalten Exp $
// -*- C++ -*-

#include <CLHEP/config/iostream.h>
#include <CLHEP/config/CLHEP.h>
#include <stdlib.h>
#ifdef DEBUG
#include <assert.h>
#endif
#include "defines.h"
#include "Plane.h"
#include "StateVector.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Matrix/Vector.h"
#include "any2str.h"


#ifdef GL
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoSelection.h>
#include <Inventor/nodes/SoMaterial.h>
#include "hepvis/SoHelicalTrack.hh"

vector <StateVector> sv_list;
extern int hlx_shows;
extern void addChild ( SoSelection *, SoNode *, string (*info)(int), int,bool);
extern void addChild ( SoSelection *, SoNode *, string (*info)(int), int );
extern void addChild ( SoSelection *root, SoNode *node );

extern int matrix_prec;
#endif /* ifdef GL */

extern void file_not_found (string filename);
inline string del_leading_spaces(string line) {
	while (line.find_first_of(" 	")==0) {
		line.erase(0,1);  // delete leading spaces and tabs
	};
	return line;
};

inline HepDouble random_interval ( HepDouble min , HepDouble max ) {
	return drand48() * (max - min) + min;
}

inline HepDouble random_interval ( interval iv ) {
	if ( iv.max > iv.min )
		return drand48() * (iv.max - iv.min) + iv.min;
	else {
		#ifdef DEBUG
		if (iv.max < iv.min ) {
			WARNING("iv.max " + any2str(iv.max) + " < iv.min " + 
					any2str(iv.min));
		};
		#endif
		return iv.min;
	};
}

HepDouble StateVector::Phi()
{
	#ifdef DEBUG
	HepDouble delta=fmod ( _Phi - _phi + _beta, M_2PI);
	if ( fabs ( delta) > .0000000001 && fabs(delta) < M_2PI-.000000000001 ) {
		ERROR ("Consistency check failed. delta="+
				any2str( delta , 11) +
				(string) " Phi="+any2str(_Phi)+
				(string) " phi="+any2str(_phi)+
				(string) " beta="+any2str(_beta));
	};
	#endif
	if (!useDelphi) CMS2Delphi();
	return _Phi;
};

HepDouble StateVector::phi()
{
	if (!useDelphi) CMS2Delphi();
	#ifdef DEBUG
	HepDouble delta=fmod ( _Phi - _phi + _beta, M_2PI);
	if ( fabs ( delta) > .0000000001 && fabs(delta) < M_2PI-.000000000001 ) {
		ERROR ("Consistency check failed: delta="+
				any2str( delta )+
				(string) " Phi="+any2str(_Phi)+
				(string) " phi="+any2str(_phi)+
				(string) " beta="+any2str(_beta)
				);
	};
	#endif
	return _phi;
};

HepDouble StateVector::beta()
{
	if (!useDelphi) CMS2Delphi();
	#ifdef DEBUG
	HepDouble delta=fmod ( _Phi - _phi + _beta, M_2PI);
	if ( fabs ( delta) > .0000000001 && fabs(delta) < M_2PI-.000000000001 ) {
		ERROR ("Consistency check failed: delta="+
				any2str(delta)+
				(string) " Phi="+any2str(_Phi)+
				(string) " phi="+any2str(_phi)+
				(string) " beta="+any2str(_beta));
	};
	#endif
	return _beta;
}

HepDouble StateVector::kappa()
{
	if (!useDelphi) CMS2Delphi();
	return _kappa;
}

void StateVector::clear ()
{
	useDelphi=false; useCMS=false;
	_Phi=0; _beta=0; _kappa=0; _z=0; _theta=0;
	_psi=0; _a0=0; _kappa_cms=0; _phi=0; _x=0; _y=0;
};
	
void StateVector::setDelphi ( const HepDouble & a, const HepDouble & b, 
    const HepDouble & c )
{
	if ( _el.Type() == DIS ) {
		// a = Phi , b = r
		_x= b * cos (a); _y= b * sin (a); _phi = c; _Phi= a;
		_beta= _phi - _Phi;
		#ifdef DEBUG
		if ( _Phi < 0 ) {
			_Phi += M_2PI;
			WARNING ( "watch your Phi: Phi e [0, 2pi]" );
		};
		#endif
	} else {
		_Phi=a; _beta=b; _kappa=c; _z=0; _theta=0;
		_phi=_beta + _Phi;
		#ifdef DEBUG
		if ( _Phi < 0 ) {
			_Phi += M_2PI;
			WARNING ( "watch your Phi: Phi e [0, 2pi]" );
		};
		#endif
	};
	if (_beta < - M_PI ) _beta += M_2PI;
	useDelphi=true; useCMS=false;
};

void StateVector::setDelphi ( const HepDouble & a, const HepDouble & b, 
    const HepDouble & c, const HepDouble & d, const HepDouble & e )
{
	if ( _el.Type() == DIS ) {
		_x=b * cos (a); _y=b * sin (a); _theta=c; _phi=d; _kappa=e;
		_Phi= a; 
		if ( _Phi < 0 ) {
			MESSAGE (10, "Phi not element [0, 2pi]. Phi="+
					any2str(_Phi)+" beta="+any2str(_beta) 
					+" phi="+any2str(_phi) );
			_Phi += M_2PI;
		};
		if ( _phi < 0 ) {
			MESSAGE (15, "phi not element [0, 2pi]. Phi="+
					any2str(_Phi)+" beta="+any2str(_beta) 
					+" phi="+any2str(_phi) );
			_phi += M_2PI;
		};
		_beta = _phi - _Phi;
	} else {
		_Phi=a; _z=b; _theta=c; _beta=d; _kappa=e;
		if ( _Phi < 0 ) {
			MESSAGE (10, "Phi not element [0, 2pi]" );
			_Phi += M_2PI;
		};
		/*
		if ( _phi < 0 ) {
			_phi += M_2PI;
			MESSAGE (10, "phi not element [0, 2pi]. Phi="+
					any2str(_Phi)+" beta="+any2str(_beta) 
					+" phi="+any2str(_phi) );
		};
		*/
		_phi=_beta + _Phi;
	};
	#ifdef DEBUG
	if (_beta < - M_PI ) {
		MESSAGE(10,"we just corrected beta.");
		_beta += M_2PI;
	};
	#endif
	useDelphi=true; useCMS=false;
};

void StateVector::setCMS ( const HepDouble & a0, const HepDouble & z, 
    const HepDouble & theta, const HepDouble & psi,const HepDouble & kappa )
{
	_psi=psi; _a0=a0; _kappa_cms=kappa;
	_theta=theta, _z=z;
	useCMS=true; useDelphi=false;
}
void StateVector::setCMS ( const HepDouble & a0, const HepDouble & psi, 
    const HepDouble & kappa )
{
	_psi=psi; _a0=a0; _kappa_cms=kappa; _z=0; _theta=0;
	useCMS=true; useDelphi=false;
}

void StateVector::Delphi2CMS ()
{
	#ifdef DEBUG
	if (!useDelphi) {
		ERROR ("cannot convert StateVector from Delphi to CMS");
		exit (-1);
	};
	#endif
	double PsiPrime; // D=beta+_Phi;
	PsiPrime=Phi()+beta();
	HepDouble Xic = -1 * sin (PsiPrime) / _kappa;
	HepDouble Yic = cos (PsiPrime) / _kappa;
	HepDouble Xi  = cos ( Phi() ) * R();
	HepDouble Yi  = sin ( Phi() ) * R();
	HepDouble u0  = Xi + Xic;
	HepDouble v0  = Yi + Yic;
	HepDouble Rh = 1 / fabs ( kappa() );
	
	HepDouble XY0_norm = norm(u0,v0);
	HepDouble nX0 = u0 / XY0_norm; // = nR0
	HepDouble nY0 = v0 / XY0_norm; // = nR0
	HepDouble RVecX = Rh *nX0;
	HepDouble RVecY = Rh *nY0;
	HepDouble XVecX = R() * cos ( Phi() ); // XXX
	HepDouble XVecY = R() * sin ( Phi() );
	
	HepDouble a0VecX = u0 - Rh * nX0;
	HepDouble a0VecY = v0 - Rh * nY0;
	
	HepDouble XVecPrimeX= XVecX - a0VecX;
	HepDouble XVecPrimeY= XVecY - a0VecY;
	HepDouble norm_PsiVec = 
		sign( XVecPrimeX * RVecY-XVecPrimeY*RVecX) / norm(RVecX, RVecY);
	HepDouble nPsiVecX = norm_PsiVec * RVecY;
	HepDouble nPsiVecY = - norm_PsiVec * RVecX;
	_psi=atan2(nPsiVecY,nPsiVecX);
	
	HepDouble sign_a0   = - sign ( a0VecX * nPsiVecY - a0VecY * nPsiVecX);
	_a0= sign_a0 * norm(a0VecX,a0VecY);
	_kappa_cms= - _kappa;
	useCMS=true;
	if (_el.Type() == DIS ) {
		HepDouble diff = fmod ( _phi - _psi + M_PI, M_2PI ) - M_PI;
/*		HepDouble diff = _phi - _psi;
		if ( diff < - M_PI ) diff += M_2PI;
		if ( diff > M_PI ) diff -= M_2PI;*/
		_z= _el.Z()+diff * tan ( M_PI_2 - theta() ) / _kappa_cms;
	};
};

void StateVector::CMS2Delphi ( HepDouble myr )
{
	#ifdef DEBUG
	if (!useCMS) {
		ERROR ("cannot convert StateVector from CMS to Delphi");
		exit (-1);
	};
	#endif
	HepDouble koeff = _a0 - 1 / _kappa_cms;
	HepDouble Rh= 1 / fabs(_kappa_cms);
	HepDouble CosDeltaPhi;
	HepDouble acos_cosdeltaphi, alpha;
	HepDouble CentrePointx = - koeff * sin( _psi );
	HepDouble CentrePointy = koeff * cos( _psi );
	CosDeltaPhi = (sqr( Rh ) + sqr(CentrePointx)+sqr(CentrePointy) \
		  	- sqr( myr ) ) / (2 * Rh * norm (CentrePointx,CentrePointy));
	if (CosDeltaPhi>=1.) { // an ugly hack: acos(1.) is out of range!
		WARNING("problem with acos(" + any2str(CosDeltaPhi)+ ")"); 
		acos_cosdeltaphi=0;
	} else {
		acos_cosdeltaphi=acos(CosDeltaPhi);
	};
	HepDouble RPrime = sqrt ( 2* sqr(Rh)*(1-CosDeltaPhi));
	HepDouble PsiPrime = _psi- sign(_kappa_cms)* acos_cosdeltaphi;
	alpha=PsiPrime;
	HepDouble PhiPrime=_psi-0.5*sign(_kappa_cms)*acos_cosdeltaphi;
	HepDouble xprime=RPrime*cos(PhiPrime);
	HepDouble yprime=RPrime*sin(PhiPrime);
	
	HepVector X = HepVector ( 2, 0);
	X(1)= -1 * _a0 * sin ( _psi )+ xprime;
	X(2)= _a0 * cos ( _psi )+ yprime;
	_x=X(1);
	_y=X(2);
	_Phi= atan2(X(2),X(1));
	_Phi = _Phi < 0 ? _Phi + M_2PI : _Phi ;
	_beta= alpha - _Phi;
//	_beta = fmod ( _beta + M_PI , M_2PI ) - M_PI;
	_beta = _beta < - M_PI ? _beta + M_2PI : _beta;
	_phi = _Phi + _beta; // the M_PI part is a lie
	_kappa= - _kappa_cms; // FIXME is this true?
	useDelphi=true;
}

void StateVector::CMS2Delphi ()
{
	#ifdef DEBUG
	if (!useCMS) {
		ERROR ("cannot convert StateVector from CMS to Delphi");
		exit (-1);
	};
	#endif
	HepDouble myr;
	HepDouble koeff = _a0 - 1 / _kappa_cms;
	if ( _el.Type() == CYL ) {
		myr=R();
	} else {
		HepDouble mphi = _psi + tan ( _theta ) * _kappa_cms * ( _z - _el.Z() );
		HepDouble mx = - koeff * sin ( _psi ) - sin ( mphi ) / _kappa_cms;
		HepDouble my = koeff * cos ( _psi ) + cos ( mphi )/ _kappa_cms;
		myr = norm ( mx , my );
	};
	CMS2Delphi ( myr );
}


void StateVector::getPlane ( Plane & pl , const DetElement & el, 
    const HepDouble & r )
{
	#ifdef DEBUG
	if ( _el.Type() == UNK )
		ERROR("5Vector uninitialized where it shouldnt be!");
	#endif
	_el=el; _r=el.R(); // _z=el.Z();
	HepDouble mz=pl.z, mtheta=pl.theta;
	setCMS ( pl.a0(), pl.psi(), pl.kappa() );
	CMS2Delphi ( r );
	HepDouble diff = fmod ( _phi - _psi + M_PI, M_2PI ) - M_PI;
	mz= mz+diff * tan ( M_PI_2 - mtheta ) / pl.kappa();
	setCMS ( pl.a0(), mz, mtheta, pl.psi(), pl.kappa() );
	CMS2Delphi();
};

void StateVector::getPlane ( Plane & pl ,const DetElement & el )
{
	#ifdef DEBUG
	if ( _el.Type() == UNK )
		ERROR("5Vector uninitialized where it shouldnt be!");
	#endif
	_el=el; _r=el.R(); // _z=el.Z();
	setCMS ( pl.a0(), pl.z, pl.theta, pl.psi(), pl.kappa() );
	CMS2Delphi();
};

HepDouble StateVector::pr()
{
	return fabs ( CNV * _el.Bz() / kappa() / sin (theta()));
};

//           ----------------< write routines >----------

void StateVector::write_raw ()
{
	cout << "SV" << _el.TypeName() << short_id << " " << DelphiValues(8);
};

void StateVector::write_raw (char *file )
{
	if (!strcmp(file,"")) {
		StateVector::write_raw(); // hm. How can we merge the two?
	} else {
		static ofstream fout(file,ios::app);
		fout << "SV" << _el.TypeName() << short_id << " " << DelphiValues(8);
	};
};

void StateVector::write()
{
	writeDelphi();
};

void StateVector::writeDelphi()
{
	cout << identifier << " (Delphi):";
	if (_el.Type()==CYL) {
		cout << "   Element #" << _el.num << " el.R()=" << R() << endl
	         << "   Phi=" << setprecision(9) << C_GREEN << Phi() << C_DEFAULT << " z=" << 
			 C_GREEN << z() << C_DEFAULT << " theta=" << C_GREEN 
			 << theta() << C_DEFAULT
	         << " beta=" << C_GREEN << beta() << C_DEFAULT << 
			 " kappa=" << C_GREEN << kappa() << C_DEFAULT << endl;
	} else {
		cout << "   Element #" << _el.num << " el.Z()=" << _el.Z() << endl
	         << setprecision(9) << "     Phi=" << C_GREEN << Phi() << C_DEFAULT << " R=" 
			 << C_GREEN << R() << C_DEFAULT << " theta=" << 
			 C_GREEN << theta() << C_DEFAULT
	         << "  phi=" <<  C_GREEN << phi() << C_DEFAULT << " kappa=" 
			 << C_GREEN << kappa() << C_DEFAULT << endl;
	#ifdef DEBUG
	if (_el.Type()==CYL) 
		WARNING ("this phi should only occur with DIS coordinates.");
	#endif
	};
};

void StateVector::writeCMS()
{
	cout << identifier << " (CMS):" << "   Element #" << _el.num;
	if ( _el.Type() == CYL ) {
		cout << " el.R()=" << _el.R() << endl;
	} else {
		cout << " el.Z()=" << _el.Z() << endl;
	};
	
	cout << "   a0=" << C_GREEN << a0() << C_DEFAULT << " z=" << C_GREEN << 
		 z() << C_DEFAULT << " theta=" << C_GREEN << theta() << C_DEFAULT << 
		    " psi=" << C_GREEN << psi() << C_DEFAULT << " kappa_cms=" << 
			C_GREEN << kappa_cms() << C_DEFAULT << endl;
};

HepVector StateVector::Delphi()
{
	HepVector d(5,0);
	if (_el.Type()==CYL) {
		d(1)=Phi();
		d(2)=z();
		d(3)=theta();
		d(4)=beta();
		d(5)=kappa();
	} else {
		d(1)=Phi();
		d(2)=R();
		d(3)=theta();
		d(4)=phi();
		d(5)=kappa();
	};
	return d;
};

// scatter deterministically, if det == true
void StateVector::Scatter ( bool det )
{
	if (!det) {
		Scatter();
		return;
	};
	HepDouble Bz=_el.Bz();
	/// \bug Again, ashamed and embarrassed I specify to B={0,0,Bz}
	HepDouble pT, p, X, sigms_sim; // , randa,randb;
	pT= CNV * Bz / _kappa;
	p= pT / sin( _theta ) ;
	X= eff_thick(); // _el.Thick() / sin ( _theta );
	sigms_sim=0.015 / fabs (p) * sqrt (X) * ( 1+.038 * log (X));
	HepDouble ran1=sigms_sim;
	HepDouble ran2=sigms_sim / sin ( _theta);
	_theta=_theta+ran1;
	if ( _el.Type() == CYL ) {
		_beta += ran2;
		_phi = _beta + _Phi;
	} else {
		_phi += ran2;
		_beta = _phi - _Phi; // unneccessary?
	};
	pT=p*sin(_theta);
	_kappa=CNV * Bz / pT;
	
};
void StateVector::Scatter ( )
{
	HepDouble Bz=_el.Bz();
	HepDouble pT, p, X, sigms_sim; // , randa,randb;
	pT= CNV * Bz / _kappa;
	p= pT / sin( _theta ) ;
	X= eff_thick(); // _el.Thick() / sin ( _theta );
	sigms_sim=0.015 / fabs ( p ) * sqrt (X) * ( 1+.038 * log (X));
	HepDouble ran1 = RandGauss::shoot(0,sigms_sim);
	HepDouble ran2 = RandGauss::shoot(0,sigms_sim / sin (_theta));
	_theta=_theta+ran1;
	if ( _el.Type() == CYL ) {
		_beta += ran2;
		_phi = _beta + _Phi;
	} else {
		_phi += ran2;
		_beta = _phi - _Phi; // unneccessary?
	};
	pT=p*sin(_theta);
	_kappa=CNV * Bz / pT;
};

StateVector StateVector::Propagate ( const DetElement & nxt , bool Drvtvs )
{
	switch (_el.Type()) {
		case CYL:
			return PropagateCyl ( nxt, Drvtvs );
		case DIS:
			return PropagateDis ( nxt, Drvtvs );
		default:
			ERROR("dunno what kinda element this is.");
	};
	exit (-1);
};

string StateVector::DelphiValues() {
	return DelphiValues ( 7 );
};

string StateVector::CMSValues() {
	return CMSValues ( 7 );
};

string StateVector::DelphiValues( char pre ) {
	char line[80];
	if ( Element().Type() == CYL ) {
		snprintf(line,80,"%.*e %.*e %.*e %.*e %.*e\n",pre,Phi()
				,pre,z(),pre, theta(),pre, beta(),pre, kappa());
	} else {
		snprintf(line,80,"%.*e %.*e %.*e %.*e %.*e\n",pre,Phi()
				,pre,R(),pre, theta(),pre, phi(),pre, kappa());
	};
	return line;
};
	
string StateVector::CMSValues( char pre ) {
	char line[80];
	snprintf(line,80,"%.*e %.*e %.*e %.*e %.*e\n",pre,a0()
			,pre,z(),pre, theta(),pre, psi(),pre, kappa_cms());
	return line;
};

string StateVector::DelphiValues( StateVector tv) {
	char line[80];
	if ( Element().Type() == CYL ) {
		snprintf(line,80,"%.7e %.7e %.7e %.7e %.7e\n",Phi()-tv.Phi()
				,z()-tv.z(), theta()-tv.theta(),
				beta()-tv.beta(), kappa()-tv.kappa()
		);
	} else {
		snprintf(line,80,"%.7e %.7e %.7e %.7e %.7e\n",Phi()-tv.Phi()
				,R()-tv.R(), theta()-tv.theta(),
				phi()-tv.phi(), kappa()-tv.kappa()
		);
	};
	return line;
};

string StateVector::CMSValues( StateVector tv) {
	char line[80];
	snprintf(line,80,"%.7e %.7e %.7e %.7e %.7e\n",a0()-tv.a0()
			,z()-tv.z(), theta()-tv.theta(),
			psi()-tv.psi(), kappa_cms()-tv.kappa_cms()
	);
	return line;
};

HepDouble StateVector::R() const 
{ 
	if ( _el.Type() == CYL ) {
		return _r; 
	} else {
		return norm(_x,_y);
	};
};

void StateVector::unwrapPhi( StateVector *prev )
{
	if ( prev->Phi() - _Phi > M_PI ) _Phi += M_2PI;
	if ( _Phi - prev->Phi() > M_PI ) _Phi -= M_2PI;
	if ( prev->phi() - _phi > M_PI ) _phi += M_2PI;
	if ( _phi - prev->phi() > M_PI ) _phi -= M_2PI;
	if ( prev->beta() - _beta > M_PI ) _beta += M_2PI;
	if ( _beta - prev->beta() > M_PI ) _beta -= M_2PI;
}

StateVector StateVector::PropagateDis ( const DetElement & nxt , bool Drvtvs )
///< Propagates a StateVector
{
	#ifdef DEBUG
	if (!useDelphi) {
		// this one is based on Delphi coordinates.
		CMS2Delphi();
		WARNING("Oops! we need to convert here! Change the code.");
	};
	#endif
	
	StateVector _neu;
	
	#ifdef DEBUG
	if (nxt.Type()==UNK || nxt.Type()==NON) { 
		// we cannot propagate any further.
		ERROR ("could not determine next detector element.\n       element #"
			  	+ any2str(_el.num) + " type " + _el.TypeName());
		return _neu;
	};
	HepDouble rmin=-200, rmax=200;
	#endif
	HepDouble dphimn=pow(10,-4);
	
	HepDouble dz=nxt.Z()-_el.Z();
	HepDouble tanth=tan(theta());
	HepDouble rdphi=dz*tanth;
	HepDouble rtrk= 1 / kappa();
	HepDouble cosf0 = cos ( phi() );
	HepDouble sinf0 = sin ( phi() );
	HepDouble xc=x()-rtrk*sinf0;
	HepDouble yc=y()+rtrk*cosf0;
	HepDouble dphi=kappa()*rdphi;
	HepDouble phi1=fmod (phi()+dphi, M_2PI );
	HepDouble cosf1=cos(phi1);
	HepDouble sinf1=sin(phi1);
	HepDouble x1=xc+rtrk*sinf1;
	HepDouble y1=yc-rtrk*cosf1;
	#ifdef DEBUG
	HepDouble r1=norm ( x1,y1);
	if (r1 < rmin || r1 > rmax) {
		ERROR("Did not intersect: " + any2str(rmin) 
				+ " " + any2str(r1) + " " + any2str(rmax));
		exit (-1);
	};
	#endif
	
	_neu=StateVector( nxt );
	_neu._x     = x1;
	_neu._y     = y1;
	_neu._theta = theta();
	_neu._phi   = phi1;
	if (phi1 < 0. ) _neu._phi += M_2PI;
	_neu._Phi   = myatan ( y1 , x1);
	if (_neu._Phi < 0 ) _neu._Phi += M_2PI;
	_neu._beta  = phi1 - _neu._Phi;
	_neu._kappa = kappa();
	if (Drvtvs) { // ---------- computation of derivatives -----------
		_neu._Deriv=HepMatrix(5,5,1);
		HepDouble ct2inv=1.+sqr(tanth);
		_neu._Deriv(1,3)=ct2inv*dz*cosf1;
		_neu._Deriv(2,3)=ct2inv*dz*sinf1;
		_neu._Deriv(4,3)=dz*kappa()*ct2inv;
		_neu._Deriv(4,5)=rdphi;
		
		if ( fabs(dphi) >= dphimn ) {
			HepDouble dcosf=cosf1-cosf0;
			HepDouble dsinf=sinf1-sinf0;
			_neu._Deriv(1,4)=rtrk*dcosf;
			_neu._Deriv(1,5)=sqr(rtrk)*(dphi*cosf1-dsinf);
			_neu._Deriv(2,4)=rtrk*dsinf;
			_neu._Deriv(2,5)=sqr(rtrk)*(dphi*sinf1+dcosf);
		} else {
			_neu._Deriv(1,4)=-rdphi*sinf0;
			_neu._Deriv(1,5)=.5*rdphi*_neu._Deriv(1,4);
			_neu._Deriv(2,4)=rdphi*cosf0;
			_neu._Deriv(2,5)=.5*rdphi*_neu._Deriv(2,4);
		};

		// transform to Phi, R.
		// FIXME optimize it?
		HepMatrix Ti ( 5,5,1), Tf (5,5,1);
		Ti(1,1)=-R() * sin (Phi());
		Ti(2,2)= sin (Phi() );
		Ti(1,2)= cos ( Phi() );
		Ti(2,1)= R() * cos ( Phi());
		Tf(1,1)=- sin (_neu.Phi()) / _neu.R();
		Tf(2,2)= sin ( _neu.Phi() );
		Tf(1,2)= cos ( _neu.Phi() ) / _neu.R();
		Tf(2,1)= cos ( _neu.Phi());
		_neu._Deriv = Tf * _neu._Deriv * Ti;
		_neu._Der=_neu._Deriv * _Der; // chain rule
	};
	
	_neu.useDelphi=true;
	_neu.useCMS=false;
	#ifdef GL
	_neu.Delphi2CMS();
	_neu.identifier="Propagated StateVector";
	_neu.short_id="PS";
	#endif
//	_neu.unwrapPhi ( this );
	return _neu;
};

StateVector StateVector::PropagateCyl ( const DetElement & nxt , bool Drvtvs )
///< Propagates a StateVector
{
	// this one is based on Delphi coordinates.
	#ifdef DEBUG
	if (!useDelphi) {
		CMS2Delphi();
		WARNING("oops! we needa convert here!");
	};
	#endif
	
	// default values
	#ifdef DEBUG
	HepDouble sinbmx=0.9, zmin= -200, zmax = 200;
	#endif
	// end default values
	//
	
	HepDouble dphimn=pow(10,-4);
	
	HepDouble cosb, sinb, cotth, rtrk, xc, yc, rc2, rrr, delt;
	HepDouble phif, zf, xf, yf, tth, alff, dr2, rcosb, radi, phii, radf;
	HepDouble dphi, alrphi, pph, sinf, cosf, sinbf, dcosf, aa, dsinf, rdphi;
	StateVector _neu;
	
	#ifdef DEBUG
	if (nxt.Type()==UNK || nxt.Type()==NON) { 
		// we cannot propagate any further.
		ERROR ("could not determine next detector element.\n       element #"
			  	+ any2str(_el.num) + " type " + _el.TypeName());
		return _neu;
	};
	#endif
	_neu=StateVector( nxt );
	
	cosb= cos(_beta);
	sinb= sin(_beta);
	cotth=(1/tan(_theta)); 
	rtrk= 1/ _kappa; // FIXME we dont need this one
	phii= _Phi; radf=_neu.R();
	radi = _el.R();
	
	xc=radi-rtrk*sinb;
	yc=rtrk*cosb;
	rc2=pow(xc,2)+pow(yc,2);
	rrr=(pow(radf,2)-pow(rtrk,2)-rc2) / ( 2 * rtrk );
	delt = rc2 - pow(rrr,2);
	#ifdef DEBUG
	if (delt <=0 ) {
		ERROR("did not intersect. delt = " + any2str(delt) + "< 0");
		exit (-1);
	};
	#endif
	delt=sqrt(delt);
	sinf=(xc*rrr+yc*delt)/rc2;
	cosf=(xc*delt-yc*rrr)/rc2;
	xf=xc+rtrk*sinf;
	yf=yc-rtrk*cosf;
	sinbf=(sinf*xf-cosf*yf)/radf;
	#ifdef DEBUG
	if ( fabs ( sinbf) > sinbmx ) { // FIXME is this ok?
		ERROR("Beta too large");
		exit (-1);
	};
	#endif
	alff=atan2(sinf,cosf);
	dphi=alff-_beta;
	alrphi=rtrk*dphi;
	if ( fabs( dphi )>=dphimn ) { // FIXME is this ok?
		zf=_z+cotth*rtrk*dphi;	
		tth= M_PI / 2 - atan (zf / radf);
		phif=atan2 (yf, xf);
		pph=phif;
		if ( pph < 0 ) {pph=pph+ M_2PI;};
		#ifdef DEBUG
		if (zf < zmin || zf > zmax ) {
			ERROR("z not in the given interval (z=" + any2str(zf) +
			     " zmin,zmax=" + any2str(zmin)+ ", " + any2str(zmax)+")");
			exit (-1);
		};
		#endif
		_neu._Phi  = fmod(phii + phif + M_2PI, M_2PI);
		_neu._z     = zf;
		_neu._theta = _theta;
		_neu._beta  = alff-phif;
		_neu._kappa = _kappa;
		_neu._phi = _neu._Phi + _neu._beta;
		if (Drvtvs) { // ---------- computation of derivatives -----------
			_neu._Deriv=HepMatrix(5,5,1);
		//	for (int i=1;i<6; i++) {_neu._Deriv(i,i)=1;};
			HepDouble cosbf = sqrt(1-sqr(sinbf));
			HepDouble ccpsi = radi- rtrk * sinb;
			HepDouble scpsi = rtrk * cosb;
			HepDouble ccpsf = radf-rtrk*sinbf;
			// HepDouble scpsf = rtrk*cosbf;
			HepDouble cpsii = rtrk- radi*sinb;
			HepDouble spsii = -radi*cosb;
			HepDouble cpsif = rtrk - radf * sinbf;
			HepDouble spsif = -radf * cosbf;
			HepDouble sdphi = sin ( dphi );
			HepDouble cdphi = cos ( dphi );
			HepDouble fact = - rtrk / spsif;
			_neu._Deriv(1,4)=sdphi*fact;
			_neu._Deriv(1,5)=fact*rtrk*(1-cdphi);
			_neu._Deriv(2,3)=-rtrk*dphi*(1+sqr(cotth));
			_neu._Deriv(2,4)=rtrk*cotth*(radf*ccpsf*spsii/spsif-radi*ccpsi)/rc2;
			_neu._Deriv(2,5)=sqr(rtrk)*cotth*(-dphi+sinbf/cosbf- \
					(radi*scpsi+radf*ccpsf*cpsii/spsif)/rc2);
			_neu._Deriv(4,4)=spsii/spsif;
			_neu._Deriv(4,5)=rtrk*(cpsif-cpsii)/spsif;
			_neu._Der=_neu._Deriv * _Der; // chain rule
		};
	} else {
		dr2=pow(radf,2)-pow(radi,2);
		rcosb=radi*cosb;
		aa=1-radi*sinb/rtrk;
		delt=pow(rcosb,2)+aa*dr2;
		#ifdef DEBUG
		if (delt<=0) {
			ERROR("rcsob=" + any2str(rcosb) + " _neu-r=" + any2str(_neu.R())+ " r=" 
				+ any2str(_r)+ " aa=" + any2str(aa) + " dr2=" + any2str(dr2)
				+ "did not intersect. delt = " + any2str(delt)+ "< 0");
			return StateVector();
		};
		#endif
		rdphi=(sqrt(delt)-rcosb)/aa;
		dphi=rdphi/rtrk;
		dcosf= - sinb - .5 *cosb*dphi;
		cosf=cosb+dcosf*dphi;
		yf=-rdphi*dcosf;
		dsinf=cosb-.5*sinb*dphi;
		sinf=sinb+dsinf*dphi;
		xf=radi+rdphi*dsinf;
		sinbf=(sinf*xf-cosf*yf)/radf;
		zf=_z+cotth*rdphi;
		tth=M_PI_2-atan(zf/radf);
		phif=atan2(yf,xf);
		pph=phif;
		if (pph < 0) { pph=pph + M_2PI;};
		#ifdef DEBUG
		if (zf < zmin || zf > zmax ) {
			ERROR("z not in the given interval (z=" + any2str(zf) +
			     " zmin,zmax=" + any2str(zmin) + ", " + any2str(zmax)+")");
			exit (-1);
		};
		#endif
		phif = atan2(yf,xf);	
		_neu._Phi   = fmod ( _Phi + phif + M_2PI, M_2PI );
		_neu._z     = zf;
		_neu._theta = _theta;
		_neu._beta  = _beta+dphi-phif;
		_neu._kappa = _kappa;
		_neu._phi   = _neu._Phi + _neu._beta;
		if (Drvtvs) { // ------------- computation of derivatives ---------
			//_neu._Deriv=HepMatrix(5,5,1); whiy not?
			for (int i=1;i<6; i++) {_neu._Deriv(i,i)=1;};
			HepDouble cosbf = sqrt(1-sqr(sinbf));
			HepDouble sphif = yf / radf;
			_neu._Deriv(1,4) = rdphi / ( radf * cosbf );
			_neu._Deriv(1,5)=.5*rdphi*_Deriv(1,4);
			_neu._Deriv(2,3)=-rdphi*(1.+sqr(cotth));
			_neu._Deriv(2,4)=radi*cotth*sphif/cosbf;
			_neu._Deriv(2,5)=.5*rdphi*_Deriv(2,4);
			_neu._Deriv(4,4)=(radi*cosb)/(radf*cosbf);
			_neu._Deriv(4,5)=.5*rdphi*(1.+_Deriv(4,4));
			_neu._Der=_neu._Deriv * _Der; // chain rule
		};
	};
	
	_neu.useDelphi=true;
	_neu.useCMS=false;
	#ifdef GL
	_neu.Delphi2CMS();
	_neu.identifier="Propagated StateVector";
	_neu.short_id="PS";
	#endif
	return _neu;
};

StateVector::StateVector()
{
	_Phi=0; _theta=0; _beta=0; _kappa=0; _z=0; _r=0; 
	_a0=0; _psi=0; _kappa_cms=0; _phi=0;
	useCMS=false;
	useDelphi=false;
	_Deriv=HepSymMatrix(5,0);
	_Der=HepSymMatrix(5,1);
	#ifdef GL
	identifier="zero StateVector";
	short_id="0";
	#endif
};

HepDouble StateVector::eff_thick()
{
	/*
	HepDouble coeff;
	if (_el.Type() == CYL ) {
		coeff = 1 / sin(theta());
	} else {
		coeff = 1 / cos(theta());
	};
	return _el.Thick() * coeff;
	*/
	if (_el.Type() == CYL ) {
		return _el.Thick() / sin(theta());
	} else {
		return _el.Thick() / cos(theta());
	};
};
		
#ifdef GL
string StateVector::getInfo ( int what )
{
	switch (what) {
		case 0: // CMS
		case 1:
		case 2:
		case 4: // all cms
			return identifierCMS();
		case 3: // delphi
			return identifierDelphi();
		case 99:
			return identifierCMS();
		default:
			return "error";
	};
};

string getsvInfo (int a)
{
	if (hlx_shows < 5)
		return sv_list[a].getInfo(hlx_shows);
};

void StateVector::drawonDis ( SoSelection *root , string idnt, HepDouble lngth ) {
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
		sv_list.push_back(*this);
		identifier=idnt+" !illegal!!";
		addChild ( root, newHelix , (getsvInfo) , sv_list.size()-1 , true );
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
		sv_list.push_back(*this);
		identifier=idnt;
		short_id=idnt;
		addChild ( root, newHelix , (getsvInfo) , sv_list.size()-1 , true );
	};
};

void StateVector::drawonCyl ( SoSelection *root , string idnt, HepDouble lngth ) {
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
		sv_list.push_back(*this);
		identifier=idnt+" !illegal!!";
		addChild ( root, newHelix , (getsvInfo) , sv_list.size()-1 , true );
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
		sv_list.push_back(*this);
		identifier=idnt;
		short_id=idnt;
		addChild ( root, newHelix , (getsvInfo) , sv_list.size()-1 , true );
	};
};

void StateVector::GLdraw ( SoSelection *root ,HepDouble r, HepDouble g,\
		HepDouble b, HepDouble lngth ) {
	GLinitMaterial ( root, r , g , b );
	if ( _el.Type() == CYL ) {
		drawonCyl ( root , "Track.", lngth );
	} else {
		drawonDis ( root , "Track.", lngth );
	};
};

void StateVector::GLdraw ( SoSelection *root , HepDouble lngth ) {
//	GLinitMaterial ( root );
	if ( _el.Type() == CYL ) {
		drawonCyl ( root , "Track.", lngth );
	} else {
		drawonDis ( root , "Track.", lngth );
	};
};

	
void StateVector::GLdraw ( SoSelection *root , const string & idn , HepDouble r, \
		HepDouble g, HepDouble b , HepDouble lngth ) {
	GLinitMaterial ( root, r , g , b );
	if ( _el.Type() == CYL ) {
		drawonCyl ( root , idn, lngth );
	} else {
		drawonDis ( root , idn, lngth );
	};
};

void StateVector::GLdraw ( SoSelection *root , const string & idn , HepDouble lngth ) {
//	GLinitMaterial ( root );
	if ( _el.Type() == CYL ) {
		drawonCyl ( root , idn, lngth );
	} else {
		drawonDis ( root , idn, lngth );
	};
};

void StateVector::GLinitMaterial ( SoSelection *root ) const {
	SoMaterial *myMaterial = new SoMaterial;
	myMaterial->diffuseColor.setValue(0.2, 0.3, 0.7);
	myMaterial->emissiveColor.setValue(0.2, 0.3, 0.7);
	addChild(root,myMaterial);
};

void StateVector::GLinitMaterial ( SoSelection *root, HepDouble red, \
		HepDouble green , HepDouble blue ) const {
	SoMaterial *myMaterial = new SoMaterial;
	myMaterial->diffuseColor.setValue ( red, green, blue );
	myMaterial->emissiveColor.setValue( red, green, blue );
	addChild(root,myMaterial);
};
#endif

StateVector::StateVector ( const DetElement & el, const HepDouble & Phi , 
    const HepDouble & z , const HepDouble & theta, const HepDouble & beta, 
    const HepDouble & kappa )
{
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
		_phi=_Phi+_beta;
	} else {
		ERROR ("not implemented");
	};
	_kappa_cms=0;
	useCMS=false;
	useDelphi=true;
	_Deriv=HepSymMatrix(5,0);
	_Der=HepSymMatrix(5,1);
	#ifdef GL
	identifier="'initialized' StateVector";
	short_id="0";
	#endif
};

StateVector::StateVector( const DetElement & elm )
{
	_el=elm;
//	_s=0;
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
	identifier="zero StateVector on DetEl";
	short_id="0";
	#endif
};

HepVector StateVector::XY0() {
	// coc: center of circle
	// wrong
	//	return R()*cos(Phi()) - 1/ kappa() * sin (beta() + Phi())
	HepDouble r=_el.R(); // FIXME
	HepVector ret= HepVector ( 2, 0);
	HepDouble Rh = 1 / fabs ( kappa() );
	HepDouble CentreX = ( - a0() + 1 / kappa() ) * sin(psi());
	HepDouble CentreY = ( a0() - 1 / kappa() ) * cos (psi());
	HepDouble CosDeltaPhi = ( sqr(Rh)+sqr(CentreX) + sqr(CentreY) - sqr(r)) / \
		(2 * Rh * norm (CentreX, CentreY));
	HepDouble RPrime = sqrt ( 2*sqr(Rh)*(1-CosDeltaPhi));
	HepDouble PhiPrime=psi()-0.5*sign(kappa()) * acos(CosDeltaPhi);
	HepDouble xprime = RPrime * cos (PhiPrime);
	HepDouble yprime = RPrime * sin (PhiPrime);
	ret(1)= -1 * a0() * sin (psi() )+ xprime;
	ret(2)= a0() * cos ( psi() )+ yprime;
	return ret;
};

HepDouble StateVector::a0()
{
	if (!(useCMS)) Delphi2CMS();
	return _a0;
};

HepDouble StateVector::psi()
{
	if (!(useCMS)) Delphi2CMS();
	return _psi;
};

HepDouble StateVector::kappa_cms()
{
	if (!(useCMS)) Delphi2CMS();
	return _kappa_cms;
};

void StateVector::writeParam() { 
	if ( _el.Type() == CYL ) {
		cout << Phi() << " " << z() << " " << theta()
			<< " " << beta() << " " << kappa() << endl;
	} else {
		cout << Phi() << " " << R() << " " << theta()
			<< " " << phi() << " " << kappa() << endl;
	};
};
HepDouble StateVector::Radius() const
{
	if (useDelphi) {
		return 1 / _kappa;
	} else if (useCMS) {
		return 1 / fabs (_kappa_cms);
	};
	ERROR("cannot determine which coordinate system is used here.");
	exit (-1);
};

void StateVector::Randomize( interval a, interval b, interval c, interval d, 
	   interval e )
{
	/*
	#ifdef DEBUG
	assert ( a.min>0 );
	#endif
	*/
	HepDouble B=_el.Bz();
	HepDouble pT = random_interval ( a ); // GeV
	#ifdef GL
	// char name[109];
	identifier="Randomized StateVector";
	short_id="RS";
	#endif
	if ( _el.Type() == DIS ) {
		/*
		 * Randomize disc
		 */
		_phi=random_interval( c );
		if (_phi < 0) _Phi+=_Phi;
		_theta=random_interval ( b );
		_x=random_interval ( d ); 
		_y=random_interval ( e ); 

		_Phi = myatan ( _y , _x );
		_kappa=CNV * B / pT; //  / sin ( _theta );
		if (_Phi < 0) _Phi+=_Phi;
		_beta= fmod ( _phi - _Phi + M_PI_2 , M_2PI ) - M_PI_2;
	
		// okay. this is a little different here. we
		// have to propagate  now to the first element.
		// (above values are for z=0)
		DetElement dummy_el ( DIS , 0 , 0 , 0 , 0 );
		StateVector sv ( dummy_el );
		sv.setDelphi ( Phi() , R() , theta() , phi() , kappa() );
		StateVector nv=sv.PropagateDis ( Element(), false );
		// now take nv values as initial values.
		
		_x= nv.x() ; _y= nv.y(); _theta=nv.theta(); _phi=nv.phi();
		_Phi = nv.Phi(); _beta= _phi - _Phi;
		_kappa=nv.kappa();
		if ( _Phi < 0 ) _Phi += M_2PI;
		if ( _Phi > M_2PI ) _Phi -= M_2PI;
		_beta= fmod ( _phi - _Phi + M_PI_2 , M_2PI ) - M_PI_2;
		_kappa=CNV * B / pT; //  / sin (_theta );
	} else {
		/*
		 * Randomize Cylinder
		 */
		// FIXME this is currently specific to Br=0, Bphi=0
		// 
		_Phi=  random_interval ( c );
		_theta=random_interval ( b );
		_z =   random_interval ( d );
		_beta= random_interval ( e ); // is this good?
		
		_kappa = CNV * B / pT; // / sin ( _theta );
		_phi=_beta+_Phi;
		// _z = 0; _Phi=3; _theta=1;
	};
	
	// HepDouble eps = pow(10,-15); // precision of the calculations
	useDelphi=true;
	useCMS=false;
	MESSAGE(10,"[random]   Random TrueVector generated.\n     '-> "+DelphiFixedString());
};

void StateVector::Randomize( const string & filename )
{
	if (filename == "") {
		MESSAGE(5,"[random]   Using the hard coded intervals. Fix this.");
		Randomize();
		return;
	};
	static bool ronce=false;
	if (!ronce) {
		MESSAGE(5,"[random]   intervals: "+filename);
		ronce=true;
	};
	ifstream fin(filename.c_str());
	if (!fin) file_not_found(filename);
	char word[MLC];
	string line;
	string nxtword, min_val, max_val;
	HepDouble mmin, mmax;
	interval a,b,c,d,e;
	a.min=0; a.max=0;
	b.min=0; b.max=0;
	c.min=0; c.max=0;
	d.min=0; d.max=0;
	e.min=0; e.max=0;
	/// get the description from the file.
	/// parse init file!!
	while (fin.getline(word,MLC)) {
		line=word;
		line=del_leading_spaces(line);
		line=line.substr(0,line.find("#")); // delete comments
		nxtword=line.substr(0,line.find(" "));
		min_val="";
		if (nxtword.size()) {
			line=del_leading_spaces(line);
			line.erase(0,line.find_first_of(" 	")+1);
			line=del_leading_spaces(line);
			min_val=line.substr(0,line.find(" "));
			max_val=min_val;
			if (min_val.size()) {
				line=del_leading_spaces(line);
				line.erase(0,line.find_first_of(" 	")+1);
				line=del_leading_spaces(line);
				max_val=line.substr(0,line.find(" "));
				max_val=del_leading_spaces(max_val);
				if (max_val=="") max_val=min_val;
			};
			mmin=atof(min_val.c_str());
			mmax=atof(max_val.c_str());
			if (nxtword == "p" ) {
					ERROR("we dont use p anymore. Use pT!");
					a.min=mmin; a.max=mmax;
			} else if (nxtword == "pT" ) {
					a.min=mmin; a.max=mmax;
			} else if (nxtword == "theta" ) {
					b.min=mmin; b.max=mmax;
			} else if (nxtword == "phi"  && Element().Type()==DIS  ) {
					c.min=mmin; c.max=mmax;
			} else if (nxtword == "x"   && Element().Type()==DIS  ) {
					d.min=mmin; d.max=mmax;
			} else if (nxtword == "y"   && Element().Type()==DIS  ) {
					e.min=mmin; e.max=mmax;
			} else if (nxtword == "Phi"   && Element().Type()==CYL  ) {
					c.min=mmin; c.max=mmax;
			} else if (nxtword == "z"   && Element().Type()==CYL  ) {
					d.min=mmin; d.max=mmax;
			} else if (nxtword == "beta" && Element().Type()==CYL  ) {
					e.min=mmin; e.max=mmax;
			} else {
					ERROR("unknown argument in init file");
			};
		};
	};
	Randomize( a, b, c, d, e);
}

void StateVector::Randomize( )
{
	interval p, theta;
	p.min=1.5;
	p.max=p.min;
	
	if ( _el.Type() == DIS ) {
		interval phi, x, y;
		phi.min=M_PI_4; phi.max=phi.min;
		theta.min=.5; theta.max=theta.min;
		x.min=.02; x.max=x.min;
		y.min=.01; y.max=y.min;
		Randomize ( p , theta, phi, x , y );
	} else {
		interval Phi, z, beta;
		Phi.min=3; Phi.max=3;
		theta.min=1; theta.max=1;
		z.min=0; z.max=0;
		beta.min = 0; beta.max = 0;
		Randomize ( p , theta, Phi, z , beta );
	};
};

HepDouble StateVector::zVth_fak() const {
	if ( Element().Type() == CYL ) {
		return pow ( sin ( theta() ) , -4);
	} else {
		return pow ( cos ( theta() ) , -4);
	};
};

string StateVector::DelphiFixedString()
{
	char out[255];
	if ( _el.Type()==CYL ) {
		snprintf(out,255,"Phi=%.4f z=%.4f theta=%.4f beta=%.4f kappa=%.4f",
				Phi(), z(), theta(), beta(), kappa() );
	} else {
		snprintf(out,255,"Phi=%.4f R=%.4f theta=%.4f phi=%.4f kappa=%.4f",
				Phi(), R(), theta(), phi(), kappa() );
	};
	return (string) out;
};
#ifdef GL
string StateVector::DelphiString()
{
	char out[255];
	if ( _el.Type()==CYL ) {
		snprintf(out,255,"Phi=%.*f z=%.*f theta=%.*f beta=%.*f kappa=%.*f",
				matrix_prec+3, Phi(), 
				matrix_prec+3, z(), 
				matrix_prec+3, theta(), 
				matrix_prec+3, beta(), 
				matrix_prec+3, kappa() );
	} else {
		snprintf(out,255,"Phi=%.*f R=%.*f theta=%.*f phi=%.*f kappa=%.*f",
				matrix_prec+3, Phi(), 
				matrix_prec+3, R(), 
				matrix_prec+3, theta(), 
				matrix_prec+3, phi(), 
				matrix_prec+3, kappa() );
	};
	return (string) out;
};

string StateVector::CMSString()
{
	char out[255];
	snprintf(out,255,"a0=%.*f z=%.*f theta=%.*f psi=%.*f kappa=%.*f",
			matrix_prec+3, a0(),
			matrix_prec+3, z(),
			matrix_prec+3, theta(),
			matrix_prec+3, psi(),
			matrix_prec+3, kappa());
	return (string) out;
};

string StateVector::DelphiShortString()
{
	char out[255];
	if ( _el.Type()==CYL ) {
		snprintf(out,255,"%.*f %.*f %.*f %.*f %.*f", 
				matrix_prec+3, Phi(), 
				matrix_prec+3, z(), 
				matrix_prec+3, theta(), 
				matrix_prec+3, beta(), 
				matrix_prec+3, kappa() );
	} else {
		snprintf(out,255,"%.*f %.*f %.*f %.*f %.*f",
				matrix_prec+3, Phi(), 
				matrix_prec+3, R(), 
				matrix_prec+3, theta(), 
				matrix_prec+3, phi(), 
				matrix_prec+3, kappa() );
	};
	return (string) out;
};


string StateVector::identifierDelphi()
{
	return identifier+" (Delphi)\n      "+DelphiString();
};

string StateVector::identifierCMS()
{
	return identifier+" (CMS)\n      "+CMSString();
};

string StateVector::identifierDeriv()
{
	return any2str ( identifier + " (Deriv, Delphi)\n      Der=", 
			_Deriv , matrix_prec+3 );
};

string StateVector::identifierDer()
{
	return any2str ( identifier + " (Der, Delphi)\n      Der=", 
			_Der, matrix_prec+3 );
};

string StateVector::CMSShortString()
{
	char out[255];
	snprintf(out,255,"%.*f %.*f %.*f %.*f %.*f",
			matrix_prec+3, a0(),
			matrix_prec+3, z(),
			matrix_prec+3, theta(),
			matrix_prec+3, psi(),
			matrix_prec+3, kappa());
	return (string) out;
};
#endif
