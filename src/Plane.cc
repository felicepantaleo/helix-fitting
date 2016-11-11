/// \file

#include <math.h>
#include "Plane.h"
#include "Detector.h"
#include "Track.h"
#include "any2str.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/DiagMatrix.h"
#include "CLHEP/config/CLHEP.h"
#include "HepMoreOperators.h"
#include <vector>
#ifdef SOLARIS
#include <ieeefp.h>
#endif


void Plane::write () const 
{ 
	cout << " c= " << c << endl << " n=" << n;
};

void RealPlane::write () const 
{ 
	Plane::write();
	cout << " Cov" << cov;
};

HepDouble Plane::Radius() const {
	return sqrt( (1-sqr(n[2])- 4*c*n[2] ) / (4*sqr(n[2]) ));
};

HepDouble Plane::X0() const {
	return -n[0]/2/n[2];
};
	
HepDouble Plane::Y0() const {
	return -n[1]/2/n[2];
}

HepVector Plane::s( Detector *det ) {
	return s( det, a0(), psi(), kappa() );
}

HepVector Plane::s( Detector *det, HepDouble la0 , HepDouble lpsi,
	HepDouble lkappa, HepDouble rref, HepDouble phiref ) {
	HepVector mys = HepVector ( det->Num , 0);
	HepVector coc = HepVector ( 2, 0); // coc=center of circle
	HepVector k1 = HepVector (2,0);
	HepVector kk = HepVector (2,0);
	coc(1) = ( - la0 + 1 / lkappa ) * sin(lpsi);
	coc(2) = ( la0 - 1 / lkappa ) * cos (lpsi);
	k1=cartXY ( det->Element[0].R(), rref, phiref ) - coc;
	HepDouble rci = 1 / fabs(lkappa);
	HepDouble rcineg2 = sqr ( lkappa );
	 
	for (int k=2; k<= (det->Num); k++ ) {
		kk=cartXY ( det->Element[k-1].R(), rref, phiref ) - coc; 
		mys(k)= acos( dot( k1, kk) * rcineg2 ) * rci ;
		#ifdef DEBUG
		if (!finite(mys(k))) 
			WARNING("Warning: arclength at " + any2str(k) + ": "\
		   + any2str (mys(k)) + " is not well defined.");
		#endif 
	}; 
	return mys;
}

HepVector Plane::s( Detector *det, HepDouble la0 , HepDouble lpsi,
	HepDouble lkappa ) {
	HepVector mys = HepVector ( det->Num , 0);
	HepVector coc = HepVector ( 2, 0); // coc=center of circle
	HepVector k1 = HepVector (2,0);
	HepVector kk = HepVector (2,0);
	coc(1) = ( - la0 + 1 / lkappa ) * sin(lpsi);
	coc(2) = ( la0 - 1 / lkappa ) * cos (lpsi);
	k1=cartXY_old ( det->Element[0].R(), la0, lpsi,lkappa) - coc;
	HepDouble rci = 1 / fabs(lkappa);
	HepDouble rcineg2 = sqr ( lkappa );
	 
	for (int k=2; k<= (det->Num); k++ ) {
		kk=cartXY_old (det->Element[k-1].R(),la0,lpsi, lkappa) - coc; 
		mys(k)= acos( dot( k1, kk) * rcineg2 ) * rci ;
		#ifdef DEBUG
		if (!finite(mys(k))) 
			WARNING("Warning: arclength at " + any2str(k) + ": "\
		   + any2str (mys(k)) + " is not well defined.");
		#endif 
	}; 
	return mys;
}

HepVector Plane::s( vector < HepDouble > rs ) {
	return s( rs, a0(), psi(), kappa() );
}


HepVector Plane::s( vector <HepDouble> rs, HepDouble rref, HepDouble phiref) {
	return s( rs, a0(), psi(), kappa(), rref, phiref );
}

HepMatrix Plane::RadVec ( vector < HepDouble > rs )
{
	return RadVec ( rs, a0(), psi(), kappa() );
};

HepMatrix Plane::RadVec ( vector <HepDouble> rs, HepDouble la0, HepDouble lpsi,
		HepDouble lkappa, HepDouble rref, HepDouble phiref )
{
	HepMatrix ret ( 2 , rs.size() );
	HepVector coc = HepVector ( 2, 0); // coc=center of circle
	HepVector k1 = HepVector (2,0);
	HepVector kk = HepVector (2,0);
	coc(1) = ( - la0 + 1 / lkappa ) * sin(lpsi);
	coc(2) = ( la0 - 1 / lkappa ) * cos (lpsi);
	for (unsigned int k=1; k<= rs.size(); k++ ) {
		ret.sub(1,k, cartXY ( rs[k-1], rref, phiref ) - coc);
		// ret.sub(1,k, cartXY_old ( rs[k-1], la0, lpsi, lkappa ) - coc);
	};
	return ret;
}	

HepMatrix Plane::RadVec ( vector < HepDouble > rs, HepDouble la0, HepDouble lpsi,
		HepDouble lkappa )
{
	HepMatrix ret ( 2 , rs.size() );
	HepVector coc = HepVector ( 2, 0); // coc=center of circle
	HepVector k1 = HepVector (2,0);
	HepVector kk = HepVector (2,0);
	coc(1) = ( - la0 + 1 / lkappa ) * sin(lpsi);
	coc(2) = ( la0 - 1 / lkappa ) * cos (lpsi);
	for (unsigned int k=1; k<= rs.size(); k++ ) {
		ret.sub(1,k, cartXY_old ( rs[k-1],la0,lpsi, lkappa ) - coc);
	};
	return ret;
}	

HepVector Plane::RadVec ( HepDouble rs )
{
	return RadVec ( rs, a0(), psi(), kappa() );
};

HepVector Plane::RadVec ( HepDouble rs, HepDouble la0, HepDouble lpsi,
		HepDouble lkappa, HepDouble rref, HepDouble phiref  )
{
	HepVector ret ( 2 , 0 );
	HepVector coc = HepVector ( 2, 0); // coc=center of circle
	HepVector k1 = HepVector (2,0);
	HepVector kk = HepVector (2,0);
	coc(1) = ( - la0 + 1 / lkappa ) * sin(lpsi);
	coc(2) = ( la0 - 1 / lkappa ) * cos (lpsi);
	ret=cartXY ( rs, rref, phiref ) - coc;
//	ret=cartXY_old ( rs, la0, lpsi, lkappa ) - coc;
	return ret;
}	

HepVector Plane::RadVec ( HepDouble rs, HepDouble la0, HepDouble lpsi,
		HepDouble lkappa )
{
	HepVector ret ( 2 , 0 );
	HepVector coc = HepVector ( 2, 0); // coc=center of circle
	HepVector k1 = HepVector (2,0);
	HepVector kk = HepVector (2,0);
	coc(1) = ( - la0 + 1 / lkappa ) * sin(lpsi);
	coc(2) = ( la0 - 1 / lkappa ) * cos (lpsi);
	ret=cartXY_old ( rs,la0,lpsi, lkappa ) - coc;
	return ret;
}	

HepVector Plane::s( vector < HepDouble > rs, HepDouble la0 , HepDouble lpsi, 
	HepDouble lkappa, HepDouble rref, HepDouble phiref ) {
	HepVector mys = HepVector ( rs.size() , 0);
	HepVector coc = HepVector ( 2, 0); // coc=center of circle
	HepVector k1 = HepVector (2,0);
	HepVector kk = HepVector (2,0);
	coc(1) = ( - la0 + 1 / lkappa ) * sin(lpsi);
	coc(2) = ( la0 - 1 / lkappa ) * cos (lpsi);
	k1=cartXY ( rs[0], rref, phiref ) - coc;
	// k1=cartXY_old ( rs[0], la0, lpsi, lkappa ) - coc;
	HepDouble rci = 1 / fabs(lkappa);
	HepDouble rcineg2 = sqr ( lkappa );
	
//	cerr << rs[0] << rs[1] << endl;
	for (unsigned int k=2; k<= (rs.size()); k++ ) {
		kk=cartXY ( rs[k-1], rref, phiref ) - coc;
		// kk=cartXY_old ( rs[k-1], la0, lpsi, lkappa ) - coc;
		mys(k)= acos( dot( k1, kk) * rcineg2 ) * rci;
		#ifdef DEBUG
		if (!(mys(k) < 200 )) {
			WARNING("Warning: arclength at " + any2str(k) + ": "\
		   + any2str (mys(k)) + " is not well defined.");
		};
		#endif
	};
	return mys;
}

HepVector Plane::s( vector < HepDouble > rs, HepDouble la0 , HepDouble lpsi, 
	HepDouble lkappa) {
	HepVector mys = HepVector ( rs.size() , 0);
	HepVector coc = HepVector ( 2, 0); // coc=center of circle
	HepVector k1 = HepVector (2,0);
	HepVector kk = HepVector (2,0);
	coc(1) = ( - la0 + 1 / lkappa ) * sin(lpsi);
	coc(2) = ( la0 - 1 / lkappa ) * cos (lpsi);
	k1=cartXY_old ( rs[0], la0, lpsi,lkappa ) - coc;
	HepDouble rci = 1 / fabs(lkappa);
	HepDouble rcineg2 = sqr ( lkappa );
	
//	cerr << rs[0] << rs[1] << endl;
	for (unsigned int k=2; k<= (rs.size()); k++ ) {
		kk=cartXY_old ( rs[k-1],la0,lpsi, lkappa ) - coc;
		mys(k)= acos( dot( k1, kk) * rcineg2 ) * rci;
		#ifdef DEBUG
		if (!(mys(k) < 200 )) {
			WARNING("Warning: arclength at " + any2str(k) + ": "\
		   + any2str (mys(k)) + " is not well defined.");
		};
		#endif
	};
	return mys;
}

void Plane::calculateCMS()
{
	HepDouble XY0_norm = norm(X0(),Y0());
	HepDouble nX0 = X0() / XY0_norm;
	HepDouble nY0 = Y0() / XY0_norm;
	HepDouble RVecX = Radius()*nX0;
	HepDouble RVecY = Radius()*nY0;
	HepDouble a0VecX = X0() - RVecX;
	HepDouble a0VecY = Y0() - RVecY;
	HepDouble XVec0X = myr * cos ( myphi );
	HepDouble XVec0Y = myr * sin ( myphi );
	HepDouble XVecPrimeX = XVec0X  - a0VecX;
	HepDouble XVecPrimeY = XVec0Y  - a0VecY;
	HepDouble norm_PsiVec = sign ( XVecPrimeX * RVecY -XVecPrimeY*RVecX) * \
		1 / norm(RVecX, RVecY);
	HepDouble nPsiVecX = norm_PsiVec * RVecY;
	HepDouble nPsiVecY = - norm_PsiVec * RVecX;
	int sign_a0 = - sign (  a0VecX * nPsiVecY - a0VecY * nPsiVecX );
	mya0=  sign_a0 *  norm(a0VecX,a0VecY);
	mypsi=atan2 ( nPsiVecY, nPsiVecX );
	int signKappa=sign(RVecX*nPsiVecY-RVecY*nPsiVecX);
	mykappa= signKappa  / Radius();
	CMScalculated=true;
};

Plane::Plane () {
	c=0; n=HepVector(3,0); z=0; theta=0;
	CMScalculated=false;
};


HepDouble Plane::a0() {
	if (!CMScalculated) {
		calculateCMS();
	};
	return mya0;
}

HepDouble Plane::psi() {
	if (!CMScalculated) {
		calculateCMS();
	};
	return mypsi;
}

HepDouble Plane::kappa() {
	if (!CMScalculated) {
		calculateCMS();
	};
	return mykappa;
};

/// new algorithm.
HepVector Plane::cartXY ( HepDouble r, HepDouble rref, 
		HepDouble phiref )
{
	HepMatrix ret = HepMatrix (2,1,0);
	HepDouble a;
	HepDouble s;
	HepDouble x=rref*cos(phiref);
	HepDouble y=rref*sin(phiref);
	a= c+ n(3)* sqr(r);
	HepDouble b=sqr(n(1))+sqr(n(2));
	s= sqrt(b*sqr(r)-sqr(a));
	HepDouble x1=(-n(1)*a+n(2)*s)/b;
	HepDouble x2=(-n(1)*a-n(2)*s)/b;
	HepDouble y1=(-n(2)*a-n(1)*s)/b;
	HepDouble y2=(-n(2)*a+n(1)*s)/b;
	signed int sign=1;
	if ( sqr ( x-x1)+sqr(y-y1) > sqr (x-x2)+sqr(y-y2))
		sign=-1;
	ret(1,1)=(-n(1)*a+sign*n(2)*s)/b;
	ret(2,1)=(-n(2)*a-sign*n(1)*s)/b;
	return ret;
};

/// new algorithm.
HepVector Plane::cartXY ( HepVector r, HepDouble rref, 
		HepDouble phiref )
{
	HepDouble sz=r.num_row();
	#ifdef DEBUG
		if (sz < 2) {
			WARNING("cartXY called with a strange HepVector");
		};
	#endif
	/*
	#ifdef SOLARIS // Arg!!
	if (sz > 20 ) {
		ERROR("solaris array problem");
		exit (-1);
	};
	HepMatrix ret = HepMatrix (2,20,0);
	HepVector a = HepVector ( 20 , 0 );
	HepVector s = HepVector ( 20 , 0 );
	#else
	*/
	HepMatrix ret = HepMatrix (2,sz,0);
	HepVector a = HepVector ( sz , 0 );
	HepVector s = HepVector ( sz , 0 );
	// #endif /* SOLARIS */
	HepDouble x=rref*cos(phiref);
	HepDouble y=rref*sin(phiref);
	a= c+ n(3)* sqr(r);
	HepDouble b=sqr(n(1))+sqr(n(2));
	s= sqrt(b*sqr(r)-sqr(a));
	HepDouble x1=(-n(1)*a(1)+n(2)*s(1))/b;
	HepDouble x2=(-n(1)*a(1)-n(2)*s(1))/b;
	HepDouble y1=(-n(2)*a(1)-n(1)*s(1))/b;
	HepDouble y2=(-n(2)*a(1)+n(1)*s(1))/b;
	// signed int sign=1;
	HepDouble signn2, signn1;
	if ( sqr ( x-x1)+sqr(y-y1) > sqr (x-x2)+sqr(y-y2)) {
			// sign=-1;
		signn2=-1*n(2)/b;
		signn1=1*n(1)/b;
	} else {
		signn2=1*n(2)/b;
		signn1=-1*n(1)/b;
	};
	for ( int i=1; i <=sz ; i++ ) {
//		ret(1,i)=(-n(1)*a(i)+sign*n(2)*s(i))/b;
//		ret(2,i)=(-n(2)*a(i)-sign*n(1)*s(i))/b;
		ret(1,i)=(-n(1)*a(i)+signn2*s(i));
		ret(2,i)=(-n(2)*a(i)+signn1*s(i));
	};
	return ret;
};

/// old algorithm.
HepVector Plane::cartXY_old ( HepDouble r, HepDouble la0, HepDouble lpsi, \
	HepDouble lkappa ) 
{
	HepVector ret= HepVector ( 2, 0);
	HepDouble Rh = 1 / fabs ( lkappa );
	HepDouble CentreX = ( - la0 + 1 / lkappa ) * sin(lpsi);
	HepDouble CentreY = ( la0 - 1 / lkappa ) * cos (lpsi);
	HepDouble CosDeltaPhi = ( sqr(Rh)+sqr(CentreX) + sqr(CentreY) - sqr(r)) / \
		(2 * Rh * norm (CentreX, CentreY));
	HepDouble RPrime = sqrt ( 2*sqr(Rh)*(1-CosDeltaPhi));
	HepDouble PhiPrime=lpsi-.5*sign(lkappa) * acos(CosDeltaPhi);
	HepDouble xprime = RPrime * cos (PhiPrime);
	HepDouble yprime = RPrime * sin (PhiPrime);
	ret(1)= -1 * la0 * sin ( lpsi )+ xprime;
	ret(2)= la0 * cos ( lpsi )+ yprime;
	return ret;
};

HepVector Plane::cartXY_old ( HepDouble r )
{
	return cartXY_old(r,a0(),psi(),kappa());
};
