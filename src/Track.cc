// $Id: Track.cc,v 1.1 2001/11/16 10:29:57 wwalten Exp $
// looks stable

#include "Track.h"
#include "Hit.h"
#include "FullVector.h"
#include "StateVector.h"
#include "DetElement.h"
#include "HepMoreOperators.h"
#include <time.h>
#include <sys/times.h>
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/DiagMatrix.h"
#include "CLHEP/config/CLHEP.h"
#include <vector>
#include <string>
#include <stdio.h>

#ifdef GL
#include <Inventor/nodes/SoSelection.h>
#endif
extern void mycout ( HepSymMatrix x );
extern clock_t Timing[3];
extern clock_t userTiming[3];
extern clock_t sysTiming[3];
#ifdef DEBUG
extern clock_t lastTiming[3];
#endif
extern struct tms mytms;
extern void WatchStop( int i );
extern void WatchStart( int i );

void Track::write ()
{
  cout << endl << "*****************************" << endl;
  cout << "* New Track" << endl;
//  for (int i=0; i<Num; i++ )
//    cout << hits[i].RPhi() << " " << hits[i].z() << endl;
  #ifdef DEBUG
  truevec.writeCMS();
  truevec.writeDelphi();
  #endif
  rt[0].writeCMS();
  rt[0].writeDelphi();
};

void Track::write_raw ( char *file) const
{
  for (int i=0; i<Num; i++ ) {
    hits[i].write_raw(file);
  };
};

HepDouble ret_el ( const HepMatrix & m, int row, int col )
{
  return m(row,col);
};

HepDouble ret_els ( const HepSymMatrix & m, int row, int col )
{
  return m(row,col);
};

#ifdef GL
void Track::GLinitMaterial ( )
{
  hits[0].GLinitMaterial ( detector->getRoot() );
}

void Track::GLdraw ( )
{
  hits[0].GLinitMaterial ( detector->getRoot() );
  for (int i=0; i< Num; i++ ) {
    hits[i].GLdraw( detector->getRoot() );
  };
  #ifdef DEBUG
  GLdrawTrueTrack ( );
  GLdrawRefTrack ( );
  #endif
};
#endif

Track::Track ( Detector *det )
{
  // Detector, and the B field
  detector=det;
  Num=detector->Num+1;

  computed_cov=-1;
  computed_fv_cov=-1;
  rt_computed=false;
  R_pred_computed=false;
  last_fit_was_special=false;
};

Track & Track::operator = ( Track &r )
{
  return r;
};
void Track::ReferenceTrack ( bool comp_der )
{
  if (!(rt_computed)) {
    RealPlane p;
    ReferenceTrack ( &p, comp_der);
  };
};

void Track::setReferenceTrack ( StateVector sv, bool comp_der )
{
  sv.identifier="ReferenceTrack.";
  sv.short_id="RT";
  rt.clear();
  rt.reserve ( Num );
//  rts = HepVector ( Num, 0 );
  rt.push_back ( sv );
//  rts = plane->s( allR(),hits[0].R(),hits[0].Phi() );
  for ( int k=2; k<= Num; k++ ) {

    sv = sv.Propagate( detector->Element[k-1] , comp_der );
    rt.push_back ( sv );
  };

  rt_computed=comp_der;
  MESSAGE(10,"[ref]      Reference Track computed. Comp_Der="+any2str(comp_der));
};

// use this one if you want to keep the Plane object.
void Track::ReferenceTrack ( Plane *p, bool comp_der )
{
  WatchStart(4);
  if (!(rt_computed)) {
    StateVector blub( detector->Element[0] );
    WatchStart(3);
    initCov();
    WatchStop(3);
    *p=RiemannFit_Plane ( false , false );
    if (detector->Element[0].Type()==CYL) {
      blub.getPlane ( *p , detector->Element[0] );
    } else {
      blub.getPlane ( *p , detector->Element[0], R_pred[0] );
    };
    setReferenceTrack( blub , comp_der );
//    rts=p->s();
  };
  WatchStop(4);
};

StateVector Track::Simulate ( char scatter, char smear, string init_val_file )
{
  int j;
  // First we create all StateVectors
  vector < StateVector > vecs;
  vecs.reserve(Num);
  hits.reserve(Num);
  Hit onehit;
  #ifdef DEBUG
  static bool once=false;
  string sm,ss;
  switch (smear) {
    case 0:
      sm="[track]    not smearing out measurements.";
      break;
    case 2:
      sm="[track]    smearing out measurements deterministically.";
      break;
    case 1:
      sm="[track]    smearing out measurements randomly.";
      break;
    default:
      WARNING("Dunno how to smear");
      sm="[track]    smearing unknown." ;
  };
  switch (scatter) {
    case 0:
      ss="[track]    not scattering the vectors.";
      break;
    case 2:
      ss="[track]    scattering the vectors deterministically.";
      break;
    case 1:
      ss="[track]    scattering the vectors randomly.";
      break;
    default:
      WARNING("Dunno how to scatter");
      ss="[track]    scattering unknown." ;
  };
  if (!once) {
    MESSAGE(5,sm);
    MESSAGE(5,ss);
    once=true;
  };
  #endif

  vecs.push_back( StateVector ( detector->Element[0] ));
  vecs[0].Randomize ( init_val_file ); // we only respect Bz :-((
  #ifdef DEBUG
  truevec=vecs[0];
  truevec.identifier="True Track.";
  truevec.short_id="TT";
  #endif
  onehit= Hit(vecs[0]);
  if (smear>0) onehit.Smear( smear == 2 );
  hits.push_back ( onehit );
  for ( j=1; j < detector->Num ; j++) {
    if (scatter>0) vecs[j-1].Scatter ( scatter == 2 );
    vecs.push_back( vecs[j-1].Propagate (detector->Element[j] , false ) );
     // now we 'translate' into 'Hits':
    onehit = Hit(vecs[j]);
    if (smear>0) onehit.Smear( smear == 2 );
    hits.push_back ( onehit );
  };
  Num=j;
  return truevec;
};

void Track::initCov()
{
  _cov=HepSymMatrix(2*Num,0);
  // _cov_ms=HepSymMatrix(2*Num,0);
  for (int i=0; i< Num; i++ ) {
    _cov(i+1,i+1) = sqr( hits[i].SigRPhi() );
  };
  for (int i=Num; i<  2*Num; i++ ) {
    if ( detector->Element[i-Num].Type() == CYL ) {
      _cov(i+1,i+1) = sqr(hits[i-Num].Sigz());
    } else {
      _cov(i+1,i+1) = sqr(hits[i-Num].SigR());
    };
  };
  #ifdef DEBUG
  MESSAGE(15,"[ccov]     Covariance Matrix initialized.");
  #endif
};

void Track::initCov_special()
{
  _cov=HepSymMatrix(2*Num,0);
  for (int i=1; i<= Num; i++ ) {
    _cov(i+Num,i+Num) = sqr(hits[i-1].SigR());
    _cov(i,i) = VarRPhi_pred(i);
  };
  #ifdef DEBUG
  MESSAGE(15,"[ccov]     Covariance Matrix initialized (initCov_special).");
  #endif
};

HepSymMatrix Track::computeFVCov ( /* HepMatrix B */ )
{
  if (computed_cov != computed_fv_cov) {
    // we wont check much here.
    HepMatrix V=HepMatrix(2*Num,2*Num,0);
    HepMatrix A(2*Num,5);
    HepMatrix swap_row(1,2*Num,0);
    HepMatrix swap_col(2*Num,1,0);
    // _cov = NO_corrected.
    // for level 0 we do not want an NO_corrected.
    // version.
    //

    if ( computed_cov == 0 ) {
      // reshape
      for (int i=1;i<=Num;i++ ) {
        V(i*2-1,i*2-1)=hits[i-1].VarRPhi();
        V(i*2,i*2)=hits[i-1].Varz();
      };
    } else if ( computed_cov == 3 ) {
      for (int i=1;i<=Num;i++ ) {
        V(i*2-1,i*2-1)=hits[i-1].VarRPhi();
        V(i*2,i*2)=hits[i-1].Varz();
        for (int j=1;j<=Num;j++ ) {
          V(i*2-1,j*2-1)+=_cov_ms(i,j);
//          if (computed_cov==3) {
            V(i*2,j*2-1)+=_cov_ms(i+Num,j);
            V(i*2-1,j*2)+=_cov_ms(i,j+Num);
//          };
          V(i*2,j*2)=_cov(i+Num,j+Num);
        };
      };
    } else {
      // reshape
      for (int i=1;i<=Num;i++ ) {
        for (int j=1;j<=Num;j++ ) {
          V(i*2-1,j*2-1)=_cov(i,j);
          V(i*2,j*2)=_cov(i+Num,j+Num);
        };
      };
    };
    int j=0;
    for (int k=1;k<=Num;k++) {
      for (int l=1; l<=5; l++) {
        A(j+1,l)=detector->Element[k-1].R()*rt[k-1].Der()(1,l);
        A(j+2,l)=rt[k-1].Der()(2,l);
      };
      j+=2;
    };
    int fail=0;
    HepSymMatrix G(2*Num,0);
    G=symmetrize(V).inverse(fail);
    #ifdef DEBUG
    if (fail) { ERROR ("we failed inverting a matrix");};
    #endif
    _fv_cov=HepSymMatrix(5,0);
    _fv_cov=G.similarity(A.T());
    _fv_cov.invert(fail);
    #ifdef DEBUG
    if (fail) { ERROR("we failed inverting a matrix"); };
    #endif
    _mB=HepMatrix(5,2*Num,0);
    _mB=_fv_cov*A.T()*G;
    computed_fv_cov=computed_cov;
    #ifdef DEBUG
    MESSAGE(15,"[ccov]    Covariance Matrix of FullVector computed.");
    #endif
  };
  return _fv_cov;
};

HepSymMatrix Track::CovPhi() const
{
  HepSymMatrix covphi ( _cov.sub(1,Num) );
  for ( int i=0; i< Num; i++ ) {
    for ( int j=i; j< Num; j++ )
      covphi[i][j] /= (rt[i].R() * rt[j].R());
    for ( int j=Num; j< 2*Num; j++ )
      covphi[i][j] /= rt[i].R() ;
  };
  return covphi; // ( 1 , Num );
};

HepSymMatrix Track::CovMSPhi() const
{
  HepSymMatrix covphi ( _cov_ms.sub(1,Num) );
  for ( int i=0; i< Num; i++ ) {
    for ( int j=i; j< Num; j++ )
      covphi[i][j] /= (rt[i].R() * rt[j].R());
    for ( int j=Num; j< 2*Num; j++ )
      covphi[i][j] /= rt[i].R() ;
  };
  return covphi; // ( 1 , Num );
};

HepSymMatrix Track::CovMSPhiR() const
{
  HepSymMatrix covphi ( _cov_ms );
  for ( int i=0; i< Num; i++ ) {
    for ( int j=i; j< Num; j++ )
      covphi[i][j] /= (rt[i].R() * rt[j].R());
    for ( int j=Num; j< 2*Num; j++ )
      covphi[i][j] /= rt[i].R() ;
  };
  return covphi;
};

void Track::force_computeCov ( int level )
{
  WatchStart(3);
  initCov();
  Plane plane;
  computed_cov=level;
  last_fit_was_special=false;
  computeCovMS( &plane, level );
  _cov+=_cov_ms;
  correct_NO();
  WatchStop(3);
}


void Track::computeCov ( int level, bool init_special ,
    bool fullvar)
{
        WatchStart(3);
  switch (init_special) {
  case true: {
    last_fit_was_special=true;
    initCov_special();
    if (computed_cov!=level) {
      computed_cov=level;
      Plane plane;
//      ReferenceTrack( &plane, true  );
      computeCovMS( &plane, level );
    };
    if (_cov_ms.num_row()==0) {
      _cov_ms=HepSymMatrix(2*Num,0); // no MS
    } else {
      _cov+=_cov_ms;
    };
    if (rt.size() == (unsigned) Num ) {
      correct_NO( fullvar );
    } else {
      WARNING("Could not correct NO intercept (dont have an rt)");
    };
    break;
  };
  default:
    last_fit_was_special=false;
    if (computed_cov!=level) {
      initCov();
      computed_cov=level;
      Plane plane;
//      ReferenceTrack( &plane, true  );
      computeCovMS( &plane, level );
      _cov+=_cov_ms;
      // if (level>0) correct_NO(); // OK???
      correct_NO();
    } else if (computed_cov==level && last_fit_was_special) {
      initCov();
      computed_cov=level;
      Plane plane;
//      ReferenceTrack( &plane, true  );
      computeCovMS( &plane, level );
      _cov+=_cov_ms;
    };
  };
  #ifdef DEBUG
  MESSAGE(15,"[ccov]     Covariance Matrix computed.");
  #endif
  WatchStop(3);
};

void Track::computeCov ( int level, bool init_special ,
    bool fullvar, HepDouble kappa_cms )
{
        WatchStart(3);
  switch (init_special) {
  case true: {
    initCov_special();
    if (computed_cov!=level) {
      computed_cov=level;
      Plane plane;
//      ReferenceTrack( &plane, true  );
      computeCovMS( &plane, level );
    };
    if (_cov_ms.num_row()==0) {
      _cov_ms=HepSymMatrix(2*Num,0); // no MS
    } else {
      _cov+=_cov_ms;
    };
    if (rt.size() == (unsigned) Num ) {
      correct_NO( fullvar, kappa_cms );
    } else {
      WARNING("Could not correct NO intercept (dont have an rt)");
    };
    break;
  };
  default:
    if (computed_cov!=level) {
      initCov();
      computed_cov=level;
      Plane plane;
//      ReferenceTrack( &plane, true  );
      computeCovMS( &plane, level );
      _cov+=_cov_ms;
      correct_NO();
    };
  };
  #ifdef DEBUG
  MESSAGE(15,"[ccov]     Covariance Matrix computed.");
  #endif
  WatchStop(3);
};

void Track::computeCovMS ( Plane *plane, int level )
{
  _cov_ms = HepSymMatrix ( 2*Num , 0 );
  if ( level < 1 ) return;
  // HepSymMatrix V_MS ( 2 * Num , 0 );
  HepVector sigms ( Num , 0 ); //, sigms2 ( Num , 0 );
  HepDouble X = rt[0].eff_thick();
  sigms(1)=.015 / rt[0].pr() * sqrt(X) * (1+.038 * log(X)) ;

  switch (level) {
  case 1: { //   ------------- level 1 ------------------------
        // for barrel only (currently)
        // calculation of the s
    HepVector mys = HepVector ( Num , 0);
    HepVector sigms2s(Num,0), sigms2z(Num,0);
    if ( detector->Element[0].Type() == CYL ) {
//      mys=plane->s( detector, rt[0].a0(), rt[0].psi(),rt[0].kappa_cms());
      mys=rts;
    } else {
      for ( int i=0; i< Num ; i++ ) {
        mys[i]= rt[i].Element().Z();
      };
    };

    HepDouble ss, ssz;

    // create a local copy of some variables in the loop.
    vector < HepDouble > zVth,sts;
    zVth.reserve ( Num );
    sts.reserve ( Num );
    for ( int i=0; i< Num ; i++ ) {
      zVth.push_back( rt[i].zVth_fak() );
      sts.push_back( sqr ( sin ( rt[i].theta()) ));
    };
    sigms2s(1)=sqr(sigms(1)) /sts[0] ;
    sigms2z(1)=sqr(sigms(1)) *zVth[0] ;

    for ( int k=2; k<= Num ; k++ ) {
      X = rt[k-1].eff_thick();
      sigms(k)=.015 / rt[k-1].pr() * sqrt(X) * (1+.038 * log(X)) ;
      sigms2s(k)=sqr( sigms(k)) / sts[k-1];
      sigms2z(k)=sqr( sigms(k)) * zVth[k-1];
      for ( int l=k; l<=Num ; l++ ) {
        ss=0; ssz=0;
        for ( int j=1; j<= k ; j++ ) {
          ss+= ( R_pred[k-1]-R_pred[j-1] ) * ( R_pred[l-1] -
            R_pred[j-1] ) * sigms2s(j); // /
//            sts[j-1];
          ssz+= (mys[k-1]-mys[j-1]) * (mys[l-1] - mys[j-1]) *
            sigms2z(j); // * zVth[j-1];
        };
        _cov_ms.fast(l,k)=ss;
        _cov_ms.fast(l+Num,k+Num)=ssz;

      };
    };
//    _cov_ms=V_MS;
    break;
  };
  /*
  case 99: { //   ------------- level 99 -- used to be level 1  -
         // new implementation is faster.
    HepVector Vy ( Num , 0 );
    HepVector Cthy ( Num , 0 );
    HepVector Vth ( Num , 0 );
    HepVector zVy ( Num , 0 );
    HepVector zCthy ( Num , 0 );
    HepVector zVth ( Num , 0 );
    HepDouble Delta, zDelta;

    Vth(1)=sqr(sigms(1)) / sqr ( sin ( rt[0].theta ()) );
    zVth(1)=sqr(sigms(1)) * rt[0].zVth_fak();

    // calculation of the s
    HepVector mys = HepVector ( Num , 0);
    if ( detector->Element[0].Type() == CYL ) {
      mys=plane->s( detector, rt[0].a0(), rt[0].psi(),rt[0].kappa_cms());
    };
    for ( int k=2; k<= Num; k++ ) {
      Delta= rt[k-1].R() - rt[k-2].R();
      if ( rt[k-1].Element().Type() == CYL ) {
        zDelta= mys[k-1] - mys[k-2];
      } else {
        zDelta= rt[k-1].Element().Z() - rt[k-2].Element().Z();
      };
      // zDelta= s(k)-s(k-1);
      X = rt[k-1].eff_thick();
      sigms(k)=.015 / rt[k-1].pr() * sqrt(X) * (1+.038 * log(X)) ;

      Vy(k) = Vy(k-1)+2*Delta*Cthy(k-1)+sqr(Delta)*Vth(k-1);
      Cthy(k) = Cthy(k-1) + Delta * Vth(k-1);
      Vth(k)=Vth(k-1)+sqr(sigms(k)) / sqr ( sin (rt[k-1].theta()) );

      zVy(k) = zVy(k-1)+2*zDelta*zCthy(k-1)+sqr(zDelta)*zVth(k-1);
      zCthy(k) = zCthy(k-1) + zDelta * zVth(k-1);
      zVth(k)=zVth(k-1)+sqr(sigms(k)) * rt[k-1].zVth_fak();

      _cov_ms.fast(k,k)=Vy(k);
      _cov_ms.fast(Num+k,Num+k)=zVy(k);

      for ( int j=k+1; j<= Num; j++ ) {
        Delta= rt[j-1].R() - rt[k-1].R();
        if ( rt[j-1].Element().Type() == CYL ) {
          zDelta= mys[j-1] - mys[k-1];
        } else {
          zDelta= rt[j-1].Element().Z() - rt[k-1].Element().Z();
        };

        _cov_ms.fast(j,k) = Vy(k)+Delta*Cthy(k);
        _cov_ms.fast(Num+j,Num+k) = zVy(k)+zDelta*zCthy(k);
      };
    };
//    _cov_ms=V_MS;
    break;
  };
  */
  case 2: {
    // ------------------------- level 2 ---------------
    HepDouble fac = -1 * CNV * detector->Element[0].Bz() / \
      rt[0].pr() * cos (rt[0].theta()) / sqr(sin (rt[0].theta()));
    #ifdef SOLARIS
    HepSymMatrix Cj[20]; // *argh*! XXX
    if (Num>20) { // FIXME
      ERROR("Argh. Quo vadis, Solaris CC?!");
      exit (-1);
    };
    #else
    HepSymMatrix Cj[Num];
    #endif
    // derivative d_kappa / d_theta

    Cj[0]=HepSymMatrix ( 5, 0);
    Cj[0].fast(3,3)=sqr(sigms(1));
    Cj[0].fast(4,4)=sqr(sigms(1)) / sqr(sin(rt[0].theta()));
    Cj[0].fast(5,3)=fac*sqr(sigms(1));
    Cj[0].fast(5,5)=sqr(fac)*sqr(sigms(1));

    HepDouble radk, radj, sigmsk2;
    for ( int k=2; k<= Num; k++ ) {
      radk=rt[k-1].R();
      Cj[k-1]=HepSymMatrix ( 5, 0);
      X = rt[k-1].eff_thick();
      sigms(k)=.015 / rt[k-1].pr() * sqrt(X) * (1+.038 * log(X)) ;
      sigmsk2=sqr(sigms(k));
      Cj[k-1]=Cj[k-2].similarity ( rt[k-1].Deriv() );
      _cov_ms.fast(k,k)=Cj[k-1](1,1)*sqr(radk);
      _cov_ms.fast(k+Num,k+Num)=Cj[k-1](2,2);
      Cj[k-1].fast(3,3)+=sigmsk2;
      Cj[k-1].fast(4,4)+=sigmsk2 / sqr (sin (rt[k-1].theta())) ;
      Cj[k-1].fast(5,3)+=fac*sigmsk2;
      Cj[k-1].fast(5,5)+=sqr(fac)*sigmsk2;
      HepMatrix Temp=Cj[k-1];

      for ( int j=k+1; j<= Num; j++ ) {
        radj=rt[j-1].R();
        Temp=rt[j-1].Deriv() * Temp;
        _cov_ms.fast(j,k)=Temp(1,1) * radk * radj;
        _cov_ms.fast(j+Num,k+Num)=Temp(2,2) ;
      };
    };
//    _cov_ms=V_MS;
    break;
  };
  case 3: {
    // ------------------------- level 3 ---------------
    HepDouble fac = -1 * CNV * detector->Element[0].Bz() / \
      rt[0].pr() * cos (rt[0].theta()) / sqr(sin (rt[0].theta()));
    #ifdef SOLARIS
    HepSymMatrix Cj[20]; // *argh*! XXX
    if (Num>20) { // FIXME
      ERROR("Argh. Quo vadis, Solaris CC?!");
      exit (-1);
    };
    #else
    HepSymMatrix Cj[Num];
    #endif
    // derivative d_kappa / d_theta

    Cj[0]=HepSymMatrix ( 5, 0);
    Cj[0].fast(3,3)=sqr(sigms(1));
    Cj[0].fast(4,4)=sqr(sigms(1)) / sqr(sin(rt[0].theta()));
    Cj[0].fast(5,3)=fac*sqr(sigms(1));
    Cj[0].fast(5,5)=sqr(fac)*sqr(sigms(1));

    HepDouble radk, radj, sigmsk2;
    for ( int k=2; k<= Num; k++ ) {
      radk=rt[k-1].R();
      Cj[k-1]=HepSymMatrix ( 5, 0);
      X = rt[k-1].eff_thick();
      sigms(k)=.015 / rt[k-1].pr() * sqrt(X) * (1+.038 * log(X)) ;
      sigmsk2=sqr(sigms(k));
      Cj[k-1]=Cj[k-2].similarity ( rt[k-1].Deriv() );
      _cov_ms.fast(k,k)=Cj[k-1](1,1)*sqr(radk);
      _cov_ms.fast(k+Num,k)=Cj[k-1](1,2)*radk;

      _cov_ms.fast(k+Num,k+Num)=Cj[k-1](2,2);
      Cj[k-1].fast(3,3)+=sigmsk2;
      Cj[k-1].fast(4,4)+=sigmsk2 / sqr (sin (rt[k-1].theta())) ;
      Cj[k-1].fast(5,3)+=fac*sigmsk2;
      Cj[k-1].fast(5,5)+=sqr(fac)*sigmsk2;
      HepMatrix Temp=Cj[k-1];

      for ( int j=k+1; j<= Num; j++ ) {
        radj=rt[j-1].R();
        Temp=rt[j-1].Deriv() * Temp;
        _cov_ms.fast(j,k)=Temp(1,1) * radk * radj;
        _cov_ms.fast(j+Num,k+Num)=Temp(2,2) ;
        _cov_ms.fast(k+Num,j)=Temp(1,2)*radj ;
        _cov_ms.fast(j+Num,k)=Temp(2,1)*radk ;
      };
    };
//    _cov_ms=V_MS;
    break;
  };
  default: {
      ERROR("unknown Cov level." + any2str(level));
  };
  };
  #ifdef DEBUG
  MESSAGE(15,"[ccov]    MS part of covariance Matrix computed.");
  #endif
};

inline void Track::correct_NO ( bool fullVar, HepDouble kappa_cms )
{
  // Correction for nonorthogonal intercept -
  // should actually be part of StateVectors?! XXX
  // fullVar = true: we use the 'full' Covariance Matrix for
  // the correction.
  HepSymMatrix covrphi=_cov.sub(1,Num);
  if ( detector->Element[0].Type() == CYL ) {
    HepDouble coeff = - 0.5 * kappa_cms;
    HepSymMatrix B( Num , 1 );
    for ( int i=0; i< Num; i++)
  //RF  B[i][i]=cos ( asin( coeff * detector->Element[i].R() ));
      B[i][i]=sqrt(1-sqr( coeff * detector->Element[i].R() ));
    _cov.sub(1,covrphi.similarity ( B ));
  } else {
    HepSymMatrix si ( Num, 0), co ( Num , 0);
    HepDouble coeff = .5 * kappa_cms;
    for ( int i=0; i < Num ; i++ ) {
      si[i][i]= coeff * R_pred[i];
      co[i][i]=sqrt(1-sqr(si[i][i]));
    };
    if ( fullVar ) {
      _cov.sub(1,covrphi.similarity(co)+
          _cov.sub(Num+1,2*Num).similarity(si));
    } else {
      _cov.sub(1,covrphi.similarity(co)+
          _cov_ms.sub(Num+1,2*Num).similarity(si));
    };
  };
}

inline void Track::correct_NO ( bool fullvar ) {
  correct_NO ( fullvar , rt[0].kappa_cms() );
};
inline void Track::correct_NO () {
  correct_NO ( false , rt[0].kappa_cms() );
};

inline void Track::force_initR_pred() {
  VarRPhi_pred=HepVector (Num,0);
  R_pred=HepVector (Num , 0);
  for ( int i=0; i< Num ; i++ ) {
    R_pred[i]=hits[i].R();
    VarRPhi_pred[i]=hits[i].VarRPhi();
  };
//  R_pred_computed=true;
};

inline void Track::initR_pred() {
  if (!R_pred_computed) {
    VarRPhi_pred=HepVector (Num,0);
    R_pred=HepVector (Num , 0);
    for ( int i=0; i< Num ; i++ ) {
      R_pred[i]=hits[i].R();
      VarRPhi_pred[i]=hits[i].VarRPhi();
      R_pred_computed=true;
    };
  };
};

FullVector Track::RiemannFit ( int level, int iterations )
{
  #ifdef DEBUG
  if (level>2) WARNING("RiemannFit level 2,3 are alike");
  #endif
  FullVector fv ( detector->Element[0]) ;
  RealPlane plane;
  bool MS = (level!=0);
  char idn[50];
  snprintf(idn,50,"Track, Riemann Fit [MS matrix level %d].",level);
  fv.identifier=idn;
  fv.short_id="R"+any2str(level);

  WatchStart( 2 );
  initR_pred();
  ReferenceTrack ( level > 1 );
  /*
  if ( level > 1 ) {
    ReferenceTrack ( true );
  } else {
    ReferenceTrack ( false );
  };
  */
  if ( detector->Element[0].Type()==DIS ) {
    computeCov ( level , true, true );
    HepDouble old_kappa=rt[0].kappa_cms();
    for ( int i=0; i <= iterations ; i++ ) {
      computeCov ( level , true, i==0, old_kappa );
      if (i==0) {
        force_initR_pred();
      } else if (i==1) {
        WatchStop (2); // i==1 ? stop the watch.
      };

      //RF plane=RiemannFit_Plane ( MS , true );
      plane=RiemannFit_Plane ( MS , i==iterations );
      fv.getPlane ( plane , detector->Element[0], R_pred(1) );
      old_kappa=fv.kappa_cms();
      if (iterations==0) WatchStop(2);
    };
  } else {
    computeCov ( level , false );
    WatchStop ( 2 );
    plane=RiemannFit_Plane ( MS , true );
    fv.getPlane ( plane , detector->Element[0] );
  };
  #ifdef GL
  fv.track_cov=_cov;
  HepSymMatrix C=HepSymMatrix(5,0); // computeFVCov();
  fv.setDelphiCov( C );
  #endif

  if ( detector->Element[0].Type()==DIS && level < 2 ) {
    // oops :-((
    rt_computed=false;
    R_pred_computed=false;
//    last_fit_was_special=false;
  };
  return fv;
};

FullVector Track::GlobalDis ( int level )
{
  FullVector fv ( detector->Element[0] );
  HepVector mr(2*Num,0);
  HepVector m(2*Num,0);
  HepVector dm(2*Num,0);
  HepSymMatrix V(2*Num,0), G(2*Num,0);
  HepMatrix A(2*Num,5,0);
  HepSymMatrix C(5,0);
  HepMatrix B(5,2*Num,0);

  int fail;
  WatchStart( 2 );
  initR_pred();
  ReferenceTrack ( true );
  force_computeCov ( level );
  WatchStop( 2 );
  for (int k=1;k<=Num;k++) {
    mr(k)=rt[k-1].Phi();
    mr(k+Num)=rt[k-1].R();
    m(k)=hits[k-1].Phi();
    m(k+Num)=hits[k-1].R();
    dm(k)=m(k)-mr(k);
    dm(k+Num)=m(k+Num)-mr(k+Num);
    if (dm(k) > M_PI) dm(k)-=M_2PI;
    if (dm(k) < -M_PI) dm(k)+=M_2PI;
    if (dm(k+Num) > M_PI) dm(k+Num)-=M_2PI;
    if (dm(k+Num) < -M_PI) dm(k+Num)+=M_2PI;
    V(k,k)=hits[k-1].VarPhi();
    V(k+Num,k+Num)=hits[k-1].VarR();
    for (int l=1; l<=5; l++) {
      A(k,l)=rt[k-1].Der()(1,l);
      A(k+Num,l)=rt[k-1].Der()(2,l);
    };
  };
  if (level) {
    V+=CovMSPhiR();
  };
  G=V.inverse ( fail );
  #ifdef DEBUG
  if (fail) {
    ERROR("we failed inverting a matrix. V(1:3,1:3)="+any2str(V.sub(1,3)));
  };
  #endif
  C=G.similarity ( A.T() );
  C.invert ( fail );
  #ifdef DEBUG
  if (fail) { ERROR("we failed inverting a matrix"); };
  #endif
  B=C*A.T()*G;
  HepMatrix dparam(1,5,0);
  dparam=( B*(dm)).T();      // b^T a^T = (ab)^T
/*  if ( dparam(1,1) < 0 ) dparam(1,1)+=M_2PI;
  }; */
  fv.setDelphi( rt[0].Phi()  + dparam(1,1) ,
            rt[0].R()     + dparam(1,2) ,
          rt[0].theta() + dparam(1,3),
          rt[0].phi()  + dparam(1,4) ,
          rt[0].kappa() + dparam(1,5)  );
  #ifdef GL
  fv.track_cov=_cov; // this is best to just copy.
  fv.setDelphiCov( C );
  #endif
  return fv;
};

FullVector Track::GlobalCyl ( int level )
{
  FullVector fv ( detector->Element[0] );
  HepVector mr(2*Num,0);
  HepVector m(2*Num,0);

  int j=0;
  WatchStart( 2 );
  initR_pred();
  ReferenceTrack ( true );
  computeCov ( level );
  WatchStop ( 2 );
  for (int k=1;k<=Num;k++) {
      mr(j+1)=detector->Element[k-1].R()*rt[k-1].Phi();
      mr(j+2)=rt[k-1].z();
      m(j+1)=hits[k-1].RPhi();
      m(j+2)=hits[k-1].z();
    j+=2;
  };
  HepSymMatrix C=computeFVCov( );
  HepMatrix dparam(1,5,0);
  dparam=( _mB*(m-mr)).T();      // b^T a^T = (ab)^T
  fv.setDelphi( rt[0].Phi()   + dparam(1,1),
            rt[0].z()     + dparam(1,2),
          rt[0].theta() + dparam(1,3),
          rt[0].beta()  + dparam(1,4),
          rt[0].kappa() + dparam(1,5) );
  #ifdef GL
  fv.track_cov=_cov; // this is best to just copy.
  fv.setDelphiCov( C );
  #endif
  return fv;
};

FullVector Track::GlobalFit ( int level )
{
  FullVector fv; //  ( detector->Element[0] );
  if ( detector->Element[0].Type() == DIS ) {
    fv=GlobalDis ( level );
  } else {
    fv=GlobalCyl ( level );
  };
  char idn[50];
  snprintf(idn,50,"Track, Global Fit [MS matrix level %d].",level);
  fv.identifier=idn;
  fv.short_id="G"+any2str(level);
  return fv;
};

FullVector Track::KalmanFilter ( int level )
{
  #ifdef DEBUG
  if (level>1) WARNING("Kalman filter level 1,2,3 are alike");
  #endif
  FullVector fv; //  ( detector->Element[0] );
  if ( detector->Element[0].Type() == DIS ) {
    fv=KalmanDis ( level );
  } else {
    fv=KalmanCyl ( level );
  };
  char idn[50];
  snprintf(idn,50,"Track, Kalman Filter [MS matrix level %d].",level);
  fv.identifier=idn;
  fv.short_id="K"+any2str(level);
  return fv;
};

FullVector Track::KalmanCyl ( int level )
{
  FullVector fv ( detector->Element[0] );
  HepVector mk(2,0);
  HepSymMatrix Vk(2,0), Ek(2,0);
  HepMatrix Kk(5,2);
  HepMatrix Hk (2,5,0);

  //we do this togenerate some initial values.
  // FIXME can we avoid calling initCov() implicitly?

  WatchStart( 2 );
  Plane plane;
  initR_pred();
  ReferenceTrack ( true );
  WatchStop( 2 );

  HepSymMatrix Ce(5,1);
  HepMatrix parame(1,5,0);
  HepMatrix paramf(1,5,0);
  HepSymMatrix Cf(5,0);
  HepMatrix Fk(5,5,0);
  HepDiagMatrix eye(5,1);
  HepDouble X;
  HepDouble sigms, sigmsk2;
  HepDouble fac;
  int fail;
  Ce*=1000;
  for (int k=Num; k>=1; k--) {
    mk(1)=hits[k-1].RPhi()-detector->Element[k-1].R()*rt[k-1].Phi();
    mk(2)=hits[k-1].z()-rt[k-1].z();
    Hk(1,1)=detector->Element[k-1].R();
    Hk(2,2)=1;
    Vk(1,1)=hits[k-1].VarRPhi();
    Vk(2,2)=hits[k-1].Varz();
    Ek=symmetrize( Vk+Hk*Ce*Hk.T());
    Kk=Ce*Hk.T()*Ek.inverse( fail );
    #ifdef DEBUG
    if (fail) { ERROR("inversion of matrix failed"); };
    #endif
    paramf=(parame.T()+Kk*(mk-Hk*parame.T())).T();
    Cf=symmetrize((eye-Kk*Hk)*Ce);
    if (k==1) break;
    Fk=rt[k-1].Deriv().inverse(fail);
    #ifdef DEBUG
    if (fail) { ERROR("inversion of matrix failed"); };
    #endif
    parame=paramf*Fk.T();
    Ce=Cf.similarity(Fk);
    if (level>0) {
      X = rt[k-2].eff_thick();
      sigms=.015 / rt[k-2].pr() * sqrt(X) * (1+.038 * log(X)) ;
      sigmsk2=sqr(sigms);
      fac = -1 * CNV * detector->Element[k-1].Bz() / \
        rt[k-1].pr() * cos(rt[k-1].theta())/sqr(sin (rt[k-1].theta()));
      Ce(3,3)+=sigmsk2;
      Ce(4,4)+=sigmsk2 / sqr (sin (rt[k-1].theta())) ;
      Ce(3,5)+=fac*sigmsk2;
      Ce(5,5)+=sqr(fac)*sigmsk2;
    };
  };
  paramf+=rt[0].Delphi().T();

  #ifdef GL
  fv.setDelphiCov( Cf );
  fv.track_cov=_cov;
  #endif
  fv.setDelphi ( paramf(1,1) ,
             paramf(1,2) ,
           paramf(1,3) ,
           paramf(1,4) ,
           paramf(1,5) );
  return fv;
};

#ifdef DEBUG
void Track::ReferenceTrackFromTrueTrack ( bool comp_der )
{
  if (!(rt_computed)) {
    StateVector blub( truevec );
    #ifdef GL
    blub.identifier="ReferenceTrack from 'truevec'";
    blub.short_id="RT";
    #endif
    setReferenceTrack( blub , false );
  };
};

#ifdef GL
void Track::GLdrawTrueTrack ()
{
  truevec.GLdraw( detector->getRoot(),
    .9 , .9 , .3 , detector->arclength() );
};

void Track::GLdrawRefTrack ()
{
  rt[0].GLdraw( detector->getRoot(),
    .4 , .9 , .4 , detector->arclength() );
};
#endif /* GL */
#endif /* DEBUG */

HepVector Track::allz() const
{
  HepVector ret(Num,0);
  for ( int i=0; i< Num ; i++ ) {
    ret[i]=hits[i].z();
  };
  return ret;
};

FullVector Track::KalmanDis ( int level )
{
  FullVector fv ( detector->Element[0] );
  HepVector mk(2,0);
  HepSymMatrix Ek(2,0);
  HepDiagMatrix Vk(2,0);
  HepMatrix Kk(5,2);
  HepMatrix Hk (2,5,0);

  WatchStart( 2 );
  initR_pred();
  ReferenceTrack( true );
//  computeCov ( 0 );
  WatchStop ( 2 );

  HepSymMatrix Ce(5,1);
  HepMatrix parame(1,5,0); //, parami(1,5,0);
  HepMatrix paramf(1,5,0);
  HepSymMatrix Cf(5,0);
  HepMatrix Fk(5,5,0);
  HepDiagMatrix eye(5,1);
  HepDouble X;
  HepDouble sigms, sigmsk2;
  HepDouble fac; // , p, pT;
  int fail;
  Ce*=1000;
  Hk(1,1)=1;
  Hk(2,2)=1;

  for (int k=Num; k>=1; k--) {
    mk(1)=hits[k-1].Phi()-rt[k-1].Phi();
    mk(2)=hits[k-1].R()-rt[k-1].R();
    if (mk(1) > M_PI) mk(1)-=M_2PI;
    if (mk(1) < -M_PI) mk(1)+=M_2PI;
    Vk(1,1)=hits[k-1].VarPhi();
    Vk(2,2)=hits[k-1].VarR();
    Ek= Vk+Ce.similarity(Hk);
    Kk=Ce*Hk.T()*Ek.inverse( fail );
    #ifdef DEBUG
    if (fail) { ERROR("inversion of matrix failed"); };
    #endif
    paramf=(parame.T()+Kk*(mk-Hk*parame.T())).T();
    Cf=symmetrize((eye-Kk*Hk)*Ce);
    if (k==1) break;
    Fk=rt[k-1].Deriv().inverse(fail);
    #ifdef DEBUG
    if (fail) { ERROR("inversion of matrix failed"); };
    #endif
    parame=paramf*Fk.T();
    Ce=Cf.similarity(Fk);
    if (level>0) {
      fac = -1 * CNV * detector->Element[k-1].Bz() / \
        rt[k-1].pr() * cos(rt[k-1].theta())/sqr(sin (rt[k-1].theta()));

      X = rt[k-2].eff_thick();
      sigms=.015 / rt[k-2].pr() * sqrt(X) * (1+.038 * log(X)) ;
      sigmsk2=sqr(sigms);
      Ce(3,3)+=sigmsk2;
      Ce(4,4)+=sigmsk2 / sqr (sin (rt[k-2].theta()));
      Ce(3,5)+=fac*sigmsk2;
      Ce(5,5)+=sqr(fac) * sigmsk2;
    };
  };
  paramf+=rt[0].Delphi().T();

  #ifdef GL
  fv.setDelphiCov( Cf );
  fv.track_cov=_cov;
  #endif
  fv.setDelphi ( paramf(1,1) ,
             paramf(1,2) ,
           paramf(1,3) ,
           paramf(1,4) ,
           paramf(1,5) );
  return fv;
};


vector < HepDouble > Track::allR() const
{
  vector <HepDouble > res;
  res.reserve (Num);
  if ( detector->Element[0].Type() == CYL ) {
    for ( unsigned int i=0; i< hits.size() ; i++ ) {
      res.push_back ( hits[i].R() );
    };
  } else {
    for ( unsigned int i=0; i< hits.size() ; i++ ) {
      res.push_back ( R_pred[i] );
    };
  };
  return res;
};

RealPlane Track::RiemannFit_Plane ( bool MS , bool comp_cov )
{
  RealPlane plane;

  HepVector Phi(Num,0); // A vector with all phis
  HepVector co(Num, 0); // cos(phi)
  HepVector si(Num, 0); // sin(phi)

  HepVector x(Num, 0);  // R*cos(phi)
  HepVector y(Num, 0);  // R*sin(phi)
  HepVector z(Num, 0);  // R^2

  HepVector w(Num, 0);  // (Covariance R * Phi)^-1

  HepVector Xg( 3 , 0);      // Center of gravity
  HepMatrix X( 3 , Num , 0); // Center of gravity
  HepSymMatrix D( 3 , 0 );
  HepMatrix V( 1 , Num , 0 );

  HepDouble sum_w=0, varc=0, k=0;
  plane.setRPhiPoint ( R_pred[0], hits[0].Phi() );

  for ( int i=0; i< Num; i++) {
    HepDouble tmp_phi=hits[i].Phi();
    HepDouble tmp_r;
//    if (detector->Element[0].Type() == CYL ) {
//      tmp_r=hits[i].R();
//    } else {
    tmp_r=R_pred[i];
//    };
    Phi[i]=tmp_phi;

    co[i]=cos(tmp_phi);
    si[i]=sin(tmp_phi);

//RF    x[i]=tmp_r * cos(tmp_phi);
//RF    y[i]=tmp_r * sin(tmp_phi);
    x[i]=tmp_r * co[i];
    y[i]=tmp_r * si[i];
    z[i]=tmp_r * tmp_r;
  };

  // used for noMS
  HepVector vRPhi = HepVector(Num, 0); // Covariance R * Phi

  // used for MS
  HepSymMatrix W(Num,0);
  HepMatrix id(Num,0);

  if (MS) {
    int fail;
    W=CovRPhi().inverse(fail);
    #ifdef DEBUG
    if (fail) {
      ERROR("inversion of covrphi matrix failed.\nCovRPhi=" +
          any2str (CovRPhi()) );
    };
    #endif
    for ( int i=1; i<=Num; i++) {
      w(i)+=W.fast(i,i);
      sum_w+=W.fast(i,i);
      for ( int j=i+1; j <=Num ; j++) {
/*        w(i)+=W(j,i);
        sum_w+=W(j,i);*/
        w(i)+=W.fast(j,i);
        w(j)+=W.fast(j,i);
        sum_w+=2*W.fast(j,i);
      };
    };
  } else {
    for ( int i=0; i< Num; i++) {
      vRPhi[i]=CovRPhi()[i][i];
      w[i]= 1 / vRPhi[i];
      sum_w+=w[i]; // probably still the fastest method to get the sum.
      W[i][i]=w[i];
    };
  };
  w/=sum_w;
  Xg(1)=dot(x,w);
  Xg(2)=dot(y,w);
  Xg(3)=dot(z,w);

  for (int i=0; i < Num; i++ ) {
    X[0][i]=(x[i]-Xg(1));
    X[1][i]=(y[i]-Xg(2));
    X[2][i]=(z[i]-Xg(3));
  };
  D=W.similarity(X);

  V=diagonalize(&D);
  k=minimum_sym (D);

  plane.n = V.sub(1,3,k,k);
  plane.c = - dot ( plane.n ,  Xg);

  // now fit in z and theta. z-fit
  HepMatrix Ar(Num, 2, 0);
  HepSymMatrix Gr(Num,0);
  HepSymMatrix Cr(2,0);
  HepVector fitz(2,0);
  HepVector my_z(Num,0);
  int fail;
//  HepVector mys = HepVector ( Num , 0);
  int a;
  rts=plane.s( allR(),hits[0].R(),hits[0].Phi() );

  if (detector->Element[0].Type() == CYL ) { // XXX
    Gr=Covz().inverse(fail);
    for (int k=1; k<= (Num); k++ ) {
      Ar(k,2)=rts(k);
      Ar(k,1)=1;
      my_z(k)=hits[k-1].z();
    };

    Cr = Gr.similarityT ( Ar );
    Cr.invert( a );
    #ifdef DEBUG
    if (a) ERROR ("failed to invert matrix");
    #endif
    fitz = Cr* Ar.T() * Gr * my_z;
    plane.z= fitz(1);
    plane.theta= M_PI_2 - atan ( fitz(2) );
  } else { // Disc
    Gr=CovR().inverse(fail);
    for (int k=1; k<= (Num); k++ ) {
      my_z(k)=hits[k-1].z();
      Ar(k,2)=my_z(k);
      Ar(k,1)=1;
    };

    Cr = Gr.similarityT ( Ar );
    Cr.invert( a );
    #ifdef DEBUG
    if (a) ERROR ("failed to invert matrix");
    #endif
    fitz = Cr* Ar.T() * Gr * rts;
    plane.theta= atan ( fitz(2) );
    // plane.theta= atan ( fitz(2) );
    plane.z= hits[0].z(); // -hits[0].z();

    HepVector s1(Num,0);
    s1= allz() * fitz(2) + fitz(1);
    HepVector Phi1(Num,0);
    HepVector rvec=plane.RadVec ( R_pred[0] );
     Phi1= atan2 ( rvec(2) , rvec(1)) + s1 * fabs ( plane.kappa() );
    HepVector x_new = plane.X0() + cos (Phi1) / fabs ( plane.kappa() );
    HepVector y_new = plane.Y0() + sin (Phi1) / fabs ( plane.kappa() );
    R_pred = norm ( x_new , y_new );
    for ( int i=0; i< Num ; i++ )
      VarRPhi_pred[i] = sqr ( R_pred[i] ) * hits[i].VarPhi(); // * VarPhi
  };

  // -------------< compute covariance matrices >---------------
  if (comp_cov) {
    HepDouble eigensI=D[k][k];
    HepMatrix Covn( 3 , 3, 0);
    HepVector corr_cn( 3 , 0 );

    HepSymMatrix covrcg( 2 , 0);

    for (int i=0; i < 3; i++ ) {
      if (i!=k) {
        Covn += eigensI * D[i][i] * pow( D[i][i] - eigensI, - 2 ) * \
                V.sub (1,3,i+1,i+1) * V.sub(1,3,i+1,i+1).T() ;
      };
      Covn/=(Num-5);
    };
    if (MS) {
      HepSymMatrix vx (Num,0);
      HepSymMatrix vy (Num,0);
      HepMatrix cxy (Num,Num,0);
      vx=symmetrize( CovRPhi()* si * si.T() );
      vy=symmetrize( CovRPhi()* co * co.T() );
      cxy=CovRPhi() * si * co.T();
      covrcg[0][0] = dot ( w , vx * w );
      covrcg[1][1] = dot ( w , vy * w );
      covrcg[1][0] = dot ( w , cxy * w );
    } else {
      HepVector vx(Num,0);
      HepVector vy(Num,0);
      HepVector cxy(Num,0);
      for (int i=0; i < Num; i++ ) {
        HepDouble swi=sqr( w[i] );
        vx[i]  =  vRPhi[i] * sqr ( si[i] );
        vy[i]  =  vRPhi[i] * sqr ( co[i] );
        cxy[i] = -vRPhi[i] * si[i] * co[i];
        covrcg[0][0] += swi /* sqr ( w[i] )*/ * vx[i];
        covrcg[1][1] += swi /* sqr ( w[i] )*/ * vy[i];
        covrcg[1][0] += swi /* sqr ( w[i] )*/ * cxy[i];
      };

    };
    varc = sqr ( plane.n[0] ) * covrcg[0][0] + // n' C n
           2 * plane.n[0] * plane.n[1] * covrcg[0][1] +
           sqr ( plane.n[1] ) * covrcg[1][1] +
           covrcg[0][0] * Covn[0][0] + // sum(sum(Cn x C ))
           covrcg[1][1] * Covn[1][1] +
           2 * covrcg[1][0] * Covn[1][0] +
           dot(Xg, (Covn * Xg)); // Xg' . Cn . Xg

    corr_cn = - Covn * Xg;

    //  Copy. It wont get any faster.
    /*
    plane.cov.fast(1,1)= varc;
    plane.cov.fast(2,2)= Covn[0][0];
    plane.cov.fast(3,2)= Covn[1][0];
    plane.cov.fast(3,3)= Covn[1][1];
    plane.cov.fast(4,2)= Covn[2][0];
    plane.cov.fast(4,3)= Covn[2][1];
    plane.cov.fast(4,4)= Covn[2][2];
    plane.cov.fast(2,1)= corr_cn[0];
    plane.cov.fast(3,1)= corr_cn[1];
    plane.cov.fast(4,1)= corr_cn[2];
    */
    plane.cov[0][0] = varc;
    plane.cov[1][1] = Covn[0][0];
    plane.cov[2][1] = Covn[1][0];
    plane.cov[2][2] = Covn[1][1];
    plane.cov[3][1] = Covn[2][0];
    plane.cov[3][2] = Covn[2][1];
    plane.cov[3][3] = Covn[2][2];
    plane.cov[1][0] = corr_cn[0];
    plane.cov[2][0] = corr_cn[1];
    plane.cov[3][0] = corr_cn[2];

  };
  return plane;
};
