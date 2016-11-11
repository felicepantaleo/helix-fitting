/*-----------------------------Hepvis----------------------------------------*/
/*                                                                           */
/* Node:             SoHelicalTrack                                          */
/* Description:      Represents a helical trajectory with axis along z       */
/* Author:           Joe Boudreau Nov 7 1996                                 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#include <assert.h>
#include <math.h>
#include <Inventor/SbBox.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoSeparator.h>
#include "hepvis/SoHelicalTrack.hh"
// This statement is required
SO_NODE_SOURCE(SoHelicalTrack)

// Constructor
SoHelicalTrack::SoHelicalTrack() {
  // This statement is required
  SO_NODE_CONSTRUCTOR(SoHelicalTrack);

  // Data fields are initialized like this:
  SO_NODE_ADD_FIELD(inverseRadius,       (0.0));
  SO_NODE_ADD_FIELD(cotTheta,            (0.0));
  SO_NODE_ADD_FIELD(phi0,                (0.0));
  SO_NODE_ADD_FIELD(d0,                  (0.0));
  SO_NODE_ADD_FIELD(z0,                  (0.0));
  SO_NODE_ADD_FIELD(s0,                  (0.0));
  SO_NODE_ADD_FIELD(s1,                  (1.0));
  SO_NODE_ADD_FIELD(alternateRep,        (NULL));
  children = new SoChildList(this);
}

// Destructor
SoHelicalTrack::~SoHelicalTrack() {
}


// initClass
void SoHelicalTrack::initClass(){
  // This statement is required.
  SO_NODE_INIT_CLASS(SoHelicalTrack,SoShape,"Shape");
}


// GLRrender
void SoHelicalTrack::GLRender(SoGLRenderAction *action) {
  if (!shouldGLRender(action)) return;
  if (children->getLength() == 0) generateChildren();
  if ((inverseRadius.getValue()!=myInverseRadius) ||
      (cotTheta.getValue()     !=myCotTheta)      ||
      (phi0.getValue()         !=myPhi0)          ||
      (d0.getValue()           !=myD0)            ||
      (z0.getValue()           !=myZ0)            ||
      (s0.getValue()           !=myS0)            ||
      (s1.getValue()           !=myS1)      ) {
    updateChildren();
  }
  children->traverse(action);  
  myInverseRadius=inverseRadius.getValue();
  myCotTheta=cotTheta.getValue();
  myPhi0=phi0.getValue();
  myD0  = d0.getValue();
  myZ0  = z0.getValue();
  myS0  = s0.getValue();
  myS1  = s1.getValue();

}

// generatePrimitives
void SoHelicalTrack::generatePrimitives(SoAction *action) {
  if (children->getLength()==0) generateChildren();
  updateChildren();
  children->traverse(action);
}

// getChildren
SoChildList *SoHelicalTrack::getChildren() const {
  return children;
}


// computeBBox
void SoHelicalTrack::computeBBox(SoAction *, SbBox3f &box, SbVec3f &center) {
  double WHOPPING=100000000.0;
  SbVec3f min,max;
  SbVec3f A(0,0,0),B(0,0,0),B0(WHOPPING,WHOPPING,WHOPPING),B1(-WHOPPING,-WHOPPING,-WHOPPING);
  for (int i=0;i<13;i++) {
    A=getTrackPosition((s1.getValue()-s0.getValue())*(i)/13.0 + s0.getValue());
    B+=A;
    for (int j=0;j<3;j++) {
      if (A[j]<B0[j]) B0[j]=A[j];
      if (A[j]<B0[j]) B0[j]=A[j];
      if (A[j]>B1[j]) B1[j]=A[j];
    }
    B= B/13.0;
  }
    
  center.setValue(B[0],B[1],B[2]);
  min.setValue(B0[0],B0[1],B0[2]);
  max.setValue(B1[0],B1[1],B1[2]);
  box.setBounds(min,max);
}

#define HTMPL 100
float points[HTMPL][3];

// updateChildren
void SoHelicalTrack::updateChildren() {


  // This routine simply draws the points along the helix.  It's invoked
  // whenever primitives are generated or a render action occurs.  

  assert(children->getLength()==1);
  SoSeparator   *sep           = (SoSeparator *)  ( *children)[0];
  SoCoordinate3 *controlPts    = (SoCoordinate3 *)( sep->getChild(0));
  SoLineSet  *curve            = (SoLineSet *)    ( sep->getChild(1));
  for (int i=0;i<HTMPL;i++) {
    SbVec3f TrackPosition= getTrackPosition(i*(s1.getValue()-s0.getValue())/HTMPL+s0.getValue());
    for (int j=0;j<3;j++) {
      points[i][j]=TrackPosition[j];
    }
  }
  // Add a "trajectory" with these points into the
  // separator:

  for (int h=0;h<HTMPL;h++) controlPts->point.set1Value(h,points[h][0],points[h][1],points[h][2]);

}

// generateChildren
void SoHelicalTrack::generateChildren() {

  // This routines creates one SoSeparator, one SoCoordinate3, and
  // one SoLineSet, and puts it in the child list.  This is done only
  // once, whereas redrawing the position of the coordinates occurs each
  // time an update is necessary, in the updateChildren routine. 

  assert(children->getLength() ==0);
  SoSeparator   *sep = new SoSeparator(); 
  SoCoordinate3 *controlPts = new SoCoordinate3;
  SoLineSet     *curve      = new SoLineSet;
  sep->addChild(controlPts);
  sep->addChild(curve);
  children->append(sep);
}

// generateAlternateRep
void SoHelicalTrack::generateAlternateRep() {

  // This routine sets the alternate representation to the child
  // list of this mode.  

  if (children->getLength() == 0) generateChildren();
  updateChildren();
  alternateRep = (SoSeparator *)  ( *children)[0];
}

// clearAlternateRep
void SoHelicalTrack::clearAlternateRep() {
  alternateRep = NULL;
}

// getTrackPosition
SbVec3f SoHelicalTrack::getTrackPosition(float s){

  // get the trajectory as a function of distance
  // this routine is for the track class.  Another class might
  // substitute it's own helper routines here.

  float W        = inverseRadius.getValue();
  float CotTheta = cotTheta.getValue();
  float Phi0     = phi0.getValue();
  float D0       = d0.getValue();
  float Z0       = z0.getValue();
  float Theta    = M_PI/2.0 - atan(CotTheta);
  float CosTheta = cos(Theta);
  float SinTheta = sin(Theta);

  if (s==0.0 || W==0.0) {
    float phi1     = Phi0+s*W*SinTheta;
    float xtan     = SinTheta*cos(phi1);
    float ytan     = SinTheta*sin(phi1);
    float ztan     = CosTheta;
    return SbVec3f(-D0*cos(Phi0+M_PI/2.0),-D0*sin(Phi0+M_PI/2),Z0)+SbVec3f(xtan,ytan,ztan)*s;
  }
  return SbVec3f(
		 (cos(Phi0)*sin(s*W*SinTheta)-sin(Phi0)*(-W*D0+1.0-cos(s*W*SinTheta)))/W,       
		 (sin(Phi0)*sin(s*W*SinTheta)+cos(Phi0)*(-W*D0+1.0-cos(s*W*SinTheta)))/W,       
		 s*CosTheta + Z0);
}
