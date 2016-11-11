/*-----------------------------Hepvis----------------------------------------*/
/*                                                                           */
/* Node:             SoG4Tubs                                                */
/* Description:      Represents the G4Tubs Geant Geometry entity             */
/* Author:           Joe Boudreau Nov 11 1996                                */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#include <assert.h>
#include <math.h>
#include <Inventor/SbBox.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoIndexedFaceSet.h>
#include <Inventor/nodes/SoNormal.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoNormalBinding.h>
#include "hepvis/SoG4Tubs.hh"

// This statement is required
SO_NODE_SOURCE(SoG4Tubs)

// Constructor
SoG4Tubs::SoG4Tubs() {
  // This statement is required
  SO_NODE_CONSTRUCTOR(SoG4Tubs);

  // Data fields are initialized like this:
  SO_NODE_ADD_FIELD(pRMin,                (0.0));
  SO_NODE_ADD_FIELD(pRMax,                (1.0));
  SO_NODE_ADD_FIELD(pDz,                 (10.0));
  SO_NODE_ADD_FIELD(pSPhi,                (0.0));
  SO_NODE_ADD_FIELD(pDPhi,             (2*M_PI));
  SO_NODE_ADD_FIELD(smoothDraw,          (TRUE));
  SO_NODE_ADD_FIELD(alternateRep,        (NULL));
  children = new SoChildList(this);
}

// Destructor
SoG4Tubs::~SoG4Tubs() {
}


// initClass
void SoG4Tubs::initClass(){
  // This statement is required.
  SO_NODE_INIT_CLASS(SoG4Tubs,SoShape,"Shape");
}

// GLRrender
void SoG4Tubs::GLRender(SoGLRenderAction *action) {
  if (!shouldGLRender(action)) return;
  if (children->getLength() == 0) generateChildren();
  if ( (prmin      != pRMin.getValue()) ||
       (prmax      != pRMax.getValue()) ||
       (pdz        != pDz.getValue())   ||
       (psphi      != pSPhi.getValue()) ||
       (pdphi      != pDPhi.getValue()) ) {
    // SmoothDraw is damned!
    //       (smoothdraw  && !smoothDraw.getValue())||
    // (!smoothdraw && smoothDraw.getValue()))    {
    updateChildren();
  }
  children->traverse(action);  
  prmin     =pRMin.getValue();
  prmax     =pRMax.getValue();
  pdz       =pDz.getValue();
  psphi     =pSPhi.getValue();
  pdphi     =pDPhi.getValue();
  // smoothDraw=smoothDraw.getValue();
 

}

// generatePrimitives
void SoG4Tubs::generatePrimitives(SoAction *action) {
  if (children->getLength()==0) generateChildren();
  updateChildren();
  children->traverse(action);
}

// getChildren
SoChildList *SoG4Tubs::getChildren() const {
  return children;
}


// computeBBox
void SoG4Tubs::computeBBox(SoAction *, SbBox3f &box, SbVec3f &center ){
  SbVec3f min(-pRMax.getValue(),-pRMax.getValue(),-pDz.getValue()), 
          max( pRMax.getValue(), pRMax.getValue(), pDz.getValue());
  center.setValue(0,0,0);
  box.setBounds(min,max);
}


// updateChildren
void SoG4Tubs::updateChildren() {

  // Redraw the G4Tubs....

  assert(children->getLength()==1);
  SoSeparator       *sep                = (SoSeparator *)  ( *children)[0];
  SoCoordinate3     *theCoordinates     = (SoCoordinate3 *)      ( sep->getChild(0));
  SoNormal          *theNormals         = (SoNormal *)           ( sep->getChild(1)); 
  SoNormalBinding   *theNormalBinding   = (SoNormalBinding *)    ( sep->getChild(2));
  SoIndexedFaceSet  *theFaceSet         = (SoIndexedFaceSet *)   ( sep->getChild(3));


  const int NPHI=24, NPOINTS=2*(2*NPHI+2), NFACES=4*NPHI+2, NINDICES = NFACES*5;
  float points[NPOINTS][3],normals[NFACES][3];
#ifdef INVENTOR2_0
  static long     indices[NINDICES];
#else
  static int32_t  indices[NINDICES];
#endif

  static int init=0;
  double phi, pp, DeltaPhi;

  // Indices need to be generated once! This is here to keep it close to the point
  // generation, since otherwise it will be confusing.

  int i;
  if (!init) {
    init = 1;
    // Outer face
    for (i = 0; i< NPHI; i++) {
      // 0 1 3 2;
      indices[5*i+0] =  2*i+0;
      indices[5*i+1] =  2*i+1;
      indices[5*i+2] =  2*i+3;
      indices[5*i+3] =  2*i+2;
      indices[5*i+4] = SO_END_FACE_INDEX;
    }
    // the inner face
    for (i=0;i<NPHI;i++) {
      indices[5*1*NPHI + 5*i+0] = 2*NPHI+2 + 2*i+0;  
      indices[5*1*NPHI + 5*i+1] = 2*NPHI+2 + 2*i+1;
      indices[5*1*NPHI + 5*i+2] = 2*NPHI+2 + 2*i+3;
      indices[5*1*NPHI + 5*i+3] = 2*NPHI+2 + 2*i+2;
      indices[5*1*NPHI + 5*i+4] = SO_END_FACE_INDEX;
    }
    // the top side
    for (i=0;i<NPHI;i++) {
      indices[5*2*NPHI + 5*i+0] = 2*i+0;
      indices[5*2*NPHI + 5*i+1] = 2*i+2;
      indices[5*2*NPHI + 5*i+2] = NPOINTS - (2*i+4);
      indices[5*2*NPHI + 5*i+3] = NPOINTS - (2*i+2);
      indices[5*2*NPHI + 5*i+4] = SO_END_FACE_INDEX;
    }
    // the bottom side
    for (i=0;i<NPHI;i++) {
      indices[5*3*NPHI + 5*i+0] = 2*i+1;
      indices[5*3*NPHI + 5*i+1] = NPOINTS - (2*i+1);
      indices[5*3*NPHI + 5*i+2] = NPOINTS - (2*i+3);
      indices[5*3*NPHI + 5*i+3] = 2*i+3;
      indices[5*3*NPHI + 5*i+4] = SO_END_FACE_INDEX;
    }
    // the odd side
    indices[5*4*NPHI +0] = 2*NPHI;
    indices[5*4*NPHI +1] = 2*NPHI+1;
    indices[5*4*NPHI +2] = 2*NPHI+3;
    indices[5*4*NPHI +3] = 2*NPHI+2;
    indices[5*4*NPHI +4] = SO_END_FACE_INDEX;
    // aother odd side
    indices[5*4*NPHI +5 +0] = 0;
    indices[5*4*NPHI +5 +1] = NPOINTS-2;
    indices[5*4*NPHI +5 +2] = NPOINTS-1;
    indices[5*4*NPHI +5 +3] = 1;
    indices[5*4*NPHI +5 +4] = SO_END_FACE_INDEX;
  }
  // Points need to be generated each time:

  // The outer surface
  DeltaPhi = pDPhi.getValue()/NPHI, phi = pSPhi.getValue();
  for (i = 0; i<=NPHI; i++) {
    points[2*i+0][0] = pRMax.getValue()*cos(phi); points[2*i+0][1]= pRMax.getValue()*sin(phi); points[2*i+0][2] = +pDz.getValue();
    points[2*i+1][0] = pRMax.getValue()*cos(phi); points[2*i+1][1]= pRMax.getValue()*sin(phi); points[2*i+1][2] = -pDz.getValue();
    pp = phi+DeltaPhi/2.0;
    if (i!=NPHI) normals[i][0] = cos(pp); normals[i][1] = sin(pp); normals[i][2] = 0;
    phi+=DeltaPhi;
  }
  // The inner surface
  phi = pSPhi.getValue() + pDPhi.getValue();
  for (i = 0; i<=NPHI; i++) {
    points[2*NPHI+2+2*i+0][0] = pRMin.getValue()*cos(phi); 
    points[2*NPHI+2+2*i+0][1] = pRMin.getValue()*sin(phi); 
    points[2*NPHI+2+2*i+0][2] = +pDz.getValue();
    points[2*NPHI+2+2*i+1][0] = pRMin.getValue()*cos(phi); 
    points[2*NPHI+2+2*i+1][1] = pRMin.getValue()*sin(phi); 
    points[2*NPHI+2+2*i+1][2] = -pDz.getValue();
    pp = phi-DeltaPhi/2.0;
    if (i!=NPHI) normals[NPHI+i][0] = -cos(pp); normals[NPHI+i][1] = -sin(pp); normals[NPHI+i][2] = 0;
    phi-=DeltaPhi;
  }
  // The top side
  for (i=0;i<NPHI;i++) {
    normals[2*NPHI+i][0]=normals[2*NPHI+i][1]=0; normals[2*NPHI+i][2]=  1.0;
  } 
  // The bottom side
  for (i=0;i<NPHI;i++) {
    normals[3*NPHI+i][0]=normals[3*NPHI+i][1]=0; normals[3*NPHI+i][2]= -1.0;
  } 
  // The odd side
  phi = pSPhi.getValue(); 
  normals[4*NPHI+0][0]=sin(phi); normals[4*NPHI+0][1]= -cos(phi); normals[4*NPHI+0][2]=0;
  
  // Another odd side
  phi = pSPhi.getValue()+pDPhi.getValue(); 
  normals[4*NPHI+1][0]=-sin(phi); normals[4*NPHI+1][1]= +cos(phi); normals[4*NPHI+1][2]=0;

  for (int np=0;np<NPOINTS; np++) theCoordinates->point.set1Value(np,points[np][0],points[np][1],points[np][2]);
  for (int ni=0;ni<NINDICES;ni++) theFaceSet->coordIndex.set1Value(ni,indices[ni]);
  if (smoothDraw.getValue()) {
    //    This Line is replaced by the next one because of an apparent Bug in Inventor (mem. leak).
    //theNormals->vector.deleteValues(0);
    for (int nf=0;nf<NFACES;nf++) theNormals->vector.set1Value(nf,normals[nf][0],normals[nf][1],normals[nf][2]);
    theNormalBinding->value=SoNormalBinding::PER_FACE;
  }
  else {
    for (int nf=0;nf<NFACES;nf++) theNormals->vector.set1Value(nf,normals[nf][0],normals[nf][1],normals[nf][2]);
    theNormalBinding->value=SoNormalBinding::PER_FACE;
  }
}

// generateChildren
void SoG4Tubs::generateChildren() {

  // This routines creates one SoSeparator, one SoCoordinate3, and
  // one SoLineSet, and puts it in the child list.  This is done only
  // once, whereas redrawing the position of the coordinates occurs each
  // time an update is necessary, in the updateChildren routine. 

  assert(children->getLength() ==0);
  SoSeparator      *sep              = new SoSeparator(); 
  SoCoordinate3    *theCoordinates   = new SoCoordinate3();
  SoNormal         *theNormals       = new SoNormal(); 
  SoNormalBinding  *theNormalBinding = new SoNormalBinding();
  SoIndexedFaceSet *theFaceSet       = new SoIndexedFaceSet();
  // 
  // This line costs some in render quality! but gives speed.
  // 
  sep->addChild(theCoordinates);
  sep->addChild(theNormals);
  sep->addChild(theNormalBinding);
  sep->addChild(theFaceSet);
  children->append(sep);
}

// generateAlternateRep
void SoG4Tubs::generateAlternateRep() {

  // This routine sets the alternate representation to the child
  // list of this mode.  

  if (children->getLength() == 0) generateChildren();
  updateChildren();
  alternateRep = (SoSeparator *)  ( *children)[0];
}

// clearAlternateRep
void SoG4Tubs::clearAlternateRep() {
  alternateRep = NULL;
}

