/*-----------------------------Hepvis----------------------------------------*/
/*                                                                           */
/* Node:             SoArrow                                                 */
/* Description:      Represents an Arrow, or Vector or something             */
/* Author:           Joe Boudreau September 1997                             */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#include <assert.h>
#include <math.h>
#include <Inventor/SbBox.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/nodes/SoRotationXYZ.h>
#include <Inventor/nodes/SoCylinder.h>
#include <Inventor/nodes/SoCone.h>
#include <Inventor/nodes/SoScale.h>
#include "SoArrow.h"

// This statement is required
SO_NODE_SOURCE(SoArrow)

// Constructor
SoArrow::SoArrow() {
  // This statement is required
  SO_NODE_CONSTRUCTOR(SoArrow);

  // Data fields are initialized like this:
  SO_NODE_ADD_FIELD(tail,                (0,0,0));
  SO_NODE_ADD_FIELD(tip,                 (0,0,1));
  SO_NODE_ADD_FIELD(alternateRep,        (NULL));
  children = new SoChildList(this);
}

// Destructor
SoArrow::~SoArrow() {
  delete children;
}


// initClass
void SoArrow::initClass(){
  // This statement is required.
  SO_NODE_INIT_CLASS(SoArrow,SoShape,"Shape");
}


// GLRrender
void SoArrow::GLRender(SoGLRenderAction *action) {
  if (!shouldGLRender(action)) return;
  if (children->getLength() == 0) generateChildren();
  updateChildren();
  children->traverse(action);  
}

// generatePrimitives
void SoArrow::generatePrimitives(SoAction *action) {
  if (children->getLength()==0) generateChildren();
  updateChildren();
  children->traverse(action);
}

// getChildren
SoChildList *SoArrow::getChildren() const {
  return children;
}


// computeBBox
void SoArrow::computeBBox(SoAction *, SbBox3f &box, SbVec3f &center ){
  SbVec3f mtip=tip.getValue(),
          mtail=tail.getValue();
  float size=(mtip-mtail).length();
  
  float xmin,ymin,zmin,xmax,ymax,zmax;

  if (mtip[0]>mtail[0]) {
    xmax =mtip [0]+0.2*size;
    xmin =mtail[0]-0.2*size;
  }
  else {
    xmax =mtail[0]+0.2*size;
    xmin =mtip [0]-0.2*size;
  }

  if (mtip[1]>mtail[1]) {
    ymax =mtip [1]+0.2*size;
    ymin =mtail[1]-0.2*size;
  }
  else {
    ymax =mtail[1]+0.2*size;
    ymin =mtip [1]-0.2*size;
  }

  if (mtip[2]>mtail[2]) {
    zmax =mtip [2]+0.2*size;
    zmin =mtail[2]-0.2*size;
  }
  else {
    zmax =mtail[2]+0.2*size;
    zmin =mtip [2]-0.2*size;
  }

  SbVec3f MMin(xmin,ymin,zmin), MMax(xmax,ymax,zmax);
  center.setValue((xmax+xmin)/2,(ymax+ymin)/2,(zmax+zmin)/2);
  box.setBounds(MMin,MMax);
}




// updateChildren
#define ARROWSIZE 0.02
void SoArrow::updateChildren() {
  assert(children->getLength()==1);
  if (cachedTip!=tip.getValue() || cachedTail!=tail.getValue()) {

    // Redraw the Arrow....
    
    SoSeparator       *axis                = (SoSeparator *)  ( *children)[0];
    SoTranslation     *finalTranslation    = (SoTranslation *)( axis->getChild(0));
    SoScale           *scale               = (SoScale *      )( axis->getChild(1));
    SoRotation        *myRotation          = (SoRotation *   )( axis->getChild(2));
    SoRotationXYZ     *rot                 = (SoRotationXYZ *)( axis->getChild(3));
    SoTranslation     *zTranslation        = (SoTranslation *)( axis->getChild(4));
    SoCylinder        *axisCyl             = (SoCylinder *   )( axis->getChild(5));
    SoTranslation     *cTranslation        = (SoTranslation *)( axis->getChild(6));
    SoCone            *axisCone            = (SoCone *       )( axis->getChild(7));
    
    
    double length = (tip.getValue()-tail.getValue()).length();
    finalTranslation->translation.setValue(tail.getValue());
    scale->scaleFactor.setValue(length/2,length/2,length/2);
    SbVec3f ax = SbVec3f(0,0,1).cross(tip.getValue()-tail.getValue());
    double  an = asin(ax.length()/length);
    myRotation->rotation.setValue(ax,an);
    zTranslation->translation.setValue(0,1,0);
    rot->axis=SoRotationXYZ::X;
    rot->angle=M_PI/2.0;
    axisCyl->radius.setValue(2*ARROWSIZE);
    cTranslation->translation.setValue(0,1,0);
    axisCone->bottomRadius.setValue(ARROWSIZE*4.0);
    axisCone->height.setValue(ARROWSIZE*8.0);
  }
  cachedTip=tip.getValue();
  cachedTail=tail.getValue();
}

// generateChildren
void SoArrow::generateChildren() {
  assert(children->getLength() ==0);
  SoSeparator   *axis=new SoSeparator();
  SoTranslation *finalTranslation=new SoTranslation;
  SoScale       *scale=new SoScale();
  SoRotation *myRotation = new SoRotation();
  SoRotationXYZ  *rot = new SoRotationXYZ;
  SoTranslation  *zTranslation=new SoTranslation;
  SoCylinder *axisCyl = new SoCylinder;
  SoTranslation  *cTranslation=new SoTranslation;
  SoCone     *axisCone=new SoCone;

  axis->addChild(finalTranslation);
  axis->addChild(scale);
  axis->addChild(myRotation);
  axis->addChild(rot);
  axis->addChild(zTranslation);
  axis->addChild(axisCyl);
  axis->addChild(cTranslation);
  axis->addChild(axisCone);
  children->append(axis);
}

// generateAlternateRep
void SoArrow::generateAlternateRep() {

  // This routine sets the alternate representation to the child
  // list of this mode.  

  if (children->getLength() == 0) generateChildren();
  updateChildren();
  alternateRep = (SoSeparator *)  ( *children)[0];
}

// clearAlternateRep
void SoArrow::clearAlternateRep() {
  alternateRep = NULL;
}

