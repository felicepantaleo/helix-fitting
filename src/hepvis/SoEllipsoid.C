/*-----------------------------Hepvis----------------------------------------*/
/*                                                                           */
/* Node:             SoEllipsoid                                             */
/* Description:      Represents an Ellipsoid                                 */
/* Author:           Joe Boudreau Nov 11 1996                                */
/* hacked by wwalten (no level of detail here)                               */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#include <assert.h>
#include <math.h>
#include <Inventor/SbBox.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoScale.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/nodes/SoTranslation.h>
#include "hepvis/SoEllipsoid.hh"

// This statement is required
SO_NODE_SOURCE(SoEllipsoid)

// Constructor
SoEllipsoid::SoEllipsoid() {
  // This statement is required
  SO_NODE_CONSTRUCTOR(SoEllipsoid);

  // Data fields are initialized like this:
  SO_NODE_ADD_FIELD(eigenvalues,              (1.0,1.0,1.0));
  SO_NODE_ADD_FIELD(rotation,                 (SbVec3f(0,0,0),0)  );
  SO_NODE_ADD_FIELD(center,                   (SbVec3f(0,0,0))  );
  SO_NODE_ADD_FIELD(alternateRep,       (NULL));
  children = new SoChildList(this);
}

// Destructor
SoEllipsoid::~SoEllipsoid() {
}


// initClass
void SoEllipsoid::initClass(){
  // This statement is required.
  SO_NODE_INIT_CLASS(SoEllipsoid,SoShape,"Shape");
}


// GLRrender
void SoEllipsoid::GLRender(SoGLRenderAction *action) {
  if (!shouldGLRender(action)) return;
  if (children->getLength() == 0) generateChildren();
  updateChildren();
  children->traverse(action);  
}

// generatePrimitives
void SoEllipsoid::generatePrimitives(SoAction *action) {
  if (children->getLength()==0) generateChildren();
  updateChildren();
  children->traverse(action);
}

// getChildren
SoChildList *SoEllipsoid::getChildren() const {
  return children;
}


// computeBBox
void SoEllipsoid::computeBBox(SoAction *, SbBox3f &box, SbVec3f &TheCenter ){
  float Extent = eigenvalues.getValue().length();
  SbVec3f min(-Extent,-Extent,-Extent), 
          max( Extent, Extent, Extent);
  min = min+center.getValue();
  max = max+center.getValue();
  TheCenter=center.getValue();
  box.setBounds(min,max);
}




// updateChildren
void SoEllipsoid::updateChildren() {

}

// generateChildren
void SoEllipsoid::generateChildren() {

  // A box consists of a set of scale factors and a 
  // cube.

  assert(children->getLength() ==0);
  SoSeparator      *sep              = new SoSeparator();
  SoTranslation    *trans            = new SoTranslation();
  SoRotation       *rot              = new SoRotation();
  SoScale          *scale            = new SoScale();
  SoSphere         *sphere           = new SoSphere();

  rot->rotation.connectFrom(&rotation);
  trans->translation.connectFrom(&center);
  scale->scaleFactor.connectFrom(&eigenvalues);
  sep->addChild(trans);
  sep->addChild(rot);
  sep->addChild(scale);
  sep->addChild(sphere);
  children->append(sep);
}

// generateAlternateRep
void SoEllipsoid::generateAlternateRep() {

  // This routine sets the alternate representation to the child
  // list of this mode.  

  if (children->getLength() == 0) generateChildren();
  updateChildren();
  alternateRep = (SoSeparator *)  ( *children)[0];
}

// clearAlternateRep
void SoEllipsoid::clearAlternateRep() {
  alternateRep = NULL;
}

