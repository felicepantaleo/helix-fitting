/*-----------------------------Hepvis----------------------------------------*/
/*                                                                           */
/* Node:             SoHelicalTrack                                          */
/* Description:      Represents a helical trajectory with axis along z       */
/* Author:           Joe Boudreau Nov 7 1996                                 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#ifndef __SoHelicalTrack_hh_
#define __SoHelicalTrack_hh__
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFNode.h>
#include <Inventor/nodes/SoShape.h>

class SbBox;
class SoSFNode;
class SoHelicalTrack:public SoShape {

  // The following is required:
  SO_NODE_HEADER(SoHelicalTrack);
 
public:

  // Always use Inventor Fields. This allows Inventor to detect a change to
  // the data field and take the appropriate action; e.g., redraw the scene.
  SoSFFloat inverseRadius;         // signed! positive if anticlockwise
  SoSFFloat cotTheta;              // theta is scattering angle
  SoSFFloat phi0;                  // particle direction at closest approach
  SoSFFloat d0;                    // signed! positive if Lz>0
  SoSFFloat z0;                    // z postition at closest approach
  SoSFFloat s0;                    // inner arclength
  SoSFFloat s1;                    // outer arclength
  SoSFNode alternateRep;           // the alternate representation, required

  // Constructor, required
  SoHelicalTrack();

  // Class Initializer, required
  static void initClass();

  // Generate AlternateRep, required.  Generating an alternate representation
  // must be done upon users request.  It allows an Inventor program to read
  // back the file without requiring *this* code to be dynamically linked. 
  // If the users expects that *this* code will be dynamically linked, he
  // need not invoke this method.  
  virtual void generateAlternateRep();

  // We better be able to clear it, too!
  virtual void clearAlternateRep();

protected:

  // Compute Bounding Box, required
  virtual void computeBBox(SoAction *action, SbBox3f & box, SbVec3f &center);

  // Generate Primitives, required
  virtual void generatePrimitives(SoAction *action);

  // GLRender, required
  virtual void GLRender(SoGLRenderAction *action); 

  // GetChildList, required whenever the class has hidden children
  virtual SoChildList *getChildren() const;


private: 

  // Destructor.  Required.  Always Private.
  virtual ~SoHelicalTrack();

  // Generate Children. Used to create the hidden children. Required whenever
  // the node has hidden children.  
  void generateChildren();  

  // Used to modify hidden children when a data field is changed. Required 
  // whenever the class has hidden children. 
  void updateChildren();

  // ChildList. Required whenever the class has hidden children.  
  SoChildList *children;

  // THE METHODS BELOW ARE PARTICULAR TO THE TRACK CLASS.  IN GENERAL
  // ONE CAN PUT ANY SORT OF HELPER FUNCTIONS IN PLACE OF THIS:

protected:

  // Retrieves to track position as a function of distance
  virtual SbVec3f getTrackPosition(float distance);

private:

  float myInverseRadius;
  float myCotTheta;     
  float myPhi0;         
  float myD0;           
  float myZ0;           
  float myS0;           
  float myS1;           
};

#endif
