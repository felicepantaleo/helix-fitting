/*-----------------------------Hepvis----------------------------------------*/
/*                                                                           */
/* Node:             SoG4Tubs                                                */
/* Description:      Represents the G4Tubs Geant Geometry entity             */
/* Author:           Joe Boudreau Nov 11 1996                                */
/*                                                                           */
/*---------------------------------------------------------------------------*/
#ifndef __SoG4Tubs_hh_
#define __SoG4Tubs_hh__
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoSFNode.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/nodes/SoShape.h>

class SoSFNode;
class SoG4Tubs:public SoShape {

  // The following is required:
  SO_NODE_HEADER(SoG4Tubs);
 
public:


  //  // The documentation from Geant says:
  //  //
  //  // $Id: SoG4Tubs.hh,v 1.1 2001/08/19 21:15:02 wwalten Exp $
  //  // class G4Tubs
  //  //
  //  // A tube or tube segment with curved sides parallel to
  //  // the z-axis. The tube has a specified half-length along
  //  // the z axis, about which it is centred, and a given
  //  // minimum and maximum radius. A minimum radius of 0
  //  // signifies a filled tube /cylinder. The tube segment is
  //  // specified by starting and delta
  //  // angles for phi, with 0 being the +x axis, PI/2
  //  // the +y axis. A delta angle of 2PI signifies a
  //  // complete, unsegmented tube/cylinder

  // Always use Inventor Fields. This allows Inventor to detect a change to
  // the data field and take the appropriate action; e.g., redraw the scene.

  SoSFFloat pRMin;   
  SoSFFloat pRMax;
  SoSFFloat pDz;
  SoSFFloat pSPhi;
  SoSFFloat pDPhi;             
  SoSFBool  smoothDraw;            // slightly better render, worse performance
  SoSFNode  alternateRep;          // the alternate representation, required

  // Constructor, required
  SoG4Tubs();

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

  // computeBBox, required
  virtual void computeBBox(SoAction *action, SbBox3f &box, SbVec3f &center );

  // Generate Primitives, required
  virtual void generatePrimitives(SoAction *action);

  // GLRender, required
  virtual void GLRender(SoGLRenderAction *action); 

  // GetChildList, required whenever the class has hidden children
  virtual SoChildList *getChildren() const;


private: 

  // Destructor.  Required.  Always Private.
  virtual ~SoG4Tubs();

  // Generate Children. Used to create the hidden children. Required whenever
  // the node has hidden children.  
  void generateChildren();  

  // Used to modify hidden children when a data field is changed. Required 
  // whenever the class has hidden children. 
  void updateChildren();

  // ChildList. Required whenever the class has hidden children.  
  SoChildList *children;

  float prmin;   
  float prmax;
  float pdz;
  float psphi;
  float pdphi;             
  SbBool smoothdraw;            // slightly better render, worse performance
};

#endif
