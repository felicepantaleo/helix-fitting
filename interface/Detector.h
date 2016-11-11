/*!
 * \file Detector.h
 * \brief An instantiation detector class represents a specific
 * physical detector geometry.
 *
 * The geometry is read from the filename that is supplied to
 * the constructor: Detector ( string filename )
 * The detectors geometric data can be given to stdout
 * via write (human readable) and write_raw.
 * write_raw also has the capability to write to a file
 * instead of stdout: write_raw( char *file)
 *
 * \author Wolfgang Waltenberger
 * \date Tue 10 Jul 2001 23:41:54 CEST
 */

#ifndef Detector_H
#define Detector_H

#include <string>
#include "DetElement.h"
#include "CLHEP/config/CLHEP.h"
#include "CLHEP/config/iostream.h"
#include "Detector.h"

#ifdef GL
#include <Inventor/nodes/SoSelection.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/Xt/SoXtGLWidget.h>
#include <Inventor/actions/SoLineHighlightRenderAction.h>
#include <Inventor/Xt/viewers/SoXtExaminerViewer.h>
#include <Xm/PushB.h>
#endif

using namespace std;

#ifdef GL
struct ToggleEntry { // an entry for a radiobox
	string text;
	string keybind;
	string accel;
	XtCallbackProc cb;
	bool start;
	int num;
};
#endif
/// The detector class.
class Detector {
private:
	string Name;
	void sort_elements();
	#ifdef GL
	void writeGeneralInfo ();
	SoLineHighlightRenderAction *myAction;
	SoXtExaminerViewer *myViewer;
	Widget myparent;
	
	void setupMenubar ( Widget ); ///< initializes the menubar
	/// add a button in a (sub)menu
	void addMenuEntry ( Widget , char *name, string key, char *mnemo, \
			XtCallbackProc );
	/// adds a radio box in a (sub)menu
	Widget addRadioBox ( Widget , char *name, char *mnemo, ToggleEntry * );
	/// adds a submenu or a subsubmenu or ...
	Widget addPulldownMenu ( Widget parent, char *name, char *mnemo );
	
	/* nobody should have to use these: they are implicitly called by 
	 * the above */
	void addToggleEntry ( Widget , ToggleEntry *entries );
	void addMenuEntry ( Widget , char *name, string key, char *mnemo, \
			XtCallbackProc, XtPointer );
	void addToggleEntry ( Widget parent, char *name, string key, char *mnemo,\
			XtCallbackProc _proc, bool is_set);
	void addToggleEntry ( Widget parent, char *name, string key, char *mnemo,\
			XtCallbackProc _proc, bool is_set, XtPointer);
	void addToggleEntry ( Widget parent, char *name, string key, char *mnemo,\
			XtCallbackProc _proc, bool is_set,  int );
	Widget addRadioBox ( Widget , char *name, char *mnemo );
	
	SoSelection *_root;
	Widget myWindow;
	Widget mainWindow;
	#endif
public:
	#ifdef GL
	void clearText();
	void addText(char *TextItem);
	void addText(string TextItem);
	bool writePickedPath (SoNode *root, const SbViewportRegion &viewport,
		bool verbose , const SbVec2s &cursorPosition);
	void writePath (int indx);
	int oldindx; // the last item picked
//	string BField() {return bfieldinfo; };
	string bfieldinfo();
	SoSelection *getRoot() { return _root; };
	Widget getGLparent() { return myparent; };
	
	void Lights ( eltype ) const;
	void Camera ( eltype ) const;
	
	void GLprologue ( SoSelection * , char *windowname );
	void GLdraw ( ) ;
	void GLepilogue ();
	SoSelection *SoSel() {return _root; };
	inline string Detector::keyName ( string  s );
	#endif /* GL */

	string filename;
	int Num;	///< Number of detector elements a detector is made of.
	double Br,Bphi,Bz; ///< B Field in cylinder coordinates [T]
	DetElement *Element;  ///< The Array of all DetElements.
	Detector (string filename);
	string Filename() const;
	void write();
	void write_raw();
	void write_raw(char *file);
	///< A good value for the arc length of a helix.
	HepDouble arclength ();
};
#endif
