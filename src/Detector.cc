/*
 * Detector.cc
 */

#include "Detector.h"
#include "defines.h"
#include "Hit.h"
#include "FullVector.h"
#include <string>
#include <vector>
#include <iostream.h>
#include <fstream>
#include <CLHEP/config/iostream.h>
#include <CLHEP/config/fstream.h>
#include <stdio.h>


#ifdef GL
#include "hepvis/SoArrow.h"
#include <X11/StringDefs.h>
#include <Xm/MainW.h>
#include <Xm/PanedW.h>
#include <X11/Shell.h>
#include <Xm/Text.h>
#include <Xm/Form.h>
#include <Xm/List.h>
#include <Xm/Label.h>
#include <Xm/TextF.h>
#include <Xm/ArrowB.h>
#include <Xm/RowColumn.h>
#include <Xm/CascadeB.h>
#include <Xm/PushB.h>
#include <Xm/ToggleB.h>
#include <Xm/ToggleBG.h>
#include <Xm/FileSB.h>
#include <Xm/MessageB.h>

#include <Inventor/SbViewportRegion.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/Xt/viewers/SoXtExaminerViewer.h>
#include <Inventor/Xt/SoXtRenderArea.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoSelection.h>

#include <Inventor/events/SoMouseButtonEvent.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/actions/SoLineHighlightRenderAction.h>
#include <Inventor/Xt/SoXt.h>
#endif

#ifndef XmSET
#define XmSET true
#endif

/*
#ifndef XmUNSET
#define XmUNSET true
#endif
*/

using namespace std; // solaris needs that.

inline string del_leading_spaces(string line) {
	while (line.find_first_of(" 	")==0) {
		line.erase(0,1);  // delete leading spaces and tabs
	};
	return line;
};
extern istream& operator>> (istream& s, DetElement& el);
extern string GetElementName(eltype);

// clicking on a helix:
// get the track or the fullvector?
// bool show_track=false;
#ifdef GL
int   hlx_shows=3; // what data does a helix show
int   hit_shows=1;   // what data does clicking on a hit show
int matrix_prec=1; // precision of the matrices.
Detector *my_det;

#endif

void file_not_found (string filename) {
	ERROR("Error: file '" + filename + "' not found ");
	exit (-1);
};

#ifdef GL
// this is in traxx.cpp, because we need to know the commandline arguments.
extern void invokeRestart(Widget, XtPointer , XmPushButtonCallbackStruct *);
static vector <string (*)(int)> GLinfos;
static vector <int> GLindex;
static vector <bool> GLsel;
static vector <int> GLselindices;
Widget _theTextArea;
Widget _dialogTextArea;
Widget dialog;


void Detector::addText(string TextItem) {
  if  (XtIsManaged(_theTextArea))
	XmTextInsert(_theTextArea,XmTextGetLastPosition(_theTextArea),
			(char *) TextItem.c_str());
};

void Detector::addText(char *TextItem) {
  if  (XtIsManaged(_theTextArea))
	XmTextInsert(_theTextArea,XmTextGetLastPosition(_theTextArea),TextItem);
};

void Detector::clearText() {
  if  (XtIsManaged(_theTextArea))
    XmTextReplace(_theTextArea,0,XmTextGetLastPosition(_theTextArea),"");
};

void Detector::writePath (int indx)
{
	if (indx<GLindex.size() && indx > 0) {
		clearText();
		char index[8];
		sprintf(index,"[%d]",indx);
		addText(index);
		addText(" Selected object is a: ");
		addText((char *) GLinfos[indx](GLindex[indx]).c_str());
		XmTextSetInsertionPosition(_theTextArea,0);
		oldindx=indx;
	
		// we select 'manually' in order to sync with the description
		SoNode *mynode= my_det->getRoot()->getChild(indx);
		if (mynode) {
			my_det->getRoot()->deselectAll(); // update selection
			my_det->getRoot()->select( mynode );
		};
	};
};

void Detector::writeGeneralInfo ()
{
	char line[255];
	clearText();
	#ifdef VERSION
	addText("This is the TinyCMS Project, version: "VERSION"\n\n");
	#endif

	
	snprintf(line,255,"Detector: %s (%s)\n",Name.c_str(),filename.c_str());
	addText(line);
	snprintf(line,255,"z scales: %f",(float) Z_SCALE);
	addText(line);
	XmTextSetInsertionPosition(_theTextArea,0);
	getRoot()->deselectAll(); // update selection
	my_det->oldindx=0;
};

bool Detector::writePickedPath (SoNode *root, const SbViewportRegion &viewport,
		bool verbose , const SbVec2s &cursorPosition) {
	
	SoRayPickAction myPickAction(viewport);
	
	// Set an 4-pixel wide region around the pixel
	myPickAction.setPoint(cursorPosition);
	myPickAction.setRadius( 4.0 );
	
	// Start a pick traversal
	myPickAction.apply(root);
	const SoPickedPoint *myPickedPoint = 
		myPickAction.getPickedPoint();
	if (myPickedPoint == NULL) {
		writeGeneralInfo();
		return FALSE;
	};
	SoPath *mypath = myPickedPoint->getPath();
	SbName myname = mypath->getName();
	int lngth=mypath->getLength();
	if (lngth==3) {
		int indx = mypath->getIndex(2);
			SoNode *mynode=mypath->getNode(2);
		writePath(indx);
	} else {
		writeGeneralInfo();
	};
};

void toggle_prec(Widget, XtPointer Client, XmPushButtonCallbackStruct *)
{
	matrix_prec=(int) Client;
	my_det->clearText();
	my_det->writePath(my_det->oldindx); // update description
};

void toggle_hit(Widget, XtPointer Client, XmPushButtonCallbackStruct *)
{
	hit_shows=(int) Client;
	my_det->clearText();
	my_det->writePath(my_det->oldindx); // update description
};

void cancelObject(Widget root, XtPointer , XmPushButtonCallbackStruct *)
{
	XtUnmanageChild(dialog);
};
	
// ok. we hit the 'ok' button. now let's get the info.
void writeObject(Widget root, XtPointer , XmPushButtonCallbackStruct *)
{
	my_det->clearText();
	string text;
	text=(string) XmTextGetString(_dialogTextArea);
	int _num=0;
	_num=atoi(text.c_str());
	if (_num && _num < GLindex.size()) {
		my_det->writePath(_num); // update description
		SoNode *mynode= my_det->getRoot()->getChild(_num);
		if (mynode) {
			my_det->getRoot()->deselectAll(); // update selection
			my_det->getRoot()->select( mynode );
		};
	};
	XtUnmanageChild(dialog);
};

void av_list_cb (Widget root, XtPointer Client, XmListCallbackStruct *blub)
{
	my_det->clearText();
	int _num=GLselindices[blub->item_position-1];
	if (_num && _num < GLindex.size()) {
		my_det->writePath(_num); // update description
		SoNode *mynode= my_det->getRoot()->getChild(_num);
		if (mynode) {
			my_det->getRoot()->deselectAll(); // update selection
			my_det->getRoot()->select( mynode );
		};
	};
};

// write Object. This is here, so we can hit return in the widget;
// we do not want to have to hit the 'ok' button.
void writeObject2(Widget r, XtPointer c, XmPushButtonCallbackStruct *s)
{
	string text;
	text=(string) XmTextGetString(_dialogTextArea);
	if (text[text.size()-1]=='\n') { writeObject (r,c,s);};
};

void delete_about()
{
	XtUnmanageChild(dialog);
};
	
void about_w(Widget root, XtPointer , XmPushButtonCallbackStruct *)
{
	// Text Area!
	Cardinal ac = 0;
	Arg args[20];
	
	XtSetArg (args[ac], XmNallowShellResize, True); ac++;
	XtSetArg (args[ac], XmNtitle, "About"); ac++;
	XtSetArg (args[ac], XmNiconName, "About"); ac++;
	XtSetArg (args[ac], XmNdeleteResponse, XmUNMAP); ac++;
	XtSetArg (args[ac], XmNwidth,380);ac++;
	XtSetArg (args[ac], XmNheight,156);ac++;
	dialog = XtCreatePopupShell( "About", topLevelShellWidgetClass,\
		   	my_det->getGLparent() , args,ac);
	XtManageChild(dialog);
	
	ac=0;
	XtSetArg (args[ac], XmNeditMode,XmMULTI_LINE_EDIT);ac++;
	XtSetArg (args[ac], XmNmarginHeight, 7); ac++;
	XtSetArg (args[ac], XmNmarginWidth, 7); ac++;
	Widget prompt_w = XmCreateText (dialog, "AboutText",args,ac);
	XmTextSetString(prompt_w,"   This is the TinyCMS project, Version "\
			VERSION".\n\n"\
			"Bug reports, feature requests, patches,\n"\
			"and free pizza can be deposited at:"\
			"\nWolfgang Waltenberger <wwalten@hephy.oeaw.ac.at>");
	XtManageChild(prompt_w);
	ac=0;
	XtSetArg (args[ac], XmNy,118);ac++;
	XtSetArg (args[ac], XmNx,3);ac++;
	XtSetArg (args[ac], XmNwidth,374);ac++;
	XtSetArg (args[ac], XmNmarginHeight, 7); ac++;
	XtSetArg (args[ac], XmNmarginWidth, 7); ac++;
	XtSetArg(args[ac],XmNlabelString,XmStringCreateSimple("OK")); ac++;
	Widget _okwidget = XmCreatePushButton( dialog, "OK" , args, ac);
	XtAddCallback( _okwidget , XmNactivateCallback, 
			(XtCallbackProc) delete_about,(XtPointer) 0 );
	XtManageChild(_okwidget);
};


void showObject_m(Widget root, XtPointer , XmPushButtonCallbackStruct *)
{
	// Text Area!
	Cardinal ac = 0;
	Arg args[20];
	
	XtSetArg (args[ac], XmNallowShellResize, True); ac++;
	XtSetArg (args[ac], XmNtitle, "Show Object"); ac++;
	XtSetArg (args[ac], XmNiconName, "Show Object"); ac++;
	XtSetArg (args[ac], XmNdeleteResponse, XmUNMAP); ac++;
	XtSetArg (args[ac], XmNmarginHeight, 10); ac++;
	XtSetArg (args[ac], XmNmarginWidth, 10); ac++;
	dialog = XtCreatePopupShell( "Dialog", topLevelShellWidgetClass, \
			my_det->getGLparent() , args,ac);
	XtManageChild(dialog);
	
	ac = 0;
	XtSetArg (args[ac], XmNresizePolicy, XmRESIZE_ANY); ac++;
	XtSetArg (args[ac], XmNmarginHeight, 10); ac++;
	XtSetArg (args[ac], XmNmarginWidth, 10); ac++;
	Widget lstform_w = XmCreateForm(dialog, "ListF", args, ac);
	XtManageChild (lstform_w);
	
	ac = 0;
	XtSetArg (args[ac], XmNrows,1);ac++;
	XtSetArg (args[ac], XmNwidth,40);ac++;
	XtSetArg (args[ac], XmNheight,30);ac++;
	XtSetArg (args[ac], XmNx,20);ac++;
	XtSetArg (args[ac], XmNy,40);ac++;
	XtSetArg (args[ac], XmNmarginHeight, 3); ac++;
	XtSetArg (args[ac], XmNmarginWidth, 3); ac++;
	// to get the return keys, we define a multiline window
	XtSetArg (args[ac], XmNeditMode,XmMULTI_LINE_EDIT);ac++;
	_dialogTextArea = XmCreateText(lstform_w, "textArea", args,ac);
	XtManageChild(_dialogTextArea);
	XtAddCallback( _dialogTextArea , XmNvalueChangedCallback, \
			(XtCallbackProc) writeObject2, (XtPointer) 0 );
	#if XmVersion >= 1001
	/* init kb focus here if possible */
	XmProcessTraversal (_dialogTextArea, XmTRAVERSE_CURRENT);
	XmProcessTraversal (_dialogTextArea, XmTRAVERSE_CURRENT);
	#endif
	
	
	ac=0;
	XtSetArg (args[ac], XmNalignment, XmALIGNMENT_BEGINNING); ac++;
	XtSetArg (args[ac], XmNx,5);ac++;
	XtSetArg (args[ac], XmNy,2);ac++;
	XtSetArg (args[ac], XmNmarginHeight, 10); ac++;
	XtSetArg (args[ac], XmNmarginWidth, 10); ac++;
	Widget prompt_w = XmCreateLabel (lstform_w, "Show Object #", args, ac);
	XtManageChild(prompt_w);

	ac=0;
	XtSetArg (args[ac], XmNx,80);ac++;
	XtSetArg (args[ac], XmNy,35);ac++;
	XtSetArg (args[ac], XmNmarginHeight, 7); ac++;
	XtSetArg (args[ac], XmNmarginWidth, 7); ac++;
	XtSetArg(args[ac],XmNlabelString,XmStringCreateSimple("OK")); ac++;
	Widget _okwidget = XmCreatePushButton( lstform_w, "OK" , args, ac);
	XtAddCallback( _okwidget , XmNactivateCallback, \
			(XtCallbackProc) writeObject, (XtPointer) 0 );
	XtManageChild(_okwidget);
	
	ac=0;
	XtSetArg (args[ac], XmNx,140);ac++;
	XtSetArg (args[ac], XmNy,35);ac++;
	XtSetArg (args[ac], XmNmarginHeight, 7); ac++;
	XtSetArg (args[ac], XmNmarginWidth, 7); ac++;
	XtSetArg(args[ac],XmNlabelString,XmStringCreateSimple("Cancel")); ac++;
	Widget _cancelwidget = XmCreatePushButton( lstform_w, "Cancel" , args, ac);
	XtAddCallback( _cancelwidget , XmNactivateCallback, \
			(XtCallbackProc) cancelObject, (XtPointer) 0 );
	
	XtManageChild(_cancelwidget);
	ac=0;
	XtSetArg (args[ac], XmNalignment, XmALIGNMENT_BEGINNING); ac++;
	XtSetArg (args[ac], XmNx,10);ac++;
	XtSetArg (args[ac], XmNy,80);ac++;
	XtSetArg (args[ac], XmNmarginHeight, 10); ac++;
	XtSetArg (args[ac], XmNmarginWidth, 10); ac++;
	Widget track_w = XmCreateLabel (lstform_w, "Show Track", args, ac);
	XtManageChild(track_w);
	ac=0;
	XtSetArg (args[ac], XmNx,10);ac++;
	XtSetArg (args[ac], XmNy,110);ac++;
	XtSetArg (args[ac], XmNmarginHeight, 5); ac++;
	XtSetArg (args[ac], XmNmarginWidth, 0); ac++;
	XtSetArg (args[ac], XmNvisibleItemCount, 10 ); ac++;
	XtSetArg (args[ac], XmNselectionPolicy, XmSINGLE_SELECT); ac++;
	Widget _listwidget = XmCreateScrolledList ( lstform_w, "List", args, ac);
	XtAddCallback (_listwidget, XmNsingleSelectionCallback, av_list_cb, 
			(XtPointer) NULL );
	char bla[9];
	for (int i=0; i< GLindex.size(); i++) {
		if (GLsel[i]) {
			snprintf(bla,9,"[%d] ",i);
			string sshort=bla;
			string::size_type pos;
			pos=GLinfos[i](GLindex[i]).find('.');
			#ifdef DEBUG
			if ( pos == string::npos ) 
				WARNING("Identifier strings need to have a '.' (dot)");
			#endif
			sshort+=GLinfos[i](GLindex[i]).substr(0,pos);
			XmListAddItem ( _listwidget, 
					XmStringCreateSimple( (char *) sshort.c_str() ) , 0);
			GLselindices.push_back ( i ); // remember the object number.
		};
	};
	
	XtManageChild(_listwidget);
};


void toggle_track(Widget root, XtPointer Client, XmPushButtonCallbackStruct *)
{
	hlx_shows=(int) Client;
	my_det->clearText();
	my_det->writePath(my_det->oldindx); // update description
};

void invokeQuit(Widget, XtPointer ClientData, XmPushButtonCallbackStruct *)
{
	exit (0);
};

string no_info_available (int)
{
	return "No information available";
}

void addChild ( SoSelection *root, SoNode *node )
{
	root->addChild (node);
	GLinfos.push_back( (no_info_available) );
	GLindex.push_back( 0 );
	GLsel.push_back( false );
};

// The 'C interface'

void addChild ( SoSelection *root, SoNode *node , string (*info)(int), int a )
{
	root->addChild (node);
	GLinfos.push_back( (info) );
	GLindex.push_back( a );
	GLsel.push_back( false );
};

/* use this method, if we want the *node to be in the 'selection list' */
void addChild ( SoSelection *r, SoNode *n, string (*info)(int), int a, bool b )
{
	r->addChild (n);
	GLinfos.push_back( (info) );
	GLindex.push_back( a );
	GLsel.push_back( b );
};
#endif

// looks awkward, but really, we consider the file to be already sorted,
// anyways, so we do permutations.
void Detector::sort_elements()
{
	bool is_sorted=1;
	DetElement swap;
	int i;
	float a,b;
	eltype ta,tb;
	for (i=0;i<Num-1;i++) {
		ta = Element[i].Type(); tb = Element[i+1].Type();
		if (ta == CYL) {
			a = Element[i].R(); 
		} else {
			a = Element[i].Z();
		};
		if (tb == CYL) {
			b = Element[i+1].R(); 
		} else {
			b = Element[i+1].Z();
		};
		if (( ta == tb && a > b ) || ( ta > tb ))
		{
			is_sorted=0;
			swap=Element[i]; Element[i]=Element[i+1]; Element[i+1]=swap;
		};
	};
	if (!is_sorted) {sort_elements();};
};

Detector::Detector (string fname)
{
	#ifdef GL
	oldindx=-1;
	my_det=this;
	#endif
	filename=fname;
	int i=0;
	size_t elsize;
	Br=0; Bphi=0; Bz=0;
	static DetElement *el;
	elsize=sizeof (el[0]);
	el=(DetElement *) malloc( elsize ); 
	ifstream fin(fname.c_str());
	if (!fin) file_not_found(fname);
	char word[MLC];
	
	/// get the description from the file.
	fin.getline(word,MLC); // get the detector name
	string sword=word;
	sword=del_leading_spaces(sword);
	sword=sword.substr(0,sword.find("#")); // delete comments
	if ( sword == "" ) {
		Name="Unnamed";
	} else {
		Name=sword;
	};
	while (fin >> el[i]) {
		el[i].nxt= NULL; // just to be sure
		if ( el[i].Type()==B ) { // The B field has seperate variables
			if (Br!=0 || Bphi!=0 || Bz!=0) {
				i--;
			} else {
				Br=el[i].R();
				Bphi=el[i].Phi();
				Bz=el[i].Z();
				i--;
			};
		};
		if ( el[i].Type()!=NON && el[i].Type()!=UNK ) {
			el[i].num=i; el[i-1].nxt=&el[i];
			i++;
			// lets not waste one byte. not one.
			el=(DetElement *) realloc( el, ( i + 1 )* elsize);
		};
	};
	Num=i;
	Element=el;
	sort_elements();
	
	// And now we push our overall constant B field in all DetElements.
	// FIXME sometime we should change this, and assume a changing B field.
	
	for (int j=0;j<Num; j++ ) {
		el[j].setBField( Br, Bphi, Bz );
	};
	MESSAGE(8,"[detector] '"+Name+"' initialized.");
};

void Detector::write_raw() {
	// we write the B field separately. Makes the output easier to parse.
	cout << "Element " << GetElementName(B) << " " << Br << " " << Bphi << " " << Bz << endl;
	for (int i=0;i<Num;i++) {
		if (Element[i].Type()==CYL ) {
			cout << "Element " << Element[i].TypeName() << " "
			     << Element[i].R() << " " << Element[i].sigRPhi() << " " 
				 << Element[i].sigZ() << endl;
		} else if (Element[i].Type()==DIS ) {
			cout << "Element " << Element[i].TypeName() << " "
			     << Element[i].Z() << " " << Element[i].sigRPhi() << " " 
				 << Element[i].sigR() << endl;
		};
	};
};


void Detector::write_raw(char *file) {
	if (!strcmp(file,"")) {
		Detector::write_raw(); // hm. How can we merge the two?
	} else {
		static ofstream fout(file);
		MESSAGE(4,"[system] writing raw data to file " + (string) file);
		// we write the B field separately. Makes the output easier to parse.
		fout << "Element " << GetElementName(B) << " " << Br << " " << Bphi << " " 
		     << Bz << endl;
		for (int i=0;i<Num;i++) {
			if (Element[i].Type()==CYL ) {
				fout << "Element " << Element[i].TypeName() << " "
				     << Element[i].R() << " " << Element[i].sigRPhi() << " " 
					 << Element[i].sigZ() << endl;
			} else if (Element[i].Type()==DIS ) {
				fout << "Element " << Element[i].TypeName() << " "
				     << Element[i].Z() << " " << Element[i].sigRPhi() << " " 
					 << Element[i].sigR() << endl;
			};
		};
	};
};
string Detector::Filename() const {
	return filename;
};
	
void Detector::write() {
	cout << "Detector: " << endl;
	for (int i=0;i<Num;i++) {
		cout.width(Num % 10);
		cout << Element[i].num << " " << Element[i].TypeName();
		switch (Element[i].Type()) {
		case B:
			cout << " Br=" << Element[i].R()
			     << " Bphi=" << Element[i].Phi()
			     << " Bz=" << Element[i].Z() << endl;
			break;
		case CYL:
			cout << " Radius=" << Element[i].R()
			     << " SigPhi=" << Element[i].sigPhi() 
			     << " Sigz=" << Element[i].sigZ() << endl;
			break;
		case DIS:
			cout << " z=" << Element[i].Z()
				 << " SigPhi=" << Element[i].sigPhi() 
			     << " Sigz=" << Element[i].sigZ() << endl;
			break;
		default:
			ERROR( "unknown Detector Element " + Element[i].SType() );
		};
	};
}

HepDouble Detector::arclength () { 
	if ( Element[0].Type() == CYL ) {
		return Element[Num-1].R()*1.5;
	} else {
		return Element[Num-1].Z()*1.5;
	};
};
#ifdef GL

string getBfieldInfo(int)
{
	return my_det->bfieldinfo();
};

string Detector::bfieldinfo()
{
	char name[55];
	snprintf(name,55,"Mean Magnetic Field.\n     r=%.2f Phi=%.2f z=%.2f [T]",\
			Br,Bphi,Bz);
	return (string) name;
};
void Detector::GLdraw ( )
{
	//  draw the magnetic field
	SoArrow *GLBfield = new SoArrow;
	SoMaterial *myMaterial = new SoMaterial;
	myMaterial->diffuseColor.setValue(0.9, 0.7, 0.23);
	myMaterial->emissiveColor.setValue(1. , 0.5 , 0.16 );
	myMaterial->transparency=.5;
	myMaterial->shininess=1;
	SoSFVec3f mytip;
	SoSFVec3f mytail;
	HepDouble bz=0, bx, by ,br=0, bphi=0;
	for (int i=0; i<=Num;i++) {
		bz   += Element[i].Bz();
		br   += Element[i].Br();
		bphi += Element[i].BPhi();
	};
	bz /=Num; br /=Num; bphi/=Num;
	bx=br * cos (bphi) / B_SCALE;
	by=br * cos (bphi) / B_SCALE;
	bz /= B_SCALE;
	mytip.setValue(bx , by , bz );
	mytail.setValue(-bx,-by, -bz );
	GLBfield->tip=mytip;
	GLBfield->tail=mytail;
	
	addChild(_root,myMaterial);
	addChild(_root,GLBfield, (getBfieldInfo) , 0 );
	
	//                draw the detector elements
	Element[0].GLinitMaterial ( _root );
	for ( int j=0; j < Num; j++ ) {
		Element[j].GLdraw( _root ); 
	};
	Element[0].GLdeinitMaterial ( _root );
};

// This routine is called for every mouse button event.
void myMousePressCB(void *userData, SoEventCallback *eventCB)
{
	SoSelection *root = (SoSelection *) userData;
	const SoEvent *event = eventCB->getEvent();
	
	// Check for mouse button being pressed
	if (SO_MOUSE_PRESS_EVENT(event, BUTTON1)) {
		const SbViewportRegion &myRegion = 
			eventCB->getAction()->getViewportRegion();
		my_det->writePickedPath(root, myRegion, true,
			event->getPosition(myRegion));
		eventCB->setHandled();
	};
	if (SO_MOUSE_PRESS_EVENT(event, BUTTON2)) {
		const SbViewportRegion &myRegion = 
			eventCB->getAction()->getViewportRegion();
		my_det->writePickedPath(root, myRegion, false,
			event->getPosition(myRegion));
		eventCB->setHandled();
	};
};

void Detector::GLepilogue ( )
{
	myViewer = new SoXtExaminerViewer( myWindow );
	myViewer->setDecoration( true );
	myViewer->setHeadlight( false );
	
	SoEventCallback *myEventCB = new SoEventCallback;
	myEventCB->setName("EventCB");
	addChild(_root,myEventCB);
	
	myViewer->setSceneGraph(_root);
	SbColor myHighCol(1,.95,.85 );
	myAction = new SoLineHighlightRenderAction();
	myAction->setColor(myHighCol);
	myAction->setLineWidth(2);
	myViewer->setGLRenderAction(myAction);
	myViewer->redrawOnSelectionChange(_root);
	
	myViewer->setTitle("TinyCMS");
	myViewer->setViewing(FALSE);
	myViewer->viewAll();
	myViewer->show();
	
	myEventCB->addEventCallback(
		SoMouseButtonEvent::getClassTypeId(),
		myMousePressCB,
		myViewer->getSceneManager()->getSceneGraph());
	
	// Text Area!
	Cardinal ac = 0;
	Arg args[20];
	XtSetArg(args[ac], XmNwidth, 600); ac++;
	XtSetArg(args[ac], XmNrows , 4); ac++;
	XtSetArg(args[ac], XmNscrollHorizontal , false); ac++;
	XtSetArg(args[ac], XmNscrollVertical , false); ac++;
	XtSetArg(args[ac], XmNeditable , false ); ac++;
	
	_theTextArea = XmCreateScrolledText(myWindow, "theTextArea", args,ac);
	XtManageChild(_theTextArea);
	
	SoXt::show( mainWindow );
	SoXt::mainLoop();
};

void Detector::addMenuEntry ( Widget _parent, char *name, string key, \
		char *mnemo, XtCallbackProc _proc )
{
	addMenuEntry ( _parent, name, key, mnemo, _proc , (XtPointer) this );
};

inline string Detector::keyName ( string  s )
{
	string key;
	key=s;
	string::size_type pos;
	pos=key.find(" <Key>");
	if ( pos != string::npos )
		key.replace(pos,6,"-");
	pos=key.find("<Key>");
	if ( pos != string::npos )
		key.replace(pos,5,"");
	return key;
};

void Detector::addMenuEntry ( Widget _parent, char *name, string key, \
		char *mnemo, XtCallbackProc _proc , XtPointer point )
{
	Arg args[20];
	Cardinal ac = 0;
	Widget _widget;
	
	string keyname=keyName( key );
	
	XmString accstr= XmStringCreateSimple((char *)keyname.c_str());
	XtSetArg(args[ac], XmNlabelString, XmStringCreateSimple(name));
	XtSetArg(args[ac], XmNacceleratorText, accstr);  ac++;
	XtSetArg(args[ac], XmNaccelerator, (XmString) key.c_str());  ac++;
	XtSetArg(args[ac], XmNmnemonic, XStringToKeysym( mnemo)); ac++;
	_widget = XmCreatePushButton(_parent, name ,args, ac);
	XtAddCallback(_widget , XmNactivateCallback, _proc , (XtPointer) point);
	XtManageChild(_widget);
};

Widget Detector::addPulldownMenu ( Widget parent, char *name, char *mnemo )
{
	Arg args[20];
	char buttonname[sizeof(name)+7];
	Cardinal ac = 0;
	Widget _widget;
	Widget _buttonwidget;
	sprintf(buttonname,"%sButton",name);
	
	XtSetArg(args[ac],XmNlabelString,XmStringCreateSimple(name)); ac++;
	XtSetArg(args[ac], XmNmnemonic, XStringToKeysym(mnemo)); ac++;
	_buttonwidget = XmCreateCascadeButton( parent, buttonname , args, ac);
	
	XtManageChild(_buttonwidget);
	
	ac = 0;
	_widget = XmCreatePulldownMenu(XtParent(_buttonwidget), name, args, ac);
	XtSetArg(args[ac], XmNsubMenuId, _widget); ac++;
	XtSetValues(_buttonwidget, args, ac);
	return _widget;
};

Widget Detector::addRadioBox ( Widget parent, char *name, char *mnemo, 
		ToggleEntry *te )
{
	Widget my_w= addRadioBox ( parent, name, mnemo );
	addToggleEntry ( my_w, te );
	return my_w;
}

Widget Detector::addRadioBox ( Widget parent, char *name, char *mnemo )
{
	Arg args[20];
	char buttonname[sizeof(name)+7];
	Cardinal ac = 0;
	Widget _widget;
	Widget _buttonwidget;
	sprintf(buttonname,"%sButton",name);
	
	XtSetArg(args[ac], XmNlabelString, "Hilfe"); ac++;
	XtSetArg(args[ac],XmNlabelString,XmStringCreateSimple(name)); ac++;
	XtSetArg(args[ac], XmNmnemonic, XStringToKeysym(mnemo)); ac++;
	XtSetArg (args[ac], XmNindicatorType, XmONE_OF_MANY); ac++;
	
	_buttonwidget = XmCreateCascadeButton( parent, buttonname , args, ac);
	
	XtManageChild(_buttonwidget);
	
	ac = 0;
	XtSetArg (args[ac], XmNradioBehavior, True); ac++;
	_widget = XmCreatePulldownMenu(XtParent(_buttonwidget), name, args, ac);
	XtSetArg(args[ac], XmNsubMenuId, _widget); ac++;
	XtSetValues(_buttonwidget, args, ac);
	return _widget;
};

void Detector::addToggleEntry ( Widget _parent, char *name, string key, \
		char *mnemo, XtCallbackProc _proc , bool is_set )
{
	addToggleEntry( _parent,name,key,mnemo, _proc , is_set, (XtPointer) this );
};

void Detector::addToggleEntry ( Widget _parent, char *name, string key, \
		char *mnemo, XtCallbackProc _proc , bool is_set , int i)
{
	addToggleEntry (_parent,name,key,mnemo, _proc , is_set, (XtPointer) i );
};

void Detector::addToggleEntry ( Widget w, ToggleEntry *entries )
{
	for ( int i=0; ; i++ ) {
		if (entries[i].text=="") break;
		addToggleEntry (w, (char *) entries[i].text.c_str(), 
		(char *) entries[i].keybind.c_str(),
			(char *) entries[i].accel.c_str(), (XtCallbackProc) entries[i].cb, 
			entries[i].start, entries[i].num );
	};
};
void Detector::addToggleEntry ( Widget _parent, char *name, string key, \
		char *mnemo, XtCallbackProc _proc , bool is_set, XtPointer me )
{
	Arg args[20];
	Cardinal ac = 0;
	Widget _widget;
	if (key!="") {
		string keyname=keyName(key);
		XmString accstr= XmStringCreateSimple((char *)keyname.c_str());
		XtSetArg(args[ac], XmNaccelerator, (XmString) key.c_str());  ac++;
		XtSetArg(args[ac], XmNacceleratorText, accstr);  ac++;
	};
	XtSetArg(args[ac], XmNlabelString, XmStringCreateSimple(name));
	XtSetArg(args[ac], XmNmnemonic, XStringToKeysym( mnemo)); ac++;
	// XmNset = XmSET/true
	if (is_set) { XtSetArg(args[ac], XmNset, XmSET ); ac++;};
// 	if (is_set) { XtSetArg(args[ac], XmNset, true ); ac++;};
	XtSetArg (args[ac], XmNvisibleWhenOff, True); ac++;
	_widget = XmCreateToggleButton(_parent, name ,args, ac);
	XtAddCallback(_widget , XmNvalueChangedCallback, (_proc), (XtPointer) me );
	XtManageChild(_widget);
};

void Detector::setupMenubar ( Widget _MainWindow )
{
	Widget menubar, options_m, file_w;
	Cardinal ac = 0;
	Arg args[20];
	
	menubar = XmCreateMenuBar(_MainWindow,
	    "menuBar", args, ac);
	XtManageChild(menubar);
	
	//            ------------- File ---------------
	file_w=addPulldownMenu ( menubar, "File", "F" );
	addMenuEntry ( file_w, "Quit", "<Key>q", "Q", \
			(XtCallbackProc) invokeQuit );
	addMenuEntry ( file_w, "Restart", "<Key>r", "R", \
			(XtCallbackProc) invokeRestart );
	
	//            ------------- Edit ---------------
	Widget edit_w=addPulldownMenu ( menubar, "Edit", "E" );
	addMenuEntry ( edit_w, "Show ...", "<Key>o", "S", \
			(XtCallbackProc) showObject_m );
	
	//            ------------- Options ------------
	options_m=addPulldownMenu ( menubar, "Options", "O" );
	
	//            ------------- Options: Helix -----
	ToggleEntry hlx_entries[]={
	  {"... Track (CMS)",         "<Key>c","T",toggle_track, hlx_shows==0,0},
	  {"... Track (Delphi)",      "<Key>d","D",toggle_track, hlx_shows==3,3},
	  {"... CovRPhi",             "<Key>p","R",toggle_track, hlx_shows==1,1},
	  {"... CovZ",                "<Key>z","Z",toggle_track, hlx_shows==2,2},
	  {"... CovRPhi^-1",                "","1",toggle_track, hlx_shows==4,4},
	  {"","","",0,0,0} };
	Widget hlx_w=addRadioBox ( options_m, "Helix shows ...", "", hlx_entries );
	//            ------------- Options: Hit -------
	ToggleEntry hentries[]={
	  {"... physical Hit",        "<Key>h","H",toggle_hit,hit_shows==0,0},
	  {"... StateVector (Delphi)","<Key>f","S",toggle_hit,hit_shows==1,1},
	  {"... StateVector (CMS)",   "<Key>v","C",toggle_hit,hit_shows==2,2},
	  {"... Deriv (Delphi)",            "","D",toggle_hit,hit_shows==3,3},
	  {"... Der (Delphi)",              "","E",toggle_hit,hit_shows==4,4},
	  {"","","",0,0,0} }; 
	Widget h_shows=addRadioBox ( options_m, "Hit shows ...", "", hentries );
	//           -------------- Options: Variable precision ------
	ToggleEntry pentries[]={
	  {" 3 decimals","<Key>3","3",toggle_prec  , matrix_prec==0, 0 },
	  {" 4 decimals","<Key>4","4",toggle_prec  , matrix_prec==1, 1 },
	  {" 5 decimals","<Key>5","5",toggle_prec  , matrix_prec==2, 2 },
	  {" 6 decimals","<Key>6","6",toggle_prec  , matrix_prec==3, 3 },
	  {" 7 decimals","<Key>7","7",toggle_prec  , matrix_prec==4, 4 },
	  {" 8 decimals","<Key>8","8",toggle_prec  , matrix_prec==5, 5 },
	  {" 9 decimals","<Key>9","9",toggle_prec  , matrix_prec==6, 6 },
	  {"10 decimals","<Key>0","0",toggle_prec  , matrix_prec==7, 7 },
	  {"11 decimals","<Key>1","1",toggle_prec  , matrix_prec==8, 8 },
	  {"12 decimals","<Key>2","2",toggle_prec  , matrix_prec==9, 9 },
	  {"","","",0,0,0}};
	Widget precis_w=addRadioBox ( 
			options_m, "Variable precision: ", "", pentries );
	
	//            ------------- Help    ------------
	ac=0;
	Widget help_bw = XmCreateCascadeButton( menubar, "Help" , args, ac);
	XtManageChild(help_bw);
	ac = 0;
	XtSetArg (args[ac], XmNradioBehavior, True); ac++;
	Widget help_w = XmCreatePulldownMenu(XtParent(help_bw), "Help", args, ac);
	XtManageChild(help_w);
	XtSetArg(args[ac], XmNsubMenuId, help_w); ac++;
	XtSetValues(help_bw, args, ac);
	ac=0;
	Widget about_bw = XmCreatePushButton( help_w, "About" , args, ac);
	XtAddCallback( about_bw , XmNactivateCallback, \
			(XtCallbackProc) about_w, (XtPointer) 0 );
	XtManageChild(about_bw);
	ac=0;
	XtSetArg(args[ac], XmNmenuHelpWidget, help_bw ); ac++;
	XtSetValues(menubar, args, ac );
};

void Detector::GLprologue ( SoSelection *root, char *win )
{
	_root=root;
	mainWindow = SoXt::init( win );
	if (mainWindow == NULL) exit(1);
	Cardinal ac = 0;
	Arg args[20];
	Widget _MainWindow;
	
	ac = 0;
	XtSetArg(args[ac], XmNx, 69); ac++;
	XtSetArg(args[ac], XmNy, 454); ac++;
	XtSetArg(args[ac], XmNwidth, 612); ac++;
	XtSetArg(args[ac], XmNheight, 510); ac++;
	
	_MainWindow = XmCreateMainWindow(mainWindow,"TinyCMSMain",args, ac);
	XtManageChild(_MainWindow);
	
	myWindow = XmCreatePanedWindow(_MainWindow,"TinyCMSpanedWindow",args, ac);
	XtManageChild(myWindow);
	myparent=myWindow;
	setupMenubar(_MainWindow);
	
	SoMouseButtonEvent myMouseEvent;
	root->ref();
	Camera( Element[0].Type() );
	Lights( Element[0].Type() );
};

void Detector::Lights ( eltype tp ) const
{
	SoSFVec3f sovec;
	SoDirectionalLight *Light = new SoDirectionalLight;
	if ( tp == DIS ) {
		SoDirectionalLight *Light2 = new SoDirectionalLight;
		Light->intensity = 1.;
		Light2->intensity = 1.;
		Light->direction.setValue(0 ,-1,- .02 );
		Light2->direction.setValue(.01,+1, -.01 );
		addChild(_root, Light2);
	} else {
		Light->intensity = .8;
		Light->direction.setValue(.2,-.35,-.9);
	};
	addChild(_root, Light);
};

void Detector::Camera ( eltype tp ) const
{
	SoPerspectiveCamera *myCamera   = new SoPerspectiveCamera;
	SoSFVec3f sovec;
	SoSFRotation rovec;
	addChild(_root,myCamera);
//	myCamera->farDistance=1000000000;
//	myCamera->nearDistance=.1;
	if (tp == DIS) {
		myCamera->orientation.setValue(0,-1,0,1);
		myCamera->position.setValue(0,1,0);
	} else {
		myCamera->orientation.setValue(0,0,1,1);
		myCamera->position.setValue(0,0,1);
	};
}

#endif /* GL */
