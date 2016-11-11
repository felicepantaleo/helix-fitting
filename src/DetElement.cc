//
// Detector Elements:
// DetElements and Slices
//

#include <CLHEP/config/iostream.h>
#include <CLHEP/config/CLHEP.h>
#include <stdio.h>
#include "DetElement.h"
#include "defines.h"
#include "any2str.h"
#include <assert.h>
#include <string>
#include <vector>

#ifdef GL

#ifdef TEXTURE
#include <Inventor/nodes/SoTexture2.h>
#endif
#include <Inventor/nodes/SoSelection.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoTransform.h>

#include "hepvis/SoG4Tubs.hh"
	
#endif /* GL */

#ifdef GL
extern void addChild ( SoSelection *root, SoNode *node , string (*info)(int),int , bool );
extern void addChild ( SoSelection *root, SoNode *node );
#endif

vector <DetElement> el_list;

string eltype2str ( eltype tp ) {
	char ret[3];
	snprintf(ret,5,"%d",tp);
	return ret;
};


/// Get the long name of element type 'tp'
string GetElementName(eltype tp) {
	switch (tp) {
		case CYL:
			return "Cyl";
		case DIS:
			return "Dis";
		case UNK:
			return "???";
		case B:
			return "BFd";
		case NON:
			return "Non";
		default:
			ERROR ("Element type " + eltype2str(tp) + " unknown");
			return "???";
	};
};

/// Get element type # for 'ElementName'
eltype GetElementType(string ename) {
	int i=UNK;
	string name;
	while (i<MAX_ELEMENT_TYPES-1 ) { // -1 because the last element is NON
		name=GetElementName((eltype) i);
		if (name==ename) {
			return ((eltype) i);
		};
		i++;
	};
	return UNK;
};

string DetElement::SType() const {
	char ret[3];
	snprintf(ret,5,"%d",Type());
	return ret;
};

DetElement::DetElement (const eltype & ityp, const HepDouble & a, 
    const HepDouble & b , const HepDouble & c, const HepDouble & d) {
	typ=ityp; num=0;
	_Br=0, _BPhi=0, _Bz=0;
	nxt=NULL;
	switch (ityp) {
		case CYL:
			#ifdef DEBUG
			assert( a >  0 );
			assert( b >= 0 );
			assert( c >= 0 );
			assert( d >= 0 );
			#endif
			r = a; sigz = c; sigphi=b/a; z=0; sigr=0;
			sigrphi=sigphi * r;
			thick = d;
			break;
		case DIS:
			#ifdef DEBUG
			assert( a >= 0 );
			assert( b >= 0 );
			assert( c >= 0 );
			assert( d >= 0 );
			#endif
			z = a; sigr=c; r=0; sigrphi=b ; sigz=0;
			sigphi = 0;
			thick = d;
			break;
		case B:
			r=a; phi=b; z=c; thick=0;
			break;
		case NON:
			r=a; phi=b; z=c; thick=0;
			break;
		default:
			ERROR("unknown Detector Element " + eltype2str(ityp) );
	};
};

DetElement::DetElement () {
	// I want a simple way of initializing an array of geometric objects
	r = 0; sigz=0; sigphi=0; z=0; num=0; _Br=0, _BPhi=0, _Bz=0; sigrphi=0;
	nxt=NULL;
}

HepDouble DetElement::sigRPhi() const {
	return sigrphi;
};

HepDouble DetElement::sigPhi() const {
	#ifdef DEBUG
	if ( Type() == DIS ) ERROR("sigphi undefined");
	#endif
	return sigphi;
};
void DetElement::setBField( HepDouble tr, HepDouble tPhi, HepDouble tz)
{
	_Br=tr; _BPhi=tPhi; _Bz=tz;
};

string DetElement::TypeName () const {
	return GetElementName(typ);
};
	
#ifdef GL
void DetElement::GLdeinitMaterial ( SoSelection *root ) const {
	#ifdef TEXTURE
	SoTexture2 *texture = new SoTexture2;
	texture->filename.setValue("");
	root->addChild(texture);
	#endif
};
	
void DetElement::GLinitMaterial ( SoSelection *root ) const {
	SoMaterial *myMaterial = new SoMaterial;
	if (typ==CYL) {
	 	myMaterial->diffuseColor.setValue(0.03, 0.4, 0.03);
	 	myMaterial->emissiveColor.setValue(0.1, 0.4, 0.1);
		myMaterial->transparency=.5;
		myMaterial->shininess=1;
	} else {
	 	myMaterial->diffuseColor.setValue(0.03, 0.4, 0.03);
	 	myMaterial->emissiveColor.setValue(0.1, 0.4, 0.1);
		myMaterial->transparency=.50;
		myMaterial->shininess=1;
	};
	addChild(root,myMaterial);

	#ifdef TEXTURE
	SoTexture2 *texture = new SoTexture2;
	addChild(root,texture);
	texture->filename.setValue(ETCDIR"/grey-circuit.rgb");
	texture->model=SoTexture2::DECAL;
	#endif
};
	
string DEgetInfo(int a)
{ 
	return el_list[a].getInfo();
};

string DetElement::getInfo () const {
	char name[99];
	if (typ==CYL) {
		snprintf(name,99,"Detector Element.\n"\
		"      Radius=%.4fm, SigRPhi=%.6f, Sigz=%.8f\n"\
		"      Thickness=%.4f",R(), sigRPhi(), sigz,thick );
	} else {
		snprintf(name,99,"Detector Element.\n"\
		"      z=%.4fm, SigRPhi=%.6f, Sigr=%.8f\n"\
		"      Thickness=%.4f",Z(), sigRPhi(), sigR(),thick );
	};
	return name;
};

void DetElement::GLdraw ( SoSelection *root ) const {
	SoG4Tubs *newCyl = new SoG4Tubs;
	if (typ==CYL) {
		newCyl->pRMin = R()-.25 * thick;
		newCyl->pRMax = R()+.25 * thick;
		newCyl->pDz = 0.08;
	} else { // typ==DIS
		newCyl->pRMin = 0;
		newCyl->pRMax = .05 + Z() / 5;
		newCyl->pDz = .1 * thick;
		SoTransform *t_a = new SoTransform;
		t_a->translation.setValue ( 0 , 0 , Z() );
		addChild ( root, t_a );
	};
	newCyl->smoothDraw=true;
	el_list.push_back(*this);
	addChild(root,newCyl, (DEgetInfo) , el_list.size()-1 , false );
	if ( typ== DIS ) {
		SoTransform *t_b = new SoTransform;
		t_b->translation.setValue ( 0 , 0 , - Z() );
		addChild ( root, t_b );
	};
		
};
#endif

void DetElement::write () const {
	switch (Type()) {
		case CYL:
			cout << TypeName() << " Radius=" << R()
	         << " SigPhi=" << sigPhi() << " Sigz=" << sigZ() 
			 << " Thick=" << Thick() << endl;
			break;
		case DIS:
			cout << TypeName() << " z=" << Z()
	         << " SigPhi=" << sigPhi() << " Sigz=" << sigZ() 
			 << " Thick=" << Thick() << endl;
			break;
		default:
			ERROR("unknown detector element " + SType() );
	};
};

// How can we read in detector elements from a file/stdin? ->
// by this kewl operator overloading technique!

#define MLC 255 // MCL Maximum Length for a Configuration line

/// deletes comments from a string|(char *), returning a string.

inline string del_leading_spaces(string line) {
	while (line.find_first_of(" 	")==0) {
		line.erase(0,1);  // delete leading spaces and tabs
	};
	return line;
};

inline eldesc parse_line(string line) {
	eldesc el; /*={NON,1,0,0};*/
	el.typ=NON;
	el.var[0]=1;
	el.var[1]=0;
	el.var[2]=0;
	el.var[3]=0;
	int i;
	string nxtword,myline;
	line=line.substr(0,line.find("#")); // delete comments
	line=del_leading_spaces(line);
	nxtword=line.substr(0,line.find(" "));
	if (!nxtword.size()) { // empty line - jump to the next line
		el.typ = NON;
		return el;
	};
	el.typ=GetElementType(nxtword);
	for (i=0;i<4;i++) {
		line=del_leading_spaces(line);
		line.erase(0,line.find_first_of(" 	")+1);
		line=del_leading_spaces(line);
		nxtword=line.substr(0,line.find(" "));
		if(!nxtword.size()) {
			clog << "Warning: Number of values is wrong in file " << 
			     " entry #" << i << endl;
		};
		el.var[i]=atof(nxtword.c_str());
	};
	return el;
}

istream& operator>> (istream& s, DetElement& el)
{
	char word[MLC]; // fixed length strings are ok here.
	eldesc eld;
	
	s.getline(word,MLC);
	eld=parse_line(word);
	if (eld.typ != NON) // ... then the line is not empty
		el=DetElement( eld.typ , eld.var[0] , eld.var[1] , 
		               eld.var[2], eld.var[3] );
	return s;
}

