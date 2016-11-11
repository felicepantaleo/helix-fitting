/*!
 * \file DetElement.h
 * \brief The class of elements a detector is made of 
 * (i.e. cylinder, discs, ...)
 *
 * \author Wolfgang Waltenberger
 * \date Tue 14 May 2002 10:08:17 CESTDST
 */

#ifndef DetElement_H
#define DetElement_H

#include <string>
#include "CLHEP/config/CLHEP.h"
#include "defines.h"

#ifdef GL
#include <Inventor/nodes/SoSelection.h>
#endif

using namespace std;

/// Currently we support five types of elements:
/// Unknown, Cylinder, Disc, B (magnetic field), and None 
/// B - the magnetic field - is only in a very restricted sense
/// a detector element ...
enum eltype { UNK, CYL, DIS, B, NON };
// typedef el_type eltype; // does this help?
#define MAX_ELEMENT_TYPES 5
// is there a (better) way to count
// the number of elements in an 'enum'?

struct eldesc {
	eltype typ;
	HepDouble var[4];
};

// Two little 'helpers'
extern string GetElementName ( eltype );
extern eltype GetElementType ( string ElementName );

/// The Class of detector elements.
class DetElement {
public:
	DetElement();
	DetElement ( const eltype & ityp , const HepDouble & , const HepDouble & , const HepDouble & , const HepDouble & );
		// CYLINDER: r, sigrphi, sigz, thickness
		// SLI: ?
	HepDouble R () const { return r; };
	HepDouble Z ()      const { return z; };
	HepDouble Phi ()    const { return phi; };
	HepDouble sigPhi () const;
	HepDouble sigRPhi() const;
	HepDouble sigZ ()   const { return sigz; };
	HepDouble sigR ()   const { return sigr; };
	HepDouble sigx ()   const { ERROR("no sigx yet"); return sigr; };
	HepDouble sigy ()   const { ERROR("no sigy yet"); return sigr; };
	HepDouble Thick ()  const { return thick; };
	HepDouble Br()      const { return _Br; };
	HepDouble Bx()      const { return _Br * cos (_BPhi) ; };
	HepDouble By()      const { return _Br * sin (_BPhi); };
	HepDouble BPhi()    const { return _BPhi; };
	HepDouble Bz()      const { return _Bz; };
	void setBField( HepDouble r, HepDouble Phi, HepDouble z);

	eltype Type () const {return typ; };
	/// Type(), given as a string.
	string SType () const;
	void write () const;
	#ifdef GL
	void  GLinitMaterial ( SoSelection *) const;
	void  GLdeinitMaterial ( SoSelection *) const;
	/// draws the Detector element.
	void  GLdraw ( SoSelection *) const;
	string getInfo () const;
	#endif

	string TypeName () const;
	/// if we are in a detector, then this is the 'num'th element.
	int num;
	DetElement *nxt; ///< a pointer to the next detector element
private:
	HepDouble _Br, _BPhi, _Bz;
	HepDouble r, sigphi, phi, z, sigrphi, sigz, sigr;
	HepDouble thick;
	eltype typ; // this eltype cause problems sometimes.
};

#endif /* DetElement_H */
