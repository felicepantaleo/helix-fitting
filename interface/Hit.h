/*!
 * \file
 * \brief The base hit class.
 *
 * The base 'Hit' class implements the concept of an 'ideal'
 *
 * \author Wolfgang Waltenberger
 * \date Thu 05 Jul 2001 15:15:46 CEST
 */

#ifndef Hit_H
#define Hit_H

#include "DetElement.h"
#include "defines.h"
#include <string>
#include "StateVector.h"
#include <CLHEP/config/iostream.h>
#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Matrix/SymMatrix.h"


#ifdef GL
#include <Inventor/nodes/SoSelection.h>
#include "hepvis/SoEllipsoid.hh"
#endif

/// The Hit base class.
class Hit {
private:
	union {
		HepDouble _RPhi;
		HepDouble _x;
	};
	union {
		HepDouble _z;
		HepDouble _y;
	};
	DetElement _el;
	HepSymMatrix _sig;
	StateVector _sv; ///< store the statevector that produced our hit.
public:
	Hit ( const DetElement & , const HepDouble & RPhi, const HepDouble & z );
	Hit ( const DetElement & el );
	Hit ( StateVector & ); ///< constructs a Hit from a StateVector
	Hit ();

	HepDouble R()   const;
	HepDouble RPhi()   const;
	HepDouble z()   const;
	HepDouble x() const;
	HepDouble y() const;
	HepDouble Phi() const;
	HepSymMatrix Cov() const { return _sig; };
	HepDouble SigRPhi() const;
	HepDouble SigR() const;
	HepDouble VarPhi() const;
	HepDouble VarRPhi() const;
	HepDouble VarR() const;
	HepDouble Varx() const;
	HepDouble Vary() const;
	HepDouble Covxy() const;
	HepDouble Sigz() const;
	HepDouble Varz() const;
	HepDouble Sigx() const;
	HepDouble Sigy() const;
	
	void setSig( const HepSymMatrix & m) { _sig=m; };
	/// Smears out the hits.
	void Smear();
	/// Smears out the hits. Deterministically if 'det' is true.
	/// Deterministically: R_meas ( n ) = R_true ( n ) + (-1)^n sigR
	void Smear( const bool det );
	/// Create a hit from the parameters of a StateVector
	void getHit( StateVector & ); 
	
	void write_raw () const;
	void write_raw ( const char * ) const;
	void write () const;
	#ifdef GL
	string Coords() const;
	void GLinitMaterial ( SoSelection * ) const;
	void GLdraw ( SoSelection * ) ;
	string getInfo () ;
	#endif
};
#endif
