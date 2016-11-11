 /*!
 * \file
 * \brief This class realizes the concept of hits and state vectors.
 *
 * This class represents 'planes', given by \f$\vec{n} . \vec{x} + c = 0\f$
 *
 * \author Wolfgang Waltenberger
 * \date Thu 05 Jul 2001 15:17:55 CEST
 */

#ifndef Plane_H
#define Plane_H

#include <string>
#include "DetElement.h"
#include "defines.h"
#include "Detector.h"
#include <vector>
#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

/// The Plane base class.
/// We define a plane by the parameters \f$c, \vec{n}\f$ , where
/// \f$ \vec{n} . \vec{x} = c \f$
class Plane {
private:
	bool CMScalculated;
	double mya0, mypsi, mykappa;
	void calculateCMS();
	HepDouble myr, myphi; // for calculateCMS.
public:
	Plane ();
	HepDouble c; // FIXME these oughta be private
	HepVector n;
	HepDouble z;
	HepDouble R;
	HepDouble Phi;
	HepDouble theta;

	void write () const;

	void setRPhiPoint ( HepDouble r, HepDouble phi ) { myr=r; myphi=phi; };
	/// Center + Radius.
	HepDouble Radius () const;
	/// Vector from Center of Circle to 'interception point'
	HepMatrix RadVec ( vector < HepDouble > rs, HepDouble la0, HepDouble lpsi, 
			HepDouble lkappa ) ;
	HepMatrix RadVec ( vector < HepDouble > rs, HepDouble la0, HepDouble lpsi, 
			HepDouble lkappa, HepDouble rref, HepDouble phiref ) ;
	HepMatrix RadVec ( vector < HepDouble > rs );
	HepVector RadVec ( HepDouble rs );
	HepVector RadVec ( HepDouble rs, HepDouble la0, HepDouble lpsi, \
		HepDouble lkappa );
	HepVector RadVec ( HepDouble rs, HepDouble la0, HepDouble lpsi, \
		HepDouble lkappa, HepDouble rref, HepDouble phiref );
	HepDouble X0 () const;
	HepDouble Y0 () const;

	/// atlas coordinates.
	HepDouble a0 ();
	HepDouble psi ();
	HepDouble kappa ();
	
	/// arc length of RPhi projection of Helix, given a specific detector
	/// or a vector of HepDoubles
	/// supplying rref and phiref results in using the improved
	/// algorithm.
	HepVector s( Detector *det ) ;
	HepVector s( Detector *det, HepDouble, HepDouble, HepDouble ) ;
	HepVector s( Detector *det, HepDouble, HepDouble, HepDouble, 
			HepDouble rref, HepDouble phiref ) ;
	HepVector s( vector < HepDouble > ) ;
	HepVector s( vector < HepDouble >, HepDouble rref, HepDouble phiref ) ;
	HepVector s( vector < HepDouble > , HepDouble, HepDouble, HepDouble );
	HepVector s( vector < HepDouble > , HepDouble, HepDouble, HepDouble,
		   HepDouble rref, HepDouble phiref	);
	
	/// Those cartesian coordinates.
	/// dunno exactly what these mean - FIXME ask fru
	/// FIXME again: cartXY ( double, double, ...) should not be here.
	/// supplying rref and phiref results in using the improved
	/// algorithm.
	HepVector cartXY_old ( HepDouble radius ) ;
	HepVector cartXY_old ( HepDouble , HepDouble, HepDouble, HepDouble ) ;
	HepVector cartXY ( HepVector radii, HepDouble rref, HepDouble phiref );
	HepVector cartXY ( HepDouble radius, HepDouble rref, HepDouble phiref );
};

/// RealPlane. Adds a covariance matrix to the Plane class.
class RealPlane : public Plane {
public:
	HepSymMatrix cov; // FIXME i want to put this into private
	HepSymMatrix Cov()  { return cov; };
	RealPlane ()  { cov = HepSymMatrix(4,0); };
	void write () const;
};

#endif /* Plane_H */
