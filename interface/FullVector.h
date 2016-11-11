/*!
 * \file FullVector.h
 * \brief FullVector = StateVector + covariance matrix
 *
 * 5x5 covariance matrix added.
 * Space for storing residuals is provided also.
 *
 * \author Wolfgang Waltenberger
 * \date Thu 05 Jul 2001 15:15:46 CEST
 */

#ifndef FullVector_H
#define FullVector_H

#include "DetElement.h"
#include "defines.h"
#include "Plane.h"
#include "StateVector.h"
#include <CLHEP/config/iostream.h>
#include "CLHEP/config/CLHEP.h"
#include "CLHEP/Matrix/SymMatrix.h"

#ifdef GL
#include <Inventor/nodes/SoSelection.h>
#endif

class FullVector : public StateVector {
private:
	HepSymMatrix _CMSCov, _DelphiCov;
	bool useCMSCov, useDelphiCov;
	void CMS2DelphiCov();
	void Delphi2CMSCov();
	#ifdef GL
	void drawonCyl( SoSelection *, const string &, HepDouble );
	void drawonDis( SoSelection *, const string &, HepDouble );
	#endif
public:
	FullVector ( const DetElement & );
	FullVector ( const DetElement & , HepDouble , HepDouble , HepDouble , \
			HepDouble , HepDouble );
	FullVector ();
	
	/// Inserts parameter from a Plane into our
	/// FullVector.
	/// \param cov get covariance matrix also?
	//
	/*
	void getPlane ( RealPlane , DetElement );
	void getPlane ( RealPlane , DetElement, HepDouble r );
	*/
	
	void setDelphiCov ( const HepSymMatrix & );
	void setCMSCov ( const HepSymMatrix & );
	HepSymMatrix DelphiCov ();
	HepSymMatrix CMSCov ();
	void writeCMS();
	void writeDelphi();

	#ifdef GL
	HepSymMatrix track_cov;
	void GLdraw( SoSelection * , HepDouble lngth );
	void GLdraw( SoSelection * , HepDouble, HepDouble, HepDouble , HepDouble );
	void GLdraw( SoSelection * , const string & , HepDouble );
	void GLdraw( SoSelection * , const string & , HepDouble, \
			HepDouble, HepDouble , HepDouble );
	string getInfo ( int ) ;
	#endif
};
#endif
