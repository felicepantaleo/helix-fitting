#include "hfile.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include HFILE"DetElement.h"
#include HFILE"Detector.h"
#include HFILE"Hit.h"
#include <vector>
#include HFILE"Plane.h"
#include HFILE"Track.h"
#include HFILE"defines.h"
#include <CLHEP/config/iostream.h>
 
int main (int argc, char **argv)
{
	bool forward=false;
	if (argc>1) if (!strcmp(argv[1],"f")) forward=true;
	string conffile;
	if (forward) {
		conffile=ETCDIR"/forward.conf"; // default config file
	} else {
		conffile=ETCDIR"/barrel.conf"; // default config file
	};

	Detector barrel( conffile );
	StateVector sv ( barrel.Element[0] );
	HepDouble d_a[5];
	// Delphi2CMS
	if (forward) {
		d_a[0]=2.77102; d_a[1]=0.0442252; d_a[2]=.3; 
		d_a[3]=2.80696; d_a[4]=1.66991;
	} else {
		d_a[0]=.054147; d_a[1]=.056305; d_a[2]=.5; 
		d_a[3]=.0000824479; d_a[4]=-.500253;
	};
	sv.setDelphi ( d_a[0], d_a[1], d_a[2], d_a[3], d_a[4] );
		
	sv.writeDelphi();
	sv.writeCMS();
	
	// CMS2Delphi
	HepDouble c_a[5];
	c_a[0]=sv.a0(); c_a[1]=sv.z(); c_a[2]=sv.theta(); 
	c_a[3]=sv.psi(); c_a[4]=sv.kappa_cms();
	sv.clear();
	sv.setCMS( c_a[0], c_a[1], c_a[2], c_a[3], c_a[4] );
	sv.writeDelphi();
	sv.writeCMS();

	// Delphi2CMS
	HepDouble d_b[5];
	if (forward) {
		d_b[0]=sv.Phi();
		d_b[1]=sv.R();
		d_b[2]=sv.theta();
		d_b[3]=sv.phi();
		d_b[4]=sv.kappa();
	} else {
		d_b[0]=sv.Phi();
		d_b[1]=sv.z();
		d_b[2]=sv.theta();
		d_b[3]=sv.beta();
		d_b[4]=sv.kappa();
	};
	sv.clear();
	sv.setDelphi(d_b[0],d_b[1],d_b[2],d_b[3],d_b[4]);
	sv.writeDelphi();
	sv.writeCMS();
	HepDouble c_b[5];
	c_b[0]=sv.a0(); c_b[1]=sv.z(); c_b[2]=sv.theta();
	c_b[3]=sv.psi(); c_b[4]=sv.kappa_cms();
	
	cout << "Result: " << endl;
	bool ok=true;
	for ( int i=0; i< 5 ; i++ ) {
		if ( c_a[i] < c_b[i]-pow(10,-9) || c_a[i] > c_b[i]+pow(10,-9) ) {
			ok=false;
			cout << "CMS coord #" << i+1 << " "
				<< c_a[i] << "<>" << c_b[i] << endl;
		};
		if ( d_a[i] < d_b[i]-pow(10,-9) || d_a[i] > d_b[i]+pow(10,-9) ) {
			ok=false;
			cout << "Delphi coord #" << i+1 << " "
				<< d_a[i] << "<>" << d_b[i] << endl;
		};
	};
	if (ok) { cout << "OK!" << endl; };
	
};
