// $Id: any2str.h,v 1.1 2001/11/16 10:29:53 wwalten Exp $

#include <stdio.h>
#include <string>
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/DiagMatrix.h"
#include "CLHEP/Matrix/Vector.h"

#ifdef SOLARIS
using namespace std;
#endif

string any2str ( char *st );

string any2str ( string st, HepMatrix x , int prec );

// put large matrices into a string
string any2str ( string st, HepSymMatrix x , int prec );

string any2str ( char *pre, HepSymMatrix x );

string any2str ( char *pre, HepSymMatrix x , int prec );

string any2str ( string pre, HepSymMatrix x );

string any2str ( char *pre, HepMatrix x );

string any2str ( HepMatrix x , int prec );

string any2str ( HepMatrix x );

string any2str ( char *pre, HepMatrix x , int prec );

string any2str ( string pre, HepMatrix x );

string any2str ( int any );
string any2str ( bool any );
string any2str ( unsigned int any );
string any2str ( long int any );

string any2str ( double any );
string any2str ( double any, int pre );
