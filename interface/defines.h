/*
 * commonly used definitions
 */

#include <iostream>
#include <stdio.h>
#include <string>
#include "any2str.h"
#include <pthread.h>
#include "CLHEP/config/CLHEP.h"

#ifndef Defines_H
#define Defines_H
//#include "CLHEP/Matrix/SymMatrix.h"

/*
#ifdef SOLARIS
typedef void *(*thread_function_t)(void *);
#else
typedef void *(*thread_function_t)(void *);
#endif
*/
#ifndef ETCDIR
#define ETCDIR "./etc"
#endif

#include "version.h"

extern char verbosity;

using namespace std;
// #define Z_SCALE .5
#define Z_SCALE 1

// how the B field in [T] is 'scaled' into a length
#define B_SCALE 14

#ifdef COLORS
#define C_GREEN "[32;11m"
#define C_RED "[31;11m"
#define C_YELLOW "[33;11m"
#define C_MAGENTA "[35;11m"
// #define C_INVERTED "8*7"
#define C_INVERTED "*"
#define C_DEFAULT "[0m"
#else
#define C_GREEN ""
#define C_RED ""
#define C_YELLOW ""
#define C_MAGENTA ""
#define C_INVERTED "*"
#define C_DEFAULT ""
#endif

#define C_UP "[A"
#define C_RIGHT "[C"

#ifdef HEPVIS
	#ifndef GL
		#define GL
	#endif
#endif
inline double sqr ( double x )
{
	return x*x;
};

inline double norm ( double x, double y )
{
	return sqrt ( sqr (x) + sqr (y) );
};

inline int sign ( double x )
{
	return x >= 0 ? 1 : -1;
};

template <class T> inline void warning( T i, string file, int line) {
	#ifndef DEBUG
	if ( line ) {
		cerr << C_YELLOW << "Warning: " << i << " at " << file \
			<< ":" << line << C_DEFAULT << endl;
	};
	#else
	static string oldstring="";
	static int count=0;
	if (oldstring==i) {
		count++;
	} else {
		if (count) {
			cerr << C_YELLOW << "(Last warning message repeated " << ntimes(count) 
				<< ".)" << C_DEFAULT << endl;
			count=0;
		};
		if ((string) i!="") {
			cerr << C_YELLOW << "Warning: " << i << " at " << file \
				<< ":" << line << C_DEFAULT << endl;
		};
		oldstring=i;
	};
	#endif
};

struct interval {
	HepDouble min;
	HepDouble max;
};

template <class T> inline void error ( T i, string file, int line ) {
	#ifndef DEBUG
	if (line ) {
		cerr << C_RED << "Error: " << i << " at " << file \
			<< ":" << line << endl;
		cerr << "Turn debugging on!" << C_DEFAULT << endl;
	};
	#else
	static string oldstring="";
	static int count=0;
	if (oldstring==i) {
		count++;
	} else {
		if (count) {
			cerr << C_RED << "(Last error message repeated " << ntimes(count) 
				<< ".)" << C_DEFAULT << endl;
			count=0;
		};
		if ((string) i!="") {
			cerr << C_RED << "Error: " << i << " at " << file \
					<< ":" << line << C_DEFAULT << endl;
		};
		oldstring=i;
	};
	#endif
};

#define WARNING(i) warning (i,__FILE__,__LINE__)
#define ERROR(i) error (i,__FILE__,__LINE__)

inline string ntimes ( int count )
{
	switch (count) {
		case 0: return "never";
		case 1: return "once";
		case 2: return "twice";
		default: return any2str(count)+" times";
	};
};

inline void MESSAGE ( char n, string as ) {
	#ifdef DEBUG
	static string oldstring="";
	static int count=0;
	if (oldstring==as) {
		count++;
	} else {
		if (count) {
			cerr << C_YELLOW << "(Last message repeated " << ntimes(count) 
				<< ".)" << C_DEFAULT << endl;
			count=0;
		};
		if (verbosity > n && (string) as!="") {
			cerr << C_YELLOW << "Message: " << as << C_DEFAULT << endl;
			oldstring=as;
		};
	};
	#endif
};

inline void FINISH_COUTS() {
	MESSAGE(0,"");
	error("",__FILE__, 0);
	warning("",__FILE__, 0);
};

inline double myatan ( double y, double x )
{
	#ifdef DEBUG
	if (y==0.) {
		WARNING ( "y=0" );
		return atan2 ( 0, 1 );
	} else {
		return atan2 ( y , x );
	};
	#else
	return atan2 ( y , x );
	#endif
};

#define MLC 255 // MCL Maximum Length for a Configuration line
#define CNV 0.299792458 // [Gev / ( c m T )] conversion

#endif /* Defines_H */
