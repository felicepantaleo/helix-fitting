#ifndef HepMoreOperators_H
#define HepMoreOperators_H

#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/DiagMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include <vector>

/// find the minimum value. Assume m to be diagonal and of dimension n x n
/// Warning: We do not check for errors!
/// Use with care!
inline int minimum_sym ( HepMatrix m )
{
	int k=1;
	for (int i=2; i<=m.num_row(); i++) {
		if ( m(i,i)<m(k,k) ) k=i;
	};
	return k;
}

inline HepVector operator +( const HepVector &v , const HepDouble &d )
{
	HepVector ret(v.num_row()) ;
	for ( int i=0; i< v.num_row() ; i++ ) ret[i]=v[i]+d;
	return ret;
};

inline HepVector operator +( const vector <HepDouble> &v , const HepDouble &d )
{
	HepVector ret(v.size()) ;
	for ( unsigned int i=0; i< v.size() ; i++ ) ret[i]=v[i]+d;
	return ret;
};

// haddamard (element wise / point wise ) product
inline HepVector operator &( const HepVector &d, const HepVector &v )
{
	HepVector ret;
	for ( int i=0; i< d.num_row() ; i++ )
		ret[i]=d[i]*v[i];
	return ret;
};

inline HepVector operator +( const HepDouble &d, const vector <HepDouble> &v )
{
	return v+d;
};

inline HepVector operator +( const HepDouble &d, const HepVector &v )
{
	return v+d;
};

inline HepVector sqrt ( const HepVector v)
{
	HepVector ret ( v.num_row() );
	for ( int i=0; i< v.num_row(); i++ )
		ret[i]=sqrt(v[i]);
	return ret;
};

inline HepVector sqr ( const HepVector v)
{
	HepVector ret ( v.num_row() );
	for ( int i=0; i< v.num_row(); i++ )
		ret[i]=sqr(v[i]);
	return ret;
};

inline vector <HepDouble> sqr ( vector <HepDouble> v)
{
	vector <HepDouble> ret ( v.size() );
	for ( unsigned int i=0; i< v.size(); i++ )
		ret[i]=sqr(v[i]);
	return ret;
};

inline vector<HepDouble> sqrt ( vector <HepDouble> v)
{
	vector <HepDouble> ret ( v.size() );
	for ( unsigned int i=0; i< v.size(); i++ )
		ret[i]=sqrt(v[i]);
	return ret;
};


inline HepVector norm ( const HepVector v , const HepVector u )
{
	HepVector ret ( v.num_row() );
	for ( int i=0; i< v.num_row(); i++ )
		ret[i]=sqrt(sqr(v[i])+sqr(u[i]));
	return ret;
};

inline HepVector cos ( const HepVector v)
{
	HepVector ret ( v.num_row() );
	for ( int i=0; i< v.num_row(); i++ )
		ret[i]=cos(v[i]);
	return ret;
};

inline HepVector sin ( const HepVector v)
{
	HepVector ret ( v.num_row() );
	for ( int i=0; i< v.num_row(); i++ )
		ret[i]=sin(v[i]);
	return ret;
};

/// make a HepSymMatrix out ouf a (hopefully symmetric) HepMatrix
/// Warning! We do not check for errors!
/// Use with care!
inline HepSymMatrix symmetrize (HepMatrix m)
{
	int i,j;
	HepSymMatrix ret(m.num_row(),0);
	for (i=0; i<m.num_row(); i++) {
		for (j=i; j<m.num_col(); j++) {
			ret[i][j]=m[i][j];
		};
	};
	return ret;
};

#endif
