// 

#include <stdio.h>
#include <string>
#include <vector>
#include "defines.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/DiagMatrix.h"
#include "CLHEP/Matrix/Vector.h"

using namespace std;

string any2str ( char *st )
{
	return (string) st;
};

// put large matrices into a string
string any2str ( string st, HepMatrix x , int prec )
{
	HepDouble min=100,max=-100,multfak=1, dx;
//	string st;
//	st=pre;
	dx=pow(10,-2*prec); // when is something identical 0?
	char tmp[79]="";
	for (int i=1; i<= x.num_row(); i++) {
		for (int j=1; j<= x.num_col(); j++) {
			if (x(i,j)<min) min=x(i,j);
			if (x(i,j)>max) max=x(i,j);
		};
	};
	if (max==0 && min==0) {
		sprintf(tmp,"%u",x.num_row());
		st+="["+(string) tmp+"x"+(string) tmp+"] zero matrix.";
		return st;
	};
	multfak=pow(10,floor(log10(1/max)));
	snprintf(tmp,79,"[%dx%d] %.1e x\n",x.num_row(),x.num_row(),1/multfak);
	st += (string) tmp;
	
	for (int i=1; i<= x.num_row(); i++) {
		for (int j=1; j<= x.num_col(); j++) {
				if (x(i,j) == 0 ) { 
					for (int i=0; i<prec+2;i++) st+=" ";
					st+="0 ";
				} else {
					HepDouble temp=x(i,j)*multfak;
					sprintf(tmp," %*.*f",prec+3,prec,temp);
					
					st+= (string) tmp;
				};
//			};
		};
		if (i<x.num_row()) st+='\n';
	};
	return st;
};
// put large matrices into a string
string any2str ( string st, HepSymMatrix x , int prec )
{
	HepDouble min=100,max=-100,multfak=1, dx;
//	string st;
//	st=pre;
	dx=pow(10,-2*prec); // when is something identical 0?
	char tmp[79]="";
	for (int i=1; i<= x.num_row(); i++) {
		for (int j=1; j<= x.num_col(); j++) {
			if (x(i,j)<min) min=x(i,j);
			if (x(i,j)>max) max=x(i,j);
		};
	};
	if (max==0 && min==0) {
		sprintf(tmp,"%u",x.num_row());
		st+="["+(string) tmp+"x"+(string) tmp+"] zero matrix.";
		return st;
	};
	multfak=pow(10,floor(log10(1/max)));
	snprintf(tmp,79,"[%dx%d] %.1e x\n",x.num_row(),x.num_row(),1/multfak);
	st += (string) tmp;
	
	for (int i=1; i<= x.num_row(); i++) {
		for (int j=1; j<= x.num_col(); j++) {
			if (j<i) { 
					for (int i=0; i<prec+4;i++) st+=" ";
			} else {
//				if (x(i,j)>= -dx && x(i,j) <= dx) { 
				if (x(i,j) == 0 ) { 
					for (int i=0; i<prec+2;i++) st+=" ";
					st+="0 ";
				} else {
					HepDouble temp=x(i,j)*multfak;
					sprintf(tmp," %*.*f",prec+3,prec,temp);
					
					st+= (string) tmp;
				};
			};
		};
		if (i<x.num_row()) st+='\n';
	};
	return st;
};

string any2str ( char *pre, HepSymMatrix x )
{
	string me=pre;
	return any2str( me, x, 4);
};

string any2str ( char *pre, HepSymMatrix x , int prec )
{
	return any2str( (string) pre, x, 4);
};

string any2str ( string pre, HepSymMatrix x )
{
	return any2str( pre, x, 4);
};

string any2str ( char *pre, HepMatrix x )
{
	return any2str( (string) pre, x, 4);
};


string any2str ( HepMatrix x , int prec )
{
	return any2str( "", x, 4);
};

string any2str ( HepMatrix x )
{
	return any2str( "", x, 4);
};

string any2str ( char *pre, HepMatrix x , int prec )
{
	return any2str( (string) pre, x, 4);
};

string any2str ( string pre, HepMatrix x )
{
	return any2str( pre, x, 4);
};

string any2str ( long int any )
{
	#ifdef SOLARIS
	if (ceil(log10(labs(any))) > 10 ) {
		ERROR("fixme. we have a problem with solaris here.");
		exit (-1);
	};
	char ret[10];
	snprintf(ret,10,"%l",any);
	#else
	int dec=(int) ceil(log10(labs(any)))+1;
	char ret[dec];
	snprintf(ret,dec,"%ld",any);
	#endif
	return ret;
};

string any2str ( bool any )
{
	return any ? "true" : "false";
};

string any2str ( int any )
{
	if ( any==1 ) { return "1"; };
	return any2str( (long int) any );
};

string any2str ( unsigned int any )
{
	return any2str( (long int) any );
};

string any2str ( double any, int pre )
{
	char ret[17];
	snprintf(ret,17,"%.*f",pre,any);
	return ret;
};

string any2str ( double any )
{
	char ret[17];
	snprintf(ret,17,"%f",any);
	return ret;
};
