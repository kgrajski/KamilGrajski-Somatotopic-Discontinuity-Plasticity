//	ran01.cpp

//	Kamil A. Grajski
//	2-November-2014
//

//
//	Random number generator.
//	From CMU MSCF Financial Computing I class.
//
#include "simmgr.h"

using namespace std;

static int idum = 0;	// static means global, but only within this file.

inline double ran0() {

	const int IA = 16807, IM = 2147483647, IQ = 127773;
	const int IR = 2836, MASK = 123459876;
	const double AM = 1.0 / double(IM);
	int k; double ans;

	idum ^= MASK;			//	XORing with MASK allwos use of zero and simple bit patterns for idum
	k = idum / IQ;

	idum = IA * ( idum - k * IQ ) - IR * k;
	if ( idum < 0 ) idum += IM;

	ans = AM * idum;
	idum ^= MASK;

	//ans = 0;

	return ans;

}; // double ran0()

double r1unif ( const double & xmin, const double & xmax ) {

	double ans = ran0();
	ans = xmin + (xmax - xmin) * ans;

	//ans = 0;
	return ans;

}; // double r1unif ( const double & xmin, const double & xmax ) {

void sran0 ( int seed ) {
	idum = seed;			//	Reset the RNG seed.
}; // void sran0 ( int seed )