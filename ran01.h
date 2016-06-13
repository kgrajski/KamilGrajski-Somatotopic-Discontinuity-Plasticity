//	ran01.h

//	Kamil A. Grajski
//	2-November-2014
//

//
//	Header file for a random number generator.
//	From CMU MSCF Financial Computing I class.
//

#ifndef RAN01_H
#define RAN01_H


double	ran0();		//	Return next uniform ( 0, 1 );

double	r1unif ( const double &, const double & );	//	Return uniform in a range.

void	sran0(int);	//	Reset the seed of the generator.

#endif
