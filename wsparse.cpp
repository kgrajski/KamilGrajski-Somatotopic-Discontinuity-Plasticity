//	wsparse.cpp

//	Kamil A. Grajski
//	2-November-2014
//

#include "wsparse.h"

	//
	//	Adapt
	//
void SparseWeightMatrix::Adapt ( const long & prevIter, const long & currentIter, const double & wTauAlpha, const double & wBeta,
									vector<double>::iterator & ix, vector<double>::iterator & iy ) {
	
	vector<long>::iterator iWhich = idListLinear.begin();
	vector<long>::iterator iWhichEnd = idListLinear.end();

	vector<double>::iterator wT = value.begin() + currentIter * numInputWeights;
	vector<double>::iterator wTMinusOne = value.begin() + prevIter * numInputWeights;

	for ( ; iWhich != iWhichEnd; ++iWhich, ++wT, ++wTMinusOne ) {

		*wT = wTauAlpha * *wTMinusOne + wBeta * *iy * *(ix + *iWhich);

	} // 	for ( ; iWhich != iWhichEnd; ++iWhich ) {

}; // 

	//
	//	AdaptRK4:
	//
void SparseWeightMatrix::AdaptRK4 ( const long & prevIter, const long & currentIter, const double & alpha, const double & wBeta, const double & h,
									vector<double>::iterator & ix, vector<double>::iterator & iy ) {
	
		// Plenty of opportunity to do some preallocations and precomputations....

	vector<long>::iterator iWhich;
	vector<double>::iterator iwnp1 = value.begin() + currentIter * numInputWeights;
	vector<double>::iterator iwn;
	vector<double>::iterator iwnStart = value.begin() + prevIter * numInputWeights;
	vector<double>::iterator iwnEnd = value.begin() + prevIter * numInputWeights + numInputWeights;

	vector<double> k1, k2, k3, k4;
	vector<double>::iterator ik1, ik2, ik3, ik4;

	double	hD2 = h / 2.0;
	double	hD6 = h / 6.0;

	k1.assign ( numInputWeights, 0 );
	k2.assign ( numInputWeights, 0 );
	k3.assign ( numInputWeights, 0 );
	k4.assign ( numInputWeights, 0 );

	for ( iwn = iwnStart, ik1 = k1.begin(); iwn != iwnEnd; ) {
		*ik1++ = alpha * *iwn++; 
	} // for ( iy = iyStart; iy != iyEnd; ) {

	for (iwn = iwnStart, ik1 = k1.begin(), ik2 = k2.begin(); iwn != iwnEnd; ) {
		*ik2++ = alpha * ( *iwn++ + hD2 * *ik1++ ); 
	} // for ( iy = iyStart; iy != iyEnd; ) {

	for ( iwn = iwnStart, ik2 = k2.begin(), ik3 = k3.begin(); iwn != iwnEnd; ) {
		*ik3++ = alpha * ( *iwn++ + hD2 * *ik2++ ); 
	} // for ( iy = iyStart; iy != iyEnd; ) {

	for ( iwn = iwnStart, ik3 = k3.begin(), ik4 = k4.begin(); iwn != iwnEnd; ) {
		*ik4++ = alpha * ( *iwn++ + h * *ik3++ ); 
	} // for ( iy = iyStart; iy != iyEnd; ) {

	for ( iwn = iwnStart, ik1 = k1.begin(), ik2 = k2.begin(), ik3 = k3.begin(), ik4 = k4.begin(), iWhich = idListLinear.begin(); iwn != iwnEnd; ) {
		*iwnp1++ = *iwn++ + hD6 * ( *ik1++ + 2.0 * *ik2++ + 2.0 * *ik3++ + *ik4++ ) + wBeta * *iy * *(ix + *iWhich++);
	} // for ( iy = iyStart; iy != iyEnd; ) {

	k1.clear();
	k2.clear();
	k3.clear();
	k4.clear();

}; // 

	//
	//	Clear variables 
	//
void SparseWeightMatrix::ClearVars ( ) {

	N = -1;
	idLinear = -1;
	numInputWeights = 0;
	maxNumInputWeights = 0;
	idListLinear.clear();
	idListRow.clear();
	idListCol.clear();
	value.clear();

}; // void SparseWeightMatrix::LoadParms ( ...

	//
	//	DotProduct (Sparse)
	//
double SparseWeightMatrix::DotProduct ( vector<double>::iterator & x, const long & offset ) {
	
	double z = 0;

	vector<long>::iterator iWhich = idListLinear.begin();
	vector<long>::iterator iWhichEnd = iWhich + numInputWeights;

	vector<double>::iterator iy = value.begin() + offset * numInputWeights;

	for ( ; iWhich != iWhichEnd; ++iWhich, ++iy) {

		z += *(x + *iWhich) * *iy;

	} // 	for ( ; iWhich != iWhichEnd; ++iWhich ) {

	return z;

}; // double SparseWeightMatrix::DotProduct ( ) {

	//
	//	Load variables 
	//
void SparseWeightMatrix::LoadVars ( const long & N_tmp, const long & cellLinearID, const long & inCellLinearID, const long & inCellRowID, const long & inCellColID,
										const long & inMaxNumInputWeights, const double & val ) {

	N = N_tmp;
	idLinear = cellLinearID;
	++numInputWeights;
	maxNumInputWeights = inMaxNumInputWeights;

	idListLinear.push_back ( inCellLinearID );
	idListRow.push_back ( inCellRowID );
	idListCol.push_back ( inCellColID );
	value.push_back ( val );

}; // void SparseWeightMatrix::LoadParms ( ...

	//
	//	Normalize the weight values.  The method operates over the incoming weights, e.g., presynaptic weight normalization.
	//
void SparseWeightMatrix::Normalize ( const long & prevIter, const long & currentIter, const double & wResource ) {

	vector<double>::iterator ix;
	vector<double>::iterator ixStart = value.begin() + prevIter * numInputWeights;
	vector<double>::iterator ixEnd = ixStart + numInputWeights;
	vector<double>::iterator iy = value.begin() + currentIter * numInputWeights;

	double dTmp = 0;
	for ( ix = ixStart; ix != ixEnd; ++ix ) {
		dTmp += *ix;
	} // 	for ( ix = ixStart; ix != ixEnd; ++ix ) {

	dTmp = ( wResource * numInputWeights ) / ( dTmp * maxNumInputWeights );
	for ( ix = ixStart; ix != ixEnd; ++ix, ++iy ) {
		*iy = *ix * dTmp;
	} // for ( ix = ixStart; ix != ixEnd; ++ix, ++iy ) {

}; // void SparseWeightMatrix::Normalize ( const long & offset, const double & wResource ) {


	//
	//	SparseWeightMatrix::WriteToFile
	//		Note: Only the values at the last iteration are written.
	//
void SparseWeightMatrix::ReadFromFile ( fstream & inFile ) {

	double tmp;

	inFile.read ( (char *) &tmp, sizeof ( double ) ); N = (long) tmp;
	inFile.read ( (char *) &tmp, sizeof ( double ) ); idLinear = (long) tmp;
	inFile.read ( (char *) &tmp, sizeof ( double ) ); numInputWeights = (long) tmp;
	inFile.read ( (char *) &tmp, sizeof ( double ) ); maxNumInputWeights = (long) tmp;

	for ( vector<long>::iterator ix = idListLinear.begin(); ix != idListLinear.end(); ++ix ) {
		inFile.read ( (char *) &tmp, sizeof ( double ) );
		*ix = (long) tmp;
	} // for ( vector<long>::iterator ix = idListLinear.begin(); ix != idListLinear.end(); ++ix ) {

	long	offsetToLastIter = ( ( (long) ( value.size() / numInputWeights ) ) - 1 ) * numInputWeights;
	vector<double>::iterator ix = value.begin() + offsetToLastIter;
	for ( ; ix != value.end(); ++ix ) {
		inFile.read ( (char *) &tmp, sizeof ( double ) );
		*ix = tmp;
	} // for ( ; ix != value.end(); ++ix ) {

}; // void WriteToFile ( const string & fName )

	//
	//	Reg to Sparse Conversion
	//
	//	Replace the values stored in the last iteration of values using the input values in wReg.
	//

void SparseWeightMatrix::RegToSparse ( const vector<double>::iterator & wReg, const long & numIters ) {

	vector<long>::iterator iWhich = idListLinear.begin();
	vector<long>::iterator iWhichEnd = iWhich + numInputWeights;

	vector<double>::iterator iy = value.begin() + (numIters - 1) * numInputWeights;

	for ( ; iWhich != iWhichEnd; ) {

		*iy++ = *(wReg + *iWhich++);

	} // 	for ( ; iWhich != iWhichEnd; ++iWhich ) {

}; // void SparseWeightMatrix::RegToSparse ( const vector<double> & ) {

	//
	//	Trial fill-in.
	//	Extend the Lists numIter times.
	//	This is for the case of weight adaptation and we want to keep the evolution at each time step.
	//
void SparseWeightMatrix::TrialFillIn ( const long & numIters ) {

	value.resize ( numIters * numInputWeights, 0 );

}; // void SparseWeightsMatrix::TrialFillIn ( const long & numIters ) {

	//
	//	SparseWeightMatrix::WriteToFile
	//		Note: Only the values at the last iteration are written.
	//
void SparseWeightMatrix::WriteToFile ( fstream & outFile ) {

	double	tmp;

	tmp = N; outFile.write ( (char *) &tmp, sizeof ( double ) );
	tmp = idLinear; outFile.write ( (char *) &tmp, sizeof ( double ) );
	tmp = numInputWeights; outFile.write ( (char *) &tmp, sizeof ( double ) );
	tmp = maxNumInputWeights; outFile.write ( (char *) &tmp, sizeof ( double ) );

	for ( vector<long>::iterator ix = idListLinear.begin(); ix != idListLinear.end(); ++ix ) {
		tmp = *ix;
		outFile.write ( (char *) &tmp, sizeof ( double ) );
	} // 	for ( vector<long>::iterator ix = idListLinear.begin(); ix != idListLinear.end(); ++ix ) {

	long offsetToLastIter = ( ( (long) ( value.size() / numInputWeights ) ) - 1 ) * numInputWeights;
	vector<double>::iterator ix = value.begin() + offsetToLastIter;
	for ( ; ix != value.end(); ++ix ) {
		tmp = *ix;
		outFile.write ( (char *) &tmp, sizeof ( double ) );
	}  // 	for ( ; ix != value.end(); ++ix ) {

}; // void SparseWeightMatrix::WriteToFile ( const fstream & outFile ) {