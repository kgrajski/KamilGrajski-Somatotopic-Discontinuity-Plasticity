//	expmgr.cpp

//	Kamil A. Grajski
//	2-November-2014
//

#include "ioutil.h"

		//
void LoadNumericVectorFromBinaryFile ( const string & fName, fstream & inFile, const long & fileManagedByCaller, vector<double> & x, const long & iLen ) {

	double tmp;
	long k = 0;

	if ( fileManagedByCaller ) {
		for ( vector<double>::iterator ix = x.begin(); k < iLen; ++ix, ++k ) {
			inFile.read ( (char *) &tmp, sizeof ( double ) );
			*ix = tmp;
		}	// for ( vector<double>::iterator ix = x.begin(); ix != x.end(); ++ix )		
	} else {
		fstream readFile;
		readFile.open( fName.c_str(), ios::in | ios::binary );
		for ( vector<double>::iterator ix = x.begin(); k < iLen; ++ix, ++k ) {
			readFile.read ( (char *) &tmp, sizeof ( double ) );
			*ix = tmp;
		}	// for ( vector<double>::iterator ix = x.begin(); ix != x.end(); ++ix )
		readFile.close();
	} // if ( fileManagedByCaller ) {
	
}; // void LoadNumericVectorFromBinaryFile ( vector<double> &, ifstream & )

void SaveNumericVectorToBinaryFile ( const string & fName, const long & createFileFlag, vector<double> & x, const long & offset, const long & iLen  ) {
	
	vector<double>::iterator ix = (x.begin() + offset);
	vector<double>::iterator ixEnd = ix + iLen;
	double tmp;
	fstream outFile;

	if ( createFileFlag ) {
		outFile.open( fName.c_str(), ios::out | ios::binary );
	}
	else {
		outFile.open( fName.c_str(), ios::out | ios::binary | ios::app );
	} // if ( createFileFlag )

	for ( ; ix != ixEnd; ++ix ) {
			tmp = *ix;
			outFile.write ( (char *) &tmp, sizeof ( double ) );
	}	// for ( vector<double>::iterator ix = x.begin(); ix != x.end(); ++ix )
	outFile.close();

}; // void saveNumericVectorToBinaryFile ( const string & fName, const long & createFileFlag, vector<double> & x, const long & offset, const long & iLen  ) {