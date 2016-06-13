//	wsparse.h

//	Kamil A. Grajski
//	19-November-2014
//

//
//	Header file for the class of sparse representation of a weight matrix.
//

#ifndef WSPARSE_H
#define WSPARSE_H

#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	//
	//  START: wsparse class declarations.
	//
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////

class SparseWeightMatrix {

public:
	
	SparseWeightMatrix () { ClearVars (); };
	~SparseWeightMatrix () { ClearVars (); };

	void			Adapt ( const long &, const long &, const double &, const double &, vector<double>::iterator &, vector<double>::iterator & );
	void			AdaptRK4 ( const long &, const long &, const double &, const double &, const double &, vector<double>::iterator &, vector<double>::iterator & );
	void			ClearVars ();
	double			DotProduct ( vector<double>::iterator &, const long & );
	void			LoadVars ( const long &, const long &,
								const long &, const long &, const long &, const long &, const double & );
	void			Normalize ( const long &, const long &, const double & );
	void			ReadFromFile ( fstream & );
	void			RegToSparse ( const vector<double>::iterator &, const long & );
	void			SetInitialConditions ( const long &, const long &, const double &, const double & );
	void			TrialFillIn ( const long & );
	void			WriteToFile ( fstream & );
			
			//	These variables are used to identify the cell in within the weight matrix.
			//	There is some redundancy to use the row & column, and the index of the equivalent linear representation.
			//	A cell may have multiple sparce weight matrices associated with it, but they are indepdendent (as implemented here).

	long			N;		//	Of the weight matrix.
	long			idLinear;

			//	These variables identify the source of the incoming weight.
			//	There is some redundancy to use the row & column, and the index of the equivalent linear representation.

	long			numInputWeights;
	long			maxNumInputWeights;	//	Needed during normalization taking into account edge effects.
	vector<long>	idListRow;
	vector<long>	idListCol;
	vector<long>	idListLinear;
	vector<double>	value;

private:

}; // class SparseWeightMatrix {

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	//
	//  End: SparseWeightMatrix class declarations.
	//
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////

#endif