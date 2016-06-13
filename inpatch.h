//	inpatch.h

//	Kamil A. Grajski
//	2-November-2014
//

//
//	Header file for the input patch class
//

#ifndef INPATCH_H
#define INPATCH_H

#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	//
	//  START: InputPatch class declarations.
	//
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////



class InputPatch {

public:
	
	InputPatch () : digit(0), nodePatch(0), columnPatchLoc(0), rowPatchLoc(0),
					widthPatch(0), lengthPatch(0), durationPatchInIters(0), offsetPatchInIters(0), ampPatch(0), inputPattern(0) {};

	~InputPatch () {};

	void	Initialize ();
	void	LoadParms ( const long &, const long &, const long &, const long &, const long &, const long &, const long &, const long &, const long &,
		const long &, const double & );

	long			digit;
	long			nodePatch;
	long			columnPatchLoc;
	long			rowPatchLoc;
	long			widthPatch;
	long			lengthPatch;
	long			durationPatchInIters;
	long			offsetPatchInIters;
	long			numValsPerIter;
	long			numValsPerTrial;
	double			ampPatch;

	vector<double>	stimCount;

	vector<double>	inputPattern;			//	Original implementation
	vector<double>	inputPatternSparse;		//	Sparse representation: the list of nodes that are ON for durationPatchInIters starting at offsetPatchInIters

private:

}; // class InputPatch


	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	//
	//  End: InputPatch class declarations.
	//
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////

#endif