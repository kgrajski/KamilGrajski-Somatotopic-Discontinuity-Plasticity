//	inpatch.cpp

//	Kamil A. Grajski
//	2-November-2014
//

#include "inpatch.h"

void InputPatch::Initialize () {

	digit = 0;
	nodePatch = 0;
	columnPatchLoc = 0;
	rowPatchLoc = 0;
	widthPatch = 0;
	lengthPatch = 0;
	durationPatchInIters = 0;
	offsetPatchInIters = 0;
	numValsPerIter = 0;
	numValsPerTrial = 0;
	ampPatch = 0;
	inputPattern.clear();
	inputPatternSparse.clear();
	stimCount.clear();

}; // void Initialize ()

void InputPatch::LoadParms ( const long & iDigit, const long & xNodePatch, const long & iCol, const long & iRow, const long & xWidthPatch,
							const long & xLengthPatch, const long & xDurationPatchInIters, const long & xOffsetPatchInIters, const long & xNumValsPerIter,
											const long & xNumValsPerTrial, const double & xAmpPatch ) {
	digit = iDigit;
	nodePatch = xNodePatch;
	columnPatchLoc = iCol;
	rowPatchLoc = iRow;
	widthPatch = xWidthPatch;
	lengthPatch = xLengthPatch;
	durationPatchInIters = xDurationPatchInIters;
	offsetPatchInIters = xOffsetPatchInIters;
	numValsPerIter = xNumValsPerIter;
	numValsPerTrial = xNumValsPerTrial;
	ampPatch = xAmpPatch;

}; // void Initialize ()
