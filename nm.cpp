//	nm.cpp

//	Kamil A. Grajski
//	2-November-2014
//

#include "nm.h"
#include "inpatch.h"

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	//
	//  START: NM class definitions.
	//
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////

	//
	//	AddNoise.
	//
void NM::AddNoise (  vector<double> & x, const long & offset, const long & iLen, const double & xMin, const double & xMax ) {
	vector<double>::iterator ix = x.begin() + offset;
	vector<double>::iterator iEnd = ix + iLen;
	for ( ; ix != iEnd; ++ix ) {
		*ix += r1unif ( xMin, xMax );
	}
}; // void NM::setInitialConditions (  const vector<double> & x, const long & iLen, const double & xMin, const double & xMax )

	//
	//	Vector Scalar Multiplication: y_bar = a * x_bar
	//
	//AddVecScalarMult ( w1E0, offsetToCurrentIterWeights, w1E0, offsetToPreviousIterWeights, wTauAlpha, numWeightsPerIter );
void NM::AddVecScalarMult ( vector<double> & y, long & yOffset, vector <double> & x, long & xOffset, double & scalarVal, long & iLen ) {

	vector<double>::iterator ix = x.begin() + xOffset;
	vector<double>::iterator iEnd = ix + iLen;
	vector<double>::iterator iy = y.begin() + yOffset;

	for ( ; ix != iEnd; ++ix, ++iy ) {
		*iy += ( *ix * scalarVal );
	}

}; // void NM::VecScalarMult ( vector<double> & y, long & yOffset, vector <double> & x, long & xOffset, double & scalarVal, long & iLen ) {

	//
	//	AddVecMatMult: y_bar += v_bar * M.  (Everything in linear representation.)
	//
void NM::AddVecMatMult ( vector<double> & y, const long & yOffset, vector<double> & v, const long & vOffset, vector<double> & M, const long & mOffset, const long & iLen ) {
	vector<double>::iterator iy = y.begin() + yOffset;
	vector<double>::iterator iyEnd = iy + iLen;
	vector<double>::iterator im = M.begin() + mOffset;

	for ( ; iy != iyEnd; ++iy ) {
		
		vector<double>::iterator iv = v.begin() + vOffset;
		vector<double>::iterator ivEnd = iv + iLen;
		double tmp = 0;

		for ( ; iv != ivEnd; ++iv, ++im ) {

				tmp += ( *iv * *im );

		} // for ( ; iv != ivEnd; ++iv )

		*iy += tmp;

	} // for ( ; iy != iyEnd; ++iy )

};  // void AddVecMatMult ( vector<double> &, const long &, vector<double> &, const long &, vector<double> &, const long &, const long & ) {

	//
	//	Baseline Refinement
	//
void NM::BaselineRefinement () {

	long iStim = 0;
	long numStim = 0;
	vector<InputPatch> inputPatches;
	vector<long> orderInputPatches;

		//	Generate the input patches.
	GenDig3RefineInputPatchList ( inputPatches );
		
		//	Test
	inputPatches.clear();
	GenDig3SyndactylyInputPatchList ( inputPatches );

		//	Test
	inputPatches.clear();
	GenDig3SyndactylyControlInputPatchList ( inputPatches );

		//	Turn ON weight adaptation.
	plasticityFlag = 1;

		//	Determine a random ordering of the inputs.
	orderInputPatches.clear();
	RandomPermutation ( orderInputPatches, 0, (long) inputPatches.size() );
	numStim = (long) inputPatches.size();

		//	Apply inputs in random order.

	for ( long iStim = 0; iStim < numStim; ++iStim ) {

		for ( long iter = 1; iter < numIterPerTrial; ++iter ) {

			long iOffset = ( iter - 1 ) * numValsPerIter;

			VecCopy ( inputPatches[ orderInputPatches[iStim] ].inputPattern, iOffset, v0, iOffset, numValsPerIter );
			AddNoise (  v0, iOffset, numValsPerIter, -noiseLevel, noiseLevel );
			Sigmoid ( v0, r0, beta, iOffset, numValsPerIter );

			if ( networkConfigType == 5 ) {
				IterateSystem5 ( iter );
			} else {
				IterateSystem4 ( iter );
			} // if ( networkConfigType == 5 ) {

		} // for ( long iter = 1; iter < numIterPerTrial; ++iter ) {
		CopyLastIterValsToZeroIter();
		CopyLastIterWeightsToZeroIter();
		
	} // for ( iStim = 0; iStim < numStim; ++iStim

		//	Turn OFF weight adaptation.
	plasticityFlag = 0;

		//	Clean up.
	CleanUpDig3RefineInputPatchList ( inputPatches );
	orderInputPatches.clear();

}; // void NM::BaselineRefinement

	//
	//	CleanUpDig3RefineInputPatchList
	//
void NM::CleanUpDig3RefineInputPatchList ( vector<InputPatch> & x ) {

	vector<InputPatch>::iterator ix = x.begin();
	vector<InputPatch>::iterator ixEnd = x.end();

	for ( ; ix != ixEnd; ++ix ) {
		ix->inputPattern.clear();
		ix->stimCount.clear();
	} // for ( ; ix != ixEnd; ++ix ) {
	
	x.clear();

}; // void NM::GenDig3RefineInputPatchList ( vector<InputPatch> & ) {

	//
	//	CleanUpDig3SyndactylyInputPatchList
	//		A bit wasteful as for now is doing exactly same thing as CleanUpDig3RefineInputPatchList,
	//		but that may change at some point.
	//
void NM::CleanUpDig3SyndactylyInputPatchList ( vector<InputPatch> & x ) {

	vector<InputPatch>::iterator ix = x.begin();
	vector<InputPatch>::iterator ixEnd = x.end();

	for ( ; ix != ixEnd; ++ix ) {
		ix->inputPattern.clear();
		ix->stimCount.clear();
	} // for ( ; ix != ixEnd; ++ix ) {
	
	x.clear();

}; // void NM::CleanUpDig3SyndactylyInputPatchList ( vector<InputPatch> & ) {

	//
	//	CompareNE
	//
void NM::CompareNE ( vector<double> & x, vector<double> & y, const double & xTestValue, const long & offset, const long & iLen ) {
	
	vector<double>::iterator ix = x.begin() + offset;
	vector<double>::iterator iEnd = ix + iLen;
	vector<double>::iterator iy = y.begin();

	for ( ; ix != iEnd; ++ix, ++iy ) {
		if ( *ix != xTestValue ) {
			*iy = 1;
		} else {
			*iy = 0;
		} // if ( *ix != xTestValue )
	}

}; // void NM::CompareNE ( vector<double> & x, vector<double> & y, const double xTestValue, const long & offset, const long & iLen )

	//
	//	CopyLastIterValsToZeroIter
	//
void NM::CopyLastIterValsToZeroIter ( ) {

	long offsetToLastIter = ( numIterPerTrial - 1 ) * numValsPerIter;
	VecCopy ( v0, offsetToLastIter, v0, 0, numValsPerIter );
	VecCopy ( r0, offsetToLastIter, r0, 0, numValsPerIter );
	VecCopy ( v1E, offsetToLastIter, v1E, 0, numValsPerIter );
	VecCopy ( r1E, offsetToLastIter, r1E, 0, numValsPerIter );
	VecCopy ( v1I, offsetToLastIter, v1I, 0, numValsPerIter );
	VecCopy ( r1I, offsetToLastIter, r1I, 0, numValsPerIter );

}; // void NM::CopyLastIterValsToZeroIter ( ) {

	//
	//	CopyLastIterWeightsToZeroIter
	//
void NM::CopyLastIterWeightsToZeroIter ( ) {

	long offsetToLastIter = ( numIterPerTrial - 1 ) * numWeightsPerIter;
	VecCopy ( w1E0, offsetToLastIter, w1E0, 0, numWeightsPerIter );
	VecCopy ( w1EI, offsetToLastIter, w1EI, 0, numWeightsPerIter );
	VecCopy ( w1IE, offsetToLastIter, w1IE, 0, numWeightsPerIter );
	VecCopy ( w1EE, offsetToLastIter, w1EE, 0, numWeightsPerIter );

}; // void NM::CopyLastIterWeightsToZeroIter ( ) {

	//
	//	GenDig3RefineInputPatchList
	//
void NM::GenDig3RefineInputPatchList ( vector<InputPatch> & inputPatchList ) {

	InputPatch tmpInputPatch;
	vector<double> tmpVec;
	vector<double> tmpStimCount;

	tmpInputPatch.Initialize();

	long numDigits = 3;
	long digitWidth = N / numDigits;
	long N2 = N * N;

		//	Fix the patch width.
	long widthPatch = refinementPatchSize;

		//	Fix the patch length.
	long lengthPatch = refinementPatchSize;

		//	Fix the amplitude - so that the resulting input vector is normalized.	
	double ampPatch = 1.0 / ((double) refinementPatchSize + 1);

		//	Fix the duration.
	long durationPatchInIters = (long) ( refinementStimDuration * oneSecondNumIter );

		//	Fix the offset.
	long offsetPatchInIters = trialLengthInIters / 5;

	long timeMaskStart = offsetPatchInIters;
	long timeMaskEnd = offsetPatchInIters + durationPatchInIters;

	long iStim = 0;

		//	Clear the vector that counts how many times each input layer node is stimulated.
	tmpStimCount.assign ( numValsPerIter, 0 );

	for ( long iDigit = 0; iDigit < numDigits; ++iDigit ) {

		long iDigColStart = iDigit * digitWidth;
		long iDigColEnd = (iDigit + 1) * digitWidth - widthPatch;

		for ( long iCol = iDigColStart; iCol < iDigColEnd; ++iCol ) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			long iDigRowStart = 0;
			long iDigRowEnd = N - lengthPatch;

			for ( long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow ) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

						//	Construct the spatio-temporal image of the input.
						//		ith row of inputPatch is a linear representation of the NxN input layer at time i
						//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign ( numValsPerTrial, 0 );

				for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

					for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

						long iiWhich = GetLin ( iiRow, iiCol, N );
						tmpStimCount[iiWhich] += 1.0;

						vector<double>::iterator ix = tmpVec.begin() + timeMaskStart * numValsPerIter + iiWhich;						
						for ( long k = timeMaskStart; k < timeMaskEnd; ++k ) {

							*ix = ampPatch;
							ix += numValsPerIter;

						} // for ( long k = timeMaskStart; k < timeMaskEnd; ++k ) {

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
				tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  offsetPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
				inputPatchList.push_back ( tmpInputPatch );

				inputPatchList[iStim].inputPattern.assign ( numValsPerTrial, 0 );
				VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPattern, 0, numValsPerTrial );

				inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
				VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} // for ( iDigit in 1:numDigits ) {

	for ( iStim = 0; iStim < ((long) inputPatchList.size()); ++iStim ) {
		if ( iStim == 0 ) {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.Baseline", 1, inputPatchList[iStim].inputPattern, 0, numValsPerTrial  );
		} else {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.Baseline", 0, inputPatchList[iStim].inputPattern, 0, numValsPerTrial  );
		} // if ( iStim == 0 ) {
		SaveNumericVectorToBinaryFile ( "CheckInputStim.Baseline", 0, inputPatchList[iStim].stimCount, 0, numValsPerIter  );
	} // for ( iStim = 0; iStim < inputPatchList.size(); ++iStim ) {
	
	tmpInputPatch.inputPattern.clear();

}; // void NM::GenDig3RefineInputPatchList ( vector<InputPatch> & ... ) {

	//
	//	GenDig3RefineInputPatchListSparse
	//
void NM::GenDig3RefineInputPatchListSparse ( vector<InputPatch> & inputPatchList ) {

	InputPatch tmpInputPatch;
	vector<double> tmpVec;
	vector<double> tmpStimCount;

	tmpInputPatch.Initialize();

	long numDigits = 3;
	long digitWidth = N / numDigits;
	long N2 = N * N;

		//	Fix the patch width.
	long widthPatch = refinementPatchSize;

		//	Fix the patch length.
	long lengthPatch = refinementPatchSize;

		//	Fix the amplitude - so that the resulting input vector is normalized.	
	double ampPatch = 1.0 / ((double) refinementPatchSize + 1);

		//	Fix the duration.
	long durationPatchInIters = (long) ( refinementStimDuration * oneSecondNumIter );
	long timeMaskStart = refinementOffsetToPatchInIters;
	long timeMaskEnd = refinementOffsetToPatchInIters + durationPatchInIters;

	long numInputPatchElements = (refinementPatchSize + 1) * ( refinementPatchSize + 1);

	long iStim = 0;

			//	Clear the vector that counts how many times each input layer node is stimulated.
	tmpStimCount.assign ( numValsPerIter, 0 );

	for ( long iDigit = 0; iDigit < numDigits; ++iDigit ) {

		long iDigColStart = iDigit * digitWidth;
		long iDigColEnd = (iDigit + 1) * digitWidth - widthPatch;

		for ( long iCol = iDigColStart; iCol < iDigColEnd; ++iCol ) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			long iDigRowStart = 0;
			long iDigRowEnd = N - lengthPatch;

			for ( long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow ) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

						//	Construct the spatio-temporal image of the input.
						//		ith row of inputPatch is a linear representation of the NxN input layer at time i
						//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign ( numInputPatchElements, 0 );
				vector<double>::iterator ix = tmpVec.begin();

				for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

					for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

						long iiWhich = GetLin ( iiRow, iiCol, N );
						tmpStimCount[iiWhich] += 1.0;
						*ix++ = iiWhich;

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
				tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  refinementOffsetToPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
				inputPatchList.push_back ( tmpInputPatch );

				inputPatchList[iStim].inputPatternSparse.assign ( numInputPatchElements, 0 );
				VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements );

				inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
				VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} // for ( iDigit in 1:numDigits ) {

	for ( iStim = 0; iStim < ((long) inputPatchList.size()); ++iStim ) {
		if ( iStim == 0 ) {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.Baseline", 1, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements  );
		} else {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.Baseline", 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements  );
		} // if ( iStim == 0 ) {
		SaveNumericVectorToBinaryFile ( "CheckInputStim.Baseline", 0, inputPatchList[iStim].stimCount, 0, numValsPerIter  );
	} // for ( iStim = 0; iStim < inputPatchList.size(); ++iStim ) {
	
}; // void NM::GenDig3RefineInputPatchListSparse ( vector<InputPatch> & ... ) {

   //
   //	GenDig3RefineInputPatchListSparse2
   //
void NM::GenDig3RefineInputPatchListSparse2 (vector<InputPatch> & inputPatchList) {

		// Difference from previous is to not normalize the input.

	InputPatch tmpInputPatch;
	vector<double> tmpVec;
	vector<double> tmpStimCount;

	tmpInputPatch.Initialize();

	long numDigits = 3;
	long digitWidth = N / numDigits;
	long N2 = N * N;

	//	Fix the patch width.
	long widthPatch = refinementPatchSize;

	//	Fix the patch length.
	long lengthPatch = refinementPatchSize;

	//	Fix the amplitude - so that the resulting input vector is normalized.	
	// double ampPatch = 1.0 / ((double)refinementPatchSize + 1);
	double ampPatch = 1.0;


	//	Fix the duration.
	long durationPatchInIters = (long)(refinementStimDuration * oneSecondNumIter);

	//	Fix the offset.
	long offsetPatchInIters = trialLengthInIters / 5;

	long timeMaskStart = offsetPatchInIters;
	long timeMaskEnd = offsetPatchInIters + durationPatchInIters;

	long numInputPatchElements = (refinementPatchSize + 1) * (refinementPatchSize + 1);

	long iStim = 0;

	//	Clear the vector that counts how many times each input layer node is stimulated.
	tmpStimCount.assign(numValsPerIter, 0);

	for (long iDigit = 0; iDigit < numDigits; ++iDigit) {

		long iDigColStart = iDigit * digitWidth;
		long iDigColEnd = (iDigit + 1) * digitWidth - widthPatch;

		for (long iCol = iDigColStart; iCol < iDigColEnd; ++iCol) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			long iDigRowStart = 0;
			long iDigRowEnd = N - lengthPatch;

			for (long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

				//	Construct the spatio-temporal image of the input.
				//		ith row of inputPatch is a linear representation of the NxN input layer at time i
				//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign(numInputPatchElements, 0);
				vector<double>::iterator ix = tmpVec.begin();

				for (long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol) {

					for (long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow) {

						long iiWhich = GetLin(iiRow, iiCol, N);
						tmpStimCount[iiWhich] += 1.0;
						*ix++ = iiWhich;

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

				  //	Package up this input patch stimulus.
				tmpInputPatch.LoadParms(iDigit, GetLin(iRow, iCol, N), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters, offsetPatchInIters, numValsPerIter, numValsPerTrial, ampPatch);
				inputPatchList.push_back(tmpInputPatch);

				inputPatchList[iStim].inputPatternSparse.assign(numInputPatchElements, 0);
				VecCopy(tmpVec, 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements);

				inputPatchList[iStim].stimCount.assign(numValsPerIter, 0);
				VecCopy(tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter);

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {

	} // for ( iDigit in 1:numDigits ) {

	for (iStim = 0; iStim < ((long)inputPatchList.size()); ++iStim) {
		if (iStim == 0) {
			SaveNumericVectorToBinaryFile("CheckInputStim.Baseline", 1, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements);
		}
		else {
			SaveNumericVectorToBinaryFile("CheckInputStim.Baseline", 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements);
		} // if ( iStim == 0 ) {
		SaveNumericVectorToBinaryFile("CheckInputStim.Baseline", 0, inputPatchList[iStim].stimCount, 0, numValsPerIter);
	} // for ( iStim = 0; iStim < inputPatchList.size(); ++iStim ) {

}; // void NM::GenDig3RefineInputPatchListSparse ( vector<InputPatch> & ... ) {

	//
	//	GenDig3RefineInputPatchListSparse3
	//
void NM::GenDig3RefineInputPatchListSparse3 ( vector<InputPatch> & inputPatchList ) {

		//	Difference from GenDig3RefineInputPatchListSparse is that there is a normalized magnitude as an input parameter.

	InputPatch tmpInputPatch;
	vector<double> tmpVec;
	vector<double> tmpStimCount;

	tmpInputPatch.Initialize();

	long numDigits = 3;
	long digitWidth = N / numDigits;
	long N2 = N * N;

		//	Fix the patch width.
	long widthPatch = refinementPatchSize;

		//	Fix the patch length.
	long lengthPatch = refinementPatchSize;

		//	Fix the amplitude - so that the resulting input vector is normalized.	
	double ampPatch = refinementPatchNormalizedMagnitude / ((double) refinementPatchSize + 1);

		//	Fix the duration.
	long durationPatchInIters = (long) ( refinementStimDuration * oneSecondNumIter );
	long timeMaskStart = refinementOffsetToPatchInIters;
	long timeMaskEnd = refinementOffsetToPatchInIters + durationPatchInIters;

	long numInputPatchElements = (refinementPatchSize + 1) * ( refinementPatchSize + 1);

	long iStim = 0;

			//	Clear the vector that counts how many times each input layer node is stimulated.
	tmpStimCount.assign ( numValsPerIter, 0 );

	for ( long iDigit = 0; iDigit < numDigits; ++iDigit ) {

		long iDigColStart = iDigit * digitWidth;
		long iDigColEnd = (iDigit + 1) * digitWidth - widthPatch;

		for ( long iCol = iDigColStart; iCol < iDigColEnd; ++iCol ) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			long iDigRowStart = 0;
			long iDigRowEnd = N - lengthPatch;

			for ( long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow ) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

						//	Construct the spatio-temporal image of the input.
						//		ith row of inputPatch is a linear representation of the NxN input layer at time i
						//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign ( numInputPatchElements, 0 );
				vector<double>::iterator ix = tmpVec.begin();

				for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

					for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

						long iiWhich = GetLin ( iiRow, iiCol, N );
						tmpStimCount[iiWhich] += 1.0;
						*ix++ = iiWhich;

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
				tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  refinementOffsetToPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
				inputPatchList.push_back ( tmpInputPatch );

				inputPatchList[iStim].inputPatternSparse.assign ( numInputPatchElements, 0 );
				VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements );

				inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
				VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} // for ( iDigit in 1:numDigits ) {

	for ( iStim = 0; iStim < ((long) inputPatchList.size()); ++iStim ) {
		if ( iStim == 0 ) {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.Baseline", 1, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements  );
		} else {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.Baseline", 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements  );
		} // if ( iStim == 0 ) {
		SaveNumericVectorToBinaryFile ( "CheckInputStim.Baseline", 0, inputPatchList[iStim].stimCount, 0, numValsPerIter  );
	} // for ( iStim = 0; iStim < inputPatchList.size(); ++iStim ) {
	
}; // void NM::GenDig3RefineInputPatchListSparse ( vector<InputPatch> & ... ) {

	//
	//	GenDig3RefineInputPatchList
	//
void NM::GenDig3SelectiveStimInputPatchList ( vector<InputPatch> & inputPatchList ) {

	InputPatch tmpInputPatch;
	vector<double> tmpVec;
	vector<double> tmpStimCount;

	tmpInputPatch.Initialize();

	long numDigits = 3;
	long digitWidth = N / numDigits;
	long N2 = N * N;
	long iDigit = 0;

		//	Fix the patch width.
	long widthPatch = refinementPatchSize;

		//	Fix the patch length.
	long lengthPatch = refinementPatchSize;

		//	Fix the amplitude - so that the resulting input vector is normalized.	
	double ampPatch = 1.0 / ((double) refinementPatchSize + 1);

		//	Fix the duration.
	long durationPatchInIters = (long) ( refinementStimDuration * oneSecondNumIter );

		//	Fix the offset.
	long offsetPatchInIters = trialLengthInIters / 5;

	long timeMaskStart = offsetPatchInIters;
	long timeMaskEnd = offsetPatchInIters + durationPatchInIters;

	long iStim = 0;

		//	Clear the vector that counts how many times each input layer node is stimulated.
	tmpStimCount.assign ( numValsPerIter, 0 );

		//
		//	Part I: Typical baseline refinement stimulation pattern
		//

	for ( iDigit = 0; iDigit < numDigits; ++iDigit ) {

		long iDigColStart = iDigit * digitWidth;
		long iDigColEnd = (iDigit + 1) * digitWidth - widthPatch;

		for ( long iCol = iDigColStart; iCol < iDigColEnd; ++iCol ) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			long iDigRowStart = 0;
			long iDigRowEnd = N - lengthPatch;

			for ( long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow ) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

						//	Construct the spatio-temporal image of the input.
						//		ith row of inputPatch is a linear representation of the NxN input layer at time i
						//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign ( numValsPerTrial, 0 );

				for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

					for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

						long iiWhich = GetLin ( iiRow, iiCol, N );
						tmpStimCount[iiWhich] += 1.0;
						vector<double>::iterator ix = tmpVec.begin() + timeMaskStart * numValsPerIter + iiWhich;						
						for ( long k = timeMaskStart; k < timeMaskEnd; ++k ) {

							*ix = ampPatch;
							ix += numValsPerIter;

						} // for ( long k = timeMaskStart; k < timeMaskEnd; ++k ) {

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
				tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  offsetPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
				inputPatchList.push_back ( tmpInputPatch );
				inputPatchList[iStim].inputPattern.assign ( numValsPerTrial, 0 );
				VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPattern, 0, numValsPerTrial );
				inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
				VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} // for ( iDigit in 1:numDigits ) {

		//
		//	Part II Option A: Selective stimulation of a given FULL SECTOR.
		//

	if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 4 || selectiveStimZoneID == 7 ) {
		iDigit = 0;
	} else if ( selectiveStimZoneID == 2 || selectiveStimZoneID == 5 || selectiveStimZoneID == 8 ) {
		iDigit = 1;
	} else {
		iDigit = 2;
	} // if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 4 || selectiveStimZoneID == 7 )

	for ( long iFactor = 0; iFactor < selectiveStimFactor; ++iFactor ) {

		long iDigColStart = iDigit * digitWidth;
		long iDigColEnd = (iDigit + 1) * digitWidth - widthPatch;

		long iDigRowStart = 0;
		long iDigRowEnd = 0;

		for ( long iCol = iDigColStart; iCol < iDigColEnd; ++iCol ) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 2 || selectiveStimZoneID == 3 ) {
				iDigRowStart = 0;
				iDigRowEnd = digitWidth - lengthPatch;
			} else if ( selectiveStimZoneID == 4 || selectiveStimZoneID == 5 || selectiveStimZoneID == 6 ) {
				iDigRowStart = digitWidth;
				iDigRowEnd = 2 * digitWidth - lengthPatch;
			} else {
				iDigRowStart = 2 * digitWidth;
				iDigRowEnd = 3 * digitWidth - lengthPatch;
			} // if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 2 || selectiveStimZoneID == 3 ) {

			for ( long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow ) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

						//	Construct the spatio-temporal image of the input.
						//		ith row of inputPatch is a linear representation of the NxN input layer at time i
						//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign ( numValsPerTrial, 0 );

				for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

					for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

						long iiWhich = GetLin ( iiRow, iiCol, N );
						tmpStimCount[iiWhich] += 1.0;
						vector<double>::iterator ix = tmpVec.begin() + timeMaskStart * numValsPerIter + iiWhich;						
						for ( long k = timeMaskStart; k < timeMaskEnd; ++k ) {

							*ix = ampPatch;
							ix += numValsPerIter;

						} // for ( long k = timeMaskStart; k < timeMaskEnd; ++k ) {

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
				tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  offsetPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
				inputPatchList.push_back ( tmpInputPatch );
				inputPatchList[iStim].inputPattern.assign ( numValsPerTrial, 0 );
				VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPattern, 0, numValsPerTrial );
				inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
				VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} // for ( long iFactor = 0; iFactor < selectiveStimFactor; ++iFactor ) {

		//	Sanity checking stimulation patterns.  (See R display routines.)
	for ( iStim = 0; iStim < ((long) inputPatchList.size()); ++iStim ) {
		if ( iStim == 0 ) {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.SelStim", 1, inputPatchList[iStim].inputPattern, 0, numValsPerTrial  );
		} else {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.SelStim", 0, inputPatchList[iStim].inputPattern, 0, numValsPerTrial  );
		} // if ( iStim == 0 ) {
		SaveNumericVectorToBinaryFile ( "CheckInputStim.SelStim", 0, inputPatchList[iStim].stimCount, 0, numValsPerIter  );
	} // for ( iStim = 0; iStim < inputPatchList.size(); ++iStim ) {
	
	tmpInputPatch.inputPattern.clear();

}; // void NM::GenDig3SelectiveStimInputPatchList ( vector<InputPatch> & ... ) {

	

	//
	//	GenDig3RefineInputPatchListSparse
	//
void NM::GenDig3SelectiveStimInputPatchListSparse ( vector<InputPatch> & inputPatchList ) {

	InputPatch tmpInputPatch;
	vector<double> tmpVec;
	vector<double> tmpStimCount;

	tmpInputPatch.Initialize();

	long numDigits = 3;
	long digitWidth = N / numDigits;
	long N2 = N * N;
	long iDigit = 0;

		//	Fix the patch width.
	long widthPatch = refinementPatchSize;

		//	Fix the patch length.
	long lengthPatch = refinementPatchSize;

		//	Fix the amplitude - so that the resulting input vector is normalized.	
	double ampPatch = 1.0 / ((double) refinementPatchSize + 1);

		//	Fix the duration.
	long durationPatchInIters = (long) ( refinementStimDuration * oneSecondNumIter );

		//	Fix the offset.
	long offsetPatchInIters = trialLengthInIters / 5;

	long timeMaskStart = offsetPatchInIters;
	long timeMaskEnd = offsetPatchInIters + durationPatchInIters;

	long numInputPatchElements = (refinementPatchSize + 1) * ( refinementPatchSize + 1);

	long iStim = 0;

	double maxCount = 0;

		//	Clear the vector that counts how many times each input layer node is stimulated.
	tmpStimCount.assign ( numValsPerIter, 0 );

		//
		//	Part I: Typical baseline refinement stimulation pattern
		//

	for ( iDigit = 0; iDigit < numDigits; ++iDigit ) {

		long iDigColStart = iDigit * digitWidth;
		long iDigColEnd = (iDigit + 1) * digitWidth - widthPatch;

		for ( long iCol = iDigColStart; iCol < iDigColEnd; ++iCol ) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			long iDigRowStart = 0;
			long iDigRowEnd = N - lengthPatch;

			for ( long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow ) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

						//	Construct the spatio-temporal image of the input.
						//		ith row of inputPatch is a linear representation of the NxN input layer at time i
						//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign ( numInputPatchElements, 0 );
				vector<double>::iterator ix = tmpVec.begin();

				for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

					for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

						long iiWhich = GetLin ( iiRow, iiCol, N );
						tmpStimCount[iiWhich] += 1.0;
						*ix++ = iiWhich;

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
				tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  offsetPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
				inputPatchList.push_back ( tmpInputPatch );

				inputPatchList[iStim].inputPatternSparse.assign ( numInputPatchElements, 0 );
				VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements );

				inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
				VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );
				
				maxCount = VecMaxValue ( tmpStimCount, 0, numValsPerIter );

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} // for ( iDigit in 1:numDigits ) {

#if 0
		//
		//	Part II Option A: Selective stimulation of a given FULL SECTOR.
		//

	if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 4 || selectiveStimZoneID == 7 ) {
		iDigit = 0;
	} else if ( selectiveStimZoneID == 2 || selectiveStimZoneID == 5 || selectiveStimZoneID == 8 ) {
		iDigit = 1;
	} else {
		iDigit = 2;
	} // if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 4 || selectiveStimZoneID == 7 )

	for ( long iFactor = 0; iFactor < selectiveStimFactor; ++iFactor ) {

		long iDigColStart = iDigit * digitWidth;
		long iDigColEnd = (iDigit + 1) * digitWidth - widthPatch;

		long iDigRowStart = 0;
		long iDigRowEnd = 0;

		for ( long iCol = iDigColStart; iCol < iDigColEnd; ++iCol ) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 2 || selectiveStimZoneID == 3 ) {
				iDigRowStart = 0;
				iDigRowEnd = digitWidth - lengthPatch;
			} else if ( selectiveStimZoneID == 4 || selectiveStimZoneID == 5 || selectiveStimZoneID == 6 ) {
				iDigRowStart = digitWidth;
				iDigRowEnd = 2 * digitWidth - lengthPatch;
			} else {
				iDigRowStart = 2 * digitWidth;
				iDigRowEnd = 3 * digitWidth - lengthPatch;
			} // if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 2 || selectiveStimZoneID == 3 ) {

			for ( long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow ) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

						//	Construct the spatio-temporal image of the input.
						//		ith row of inputPatch is a linear representation of the NxN input layer at time i
						//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign ( numInputPatchElements, 0 );
				vector<double>::iterator ix = tmpVec.begin();

				for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

					for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

						long iiWhich = GetLin ( iiRow, iiCol, N );
						tmpStimCount[iiWhich] += 1.0;
						*ix++ = iiWhich;

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
				tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  offsetPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
				inputPatchList.push_back ( tmpInputPatch );

				inputPatchList[iStim].inputPatternSparse.assign ( numInputPatchElements, 0 );
				VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements );

				inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
				VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );
				
				maxCount = VecMaxValue ( tmpStimCount, 0, numValsPerIter );

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} // for ( long iFactor = 0; iFactor < selectiveStimFactor; ++iFactor ) {
#endif


		//
		//	Part II Option B: Selective stimulation of a single patch.
		//

	if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 4 || selectiveStimZoneID == 7 ) {
		iDigit = 0;
	} else if ( selectiveStimZoneID == 2 || selectiveStimZoneID == 5 || selectiveStimZoneID == 8 ) {
		iDigit = 1;
	} else {
		iDigit = 2;
	} // if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 4 || selectiveStimZoneID == 7 )

	for ( long iFactor = 0; iFactor < selectiveStimFactor; ++iFactor ) {

		long iSubSectorOffset = (long) floor ( (double) digitWidth / 2 );	// Selects sub-sector origin.
		long iCol = iDigit * digitWidth + iSubSectorOffset;
		long iRow = 0;

		long colMaskStart = iCol;
		long colMaskEnd = iCol + refinementPatchSize;

		if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 2 || selectiveStimZoneID == 3 ) {
			iRow = iSubSectorOffset;
		} else if ( selectiveStimZoneID == 4 || selectiveStimZoneID == 5 || selectiveStimZoneID == 6 ) {
			iRow = digitWidth + iSubSectorOffset;
		} else {
			iRow = 2 * digitWidth + iSubSectorOffset;
		} // if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 2 || selectiveStimZoneID == 3 ) {
			
		long rowMaskStart = iRow;
		long rowMaskEnd = iRow + lengthPatch;

					//	Construct the spatio-temporal image of the input.
					//		ith row of inputPatch is a linear representation of the NxN input layer at time i
					//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

		tmpVec.assign ( numInputPatchElements, 0 );
		vector<double>::iterator ix = tmpVec.begin();

		for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

			for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

				long iiWhich = GetLin ( iiRow, iiCol, N );
				tmpStimCount[iiWhich] += 1.0;
				*ix++ = iiWhich;

			} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

		} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
		tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  offsetPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
		inputPatchList.push_back ( tmpInputPatch );

		inputPatchList[iStim].inputPatternSparse.assign ( numInputPatchElements, 0 );
		VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements );

		inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
		VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );
				
		maxCount = VecMaxValue ( tmpStimCount, 0, numValsPerIter );

		++iStim;
	
	} // for ( long iFactor = 0; iFactor < selectiveStimFactor; ++iFactor ) {

#if 0
		//
		//	Part II Option C: Selective stimulation of a given SUB-SECTOR.
		//

	if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 4 || selectiveStimZoneID == 7 ) {
		iDigit = 0;
	} else if ( selectiveStimZoneID == 2 || selectiveStimZoneID == 5 || selectiveStimZoneID == 8 ) {
		iDigit = 1;
	} else {
		iDigit = 2;
	} // if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 4 || selectiveStimZoneID == 7 )

	for ( long iFactor = 0; iFactor < selectiveStimFactor; ++iFactor ) {

		long iSubSectorOffset = floor ( (double) digitWidth / 2 );	// Selects sub-sector origin.
		long iDigColStart = iDigit * digitWidth + iSubSectorOffset;
		long iDigColEnd = (iDigit + 1) * digitWidth - widthPatch;

		long iDigRowStart = 0;
		long iDigRowEnd = 0;

		for ( long iCol = iDigColStart; iCol < iDigColEnd; ++iCol ) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 2 || selectiveStimZoneID == 3 ) {
				iDigRowStart = iSubSectorOffset;
				iDigRowEnd = digitWidth - lengthPatch;
			} else if ( selectiveStimZoneID == 4 || selectiveStimZoneID == 5 || selectiveStimZoneID == 6 ) {
				iDigRowStart = digitWidth + iSubSectorOffset;
				iDigRowEnd = 2 * digitWidth - lengthPatch;
			} else {
				iDigRowStart = 2 * digitWidth + iSubSectorOffset;
				iDigRowEnd = 3 * digitWidth - lengthPatch;
			} // if ( selectiveStimZoneID == 1 || selectiveStimZoneID == 2 || selectiveStimZoneID == 3 ) {

			for ( long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow ) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

						//	Construct the spatio-temporal image of the input.
						//		ith row of inputPatch is a linear representation of the NxN input layer at time i
						//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign ( numInputPatchElements, 0 );
				vector<double>::iterator ix = tmpVec.begin();

				for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

					for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

						long iiWhich = GetLin ( iiRow, iiCol, N );
						tmpStimCount[iiWhich] += 1.0;
						*ix++ = iiWhich;

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
				tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  offsetPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
				inputPatchList.push_back ( tmpInputPatch );

				inputPatchList[iStim].inputPatternSparse.assign ( numInputPatchElements, 0 );
				VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements );

				inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
				VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );
				
				maxCount = VecMaxValue ( tmpStimCount, 0, numValsPerIter );

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} // for ( long iFactor = 0; iFactor < selectiveStimFactor; ++iFactor ) {

#endif

		//	Sanity checking stimulation patterns.  (See R display routines.)
	for ( iStim = 0; iStim < ((long) inputPatchList.size()); ++iStim ) {
		if ( iStim == 0 ) {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.SelStim", 1, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements  );
		} else {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.SelStim", 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements  );
		} // if ( iStim == 0 ) {
		SaveNumericVectorToBinaryFile ( "CheckInputStim.SelStim", 0, inputPatchList[iStim].stimCount, 0, numValsPerIter  );
	} // for ( iStim = 0; iStim < inputPatchList.size(); ++iStim ) {
	
	tmpInputPatch.inputPattern.clear();

}; // void NM::GenDig3SelectiveStimInputPatchListSparse ( vector<InputPatch> & ... ) {


	//
	//	GenDig3RefineInputPatchList.
	//
void NM::GenDig3SyndactylyControlInputPatchList ( vector<InputPatch> & inputPatchList ) {

	InputPatch tmpInputPatch;
	vector<double> tmpVec;
	vector<double> tmpStimCount;

	tmpInputPatch.Initialize();

	long numDigits = 3;
	long digitWidth = N / numDigits;
	long N2 = N * N;

		//	Fix the patch width.
	long widthPatch = refinementPatchSize;

		//	Fix the patch length.
	long lengthPatch = refinementPatchSize;

		//	Fix the amplitude - so that the resulting input vector is normalized.	
	double ampPatch = 1.0 / ((double) refinementPatchSize + 1);

		//	Fix the duration.
	long durationPatchInIters = (long) ( refinementStimDuration * oneSecondNumIter );

		//	Fix the offset.
	long offsetPatchInIters = trialLengthInIters / 5;

	long timeMaskStart = offsetPatchInIters;
	long timeMaskEnd = offsetPatchInIters + durationPatchInIters;

	long iStim = 0;

		//	Clear the vector that counts how many times each input layer node is stimulated.
	tmpStimCount.assign ( numValsPerIter, 0 );

	for ( long iDigit = 0; iDigit < numDigits; ++iDigit ) {

		long iDigColStart = iDigit * digitWidth;
		long iDigColEnd = (iDigit + 1) * digitWidth - widthPatch;

		for ( long iCol = iDigColStart; iCol < iDigColEnd; ++iCol ) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			long iDigRowStart = 0;
			long iDigRowEnd = N - lengthPatch;

			for ( long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow ) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

						//	Construct the spatio-temporal image of the input.
						//		ith row of inputPatch is a linear representation of the NxN input layer at time i
						//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign ( numValsPerTrial, 0 );

				for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

					for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

						long iiWhich = GetLin ( iiRow, iiCol, N );
						tmpStimCount[iiWhich] += 2.0;		//	Because we're du
						vector<double>::iterator ix = tmpVec.begin() + timeMaskStart * numValsPerIter + iiWhich;						
						for ( long k = timeMaskStart; k < timeMaskEnd; ++k ) {

							*ix = ampPatch;
							ix += numValsPerIter;

						} // for ( long k = timeMaskStart; k < timeMaskEnd; ++k ) {

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
				tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  offsetPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
				inputPatchList.push_back ( tmpInputPatch );
				inputPatchList[iStim].inputPattern.assign ( numValsPerTrial, 0 );
				VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPattern, 0, numValsPerTrial );
				inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
				VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} // for ( iDigit in 1:numDigits ) {

	for ( iStim = 0; iStim < ((long) inputPatchList.size()); ++iStim ) {
		if ( iStim == 0 ) {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.SyndactCtl", 1, inputPatchList[iStim].inputPattern, 0, numValsPerTrial  );
		} else {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.SyndactCtl", 0, inputPatchList[iStim].inputPattern, 0, numValsPerTrial  );
		} // if ( iStim == 0 ) {
		SaveNumericVectorToBinaryFile ( "CheckInputStim.SyndactCtl", 0, inputPatchList[iStim].stimCount, 0, numValsPerIter  );
	} // for ( iStim = 0; iStim < inputPatchList.size(); ++iStim ) {
	
	tmpInputPatch.inputPattern.clear();

}; // void NM::GenDig3SyndactylyControlInputPatchList ( vector<InputPatch> & ... ) {

	//
	//	Generate Inputs for a Syndactyly Experiment
	//
	//		Just assume (for now) that digits 1 and 2 are fused.
	//		The way to think about it is that there are really two digits. A  super digit and regular digit.
	//
void NM::GenDig3SyndactylyInputPatchList ( vector<InputPatch> & inputPatchList ) {

	InputPatch tmpInputPatch;
	vector<double> tmpVec;
	vector<double> tmpStimCount;

	tmpInputPatch.Initialize();

	long numDigits = 3;
	long digitWidth = N / numDigits;
	long numEffectiveDigits = 2;

		//	Fix the patch width.
	long widthPatch = refinementPatchSize;

		//	Fix the patch length.
	long lengthPatch = refinementPatchSize;

		//	Fix the amplitude - so that the resulting input vector is normalized.	
	double ampPatch = 1.0 / ((double) refinementPatchSize + 1);

		//	Fix the duration.
	long durationPatchInIters = (long) ( refinementStimDuration * oneSecondNumIter );

		//	Fix the offset.
	long offsetPatchInIters = trialLengthInIters / 5;

	long timeMaskStart = offsetPatchInIters;
	long timeMaskEnd = offsetPatchInIters + durationPatchInIters;

	long iStim = 0;

		//	Clear the vector that counts how many times each input layer node is stimulated.
	tmpStimCount.assign ( numValsPerIter, 0 );

	//for ( long iDigit = 0; iDigit <= numEffectiveDigits; ++iDigit ) {
	for ( long iDigit = 0; iDigit < numEffectiveDigits; ++iDigit ) {

		long iDigColStart = 0;
		long iDigColEnd = 0;

		if ( iDigit == 0 ) {
			iDigColStart = 0;
			iDigColEnd = 2 * digitWidth - widthPatch;
		} else {
			iDigColStart = 2 * digitWidth;
			iDigColEnd = N - widthPatch;
		} // if ( iDigit == 0 )

		for ( long iCol = iDigColStart; iCol < iDigColEnd; ++iCol ) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			long iDigRowStart = 0;
			long iDigRowEnd = N - lengthPatch;

			for ( long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow ) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

						//	Construct the spatio-temporal image of the input.
						//		ith row of inputPatch is a linear representation of the NxN input layer at time i
						//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign ( numValsPerTrial, 0 );

				for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

					for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

						long iiWhich = GetLin ( iiRow, iiCol, N );
						++tmpStimCount[iiWhich];
						vector<double>::iterator ix = tmpVec.begin() + timeMaskStart * numValsPerIter + iiWhich;						
						for ( long k = timeMaskStart; k < timeMaskEnd; ++k ) {

							*ix = ampPatch;
							ix += numValsPerIter;

						} // for ( long k = timeMaskStart; k < timeMaskEnd; ++k ) {

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
				tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  offsetPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
				inputPatchList.push_back ( tmpInputPatch );
				inputPatchList[iStim].inputPattern.assign ( numValsPerTrial, 0 );
				VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPattern, 0, numValsPerTrial );
				inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
				VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} // 	for ( long iDigit = 0; iDigit < numEffectiveDigits; ++iDigit ) {

	for ( iStim = 0; iStim < ((long) inputPatchList.size()); ++iStim ) {
		if ( iStim == 0 ) {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.Syndact", 1, inputPatchList[iStim].inputPattern, 0, numValsPerTrial  );
		} else {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.Syndact", 0, inputPatchList[iStim].inputPattern, 0, numValsPerTrial  );
		} // if ( iStim == 0 ) {
		SaveNumericVectorToBinaryFile ( "CheckInputStim.Syndact", 0, inputPatchList[iStim].stimCount, 0, numValsPerIter  );
	} // for ( iStim = 0; iStim < inputPatchList.size(); ++iStim ) {

	tmpInputPatch.inputPattern.clear();

}; // void NM:: GenDig3SyndactylyInputPatchList ( vector<InputPatch> & inputPatchList ) {

	//
	//	Generate Inputs for a Syndactyly Experiment (Sparse Representation)
	//
	//		Just assume (for now) that digits 1 and 2 are fused.
	//		The way to think about it is that there are really two digits. A  super digit and regular digit.
	//
void NM::GenDig3SyndactylyInputPatchListSparse ( vector<InputPatch> & inputPatchList ) {

	InputPatch tmpInputPatch;
	vector<double> tmpVec;
	vector<double> tmpStimCount;

	tmpInputPatch.Initialize();

	long numDigits = 3;
	long digitWidth = N / numDigits;
	long numEffectiveDigits = 2;

		//	Fix the patch width.
	long widthPatch = refinementPatchSize;

		//	Fix the patch length.
	long lengthPatch = refinementPatchSize;

		//	Fix the amplitude - so that the resulting input vector is normalized.	
	double ampPatch = syndactPatchNormalizedMagnitude / ((double) refinementPatchSize + 1);

		//	Fix the duration.
	long durationPatchInIters = (long) ( refinementStimDuration * oneSecondNumIter );

		//	Fix the offset.
	long timeMaskStart = syndactOffsetToPatchInIters;
	long timeMaskEnd = syndactOffsetToPatchInIters + durationPatchInIters;
	long numInputPatchElements = (refinementPatchSize + 1) * ( refinementPatchSize + 1);
	long iStim = 0;

		//	Clear the vector that counts how many times each input layer node is stimulated.
	tmpStimCount.assign ( numValsPerIter, 0 );

	//for ( long iDigit = 0; iDigit <= numEffectiveDigits; ++iDigit ) {
	for ( long iDigit = 0; iDigit < numEffectiveDigits; ++iDigit ) {

		long iDigColStart = 0;
		long iDigColEnd = 0;

		if ( iDigit == 0 ) {
			iDigColStart = 0;
			iDigColEnd = 2 * digitWidth - widthPatch;
		} else {
			iDigColStart = 2 * digitWidth;
			iDigColEnd = N - widthPatch;
		} // if ( iDigit == 0 )

		for ( long iCol = iDigColStart; iCol < iDigColEnd; ++iCol ) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			long iDigRowStart = 0;
			long iDigRowEnd = N - lengthPatch;

			for ( long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow ) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

						//	Construct the spatio-temporal image of the input.
						//		ith row of inputPatch is a linear representation of the NxN input layer at time i
						//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign ( numInputPatchElements, 0 );
				vector<double>::iterator ix = tmpVec.begin();

				for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

					for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

						long iiWhich = GetLin ( iiRow, iiCol, N );
						++tmpStimCount[iiWhich];
						*ix++ = iiWhich;

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
				tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  syndactOffsetToPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
				inputPatchList.push_back ( tmpInputPatch );

				inputPatchList[iStim].inputPatternSparse.assign ( numInputPatchElements, 0 );
				VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements );

				inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
				VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} // 	for ( long iDigit = 0; iDigit < numEffectiveDigits; ++iDigit ) {

	for ( iStim = 0; iStim < ((long) inputPatchList.size()); ++iStim ) {
		if ( iStim == 0 ) {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.Syndact", 1, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements  );
		} else {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.Syndact", 0, inputPatchList[iStim].inputPatternSparse, 0, numInputPatchElements  );
		} // if ( iStim == 0 ) {
		SaveNumericVectorToBinaryFile ( "CheckInputStim.Syndact", 0, inputPatchList[iStim].stimCount, 0, numValsPerIter  );
	} // for ( iStim = 0; iStim < inputPatchList.size(); ++iStim ) {

	tmpInputPatch.inputPattern.clear();

}; // void NM:: GenDig3SyndactylyInputPatchListSparse ( vector<InputPatch> & inputPatchList ) {


	//
	//	GenUniformRefineInputPatchList
	//
void NM::GenUniformRefineInputPatchList ( vector<InputPatch> & inputPatchList ) {

	InputPatch tmpInputPatch;
	vector<double> tmpVec;
	vector<double> tmpStimCount;

	tmpInputPatch.Initialize();

	long numDigits = 1;
	long digitWidth = N / numDigits;
	long N2 = N * N;

		//	Fix the patch width.
	long widthPatch = refinementPatchSize;

		//	Fix the patch length.
	long lengthPatch = refinementPatchSize;

		//	Fix the amplitude - so that the resulting input vector is normalized.	
	double ampPatch = 1.0 / ((double) refinementPatchSize + 1);

		//	Fix the duration.
	long durationPatchInIters = (long) ( refinementStimDuration * oneSecondNumIter );

		//	Fix the offset.
	long offsetPatchInIters = trialLengthInIters / 5;

	long timeMaskStart = offsetPatchInIters;
	long timeMaskEnd = offsetPatchInIters + durationPatchInIters;

	long iStim = 0;

		//	Clear the vector that counts how many times each input layer node is stimulated.
	tmpStimCount.assign ( numValsPerIter, 0 );

	for ( long iDigit = 0; iDigit < numDigits; ++iDigit ) {

		long iDigColStart = iDigit * digitWidth;
		long iDigColEnd = (iDigit + 1) * digitWidth - widthPatch;

		for ( long iCol = iDigColStart; iCol < iDigColEnd; ++iCol ) {

			long colMaskStart = iCol;
			long colMaskEnd = iCol + refinementPatchSize;

			long iDigRowStart = 0;
			long iDigRowEnd = N - lengthPatch;

			for ( long iRow = iDigRowStart; iRow < iDigRowEnd; ++iRow ) {

				long rowMaskStart = iRow;
				long rowMaskEnd = iRow + lengthPatch;

						//	Construct the spatio-temporal image of the input.
						//		ith row of inputPatch is a linear representation of the NxN input layer at time i
						//		jth col of inputPatch is the time evolution of cell j on the NxN input layer

				tmpVec.assign ( numValsPerTrial, 0 );

				for ( long iiCol = colMaskStart; iiCol <= colMaskEnd; ++iiCol ) {

					for ( long iiRow = rowMaskStart; iiRow <= rowMaskEnd; ++iiRow ) {

						long iiWhich = GetLin ( iiRow, iiCol, N );
						tmpStimCount[iiWhich] += 1.0;
						vector<double>::iterator ix = tmpVec.begin() + timeMaskStart * numValsPerIter + iiWhich;						
						for ( long k = timeMaskStart; k < timeMaskEnd; ++k ) {

							*ix = ampPatch;
							ix += numValsPerIter;

						} // for ( long k = timeMaskStart; k < timeMaskEnd; ++k ) {

					} // for ( iiRow in (rowMask[1]:rowMask[2]) ) {

				} // for ( iiCol in (colMask[1]:colMask[2]) ) {

					//	Package up this input patch stimulus.
				tmpInputPatch.LoadParms ( iDigit, GetLin ( iRow, iCol, N ), iCol, iRow, widthPatch, lengthPatch, durationPatchInIters,  offsetPatchInIters, numValsPerIter, numValsPerTrial, ampPatch );
				inputPatchList.push_back ( tmpInputPatch );
				inputPatchList[iStim].inputPattern.assign ( numValsPerTrial, 0 );
				VecCopy ( tmpVec, 0, inputPatchList[iStim].inputPattern, 0, numValsPerTrial );
				inputPatchList[iStim].stimCount.assign ( numValsPerIter, 0 );
				VecCopy ( tmpStimCount, 0, inputPatchList[iStim].stimCount, 0, numValsPerIter );

				++iStim;

			} // for ( iRow in minLegalRow:maxLegalRow ) {

		} // for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} // for ( iDigit in 1:numDigits ) {

	for ( iStim = 0; iStim < ((long) inputPatchList.size()); ++iStim ) {
		if ( iStim == 0 ) {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.Baseline", 1, inputPatchList[iStim].inputPattern, 0, numValsPerTrial  );
		} else {
			SaveNumericVectorToBinaryFile ( "CheckInputStim.Baseline", 0, inputPatchList[iStim].inputPattern, 0, numValsPerTrial  );
		} // if ( iStim == 0 ) {
		SaveNumericVectorToBinaryFile ( "CheckInputStim.Baseline", 0, inputPatchList[iStim].stimCount, 0, numValsPerIter  );
	} // for ( iStim = 0; iStim < inputPatchList.size(); ++iStim ) {
	
	tmpInputPatch.inputPattern.clear();

}; // void NM::GenUniformRefineInputPatchList ( vector<InputPatch> & ... ) {


	//
	//	GetCol
	//
long NM::GetCol ( const long & k, const long & N ) {
	return (long) floor( ( double) k / N );
}; // long NM::GetCol ( const long & k, const long & N )

	//
	//	GetLegal
	//
long NM::GetLegal ( const double & x, const double & minX, const double & maxX ) {
	long tmp = 0;
	if ( (x >= minX) && (x <= maxX) ) {
		tmp = 1;
	}
	return (long) tmp;
}; // long NM::GetLegal ( const double & x, const double & minX, const double & maxX ) {

	//
	//	GetLin
	//
long NM::GetLin  ( const long & rowVal, const long & colVal, const long & N ) {
	return (long)  ( colVal * N + rowVal );
}; // long NM::GetLin  ( const long & rowVal, const long & colVal, const long & N ) {

	//
	//	GetRow
	//
long NM::GetRow ( const long & k, const long & N ) {
	return (long) ( k - N * GetCol ( k, N )	);
}; // long NM::GetRow ( const long & k, const long & N )

	//
	//	Initialize Network and Workspace
	//
void NM::InitializeNetworkAndWorkspace () {

	sran0 ( 1 );		// Set the random number seed.
	plasticityFlag = 0;	// No weight plasticity until instructed.

		///////////////////////////////////////////////////////////////////
		//
		//	Initial Network & Workspaces.
		//	Set up the initial random network.
		//	For each variable allocate memory corresponding to a single trial as circular
		//	buffering is used to simulate time evolution of arbitrary length.
		//	Inject noise (and normalize as appropriate) just for the first iteration.
		//
		///////////////////////////////////////////////////////////////////

	v0.assign ( numValsPerTrial, 0 );
	v0CBuff.assign ( numValsPerIter, 0 );
	r0.assign ( numValsPerTrial, 0 );
	r0CBuff.assign ( numValsPerIter, 0 );
	AddNoise ( v0, 0, numValsPerIter, -noiseLevel, noiseLevel  );
	Sigmoid ( v0, r0, beta, 0, numValsPerIter );

	v1E.assign ( numValsPerTrial, 0 );
	v1ECBuff.assign ( numValsPerIter, 0 );
	r1E.assign ( numValsPerTrial, 0 );
	r1ECBuff.assign ( numValsPerIter, 0 );
	r1ERFMapRaw.assign ( numValsPerRFMapExp, 0 );
	r1ERFMapExpRes.assign ( numValsPerRFMapRes, 0 );
	AddNoise ( v1E, 0, numValsPerIter, -noiseLevel, noiseLevel  );
	Sigmoid ( v1E, r1E, beta, 0, numValsPerIter );

	v1I.assign ( numValsPerTrial, 0 );
	v1ICBuff.assign ( numValsPerIter, 0 );
	r1I.assign ( numValsPerTrial, 0 );
	r1ICBuff.assign ( numValsPerIter, 0 );
	r1IRFMapRaw.assign ( numValsPerRFMapExp, 0 );
	r1IRFMapExpRes.assign ( numValsPerRFMapRes, 0 );
	AddNoise ( v1I, 0, numValsPerIter, -noiseLevel, noiseLevel  );
	Sigmoid ( v1I, r1I, beta, 0, numValsPerIter );

	w1E0.assign ( numWeightsPerTrial, 0 );
	w1E0CBuff.assign ( numWeightsPerIter, 0 );
	InitPreSynWeights ( w1E0, g0, wghtMinValue, wghtMaxValue, 0, numWeightsPerIter  );
	NormalizePreSynWeights ( w1E0, eWResource, 0, numWeightsPerIter );

	w1IE.assign ( numWeightsPerTrial, 0 );
	w1IECBuff.assign ( numWeightsPerIter, 0 );
	InitPreSynWeights ( w1IE, g0, wghtMinValue, wghtMaxValue, 0, numWeightsPerIter  );
	NormalizePreSynWeights ( w1IE, iWResource, 0, numWeightsPerIter );

	w1EE.assign ( numWeightsPerTrial, 0 );
	w1EECBuff.assign ( numWeightsPerIter, 0 );
	InitPreSynWeights ( w1EE, g0, wghtMinValue, wghtMaxValue, 0, numWeightsPerIter  );
	NormalizePreSynWeights ( w1EE, eWResource, 0, numWeightsPerIter );

	w1EI.assign ( numWeightsPerTrial, 0 );
	w1EICBuff.assign ( numWeightsPerIter, 0 );
	InitPreSynWeights ( w1EI, g0, wghtMinValue, wghtMaxValue, 0, numWeightsPerIter  );
	NormalizePreSynWeights ( w1EI, eWResource, 0, numWeightsPerIter );

		//	For test purposes, read in known good weight matrices.

	//LoadNumericVectorFromBinaryFile ( "w1e0.dat1", w1E0, numWeightsPerIter );
	//LoadNumericVectorFromBinaryFile ( "w1ee.dat1", w1EE, numWeightsPerIter );
	//LoadNumericVectorFromBinaryFile ( "w1ei.dat1", w1EI, numWeightsPerIter );
	//LoadNumericVectorFromBinaryFile ( "w1ie.dat1", w1IE, numWeightsPerIter );

}; // void NM::InitializeNetworkAndWorkspace () {

	//
	//	InitializeWeights
	//
void NM::InitPreSynWeights ( vector<double> & x, const long & gVal, const double & xMin, const double & xMax, const long & offset, const long & iLen  ) {

	long cellID = 0;
	long iCellID = 0;
	long jCellID = 0;
	long N2 = (long) sqrt ( (double) iLen );
	long N = (long ) sqrt ( (double) N2 );
	long gridLength = 2 * gVal + 1;

	for ( ; cellID < N2; ++cellID ) {

		iCellID = GetRow ( cellID, N ) - gVal;
		jCellID = GetCol ( cellID, N ) - gVal;

		for ( long j = 0; j < gridLength; ++j ) {

			if ( GetLegal ( (jCellID + j), 0, (N-1) ) ) {

				for ( long  i = 0; i < gridLength; ++i ) {

					if ( GetLegal ( iCellID + i, 0, (N-1) )  ) {

						long itmp = GetLin ( (iCellID + i), (jCellID + j), N );
						x.at(offset + itmp*N2 + cellID ) = r1unif ( xMin, xMax );

					} // if ( GetLegal ( iCellID, 0, (N-1) ) & Legal ( jCellID, 0, (N-1) ) ) {

				} // for ( long  i = 0; i < gridLength; ++i ) {

			} // if ( GetLegal ( (jCellID + j), 0, (N-1) ) ) {

		} // for ( long j = 0; j < gridLength; ++j ) {

	} // for ( ; cellID < N2; ++cellID
		
}; // void NM::InitPreSynWeights ( vector<double> & w, const long & wGrid, const double & wMin, const double & wMax, const long & offset, const long & iLen  ) {

		//
		//	Iterate the system one time step forward.
		//
void NM::IterateSystem4 ( const long & iter ) {

	long	offsetToCurrentIter = iter * numValsPerIter;
	long	offsetToPreviousIter = ( iter - 1 ) * numValsPerIter;

	long	offsetToCurrentIterWeights = 0;
	long	offsetToPreviousIterWeights = 0;

	if ( plasticityFlag ) {
		offsetToCurrentIterWeights = iter * numWeightsPerIter;
		offsetToPreviousIterWeights = ( iter - 1 ) * numWeightsPerIter;
	} // if ( plasticityFlag )

		//	Update v0 term.
	VecScalarMult ( v0, offsetToCurrentIter, v0, offsetToPreviousIter, v0Alpha, numValsPerIter );
	AddNoise (  v0, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	
		//	Update v1E term.
	VecScalarMult ( v1E, offsetToCurrentIter, v1E, offsetToPreviousIter, v1EAlpha, numValsPerIter );
	AddNoise (  v1E, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	AddVecMatMult ( v1E, offsetToCurrentIter, r0, offsetToPreviousIter, w1E0, offsetToPreviousIterWeights, numValsPerIter );
	//AddVecMatMult ( v1E, offsetToCurrentIter, r1E, offsetToPreviousIter, w1EE, offsetToPreviousIterWeights, numValsPerIter );
	SubtractVecMatMult ( v1E, offsetToCurrentIter, r1I, offsetToPreviousIter, w1EI, offsetToPreviousIterWeights, numValsPerIter );

		//	Update v1I term.
	VecScalarMult ( v1I, offsetToCurrentIter, v1I, offsetToPreviousIter, v1IAlpha, numValsPerIter );
	AddNoise (  v1I, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	AddVecMatMult ( v1I, offsetToCurrentIter, r1E, offsetToPreviousIter, w1IE, offsetToPreviousIterWeights, numValsPerIter );

		//	Run Sigmoids.
	Sigmoid ( v0, r0, beta, offsetToCurrentIter, numValsPerIter );
	Sigmoid ( v1E, r1E, beta, offsetToCurrentIter, numValsPerIter );
	Sigmoid ( v1I, r1I, beta, offsetToCurrentIter, numValsPerIter );

		//	Adapt or simply propogate forward weights.
	if ( plasticityFlag ) {

			//	Adapt w1E0
		VecVecMultWeights ( r0, offsetToPreviousIter,  r1E, offsetToPreviousIter, w1E0, offsetToPreviousIterWeights, offsetToCurrentIterWeights, wBeta, numValsPerIter );
		AddVecScalarMult ( w1E0, offsetToCurrentIterWeights, w1E0, offsetToPreviousIterWeights, wTauAlpha, numWeightsPerIter );

			//	Adapt w1EI
		VecVecMultWeights ( r1I, offsetToPreviousIter,  r1E, offsetToPreviousIter, w1EI, offsetToPreviousIterWeights, offsetToCurrentIterWeights, wBeta, numValsPerIter );
		AddVecScalarMult ( w1EI, offsetToCurrentIterWeights, w1EI, offsetToPreviousIterWeights, wTauAlpha, numWeightsPerIter );
		
			//	Adapt w1IE
		VecVecMultWeights ( r1E, offsetToPreviousIter,  r1I, offsetToPreviousIter, w1IE, offsetToPreviousIterWeights, offsetToCurrentIterWeights, wBeta, numValsPerIter );
		AddVecScalarMult ( w1IE, offsetToCurrentIterWeights, w1IE, offsetToPreviousIterWeights, wTauAlpha, numWeightsPerIter );
		
			//	Adapt w1EE
		//VecVecMultWeights ( r1E, offsetToPreviousIter,  r1E, offsetToPreviousIter, w1EE, offsetToPreviousIterWeights, offsetToCurrentIterWeights, wBeta, numValsPerIter );
		//AddVecScalarMult ( w1EE, offsetToCurrentIterWeights, w1EE, offsetToPreviousIterWeights, wTauAlpha, numWeightsPerIter );

			//	Normalize the updated weights.
		NormalizePreSynWeights ( w1E0, eWResource, offsetToCurrentIterWeights, numWeightsPerIter );
		NormalizePreSynWeights ( w1IE, iWResource, offsetToCurrentIterWeights, numWeightsPerIter );
		//NormalizePreSynWeights ( w1EE, eWResource, offsetToCurrentIterWeights, numWeightsPerIter );
		NormalizePreSynWeights ( w1EI, eWResource, offsetToCurrentIterWeights, numWeightsPerIter );

	} // if ( plasticityFlag )

};  // void NM::IterateSystem4 ( const long & iter )

		//
		//	Iterate the system one time step forward.
		//
void NM::IterateSystem5 ( const long & iter ) {

	long	offsetToCurrentIter = iter * numValsPerIter;
	long	offsetToPreviousIter = ( iter - 1 ) * numValsPerIter;

	long	offsetToCurrentIterWeights = 0;
	long	offsetToPreviousIterWeights = 0;

	if ( plasticityFlag ) {
		offsetToCurrentIterWeights = iter * numWeightsPerIter;
		offsetToPreviousIterWeights = ( iter - 1 ) * numWeightsPerIter;
	} // if ( plasticityFlag )

		//	Update v0 term.
	VecScalarMult ( v0, offsetToCurrentIter, v0, offsetToPreviousIter, v0Alpha, numValsPerIter );
	AddNoise (  v0, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	
		//	Update v1E term.
	VecScalarMult ( v1E, offsetToCurrentIter, v1E, offsetToPreviousIter, v1EAlpha, numValsPerIter );
	AddNoise (  v1E, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	AddVecMatMult ( v1E, offsetToCurrentIter, r0, offsetToPreviousIter, w1E0, offsetToPreviousIterWeights, numValsPerIter );
	AddVecMatMult ( v1E, offsetToCurrentIter, r1E, offsetToPreviousIter, w1EE, offsetToPreviousIterWeights, numValsPerIter );
	SubtractVecMatMult ( v1E, offsetToCurrentIter, r1I, offsetToPreviousIter, w1EI, offsetToPreviousIterWeights, numValsPerIter );

		//	Update v1I term.
	VecScalarMult ( v1I, offsetToCurrentIter, v1I, offsetToPreviousIter, v1IAlpha, numValsPerIter );
	AddNoise (  v1I, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	AddVecMatMult ( v1I, offsetToCurrentIter, r1E, offsetToPreviousIter, w1IE, offsetToPreviousIterWeights, numValsPerIter );

		//	Run Sigmoids.
	Sigmoid ( v0, r0, beta, offsetToCurrentIter, numValsPerIter );
	Sigmoid ( v1E, r1E, beta, offsetToCurrentIter, numValsPerIter );
	Sigmoid ( v1I, r1I, beta, offsetToCurrentIter, numValsPerIter );

		//	Adapt or simply propogate forward weights.
	if ( plasticityFlag ) {

			//	Adapt w1E0
		VecVecMultWeights ( r0, offsetToPreviousIter,  r1E, offsetToPreviousIter, w1E0, offsetToPreviousIterWeights, offsetToCurrentIterWeights, wBeta, numValsPerIter );
		AddVecScalarMult ( w1E0, offsetToCurrentIterWeights, w1E0, offsetToPreviousIterWeights, wTauAlpha, numWeightsPerIter );

			//	Adapt w1EI
		VecVecMultWeights ( r1I, offsetToPreviousIter,  r1E, offsetToPreviousIter, w1EI, offsetToPreviousIterWeights, offsetToCurrentIterWeights, wBeta, numValsPerIter );
		AddVecScalarMult ( w1EI, offsetToCurrentIterWeights, w1EI, offsetToPreviousIterWeights, wTauAlpha, numWeightsPerIter );
		
			//	Adapt w1IE
		VecVecMultWeights ( r1E, offsetToPreviousIter,  r1I, offsetToPreviousIter, w1IE, offsetToPreviousIterWeights, offsetToCurrentIterWeights, wBeta, numValsPerIter );
		AddVecScalarMult ( w1IE, offsetToCurrentIterWeights, w1IE, offsetToPreviousIterWeights, wTauAlpha, numWeightsPerIter );
		
			//	Adapt w1EE
		VecVecMultWeights ( r1E, offsetToPreviousIter,  r1E, offsetToPreviousIter, w1EE, offsetToPreviousIterWeights, offsetToCurrentIterWeights, wBeta, numValsPerIter );
		AddVecScalarMult ( w1EE, offsetToCurrentIterWeights, w1EE, offsetToPreviousIterWeights, wTauAlpha, numWeightsPerIter );

			//	Normalize the updated weights.
		NormalizePreSynWeights ( w1E0, eWResource, offsetToCurrentIterWeights, numWeightsPerIter );
		NormalizePreSynWeights ( w1IE, iWResource, offsetToCurrentIterWeights, numWeightsPerIter );
		NormalizePreSynWeights ( w1EE, eWResource, offsetToCurrentIterWeights, numWeightsPerIter );
		NormalizePreSynWeights ( w1EI, eWResource, offsetToCurrentIterWeights, numWeightsPerIter );

	} // if ( plasticityFlag )

};  // void NM::IterateSystem5 ( const long & iter )

		//
		//	Read the simulation parameters.
		//
void NM::LoadExpParms ( const string & expParmFileName ) {
	
	ifstream inStream ( expParmFileName.c_str() );

	inStream >>		networkConfigType;

	inStream >>		experimentType;
	if ( experimentType > 0 ) {
			inStream >>	baselineFileName;
	}; // if ( experimentType > 0 )

	inStream >>		dumpRefineStepsRawToFile;
	inStream >>		dumpRefineStepsRawToFileStepSize;

	inStream >>		dumpRFMapRawToFile;
	inStream >>		doRFMapStepSize;

	inStream >>		N;		
	inStream >>		g0;
	inStream >>		g1EE;
	inStream >>		g1EI;
	inStream >>		g1IE;

	inStream >>		noiseLevel;
	inStream >>		beta;
	inStream >>		tau;
	inStream >>		tauDivisor;

	inStream >>		wghtNormalizeFlag;
	inStream >>		eWResource;
	inStream >>		iWResource;
	inStream >>		eeWResource;
	inStream >>		e0WResource;
	inStream >>		eiWResource;
	inStream >>		ieWResource;
	inStream >>		wghtMinValue;
	inStream >>		wghtMaxValue;
	inStream >>		wBeta;
	inStream >>		wBetaDecay;
	inStream >>		wBetaE0;
	inStream >>		wBetaEE;
	inStream >>		wBetaEI;
	inStream >>		wBetaIE;

	inStream >>		numSpinUpTrials;

	inStream >>		magRFProbe;
	inStream >>		durRFProbe;
	inStream >>		trialDurRFProbe;
	inStream >>		offsetToRFProbe;
	inStream >>		kRFPeakToEdgeDetect;

	inStream >>		numBaselineCycles;
	inStream >>		refinementPatchSize;
	inStream >>		refinementPatchNormalizedMagnitude;
	inStream >>		refinementOffsetToPatchStart;
	inStream >>		refinementStimDuration;

	inStream >>		numSyndactCycles;
	inStream >>		syndactPatchSize;
	inStream >>		syndactPatchNormalizedMagnitude;
	inStream >>		syndactOffsetToPatchStart;
	inStream >>		syndactStimDuration;
	inStream >>		syndactStimZoneID;

	inStream >>		numSelStimCycles;
	inStream >>		selectiveStimPatchSize;
	inStream >>		selectiveStimDuration;
	inStream >>		selectiveStimFactor;
	inStream >>		selectiveStimZoneID;
	
	inStream >>		numAmputCycles;
	inStream >>		refinementAmputPatchSize;
	inStream >>		refinementAmputStimDuration;
	inStream >>		kAmputDigit;
	inStream >>		kAmputZoneMax;

	inStream >>		numCLesionCycles;
	inStream >>		refinementCLesionPatchSize;
	inStream >>		refinementCLesionStimDuration;
	inStream >>		kCLesionDigit;
	inStream >>		kCLesionZoneMaxList;

	N2 = N * N;
	deltaT = tau/tauDivisor;
	deltaTD2 = deltaT / 2.0;							// Precomputed for use in RK4.
	deltaTD6 = deltaT / 6.0;							// Precomputed for use in RK4.
	minusOneOverTau = -1.0 / tau;						// Precomputed for use in RK4.

	noiseLevelDotDeltaT = noiseLevel * deltaT;			// Computation saver.
	oneSecondNumIter = (long) (1.0 / deltaT);
	v0Tau = tau;										// Input layer activity fade out is relatively fast.
	v0Alpha = 1.0-(deltaT/v0Tau);						// Save some computations and make differential equations more readable.
	vETau = tau;										// Time constant for e cells membrane potential relative to baseline.
	v1EAlpha = 1.0-(deltaT/vETau); 						// Save some computations and make differential equations more readable.
	vITau = tau;										// Time constant for i cells membrane potential relative to baseline.
	v1IAlpha = 1.0-(deltaT/vITau); 						// Save some computations and make differential equations more readable.
	wTau = tau * 100.0;									// Time constant for synaptic weight decay relative to baseline.
	wTauAlpha = 1.0-(deltaT/wTau);						// Save some computation and make differential equations more readable.
	minusOneOverWTau = -1.0 / wTau;						// For RK4.

	wBetaCopy = wBeta;									// Normalilze wBeta to deltaT.
	wBetaE0Copy = wBetaE0;
	wBetaEECopy = wBetaEE;
	wBetaEICopy = wBetaEI;
	wBetaIECopy = wBetaIE;

	numIterPerTrial = (long) (trialDurRFProbe * oneSecondNumIter);			// Default number of iterations max (when needed).  Some experiments cap iteration max count differently.
	numIterPerTrialMinusOne = numIterPerTrial - 1;
	trialDurRFProbeInIters = trialLengthInIters = numIterPerTrial;			// The basic experimental epoch duration in this simulation in units of iterations.
	offsetToRFProbeInIters = (long) (offsetToRFProbe  * oneSecondNumIter);	// Stimulus ONSET time relative to trial start in units of iterations.
	refinementOffsetToPatchInIters = (long) (refinementOffsetToPatchStart  * oneSecondNumIter);
	syndactOffsetToPatchInIters = (long) (syndactOffsetToPatchStart  * oneSecondNumIter);
	durRFProbeInIters = (long) (durRFProbe * oneSecondNumIter);				// Stimulus ON time in units of iterations.
	timeOffRFProbeInIters = offsetToRFProbeInIters + durRFProbeInIters;		// Stimulus OFFSET time relative to trial start in units of iterations.

	numValsPerIter = N2;
	numValsPerTrial = numIterPerTrial * numValsPerIter;
	numValsPerRFMapExp = N2 * numValsPerTrial;
	numValsPerRFMapRes = N2 * N2;
	
	numWeightsPerIter = N2 * N2;
	numWeightsPerTrial = N2 * N2 * numIterPerTrial;

	plasticityFlag = 0;

		//	Print to screen all the key parameters.

	cout << "N2: " << N2 << endl;
	cout << "deltaT: " << deltaT << endl;
	cout << "oneSecondNumIter: " << oneSecondNumIter << endl;
	cout << "v0Tau: " << v0Tau << endl;
	cout << "v0Alpha: " << v0Alpha << endl;
	cout << "v1EAlpha: " << v1EAlpha << endl;
	cout << "v1IAlpha: " << v1IAlpha << endl;
	cout << "wTauAlpha: " << wTauAlpha << endl;
	cout << "numIterPerTrial: " << numIterPerTrial << endl;
	cout << "trialDurRFProbeInIters: " << trialDurRFProbeInIters << endl;
	cout << "offsetToRFProbeInIters: " << offsetToRFProbeInIters << endl;
	cout << "durRFProbeInIters: " << durRFProbeInIters << endl;
	cout << "timeOffRFProbeInIters: " << timeOffRFProbeInIters << endl;
	cout << "refinementOffsetToPatchInIters: " << refinementOffsetToPatchInIters << endl;
	cout << "syndactOffsetToPatchInIters: " << syndactOffsetToPatchInIters << endl;

	cout << "eeWResource: " << eeWResource << endl;
	cout << "e0WResource: " << e0WResource << endl;
	cout << "eiWResource: " << eiWResource << endl;
	cout << "ieWResource: " << ieWResource << endl;

	cout << "magRFProbe: " << magRFProbe << endl;
	cout << "durRFProbe: " << durRFProbe << endl;
	cout << "trialDurRFProbe: " << trialDurRFProbe << endl;
	cout << "offsetToRFProbe: " << offsetToRFProbe << endl;
	cout << "kRFPeakToEdgeDetect: " << kRFPeakToEdgeDetect << endl;

	cout << "numBaselineCycles: " << numBaselineCycles << endl;
	cout << "refinementPatchSize: " << refinementPatchSize << endl;
	cout << "refinementPatchNormalizedMagnitude: " << refinementPatchNormalizedMagnitude << endl;
	cout << "refinementOffsetToPatchStart: " << refinementOffsetToPatchStart << endl;
	cout << "refinementStimDuration: " << refinementStimDuration << endl;

	cout << "numSyndactCycles: " << numSyndactCycles << endl;
	cout << "syndactPatchSize: " << syndactPatchSize << endl;
	cout << "syndactPatchNormalizedMagnitude: " << syndactPatchNormalizedMagnitude << endl;
	cout << "syndactOffsetToPatchStart: " << syndactOffsetToPatchStart << endl;
	cout << "syndactStimDuration: " << syndactStimDuration << endl;
	cout << "syndactStimZoneID: " << syndactStimZoneID << endl;
	
}; // void loadDataBase ( const string & )

	//
	//	Load the Network from a binary file.
	//
void NM::LoadNetworkFromFile ( const string & fBaseName, const long & index1, const long & index2 ) {

	fstream readFile;

	string indexString1 = static_cast<ostringstream*>( &(ostringstream() << index1) )->str();
	string indexString2 = static_cast<ostringstream*>( &(ostringstream() << index2) )->str();
	string fName = fBaseName + "." + indexString1 + "." + indexString2 + ".bin";

	readFile.open( fName.c_str(), ios::in | ios::binary );

	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, v0, numValsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, r0, numValsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, v1E, numValsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, r1E, numValsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, v1I, numValsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, r1I, numValsPerTrial  );

	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, w1E0, numWeightsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, w1EE, numWeightsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, w1EI, numWeightsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, w1IE, numWeightsPerTrial  );

	readFile.close();
	
}; // void NM::SaveNetworkToFile ( const string &, const long &, const long & ) {


	//
	//	Normalize the incoming (presynaptic) weights.
	//
void NM::NormalizePreSynWeights ( vector<double> & x, const double & xResource, const long & offset, const long & iLen ) {
	
	vector<double>::iterator ix;
	vector<double>::iterator iy;
	vector<double>::iterator ixEnd;

	vector<double> vTmp;
	vector<double> vNonZeroCount;

	long iCol = 0;

	double maxCount = 0;
	long N2 = (long) sqrt ( (double) iLen );

		//	Inputs for cell j are stored in jth column of x.  (Recall x itself is a giant vector linear representation.
		//	Count the number of non-zero inputs for each receiving cell.	
	vTmp.assign ( iLen, 0 );
	CompareNE ( x, vTmp, 0, offset, iLen );
	vNonZeroCount.assign ( N2, 0 );
	SumByCol ( vTmp, vNonZeroCount, 0, N2, N2 );

		//	For those counts that are zero, replace with a 1.0.  Only effect is to prevent subsequent divide-by-zero.
	for ( ix = vNonZeroCount.begin(), ixEnd = vNonZeroCount.end(); ix != ixEnd; ++ix ) {
		if ( *ix == 0 ) { *ix = 1.0; }
	} // 	for ( ix = vNonZeroCount.begin(), ixEnd = vNonZeroCount.end(); ix != ixEnd; ++ix ) {

		//	Get the maximum count of incoming weights across all the cells.
		//	Generate the normalization factor.
	maxCount = VecMaxValue ( vNonZeroCount, 0, N2 );
	vTmp.assign ( N2, 0 );
	SumByCol ( x, vTmp, offset, N2, N2 );
	for ( ix = vTmp.begin(), ixEnd = vTmp.end(), iy = vNonZeroCount.begin(); ix != ixEnd; ++ix, ++iy ) {
		*ix = *ix / ( ( xResource * *iy ) / maxCount );
	} // 	for ( ix = vNonZeroCount.begin(), ixEnd = vNonZeroCount.end(); ix != ixEnd; ++ix ) {
	//tmp2 = ( apply ( w, 2, sum ) / ( (wResource * tmp) / max ( tmp ) ) );

		//	Now apply the normalization factor to each weight.	
	for ( iy = vTmp.begin(), iCol = 0; iCol < N2; ++iy, ++iCol ) {
		for ( ix = ( x.begin() + offset + iCol * N2 ), ixEnd = ( ix + N2 ); ix != ixEnd; ++ix ) {
			*ix /= *iy;
		} // for ( ixEnd = ix + N2; ix != ixEnd; ++ix ) {
	} // for ( ix = x.begin(), iy = vTmp.begin(), iCol = 0; iCol < N2; ++iy ) {

		//	Sanity Check
	// vTmp.assign ( N2, 0 ); SumByCol ( x, vTmp, offset, N2, N2 );

		//	Clean up.
	vTmp.clear();
	vNonZeroCount.clear();

};  // void NM::NormalizePresynapticWeights ( vector<double> &, const long &, const double & ) {

	//
	//	Take the last Iter and park it.  See also UnParkLastIter.
	//
void NM::ParkLastIter () {

	long offsetToLastIter = ( numIterPerTrial - 1 ) * numValsPerIter;
	VecCopy ( v0, offsetToLastIter, v0CBuff, 0, numValsPerIter );
	VecCopy ( r0, offsetToLastIter, r0CBuff, 0, numValsPerIter );
	VecCopy ( v1E, offsetToLastIter, v1ECBuff, 0, numValsPerIter );
	VecCopy ( r1E, offsetToLastIter, r1ECBuff, 0, numValsPerIter );
	VecCopy ( v1I, offsetToLastIter, v1ICBuff, 0, numValsPerIter );
	VecCopy ( r1I, offsetToLastIter, r1ICBuff, 0, numValsPerIter );

}; // void NM::ParkLastIter ( ) {

	//
	//	Random permutation.  Used to shuffle order of input stimulation.
	//
void NM::RandomPermutation ( vector<long> & rOrder, const long & min, const long & max ) {

	double	xmin = min;
	double	xmax = max;

	long	numItems = max - min;

	vector<long> usedFlag;
	vector<long> referenceCandidates;
	vector<long> remainingCandidates;

	usedFlag.assign ( numItems, 0 );
	referenceCandidates.clear();
	remainingCandidates.clear();
	
		//  Build the initial list of candidate labels.
	for ( long i = 0; i < numItems; ++ i) {
		remainingCandidates.push_back ( i );
		referenceCandidates.push_back ( min + i );
	} // for ( vector<long>::iterator ix = candidateList.begin(), long i = min; i < max; ++i ) 

	while ( remainingCandidates.size() > 1 ) {

			//	Choose one at random.
		long iCandidate = (long) floor ( r1unif ( 0, (double) remainingCandidates.size() ) );
		rOrder.push_back ( referenceCandidates[ remainingCandidates[iCandidate] ] );
		usedFlag[ remainingCandidates[iCandidate] ] = 1;

			//	Build the list of remaining candidate labels.
		remainingCandidates.clear();
		for ( long i = 0; i < numItems; ++ i ) {
			if ( !usedFlag[i] ) {
				remainingCandidates.push_back ( i );
			} // if ( !usedFlag[i] ) {
		} // for ( long i = 0; i < numItems; ++i )

	} // while ( remainingCandidates.size() > 1 )

		//	Pick up the last on.
	rOrder.push_back ( referenceCandidates[ remainingCandidates[0] ] );

		//	Clean up.
	referenceCandidates.clear();
	remainingCandidates.clear();
	usedFlag.clear();
			
}; // void RandomPermutation ( vector<long> &, const long &, const long & ) {

	//
	//	Receptive Field Map.
	//
void NM::ReceptiveFieldMap () {

	r1ERFMapRaw.assign ( numValsPerRFMapExp, 0 );
	r1IRFMapRaw.assign ( numValsPerRFMapExp, 0 );

	ParkLastIter ();
	UnParkLastIterToFirstIter ();

		//
		//	Do the receptive field mapping experiment.
		//	Iterate over each S layer node with an RF Probe and store the network response.
		//
	for ( long iCell = 0; iCell < N2; ++iCell ) {

		long offsetToRFMapRaw = 0;
		vector<double>::iterator iv0 = v0.begin() + iCell;
		vector<double>::iterator ir0 = r0.begin() + iCell;
				
		for ( long iter = 1; iter < numIterPerTrial; ++iter ) {

				//	Stimulus ON/OFF.
			if ( (iter > offsetToRFProbeInIters ) && (iter <= timeOffRFProbeInIters ) ) {
				long itmp = ( iter - 1 ) * numValsPerIter;
				*( iv0 + itmp ) = magRFProbe + r1unif ( -noiseLevel, noiseLevel );
				*( ir0 + itmp ) = Sigmoid1 ( *( iv0 + itmp ), beta );
			} // 			if ( (iter > offsetToRFProbeInIters ) && (iter <= timeOffRFProbeInIters ) ) {				

			if ( networkConfigType == 5 ) {
				IterateSystem5 ( iter );
			} else {
				IterateSystem4 ( iter );
			} // if ( networkConfigType == 5 ) {

		} // for ( long iter = 1; iter < numIterPerTrial; ++iter ) {

			//	Store the raw r1E and r1I data for later analysis.
		offsetToRFMapRaw = iCell * numValsPerTrial;
		VecCopy ( r1E, 0, r1ERFMapRaw, offsetToRFMapRaw, numValsPerTrial );
		VecCopy ( r1I, 0, r1IRFMapRaw, offsetToRFMapRaw, numValsPerTrial );

			//	Restore the initial conditions to be the same as upon entry.
		UnParkLastIterToFirstIter ();

	}  // for ( long iCell = 0; iCell < N2; ++iCell )

	ReceptiveFieldProbeAnalysis ();	//	Process the raw RF probe data to generate Experimental Results (an N^2 x N^2 matrix).

}; // void NM::ReceptiveFieldMap

	//
	//	ReceptiveFieldProbeAnalysis
	//
void NM::ReceptiveFieldProbeAnalysis () {

	vector<double>::iterator ixR1EResults = r1ERFMapExpRes.begin();
	vector<double>::iterator ixR1IResults = r1IRFMapExpRes.begin();

	double meanCorrection = (double) offsetToRFProbeInIters / durRFProbeInIters;

		//	iCell indexes the cells on Layer S.
	for ( long iCell = 0; iCell < N2; ++iCell ) {

			//	jCell indexes the cells on Layer C.
		for ( long jCell = 0; jCell < N2; ++jCell ) {

				//	Determine the average activity in the pre-stimulus time period.

			vector<double>::iterator ixR1E = r1ERFMapRaw.begin() + iCell * numValsPerTrial + jCell;		
			vector<double>::iterator ixR1I = r1IRFMapRaw.begin() + iCell * numValsPerTrial + jCell;
			double r1EPreTmp = 0;
			double r1IPreTmp = 0;

			for ( long k = 0; k < offsetToRFProbeInIters; ++k ) {

				r1EPreTmp += *ixR1E;
				r1IPreTmp += *ixR1I;

				ixR1E += numValsPerIter;		//	Same cell, but next iteration.
				ixR1I += numValsPerIter;		//	Same cell, but next iteration.

			} // for ( long k = 0; k < offsetToRFProbeInIters; ++k ) {		

				//	Determine the average activity during the stimulus time period.

			ixR1E = r1ERFMapRaw.begin() + iCell * numValsPerTrial + (offsetToRFProbeInIters + 1) * numValsPerIter + jCell;		
			ixR1I = r1IRFMapRaw.begin() + iCell * numValsPerTrial + (offsetToRFProbeInIters + 1) * numValsPerIter + jCell;
			double r1EStimTmp = 0;
			double r1IStimTmp = 0;
			for ( long k = 0; k < durRFProbeInIters; ++k ) {

				r1EStimTmp += *ixR1E;
				r1IStimTmp += *ixR1I;

				ixR1E += numValsPerIter;		//	Same cell, but next iteration.
				ixR1I += numValsPerIter;		//	Same cell, but next iteration.

			} // for ( long k = 0; k < durRFProbeInIters; ++k ) {

				//	Store the ratio average responses.

			*ixR1EResults = ( r1EStimTmp / r1EPreTmp ) * meanCorrection;	
			*ixR1IResults = ( r1IStimTmp / r1IPreTmp ) * meanCorrection;

			++ixR1EResults;
			++ixR1IResults;

		} // 		for ( long jCell = 0; jCell < N2; ++jCell ) {

	}  // for ( long iCell = 0; iCell < N2; ++iCell )
	
}; // void NM::ReceptiveFieldProbeAnalysis () {

	//
	//	Save the Network to a binary file.
	//
void NM::SaveNetworkToFile ( const string & fBaseName, const long & index1, const long & index2 ) {

	string indexString1 = static_cast<ostringstream*>( &(ostringstream() << index1) )->str();
	string indexString2 = static_cast<ostringstream*>( &(ostringstream() << index2) )->str();
	string fName = fBaseName + "." + indexString1 + "." + indexString2 + ".bin";

	SaveNumericVectorToBinaryFile ( fName, 1, v0, 0, numValsPerTrial  );
	SaveNumericVectorToBinaryFile ( fName, 0, r0, 0, numValsPerTrial  );
	SaveNumericVectorToBinaryFile ( fName, 0, v1E, 0, numValsPerTrial  );
	SaveNumericVectorToBinaryFile ( fName, 0, r1E, 0, numValsPerTrial  );
	SaveNumericVectorToBinaryFile ( fName, 0, v1I, 0, numValsPerTrial  );
	SaveNumericVectorToBinaryFile ( fName, 0, r1I, 0, numValsPerTrial  );

	SaveNumericVectorToBinaryFile ( fName, 0, w1E0, 0, numWeightsPerTrial  );
	SaveNumericVectorToBinaryFile ( fName, 0, w1EE, 0, numWeightsPerTrial  );
	SaveNumericVectorToBinaryFile ( fName, 0, w1EI, 0, numWeightsPerTrial  );
	SaveNumericVectorToBinaryFile ( fName, 0, w1IE, 0, numWeightsPerTrial  );
	
}; // void NM::SaveNetworkToFile ( const string &, const long &, const long & ) {

	//
	//	Save the RF Exp Data to a binary file.
	//
void NM::SaveRFExpResToFile ( const string & fBaseName, const long & index1, const long & index2 ) {

	string indexString1 = static_cast<ostringstream*>( &(ostringstream() << index1) )->str();
	string indexString2 = static_cast<ostringstream*>( &(ostringstream() << index2) )->str();
	string fName = fBaseName + "." + indexString1 + "." + indexString2 + ".bin";

	SaveNumericVectorToBinaryFile ( fName, 1, r1ERFMapExpRes, 0, numValsPerRFMapRes  );
	SaveNumericVectorToBinaryFile ( fName, 0, r1IRFMapExpRes, 0, numValsPerRFMapRes  );
	
}; // void NM::SaveNetworkToFile ( const string &, const long &, const long & ) {

	//
	//	Save the RF Exp Data to a binary file - specialized for Focal Stim experiments.
	//
void NM::SaveRFExpResToFileFocalStim ( const string & fBaseName, const long & index1, const long & index2, const long & index3, const long & index4 ) {

	string indexString1 = static_cast<ostringstream*>( &(ostringstream() << index1) )->str();
	string indexString2 = static_cast<ostringstream*>( &(ostringstream() << index2) )->str();
	string indexString3 = static_cast<ostringstream*>( &(ostringstream() << index3) )->str();
	string indexString4 = static_cast<ostringstream*>( &(ostringstream() << index4) )->str();
	string fName = fBaseName + "." + indexString1 + "." + indexString2 + "." + indexString3 + "." + indexString4 +".bin";

	SaveNumericVectorToBinaryFile ( fName, 1, r1ERFMapExpRes, 0, numValsPerRFMapRes  );
	SaveNumericVectorToBinaryFile ( fName, 0, r1IRFMapExpRes, 0, numValsPerRFMapRes  );
	
}; // void NM::SaveRFExpResToFile ( const string & fBaseName, const long & index1, const long & index2, const long & index3, const long & index4 ) {

	//
	//	SetNoise.
	//
void NM::SetNoise (  vector<double> & x, const long & offset, const long & iLen, const double & xMin, const double & xMax ) {
	vector<double>::iterator ix = x.begin() + offset;
	vector<double>::iterator iEnd = ix + iLen;
	for ( ; ix != iEnd; ++ix ) {
		*ix = r1unif ( xMin, xMax );
	}
}; // void NM::SetNoise (  vector<double> & x, const long & offset, const long & iLen, const double & xMin, const double & xMax ) {

   //
   //	SetVecMatMult: y_bar = v_bar * M.  (Everything in linear representation.)
   //
void NM::SetVecMatMult(vector<double> & y, const long & yOffset, vector<double> & v, const long & vOffset, vector<double> & M, const long & mOffset, const long & iLen) {
	vector<double>::iterator iy = y.begin() + yOffset;
	vector<double>::iterator iyEnd = iy + iLen;
	vector<double>::iterator im = M.begin() + mOffset;

	for (; iy != iyEnd; ++iy) {

		vector<double>::iterator iv = v.begin() + vOffset;
		vector<double>::iterator ivEnd = iv + iLen;
		double tmp = 0;

		for (; iv != ivEnd; ++iv, ++im) {

			tmp += (*iv * *im);

		} // for ( ; iv != ivEnd; ++iv )

		*iy = tmp;

	} // for ( ; iy != iyEnd; ++iy )

};  // void SetVecMatMult ( vector<double> &, const long &, vector<double> &, const long &, vector<double> &, const long &, const long & ) {

	//
	//	Sigmoid function.
	//
void NM::Sigmoid ( vector<double> & x, vector<double> & y, const double & beta, const long & offset, const long & iLen ) {
	vector<double>::iterator ix = x.begin() + offset;
	vector<double>::iterator iy = y.begin() + offset;
	vector<double>::iterator iEnd = ix + iLen;

	for ( ; ix != iEnd; ++ix, ++iy ) {
		*iy = 0.5 * (1.0 + tanh(beta*(*ix - 0.5)));
		//*iy = 0.5 * (1.0 + tanh(beta*(*ix)));
	}

}; // void NM::sigmoid ( vector<double> & x, const double beta, const long iLen ) {

	//
	//	Sigmoid1 function.
	//
double NM::Sigmoid1 ( const double & x, const double & beta ) {

	double ans = 0.5 * (1.0 + tanh(beta*(x - 0.5)));

	return ans;

}; // void NM::sigmoid ( vector<double> & x, const double beta, const long iLen ) {

	//
	//	SubtractVecMatMult: y_bar -= v_bar * M.  (Everything in linear representation.)
	//
void NM::SubtractVecMatMult ( vector<double> & y, const long & yOffset, vector<double> & v, const long & vOffset, vector<double> & M, const long & mOffset, const long & iLen ) {
	vector<double>::iterator iy = y.begin() + yOffset;
	vector<double>::iterator iyEnd = iy + iLen;
	vector<double>::iterator im = M.begin() + mOffset;

	for ( ; iy != iyEnd; ++iy ) {
		
		vector<double>::iterator iv = v.begin() + vOffset;
		vector<double>::iterator ivEnd = iv + iLen;
		double tmp = 0;

		for ( ; iv != ivEnd; ++iv, ++im ) {

				tmp += *iv * *im;

		} // for ( ; iv != ivEnd; ++iv )

		*iy -= tmp;

	} // for ( ; iy != iyEnd; ++iy )

};  // void AddVecMatMult ( vector<double> &, const long &, vector<double> &, const long &, vector<double> &, const long &, const long & ) {


	//
	//	SumByCol.  The usage here can be confusing.  2D matrices are represented as linear vectors.
	//	We may have cause to sum by column for a whole collection of such 2D matrices.
	//
void NM::SumByCol ( vector<double> & x, vector<double> & y, const long & offset, const long & numCols, const long & lenCol ) {

	vector<double>::iterator ix = x.begin() + offset;
	vector<double>::iterator iy = y.begin();

	for ( long iCol = 0; iCol < numCols; ++iCol, ++iy ) {
		for ( vector<double>::iterator ixEnd = ix + lenCol; ix != ixEnd; ++ix ) {
			*iy += *ix;
		} // for ( vector<double>::iterator ixEnd = ix + lenCol; ix != ixEnd; ++ix ) {
	} // for ( long i = 0; i < numCols; ++i )

}; // void NM::SumByCol ( vector<double> & x, vector<double> & y, const long & offset, const long & iLen )

	//
	//	SumByCol.  The usage here can be confusing.  2D matrices are represented as linear vectors.
	//	We may have cause to sum by column for a whole collection of such 2D matrices.
	//
void NM::SumByRow ( vector<double> & x, vector<double> & y, const long & offset, const long & numRows, const long & numCols, const long & lenCol ) {
	
	vector<double>::iterator iy = y.begin();
	long iRow = 0;

	for ( ; iRow < numRows; ++iRow, ++iy ) {

		vector<double>::iterator ix = (x.begin() + offset + iRow);
		vector<double>::iterator ixEnd = ix + lenCol * ( numCols - 1 );

		for ( ; ix != ixEnd; ix += lenCol ) {
			*iy += *ix;
		} // for ( vector<double>::iterator ix = x.begin() + offset + iRow, ...

	} // for ( long i = 0; i < numCols; ++i )

}; // void NM::SumByRow ( vector<double> & x, vector<double> & y, const long & offset, const long & numRows, const long & numCols, const long & lenCol ) {

	//
	//	Syndactyly
	//		Shares many flow and features of Baseline Refinement
	//
void NM::Syndactyly () {

	long iStim = 0;
	long numStim = 0;
	vector<InputPatch> inputPatches;
	vector<long> orderInputPatches;

		//	Generate the input patches.
	GenDig3SyndactylyInputPatchList ( inputPatches );

		//	Turn ON weight adaptation.
	plasticityFlag = 1;

		//	Determine a random ordering of the inputs.
	orderInputPatches.clear();
	RandomPermutation ( orderInputPatches, 0, (long) inputPatches.size() );
	numStim = (long) inputPatches.size();

		//	Apply inputs in random order.

	for ( long iStim = 0; iStim < numStim; ++iStim ) {

		for ( long iter = 1; iter < numIterPerTrial; ++iter ) {

			long iOffset = ( iter - 1 ) * numValsPerIter;

			VecCopy ( inputPatches[ orderInputPatches[iStim] ].inputPattern, iOffset, v0, iOffset, numValsPerIter );
			AddNoise (  v0, iOffset, numValsPerIter, -noiseLevel, noiseLevel );
			Sigmoid ( v0, r0, beta, iOffset, numValsPerIter );

			if ( networkConfigType == 5 ) {
				IterateSystem5 ( iter );
			} else {
				IterateSystem4 ( iter );
			} // if ( networkConfigType == 5 ) {

		} // for ( long iter = 1; iter < numIterPerTrial; ++iter ) {
		CopyLastIterValsToZeroIter();
		CopyLastIterWeightsToZeroIter();
		
	} // for ( iStim = 0; iStim < numStim; ++iStim

		//	Turn OFF weight adaptation.
	plasticityFlag = 0;

		//	Clean up.
	CleanUpDig3SyndactylyInputPatchList ( inputPatches );
	orderInputPatches.clear();

}; // void NM::Syndactyly

	//
	//	SyndactylyControl
	//		Shares many flow and features of Baseline Refinement
	//
void NM::SyndactylyControl () {

	long iStim = 0;
	long numStim = 0;
	vector<InputPatch> inputPatches;
	vector<long> orderInputPatches;

		//	Generate the input patches.
	GenDig3SyndactylyInputPatchList ( inputPatches );

		//	Turn ON weight adaptation.
	plasticityFlag = 1;

		//	Determine a random ordering of the inputs.
	orderInputPatches.clear();
	RandomPermutation ( orderInputPatches, 0, (long) inputPatches.size() );
	numStim = (long) inputPatches.size();

		//	Apply inputs in random order.

	for ( long iStim = 0; iStim < numStim; ++iStim ) {

		for ( long iter = 1; iter < numIterPerTrial; ++iter ) {

			long iOffset = ( iter - 1 ) * numValsPerIter;

			VecCopy ( inputPatches[ orderInputPatches[iStim] ].inputPattern, iOffset, v0, iOffset, numValsPerIter );
			AddNoise (  v0, iOffset, numValsPerIter, -noiseLevel, noiseLevel );
			Sigmoid ( v0, r0, beta, iOffset, numValsPerIter );

			if ( networkConfigType == 5 ) {
				IterateSystem5 ( iter );
			} else {
				IterateSystem4 ( iter );
			} // if ( networkConfigType == 5 ) {

		} // for ( long iter = 1; iter < numIterPerTrial; ++iter ) {
		CopyLastIterValsToZeroIter();
		CopyLastIterWeightsToZeroIter();
		
	} // for ( iStim = 0; iStim < numStim; ++iStim

		//	Turn OFF weight adaptation.
	plasticityFlag = 0;

		//	Clean up.
	CleanUpDig3SyndactylyInputPatchList ( inputPatches );
	orderInputPatches.clear();

}; // void NM::SyndactylyControl



	//
	//	Spin up the system with a few trial following initial conditions.  No plasticity during this time.
	//
void NM::SystemSpinUp () {

	for ( long iTrials = 0; iTrials < numSpinUpTrials; ++iTrials ) {
			
		for ( long iter = 1; iter < numIterPerTrial; ++iter ) {
			
			if ( networkConfigType == 5 ) {
				IterateSystem5 ( iter );
			} else {
				IterateSystem4 ( iter );
			} // if ( networkConfigType == 5 ) {

		} // for ( long iter = 1; iter < numIterPerTrial; ++iter ) {
		CopyLastIterValsToZeroIter();
		if ( plasticityFlag ) {
			CopyLastIterWeightsToZeroIter();
		} // if ( plasticityFlag )
		
	} // for ( long iTrials = 0; iTrials < nm.numSpinUpTrials; ++iTrials ) {

}; // void NM::SystemSpinUp () {

	//
	//	Take the last Iter out of parking and put it in first iter.  See also ParkLastIter.
	//
void NM::UnParkLastIterToFirstIter ( ) {

	VecCopy ( v0CBuff, 0, v0, 0, numValsPerIter );
	VecCopy ( r0CBuff, 0, r0, 0, numValsPerIter );
	VecCopy ( v1ECBuff, 0, v1E, 0, numValsPerIter );
	VecCopy ( r1ECBuff, 0, r1E, 0, numValsPerIter );
	VecCopy ( v1ICBuff, 0, v1I, 0, numValsPerIter );
	VecCopy ( r1ICBuff, 0, r1I, 0, numValsPerIter );

}; // void NM::UnParkLastIter ( ) {

   //
   //	Add vectors: z = x + y;
   //
void NM::VecAdd(vector<double> & z, const long & zOffset, vector<double> & x, const long & xOffset, vector<double> & y, const long & yOffset, const long & iLen) {

	vector<double>::iterator ix = x.begin() + xOffset;
	vector<double>::iterator iy = y.begin() + yOffset;
	vector<double>::iterator iz = z.begin() + zOffset;
	vector<double>::iterator iEnd = ix + iLen;

	for (; ix != iEnd; ++ix, ++iy, ++iz) {
		*iz = *ix + *iy;
	}

}; // void NM::VecCopy ( vector<double> & x, const long & xOffset, vector<double> & y, const long & yOffset, const long & iLen ) {


	//
	//	Copy vector x to vector y.
	//
void NM::VecCopy ( vector<double> & x, const long & xOffset, vector<double> & y, const long & yOffset, const long & iLen ) {
	
	vector<double>::iterator ix = x.begin() + xOffset;
	vector<double>::iterator iy = y.begin() + yOffset;
	vector<double>::iterator iEnd = ix + iLen;

	for ( ; ix != iEnd; ++ix, ++iy ) {
		*iy = *ix;
	} 

}; // void NM::VecCopy ( vector<double> & x, const long & xOffset, vector<double> & y, const long & yOffset, const long & iLen ) {

	//
	//	Vector maximum value.
	//
double NM::VecMaxValue ( vector<double> & x, const long & offset, const long & iLen ) {
	
	vector<double>::iterator ix = x.begin() + offset;
	vector<double>::iterator iEnd = ix + iLen;
	double ans;

	ans = *ix++;
	for ( ; ix != iEnd; ++ix ) {
		if ( *ix > ans ) { ans = *ix; }
	}

	return ans;
}; // void NM::sigmoid ( vector<double> & x, const double beta, const long iLen ) {
	
	//
	//	Vector Scalar Multiplication: y_bar = a * x_bar
	//
void NM::VecScalarMult ( vector<double> & y, long & yOffset, vector <double> & x, long & xOffset, double & scalarVal, long & iLen ) {

	vector<double>::iterator ix = x.begin() + xOffset;
	vector<double>::iterator iEnd = ix + iLen;
	vector<double>::iterator iy = y.begin() + yOffset;

	for ( ; ix != iEnd; ++ix, ++iy ) {
		*iy = *ix * scalarVal;
	}

}; // void NM::VecScalarMult ( vector<double> & y, long & yOffset, vector <double> & x, long & xOffset, double & scalarVal, long & iLen ) {

	//
	//	RK4 Iteration
	//
void NM::RK4 ( vector<double> & ynp1, long & ynp1Offset, vector <double> & y, long & yOffset, double & alpha, long & iLen  ) {

	vector<double>::iterator iynp1 = ynp1.begin() + ynp1Offset;		// This is the output vector at time n+1.
	vector<double>::iterator iyStart = y.begin() + yOffset;			// This is the input vector at time n.
	vector<double>::iterator iyEnd = iyStart + iLen;
	vector<double>::iterator iy, ik1, ik2, ik3, ik4;
		
	for ( iy = iyStart, ik1 = k1.begin(); iy != iyEnd; ) {
		*ik1++ = alpha * *iy++; 
	} // for ( iy = iyStart; iy != iyEnd; ) {

	for ( iy = iyStart, ik1 = k1.begin(), ik2 = k2.begin(); iy != iyEnd; ) {
		*ik2++ = alpha * ( *iy++ + deltaTD2 * *ik1++ ); 
	} // for ( iy = iyStart; iy != iyEnd; ) {

	for ( iy = iyStart, ik2 = k2.begin(), ik3 = k3.begin(); iy != iyEnd; ) {
		*ik3++ = alpha * ( *iy++ + deltaTD2 * *ik2++ ); 
	} // for ( iy = iyStart; iy != iyEnd; ) {

	for ( iy = iyStart, ik3 = k3.begin(), ik4 = k4.begin(); iy != iyEnd; ) {
		*ik4++ = alpha * ( *iy++ + deltaT * *ik3++ ); 
	} // for ( iy = iyStart; iy != iyEnd; ) {

	for ( iy = iyStart, ik1 = k1.begin(), ik2 = k2.begin(), ik3 = k3.begin(), ik4 = k4.begin(); iy != iyEnd; ) {
		*iynp1++ = *iy++ + deltaTD6 * ( *ik1++ + 2.0 * *ik2++ + 2.0 * *ik3++ + *ik4++ );
	} // for ( iy = iyStart; iy != iyEnd; ) {
	
}; // void RK4 ( vector<double> & ynp1, long & ynp1Offset, vector <double> & y, long & yOffset, double & alpha, long & iLen );

	//
	//	Vector sum
	//
double NM::VecSum ( vector<double> & x ) {

	double ans = 0;
	vector<double>::iterator ix = x.begin();
	vector<double>::iterator ixEnd = x.end();
	for ( ; ix != ixEnd; ++ix ) {
		ans += *ix;
	} 

	return ( ans );

}; // void NM::VecSum ( vector<double> & x )

	//
	//	Vector * Matrix Multiplication: Z (N^2 x N^2) = scalar * X_bar (N^2 x 1) * Y_bar (1 x N^2)
	//
void NM::VecVecMult ( vector<double> & x, long & xOffset, vector<double> & y, long & yOffset, vector<double> & z, long & zOffset, const double & aCoeff, const long & iLen ) {

	vector<double>::iterator ixStart = x.begin() + xOffset;
	vector<double>::iterator ixEnd = ixStart + iLen;
	vector<double>::iterator iy = y.begin() + yOffset;
	vector<double>::iterator iyEnd = iy + iLen;
	vector<double>::iterator iz = z.begin() + zOffset;

	for ( ; iy != iyEnd; ++iy  ) {

		for ( vector<double>::iterator ix = ixStart; ix != ixEnd; ++ix ) {

			*iz++ = aCoeff * *iy * *ix;
		} // 

	} // 	for ( ; iy != iyEnd; ++iy  ) {

}; // void NM::VecVecMult ( vector<double> & x, long & xOffset, vector<double> & y, long & yOffset, vector<double> & z, long & zOffset, const double & aCoeff, const long & iLen ) {

	//
	//	Vector * Matrix Multiplication: Z (N^2 x N^2) = scalar * X_bar (N^2 x 1) * Y_bar (1 x N^2)
	//
	//	Take extra steps to make sure that if a weight is zero going into the adaptation step, it remains zero coming out.
	//
	//	VecVecMultWeights ( r0, offsetToPreviousIter,  r1E, offsetToPreviousIter, w1E0, offsetToPreviousIterWeights, offsetToCurrentIterWeights, wBeta, numValsPerIter );
void NM::VecVecMultWeights ( vector<double> & x, long & xOffset, vector<double> & y, long & yOffset, vector<double> & z, long & chkForZeroOffset, long & zOffset, const double & aCoeff, const long & iLen ) {

	vector<double>::iterator ix = x.begin() + xOffset;
	vector<double>::iterator iy = y.begin() + yOffset;
	vector<double>::iterator iyEnd = iy + iLen;
	vector<double>::iterator iz = z.begin() + zOffset;
	vector<double>::iterator iChk = z.begin() + chkForZeroOffset;
	vector<double>::iterator iTmp;
	vector<double>::iterator iTmpEnd;
	vector<double> tmpVec;

		//	Do the scalar multiplication once.
	tmpVec.assign ( iLen, 0 );
	for ( iTmp = tmpVec.begin(), iTmpEnd = tmpVec.end(); iTmp != iTmpEnd; ) {
		*iTmp++ = aCoeff * *ix++;
	} // for ( iTmp = tmpVec.begin(), iTmpEnd = tmpVec.end(); iTmp != iTmpEnd; ++iTmp, ++ix ) {

	for ( ; iy != iyEnd; ++iy  ) {

		for ( iTmp = tmpVec.begin(), iTmpEnd = tmpVec.end(); iTmp != iTmpEnd; ++iTmp, ++iz, ++iChk ) {

			*iz = (*iChk != 0) ? (*iy * *iTmp) : (0.0);

		} // for ( vector<double>::iterator ix = ixStart; ix != ixEnd; ++ix ) {

	} // 	for ( ; iy != iyEnd; ++iy  ) {
	tmpVec.clear();

}; // void NM::VecScalarMultWeights ( vector<double> & y, long & yOffset, vector <double> & x, long & xOffset, double & scalarVal, long & iLen ) {

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	//
	//  END: NM class definitions.
	//
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////





	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	//
	//  START: NM class definitions that leverage sparse matrices.
	//
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////

	//
	//		Adapt Sparse Weight Matrix
	//
void NM::AdaptSparseWeightMatrix ( vector<SparseWeightMatrix> & w, const long & prevIter, const long & currentIter, const double & wTauAlpha, const double & wBeta,
									vector<double>::iterator & ix, vector<double>::iterator & iy ) {

	for ( long cellID = 0; cellID < N2; ++cellID ) {
		w[cellID].Adapt ( prevIter, currentIter, wTauAlpha, wBeta, ix, ( iy + cellID ) );
	} // for ( long cellID = 0; cellID < N2; ++cellID ) {

}; // void AdaptSparseWeightMatrix ( w1E0_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBeta, r0, offsetToPreviousIter, r1E, offsetToPreviousIter );

	//
	//		Adapt Sparse Weight Matrix
	//
void NM::AdaptSparseWeightMatrixRK4 ( vector<SparseWeightMatrix> & w, const long & prevIter, const long & currentIter, const double & alpha, const double & wBeta,
									vector<double>::iterator & ix, vector<double>::iterator & iy ) {

	for ( long cellID = 0; cellID < N2; ++cellID ) {
		w[cellID].AdaptRK4 ( prevIter, currentIter, alpha, wBeta, deltaT, ix, ( iy + cellID ) );
	} // for ( long cellID = 0; cellID < N2; ++cellID ) {

}; // void AdaptSparseWeightMatrixRK4 ( w1E0_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBeta, r0, offsetToPreviousIter, r1E, offsetToPreviousIter );
	
	//
	//	AddInputPatchSparse
	//
void NM::AddInputPatchSparse ( long & iter, InputPatch & inputPatch, vector<double>::iterator & xStart ) {

	if ( (iter >= inputPatch.offsetPatchInIters) && (iter < inputPatch.offsetPatchInIters + inputPatch.durationPatchInIters) ) {
		
		vector<double>::iterator tmpPtr;
		vector<double>::iterator tmpPtrEnd;

		for ( tmpPtr = inputPatch.inputPatternSparse.begin(), tmpPtrEnd = tmpPtr + inputPatch.inputPatternSparse.size(); tmpPtr != tmpPtrEnd; ++tmpPtr ) {

			*(xStart + ((long) *tmpPtr)) += inputPatch.ampPatch;

		} // for ( tmpPtr = inputPatch.inputPatternSparse.begin(), tmpPtrEnd = tmpPtr + inputPatch.inputPatternSparse.size(); tmpPtr != tmpPtrEnd; ++tmpPtr ) {

	} // if ( (iter >= inputPatch.offsetPatchInIters) && (iter < inputPatch.offsetPatchInIters + inputPatch.durationPatchInIters) ) {

}; // void NM::LoadInputPatchSparse ( long &, vector<InputPatch> &, vector<double> &, long & ) {


	//
	//	AddVecMatMultSparse: y_bar += v_bar * M.  (M is a sparse matrix.)
	//
void NM::AddVecMatMultSparse ( vector<double> & y, const long & yOffset, vector<double> & v, const long & vOffset, vector<SparseWeightMatrix> & M, 
								const long & mOffset, const long & iLen ) {
	vector<double>::iterator iy = y.begin() + yOffset;
	vector<double>::iterator iyEnd = iy + iLen;
	vector<SparseWeightMatrix>::iterator im = M.begin();

	for ( ; iy != iyEnd; ++iy, ++im ) {
		
		*iy += im->DotProduct ( (v.begin() + vOffset), mOffset );

	} // for ( ; iy != iyEnd; ++iy )

};  // void AddVecMatMultSparse ( vector<double> &, const long &, vector<double> &, const long &, vector<double> &, const long &, const long & ) {

	//
	//	Baseline Refinement Sparse
	//
void NM::BaselineRefinementSparse () {

	long iStim = 0;
	long numStim = 0;
	vector<InputPatch> inputPatches;
	vector<long> orderInputPatches;

		//	Generate the input patches.
	//GenDig3RefineInputPatchList ( inputPatches );
	GenDig3RefineInputPatchListSparse ( inputPatches );

		//	Turn ON weight adaptation.
	plasticityFlag = 1;

		//	Determine a random ordering of the inputs.
	orderInputPatches.clear();
	RandomPermutation ( orderInputPatches, 0, (long) inputPatches.size() );
	numStim = (long) inputPatches.size();

		//	Apply inputs in random order.

	for ( long iStim = 0; iStim < numStim; ++iStim ) {

		for ( long iter = 0; iter < numIterPerTrial; ++iter ) {

			long iOffset = ( ( iter == 0 ) ? ( numIterPerTrialMinusOne ) : ( iter - 1 ) ) * numValsPerIter;

			SetNoise (  v0, iOffset, numValsPerIter, -noiseLevel, noiseLevel );
			AddInputPatchSparse ( iter, inputPatches [ orderInputPatches[iStim] ], (v0.begin() + iOffset) );
			Sigmoid ( v0, r0, beta, iOffset, numValsPerIter );
					
			IterateSystem5SparseNewNorm ( iter );
			//IterateSystem5SparseUpdateECellWeightsOnly ( iter );

		} // for ( long iter = 0; iter < numIterPerTrial; ++iter ) {
		
	} // for ( iStim = 0; iStim < numStim; ++iStim

		//	Turn OFF weight adaptation.
	plasticityFlag = 0;

		//	Clean up.
	CleanUpDig3RefineInputPatchList ( inputPatches );
	orderInputPatches.clear();

}; // void NM::BaselineRefinementSparse ()

	//
	//	Baseline Refinement Sparse 2.
	//
void NM::BaselineRefinementSparse2 ( const string & fBaseName, const long & numCycles, const long & iCycle ) {

	// Difference from previous is feature to dump all of the
	//	time series data for each input patch.  Aids debugging and movie making.

	string indexString1 = static_cast<ostringstream*>(&(ostringstream() << numCycles))->str();
	string indexString2 = static_cast<ostringstream*>(&(ostringstream() << iCycle))->str();
	string fName = fBaseName + "." + indexString1 + "." + indexString2 + ".bin";

		// Check whether to dump raw data to file.
	long fileOpenFlag = 1;
	long dumpRawToFile = ( ( dumpRefineStepsRawToFile ) && ( ( iCycle % dumpRefineStepsRawToFileStepSize == 0 ) || ( iCycle == 1 ) ) ) ? (1) : (0);
	
	long iStim = 0;
	long numStim = 0;
	vector<InputPatch> inputPatches;
	vector<long> orderInputPatches;

		//	Generate the input patches.
	GenDig3RefineInputPatchListSparse3 ( inputPatches );

		//	Turn ON weight adaptation.
	plasticityFlag = 1;

		//	Determine a random ordering of the inputs.
	orderInputPatches.clear();
	RandomPermutation ( orderInputPatches, 0, (long) inputPatches.size() );
	numStim = (long) inputPatches.size();

		//	Apply inputs in random order.

	for ( long iStim = 0; iStim < numStim; ++iStim ) {

		for ( long iter = 0; iter < numIterPerTrial; ++iter ) {

			long iOffset = ( ( iter == 0 ) ? ( numIterPerTrialMinusOne ) : ( iter - 1 ) ) * numValsPerIter;

			SetNoise (  v0, iOffset, numValsPerIter, -noiseLevel, noiseLevel );
			AddInputPatchSparse ( iter, inputPatches [ orderInputPatches[iStim] ], (v0.begin() + iOffset) );
			Sigmoid ( v0, r0, beta, iOffset, numValsPerIter );					
			
			IterateSystem5SparseNewNormRK4 ( iter );

		} // for ( long iter = 0; iter < numIterPerTrial; ++iter ) {

		if ( dumpRawToFile ) {
			SaveNumericVectorToBinaryFile(fName, fileOpenFlag, v1E, 0, numValsPerTrial);
			fileOpenFlag = 0;	// This is to make sure all subsequent writes are appends.
 			SaveNumericVectorToBinaryFile(fName, fileOpenFlag, v1I, 0, numValsPerTrial);
			SaveNumericVectorToBinaryFile(fName, fileOpenFlag, v0, 0, numValsPerTrial);
		} // if ( dumpRawToFile ) {
		
	} // for ( iStim = 0; iStim < numStim; ++iStim

		//	Turn OFF weight adaptation.
	plasticityFlag = 0;

		//	Clean up.
	CleanUpDig3RefineInputPatchList ( inputPatches );
	orderInputPatches.clear();

}; // void NM::BaselineRefinementSparse2 ( const string & fBaseName, const long & numCycles, const long & iCycle ) {

	//
	//	ConvertToSparse
	//
void NM::ConvertToSparse ( vector<SparseWeightMatrix> & wSparse, vector<double> & wReg ) {
	
	for ( long cellID = 0; cellID < N2; ++cellID ) {
		wSparse[cellID].RegToSparse ( wReg.begin() + (cellID * N2), numIterPerTrial );
	} // for ( long cellID = 0; cellID < N2; ++cellID ) {

}; // void NM::ConvertToSparse ( vector<SparseWeightMatrix> &, vector<double> & ) {

	//
	//	Initialize Network and Workspace
	//	Use Sparse Representation on the Weight Matrices
	//
void NM::InitializeNetworkAndWorkspaceSparse () {

	sran0 ( 1 );		// Set the random number seed.
	plasticityFlag = 0;	// No weight plasticity until instructed.

		///////////////////////////////////////////////////////////////////
		//
		//	Initial Network & Workspaces.
		//	Set up the initial random network.
		//	For each variable allocate memory corresponding to a single trial as circular
		//	buffering is used to simulate time evolution of arbitrary length.
		//	Inject noise (and normalize as appropriate) just for the first iteration.
		//
		///////////////////////////////////////////////////////////////////
	
	long offsetToLastIter = numIterPerTrialMinusOne * numValsPerIter;

	vWork.assign(numValsPerIter, 0);
	v0.assign ( numValsPerTrial, 0 );
	r0.assign ( numValsPerTrial, 0 );
	AddNoise ( v0, offsetToLastIter, numValsPerIter, -noiseLevel, noiseLevel  );
	Sigmoid ( v0, r0, beta, offsetToLastIter, numValsPerIter );

	v1E.assign ( numValsPerTrial, 0 );
	r1E.assign ( numValsPerTrial, 0 );
	//r1ERFMapRaw.assign ( numValsPerRFMapExp, 0 );
	r1ERFMapExpRes.assign ( numValsPerRFMapRes, 0 );
	AddNoise ( v1E, offsetToLastIter, numValsPerIter, -noiseLevel, noiseLevel  );
	Sigmoid ( v1E, r1E, beta, offsetToLastIter, numValsPerIter );

	v1I.assign ( numValsPerTrial, 0 );
	r1I.assign ( numValsPerTrial, 0 );
	//r1IRFMapRaw.assign ( numValsPerRFMapExp, 0 );
	r1IRFMapExpRes.assign ( numValsPerRFMapRes, 0 );
	AddNoise ( v1I, offsetToLastIter, numValsPerIter, -noiseLevel, noiseLevel  );
	Sigmoid ( v1I, r1I, beta, offsetToLastIter, numValsPerIter );

		//	Initialize and Normalize sparse weight matrices.
	InitializeSparseWeightMatrix ( w1E0_Sparse, g0, numIterPerTrial );						// Initial values placed into the 0th iteration position.
	NormalizeSparseWeightMatrix ( w1E0_Sparse, 0, numIterPerTrialMinusOne, e0WResource );	// Normalized value to last iteration position; seeds the circular buffer.

	InitializeSparseWeightMatrix ( w1EE_Sparse, g0, numIterPerTrial );
	NormalizeSparseWeightMatrix ( w1EE_Sparse, 0, numIterPerTrialMinusOne, eeWResource );

	InitializeSparseWeightMatrix ( w1EI_Sparse, g0, numIterPerTrial );
	NormalizeSparseWeightMatrix ( w1EI_Sparse, 0, numIterPerTrialMinusOne, eiWResource );

	InitializeSparseWeightMatrix ( w1IE_Sparse, g0, numIterPerTrial );
	NormalizeSparseWeightMatrix ( w1IE_Sparse, 0, numIterPerTrialMinusOne, ieWResource );

		//	Preallocate memory for use in RK4.
	k1.assign ( numValsPerIter, 0 );
	k2.assign ( numValsPerIter, 0 );
	k3.assign ( numValsPerIter, 0 );
	k4.assign ( numValsPerIter, 0 );

		//	When debugging, read in known good weight matrices (e.g., from R).  Convert to sparse representation.
#if 0
	fstream tmpFile;
	w1E0.assign ( numWeightsPerTrial, 0 );
	LoadNumericVectorFromBinaryFile ( "w1e0.dat1", tmpFile, 0, w1E0, numWeightsPerIter );	
	ConvertToSparse ( w1E0_Sparse, w1E0 );

	w1EE.assign ( numWeightsPerTrial, 0 );
	LoadNumericVectorFromBinaryFile ( "w1ee.dat1", tmpFile, 0, w1EE, numWeightsPerIter );
	ConvertToSparse ( w1EE_Sparse, w1EE );

	w1EI.assign ( numWeightsPerTrial, 0 );
	LoadNumericVectorFromBinaryFile ( "w1ei.dat1", tmpFile, 0, w1EI, numWeightsPerIter );
	ConvertToSparse ( w1EI_Sparse, w1EI );

	w1IE.assign ( numWeightsPerTrial, 0 );
	LoadNumericVectorFromBinaryFile ( "w1ie.dat1", tmpFile, 0, w1IE, numWeightsPerIter );
	ConvertToSparse ( w1IE_Sparse, w1IE );
#endif
	
}; // void NM::InitializeNetworkAndWorkspace () {

	//
	//	Set initial conditions in the sparse weight matrices.
	//		Values deposited into the 0th iteration positions.
	//
void NM::InitializeSparseWeightMatrix ( vector<SparseWeightMatrix> & w, const long & gVal, const long & numIters ) {

	long cellID = 0;
	long iCellID = 0;
	long jCellID = 0;
	long gridLength = 2 * gVal + 1;
	long gridLengthSquared = gridLength * gridLength;

	SparseWeightMatrix wTmp;

	for ( ; cellID < N2; ++cellID ) {		//	Iterate over each cell that will receive inputs, e.g., Cortical Layer E cells.

		iCellID = GetRow ( cellID, N ) - gVal;		//	Figure out the (i,j) matrix position.
		jCellID = GetCol ( cellID, N ) - gVal;
		
		for ( long j = 0; j < gridLength; ++j ) {	//	Overlay the patch on (i,j) and keep those that map onto the input layer.

			if ( GetLegal ( (jCellID + j), 0, (N-1) ) ) {

				for ( long  i = 0; i < gridLength; ++i ) {

					if ( GetLegal ( iCellID + i, 0, (N-1) )  ) {	//	Sets only the 0th iteration, or initial conditions.

						wTmp.LoadVars ( N, cellID, GetLin ( (iCellID + i), (jCellID + j), N ), (iCellID + i), (jCellID + j), gridLengthSquared,
											r1unif ( wghtMinValue, wghtMaxValue ) );

					} // if ( GetLegal ( iCellID, 0, (N-1) ) & Legal ( jCellID, 0, (N-1) ) ) {

				} // for ( long  i = 0; i < gridLength; ++i ) {

			} // if ( GetLegal ( (jCellID + j), 0, (N-1) ) ) {

		} // for ( long j = 0; j < gridLength; ++j ) {

		wTmp.TrialFillIn ( numIterPerTrial );	//	Set the weight values (to zero) for the remainder of the trial.
		w.push_back ( wTmp );	//	This completes the weight initialization for this connection type for this cell.  Pop it onto the list.
		wTmp.ClearVars();		//	Clear the temp variable before processing the next cell.

	} // for ( ; cellID < N2; ++cellID ) {

}; // void NM::InitializeSparseWeightMatrix ( vector<SparseWeightMatrix> & w, const long & gVal ) {

		//
		//	Iterate the system one time step forward.
		//
void NM::IterateSystem5Sparse ( const long & iter ) {

	long	jter = ( iter == 0 ) ? ( numIterPerTrialMinusOne ) : ( iter - 1 );

	long	offsetToCurrentIter = iter * numValsPerIter;
	long	offsetToPreviousIter = jter * numValsPerIter;

	long	wCurrentIter = numIterPerTrialMinusOne;		// Recall that initial wghts are written into the buffer zone for the last iteration.  This seeds the circular buffer.
	long	wPreviousIter = numIterPerTrialMinusOne;

	if ( plasticityFlag ) {
		wCurrentIter = iter;
		wPreviousIter = jter;
	} // if ( plasticityFlag )

		//	Update v0 term.
	VecScalarMult ( v0, offsetToCurrentIter, v0, offsetToPreviousIter, v0Alpha, numValsPerIter );
	AddNoise (  v0, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	
		//	Update v1E term.
	VecScalarMult ( v1E, offsetToCurrentIter, v1E, offsetToPreviousIter, v1EAlpha, numValsPerIter );
	AddNoise (  v1E, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );

	AddVecMatMultSparse ( v1E, offsetToCurrentIter, r0, offsetToPreviousIter, w1E0_Sparse, wPreviousIter, numValsPerIter );
	AddVecMatMultSparse ( v1E, offsetToCurrentIter, r1E, offsetToPreviousIter, w1EE_Sparse, wPreviousIter, numValsPerIter );
	SubtractVecMatMultSparse ( v1E, offsetToCurrentIter, r1I, offsetToPreviousIter, w1EI_Sparse, wPreviousIter, numValsPerIter );

		//	Update v1I term.
	VecScalarMult ( v1I, offsetToCurrentIter, v1I, offsetToPreviousIter, v1IAlpha, numValsPerIter );
	AddNoise (  v1I, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );

	AddVecMatMultSparse ( v1I, offsetToCurrentIter, r1E, offsetToPreviousIter, w1IE_Sparse, wPreviousIter, numValsPerIter );

		//	Run Sigmoids.
	Sigmoid ( v0, r0, beta, offsetToCurrentIter, numValsPerIter );
	Sigmoid ( v1E, r1E, beta, offsetToCurrentIter, numValsPerIter );
	Sigmoid ( v1I, r1I, beta, offsetToCurrentIter, numValsPerIter );

	if ( plasticityFlag ) {

		vector<double>::iterator r1E_Ptr = r1E.begin() + offsetToPreviousIter;
		vector<double>::iterator r1I_Ptr = r1I.begin() + offsetToPreviousIter;

		AdaptSparseWeightMatrix ( w1E0_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBeta, (r0.begin() + offsetToPreviousIter), r1E_Ptr );
		AdaptSparseWeightMatrix ( w1EI_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBeta, r1I_Ptr, r1E_Ptr );
		AdaptSparseWeightMatrix ( w1IE_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBeta, r1E_Ptr, r1I_Ptr );
		AdaptSparseWeightMatrix ( w1EE_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBeta, r1E_Ptr, r1E_Ptr );

		NormalizeSparseWeightMatrix ( w1E0_Sparse, wCurrentIter, wCurrentIter, eWResource );
		NormalizeSparseWeightMatrix ( w1EE_Sparse, wCurrentIter, wCurrentIter, eWResource );
		NormalizeSparseWeightMatrix ( w1EI_Sparse, wCurrentIter, wCurrentIter, eWResource );
		NormalizeSparseWeightMatrix ( w1IE_Sparse, wCurrentIter, wCurrentIter, iWResource );

	} // if ( plasticityFlag )

};  // void NM::IterateSystem5Sparse ( const long & iter )

		//
		//	Iterate the system one time step forward.  New Normalization strategy.
		//
void NM::IterateSystem5SparseNewNorm ( const long & iter ) {

	long	jter = ( iter == 0 ) ? ( numIterPerTrialMinusOne ) : ( iter - 1 );

	long	offsetToCurrentIter = iter * numValsPerIter;
	long	offsetToPreviousIter = jter * numValsPerIter;

	long	wCurrentIter = numIterPerTrialMinusOne;		// Recall that initial wghts are written into the buffer zone for the last iteration.  This seeds the circular buffer.
	long	wPreviousIter = numIterPerTrialMinusOne;

	if ( plasticityFlag ) {
		wCurrentIter = iter;
		wPreviousIter = jter;
	} // if ( plasticityFlag )

		//	Update v0 term.
	VecScalarMult ( v0, offsetToCurrentIter, v0, offsetToPreviousIter, v0Alpha, numValsPerIter );
	AddNoise (  v0, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	
		//	Update v1E term.
	VecScalarMult ( v1E, offsetToCurrentIter, v1E, offsetToPreviousIter, v1EAlpha, numValsPerIter );
	AddNoise (  v1E, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	AddVecMatMultSparse ( v1E, offsetToCurrentIter, r0, offsetToPreviousIter, w1E0_Sparse, wPreviousIter, numValsPerIter );
	AddVecMatMultSparse ( v1E, offsetToCurrentIter, r1E, offsetToPreviousIter, w1EE_Sparse, wPreviousIter, numValsPerIter );
	SubtractVecMatMultSparse ( v1E, offsetToCurrentIter, r1I, offsetToPreviousIter, w1EI_Sparse, wPreviousIter, numValsPerIter );

		//	Update v1I term.
	VecScalarMult ( v1I, offsetToCurrentIter, v1I, offsetToPreviousIter, v1IAlpha, numValsPerIter );
	AddNoise (  v1I, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	AddVecMatMultSparse ( v1I, offsetToCurrentIter, r1E, offsetToPreviousIter, w1IE_Sparse, wPreviousIter, numValsPerIter );

		//	Run Sigmoids.
	Sigmoid ( v0, r0, beta, offsetToCurrentIter, numValsPerIter );
	Sigmoid ( v1E, r1E, beta, offsetToCurrentIter, numValsPerIter );
	Sigmoid ( v1I, r1I, beta, offsetToCurrentIter, numValsPerIter );

	if ( plasticityFlag ) {

		vector<double>::iterator r1E_Ptr = r1E.begin() + offsetToPreviousIter;
		vector<double>::iterator r1I_Ptr = r1I.begin() + offsetToPreviousIter;

		AdaptSparseWeightMatrix ( w1E0_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBetaE0, ( r0.begin() + offsetToPreviousIter ), r1E_Ptr );
		AdaptSparseWeightMatrix ( w1EI_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBetaEI, r1I_Ptr, r1E_Ptr );
		AdaptSparseWeightMatrix ( w1IE_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBetaIE, r1E_Ptr, r1I_Ptr );
		AdaptSparseWeightMatrix ( w1EE_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBetaEE, r1E_Ptr, r1E_Ptr );

		if ( iter == numIterPerTrialMinusOne ) {		//	Normalize at end of the trial.
			NormalizeSparseWeightMatrix ( w1E0_Sparse, wCurrentIter, wCurrentIter, eWResource );
			NormalizeSparseWeightMatrix ( w1EE_Sparse, wCurrentIter, wCurrentIter, eWResource );
			NormalizeSparseWeightMatrix ( w1EI_Sparse, wCurrentIter, wCurrentIter, eWResource );
			NormalizeSparseWeightMatrix ( w1IE_Sparse, wCurrentIter, wCurrentIter, iWResource );
		} // if ( iter == numIterPerTrialMinusOne ) {

	} // if ( plasticityFlag )

};  // void NM::IterateSystem5SparseNewNorm ( const long & iter )

	//
	//	Iterate the system one time step forward.  New Normalization strategy.
	//
void NM::IterateSystem5ASparseNewNorm (const long & iter) {

	long	jter = (iter == 0) ? (numIterPerTrialMinusOne) : (iter - 1);

	long	offsetToCurrentIter = iter * numValsPerIter;
	long	offsetToPreviousIter = jter * numValsPerIter;

	long	wCurrentIter = numIterPerTrialMinusOne;		// Recall that initial wghts are written into the buffer zone for the last iteration.  This seeds the circular buffer.
	long	wPreviousIter = numIterPerTrialMinusOne;

	long	lZero = 0;

	double	dOne = 1.0;

	if (plasticityFlag) {
		wCurrentIter = iter;
		wPreviousIter = jter;
	} // if ( plasticityFlag )

		//	Update v0 term.
	VecScalarMult(v0, offsetToCurrentIter, v0, offsetToPreviousIter, v0Alpha, numValsPerIter);
	AddNoise(v0, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel);

		//	Update v1E term.
	VecScalarMult(v1E, offsetToCurrentIter, v1E, offsetToPreviousIter, v1EAlpha, numValsPerIter);
	SetNoise(vWork, lZero, numValsPerIter, -noiseLevel, noiseLevel);
	AddVecMatMultSparse(vWork, lZero, r0, offsetToPreviousIter, w1E0_Sparse, wPreviousIter, numValsPerIter);
	AddVecMatMultSparse(vWork, lZero, r1E, offsetToPreviousIter, w1EE_Sparse, wPreviousIter, numValsPerIter);
	SubtractVecMatMultSparse(vWork, lZero, r1I, offsetToPreviousIter, w1EI_Sparse, wPreviousIter, numValsPerIter);
	VecScalarMult(vWork, lZero, vWork, lZero, dOne, numValsPerIter);
	VecAdd(v1E, offsetToCurrentIter, v1E, offsetToCurrentIter, vWork, lZero, numValsPerIter);

		//	Update v1I term.
	VecScalarMult(v1I, offsetToCurrentIter, v1I, offsetToPreviousIter, v1IAlpha, numValsPerIter);
	SetNoise(vWork, lZero, numValsPerIter, -noiseLevel, noiseLevel);
	AddVecMatMultSparse(vWork, lZero, r1E, offsetToPreviousIter, w1IE_Sparse, wPreviousIter, numValsPerIter);
	VecScalarMult(vWork, lZero, vWork, lZero, dOne, numValsPerIter);
	VecAdd(v1I, offsetToCurrentIter, v1I, offsetToCurrentIter, vWork, lZero, numValsPerIter);

		//	Run Sigmoids.
	Sigmoid(v0, r0, beta, offsetToCurrentIter, numValsPerIter);
	Sigmoid(v1E, r1E, beta, offsetToCurrentIter, numValsPerIter);
	Sigmoid(v1I, r1I, beta, offsetToCurrentIter, numValsPerIter);

	if (plasticityFlag) {

		vector<double>::iterator r1E_Ptr = r1E.begin() + offsetToPreviousIter;
		vector<double>::iterator r1I_Ptr = r1I.begin() + offsetToPreviousIter;

		AdaptSparseWeightMatrix(w1E0_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBetaE0, (r0.begin() + offsetToPreviousIter), r1E_Ptr);
		AdaptSparseWeightMatrix(w1EI_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBetaEI, r1I_Ptr, r1E_Ptr);
		AdaptSparseWeightMatrix(w1IE_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBetaIE, r1E_Ptr, r1I_Ptr);
		AdaptSparseWeightMatrix(w1EE_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBetaEE, r1E_Ptr, r1E_Ptr);

		if (iter == (numIterPerTrial - 1)) {	//	Normalize at end of the trial.
			NormalizeSparseWeightMatrix(w1E0_Sparse, wCurrentIter, wCurrentIter, eWResource);
			NormalizeSparseWeightMatrix(w1EE_Sparse, wCurrentIter, wCurrentIter, eWResource);
			NormalizeSparseWeightMatrix(w1EI_Sparse, wCurrentIter, wCurrentIter, eWResource);
			NormalizeSparseWeightMatrix(w1IE_Sparse, wCurrentIter, wCurrentIter, iWResource);
		} // if ( iter == (numIterPerTrial - 1) ) {

	} // if ( plasticityFlag )

};  // void NM::IterateSystem5ASparseNewNorm ( const long & iter )

		//
		//	Iterate the system one time step forward.  New Normalization strategy.
		//
void NM::IterateSystem5SparseNewNormRK4 ( const long & iter ) {

	long	jter = ( iter == 0 ) ? ( numIterPerTrialMinusOne ) : ( iter - 1 );

	long	offsetToCurrentIter = iter * numValsPerIter;
	long	offsetToPreviousIter = jter * numValsPerIter;

	long	wCurrentIter = numIterPerTrialMinusOne;		// Recall that initial wghts are written into the buffer zone for the last iteration.  This seeds the circular buffer.
	long	wPreviousIter = numIterPerTrialMinusOne;

	if ( plasticityFlag ) {
		wCurrentIter = iter;
		wPreviousIter = jter;
	} // if ( plasticityFlag )

		//	Update v0 term.
	//VecScalarMult ( v0, offsetToCurrentIter, v0, offsetToPreviousIter, v0Alpha, numValsPerIter );
	RK4 ( v0, offsetToCurrentIter, v0, offsetToPreviousIter, minusOneOverTau, numValsPerIter );
	AddNoise (  v0, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	
		//	Update v1E term.
	//VecScalarMult ( v1E, offsetToCurrentIter, v1E, offsetToPreviousIter, v1EAlpha, numValsPerIter );
	RK4 ( v1E, offsetToCurrentIter, v1E, offsetToPreviousIter, minusOneOverTau, numValsPerIter );
	AddNoise (  v1E, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	AddVecMatMultSparse ( v1E, offsetToCurrentIter, r0, offsetToPreviousIter, w1E0_Sparse, wPreviousIter, numValsPerIter );
	AddVecMatMultSparse ( v1E, offsetToCurrentIter, r1E, offsetToPreviousIter, w1EE_Sparse, wPreviousIter, numValsPerIter );
	SubtractVecMatMultSparse ( v1E, offsetToCurrentIter, r1I, offsetToPreviousIter, w1EI_Sparse, wPreviousIter, numValsPerIter );

		//	Update v1I term.
	//VecScalarMult ( v1I, offsetToCurrentIter, v1I, offsetToPreviousIter, v1IAlpha, numValsPerIter );
	RK4 ( v1I, offsetToCurrentIter, v1I, offsetToPreviousIter, minusOneOverTau, numValsPerIter );
	AddNoise (  v1I, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	AddVecMatMultSparse ( v1I, offsetToCurrentIter, r1E, offsetToPreviousIter, w1IE_Sparse, wPreviousIter, numValsPerIter );

		//	Run Sigmoids.
	Sigmoid ( v0, r0, beta, offsetToCurrentIter, numValsPerIter );
	Sigmoid ( v1E, r1E, beta, offsetToCurrentIter, numValsPerIter );
	Sigmoid ( v1I, r1I, beta, offsetToCurrentIter, numValsPerIter );

	if ( plasticityFlag ) {

		vector<double>::iterator r1E_Ptr = r1E.begin() + offsetToPreviousIter;
		vector<double>::iterator r1I_Ptr = r1I.begin() + offsetToPreviousIter;

		AdaptSparseWeightMatrixRK4 ( w1E0_Sparse, wPreviousIter, wCurrentIter, minusOneOverWTau, wBetaE0, ( r0.begin() + offsetToPreviousIter ), r1E_Ptr );
		AdaptSparseWeightMatrixRK4 ( w1EI_Sparse, wPreviousIter, wCurrentIter, minusOneOverWTau, wBetaEI, r1I_Ptr, r1E_Ptr );
		AdaptSparseWeightMatrixRK4 ( w1IE_Sparse, wPreviousIter, wCurrentIter, minusOneOverWTau, wBetaIE, r1E_Ptr, r1I_Ptr );
		AdaptSparseWeightMatrixRK4 ( w1EE_Sparse, wPreviousIter, wCurrentIter, minusOneOverWTau, wBetaEE, r1E_Ptr, r1E_Ptr );

		if ( iter == numIterPerTrialMinusOne ) {		//	Normalize at end of the trial.
			NormalizeSparseWeightMatrix ( w1E0_Sparse, wCurrentIter, wCurrentIter, e0WResource );
			NormalizeSparseWeightMatrix ( w1EE_Sparse, wCurrentIter, wCurrentIter, eeWResource );
			NormalizeSparseWeightMatrix ( w1EI_Sparse, wCurrentIter, wCurrentIter, eiWResource );
			NormalizeSparseWeightMatrix ( w1IE_Sparse, wCurrentIter, wCurrentIter, ieWResource );
		} // if ( iter == numIterPerTrialMinusOne ) {

	} // if ( plasticityFlag )

};  // void NM::IterateSystem5SparseNewNormRK4 ( const long & iter )

		//
		//	Iterate the system one time step forward.
		//
void NM::IterateSystem5SparseUpdateECellWeightsOnly ( const long & iter ) {

	long	jter = ( iter == 0 ) ? ( numIterPerTrialMinusOne ) : ( iter - 1 );

	long	offsetToCurrentIter = iter * numValsPerIter;
	long	offsetToPreviousIter = jter * numValsPerIter;

	long	wCurrentIter = numIterPerTrialMinusOne;		// Recall that initial wghts are written into the buffer zone for the last iteration.  This seeds the circular buffer.
	long	wPreviousIter = numIterPerTrialMinusOne;

	if ( plasticityFlag ) {
		wCurrentIter = iter;
		wPreviousIter = jter;
	} // if ( plasticityFlag )

		//	Update v0 term.
	VecScalarMult ( v0, offsetToCurrentIter, v0, offsetToPreviousIter, v0Alpha, numValsPerIter );
	AddNoise (  v0, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );
	
		//	Update v1E term.
	VecScalarMult ( v1E, offsetToCurrentIter, v1E, offsetToPreviousIter, v1EAlpha, numValsPerIter );
	AddNoise (  v1E, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );

	AddVecMatMultSparse ( v1E, offsetToCurrentIter, r0, offsetToPreviousIter, w1E0_Sparse, wPreviousIter, numValsPerIter );
	AddVecMatMultSparse ( v1E, offsetToCurrentIter, r1E, offsetToPreviousIter, w1EE_Sparse, wPreviousIter, numValsPerIter );
	SubtractVecMatMultSparse ( v1E, offsetToCurrentIter, r1I, offsetToPreviousIter, w1EI_Sparse, wPreviousIter, numValsPerIter );

		//	Update v1I term.
	VecScalarMult ( v1I, offsetToCurrentIter, v1I, offsetToPreviousIter, v1IAlpha, numValsPerIter );
	AddNoise (  v1I, offsetToCurrentIter, numValsPerIter, -noiseLevel, noiseLevel );

	AddVecMatMultSparse ( v1I, offsetToCurrentIter, r1E, offsetToPreviousIter, w1IE_Sparse, wPreviousIter, numValsPerIter );

		//	Run Sigmoids.
	Sigmoid ( v0, r0, beta, offsetToCurrentIter, numValsPerIter );
	Sigmoid ( v1E, r1E, beta, offsetToCurrentIter, numValsPerIter );
	Sigmoid ( v1I, r1I, beta, offsetToCurrentIter, numValsPerIter );

	if ( plasticityFlag ) {

		vector<double>::iterator r1E_Ptr = r1E.begin() + offsetToPreviousIter;
		vector<double>::iterator r1I_Ptr = r1I.begin() + offsetToPreviousIter;

		AdaptSparseWeightMatrix ( w1E0_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBeta, (r0.begin() + offsetToPreviousIter), r1E_Ptr );
		AdaptSparseWeightMatrix ( w1EI_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBeta, r1I_Ptr, r1E_Ptr );
		//AdaptSparseWeightMatrix ( w1IE_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBeta, r1E_Ptr, r1I_Ptr );
		//AdaptSparseWeightMatrix ( w1EE_Sparse, wPreviousIter, wCurrentIter, wTauAlpha, wBeta, r1E_Ptr, r1E_Ptr );

		AdaptSparseWeightMatrix ( w1IE_Sparse, wPreviousIter, wCurrentIter, 1.0, 0.0, r1E_Ptr, r1I_Ptr );
		AdaptSparseWeightMatrix ( w1EE_Sparse, wPreviousIter, wCurrentIter, 1.0, 0.0, r1E_Ptr, r1E_Ptr );

		NormalizeSparseWeightMatrix ( w1E0_Sparse, wCurrentIter, wCurrentIter, eWResource );
		NormalizeSparseWeightMatrix ( w1EE_Sparse, wCurrentIter, wCurrentIter, eWResource );
		NormalizeSparseWeightMatrix ( w1EI_Sparse, wCurrentIter, wCurrentIter, eWResource );
		NormalizeSparseWeightMatrix ( w1IE_Sparse, wCurrentIter, wCurrentIter, iWResource );

	} // if ( plasticityFlag )

};  // void NM::IterateSystem5SparseUpdateECellWeightsOnly ( const long & iter )

	//
	//	Load the Network from a binary file.
	//
void NM::LoadNetworkFromFileSparse ( const string & fBaseName, const long & index1, const long & index2 ) {

	fstream readFile;

	string indexString1 = static_cast<ostringstream*>( &(ostringstream() << index1) )->str();
	string indexString2 = static_cast<ostringstream*>( &(ostringstream() << index2) )->str();
	string fName = fBaseName + "." + indexString1 + "." + indexString2 + ".bin";

	readFile.open( fName.c_str(), ios::in | ios::binary );

	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, v0, numValsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, r0, numValsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, v1E, numValsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, r1E, numValsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, v1I, numValsPerTrial  );
	LoadNumericVectorFromBinaryFile ( fName, readFile, 1, r1I, numValsPerTrial  );

	readFile.close();

	fName = fBaseName + "." + indexString1 + "." + indexString2 + ".w1E0";
	ReadSparseWeightMatrixFromFile ( fName, w1E0_Sparse.begin(), w1E0_Sparse.end() );

	fName = fBaseName + "." + indexString1 + "." + indexString2 + ".w1EE";
	ReadSparseWeightMatrixFromFile ( fName, w1EE_Sparse.begin(), w1EE_Sparse.end()  );

	fName = fBaseName + "." + indexString1 + "." + indexString2 + ".w1IE";
	ReadSparseWeightMatrixFromFile ( fName, w1IE_Sparse.begin(), w1IE_Sparse.end()  );

	fName = fBaseName + "." + indexString1 + "." + indexString2 + ".w1EI";
	ReadSparseWeightMatrixFromFile ( fName, w1EI_Sparse.begin(), w1EI_Sparse.end()  );
	
}; // void NM::SaveNetworkToFile ( const string &, const long &, const long & ) {

	//
	//	Normalize Sparse Weight Matrix
	//
void NM::NormalizeSparseWeightMatrix ( vector<SparseWeightMatrix> & w, const long & prevIter, const long & currentIter, const double & wResource ) {

	for ( long cellID = 0; cellID < N2; ++cellID ) {
		w[cellID].Normalize ( prevIter, currentIter, wResource );
	} // for ( long cellID = 0; cellID < N2; ++cellID ) {

}; // void SparseWeightMatrix::NormalizeSparseWeightMatrix ( vector<SparseWeightMatrix> & w, const long & offste, const double & eWResource ) {

	//
	//	ReadSpareseWeightMatrixToFile
	//		Only the weight values at the last iteration in the trial are read.

void NM::ReadSparseWeightMatrixFromFile ( const string & fName, vector<SparseWeightMatrix>::iterator & wPtrStart, vector<SparseWeightMatrix>::iterator & wPtrEnd ) {
	
	fstream inFile;
	inFile.open( fName.c_str(), ios::in | ios::binary );
	for ( vector<SparseWeightMatrix>::iterator wPtr = wPtrStart; wPtr != wPtrEnd; ++wPtr ) {
		wPtr->ReadFromFile ( inFile );
	} // for ( vector<SparseWeightMatrix>::iterator wPtr = wPtrStart; wPtr != wPtrEnd; ++wPtr ) {
	inFile.close();

}; // void NM::WriteSparseWeightMatrixToFile ( const string & fName, vector<SparseWeightMatrix>::iterator & wPtr, vector<SparseWeightMatrix>::iterator & wPtrEnd ) { 

	//
	//	Receptive Field Map Sparse.
	//
void NM::ReceptiveFieldMapSparse () {

	r1ERFMapRaw.assign ( numValsPerRFMapExp, 0 );
	r1IRFMapRaw.assign ( numValsPerRFMapExp, 0 );

		//
		//	Do the receptive field mapping experiment.
		//	Iterate over each S layer node with an RF Probe and store the network response.
		//
	for ( long iCell = 0; iCell < N2; ++iCell ) {

		long offsetToRFMapRaw = 0;
		vector<double>::iterator iv0 = v0.begin() + iCell;
		vector<double>::iterator ir0 = r0.begin() + iCell;
		
		for ( long iter = 1; iter < numIterPerTrial; ++iter ) {

				//	Stimulus ON/OFF.
			if ( (iter > offsetToRFProbeInIters ) && (iter <= timeOffRFProbeInIters ) ) {
				long itmp = ( iter - 1 ) * numValsPerIter;
				*( iv0 + itmp ) = magRFProbe + r1unif ( -noiseLevel, noiseLevel );
				*( ir0 + itmp ) = Sigmoid1 ( *( iv0 + itmp ), beta );
			} // if ( (iter > offsetToRFProbeInIters ) && (iter <= timeOffRFProbeInIters ) ) {				

			IterateSystem5Sparse ( iter );
			//IterateSystem5SparseUpdateECellWeightsOnly ( iter );

		} // for ( long iter = 1; iter < numIterPerTrial; ++iter ) {

			//	Store the raw r1E and r1I data for later analysis.
		offsetToRFMapRaw = iCell * numValsPerTrial;
		VecCopy ( r1E, 0, r1ERFMapRaw, offsetToRFMapRaw, numValsPerTrial );
		VecCopy ( r1I, 0, r1IRFMapRaw, offsetToRFMapRaw, numValsPerTrial );

	}  // for ( long iCell = 0; iCell < N2; ++iCell )

	ReceptiveFieldProbeAnalysis ();	//	Process the raw RF probe data to generate Experimental Results (an N^2 x N^2 matrix).

}; // void NM::ReceptiveFieldMapSparse ()

  //
  //	Receptive Field Map Sparse Streamed.
  //
void NM::ReceptiveFieldMapSparseStreamed ( const string & fBaseName, const long & numCycles, const long & iCycle ) {

	//
	//	Difference from ReceptiveFieldMapSparse is that this routine does not operate in batch mode.
	//	Iterate over each S layer node with an RF Probe calculate the cortical layer response maps, but dump the trial to a file.
	//

	string indexString1 = static_cast<ostringstream*>(&(ostringstream() << numCycles))->str();
	string indexString2 = static_cast<ostringstream*>(&(ostringstream() << iCycle))->str();
	string fName = fBaseName + "." + indexString1 + "." + indexString2 + ".bin";

	vector<double>::iterator ixR1E_Results = r1ERFMapExpRes.begin();
	vector<double>::iterator ixR1I_Results = r1IRFMapExpRes.begin();

		// Check whether to dump raw data to file.
	long fileOpenFlag = 1;
	long dumpRawToFile = 0;
	if (dumpRFMapRawToFile) {
		dumpRawToFile = ( ( ( (iCycle) % dumpRFMapRawToFile) == 0 ) || ( iCycle == 1 ) ) ? 1 : 0;
	}

	for (long iCell = 0; iCell < N2; ++iCell) {

		vector<double>::iterator iv0 = v0.begin() + iCell;
		vector<double>::iterator ir0 = r0.begin() + iCell;

		for (long iter = 1; iter < numIterPerTrial; ++iter) {
		//for (long iter = 1; iter <= timeOffRFProbeInIters; ++iter) {
				//	Stimulus ON/OFF.
			if ((iter > offsetToRFProbeInIters) && (iter <= timeOffRFProbeInIters)) {
				long itmp = (iter - 1) * numValsPerIter;
				*(iv0 + itmp) = magRFProbe + r1unif(-noiseLevel, noiseLevel);
				*(ir0 + itmp) = Sigmoid1(*(iv0 + itmp), beta);
			} // if ( (iter > offsetToRFProbeInIters ) && (iter <= timeOffRFProbeInIters ) ) {				

			IterateSystem5SparseNewNormRK4 ( iter );

		} // for ( long iter = 1; iter < numIterPerTrial; ++iter ) {

		//StimRespMag2(r1E, ( r1ERFMapExpRes.begin() + iCell * N2 ) );
		//StimRespMag2(r1I, ( r1IRFMapExpRes.begin() + iCell * N2 ) );

		StimRespMag(r1E, ( r1ERFMapExpRes.begin() + iCell * N2 ) );
		StimRespMag(r1I, ( r1IRFMapExpRes.begin() + iCell * N2 ) );

		if ( 0 ) {
			SaveNumericVectorToBinaryFile(fName, fileOpenFlag, v1E, 0, numValsPerTrial);
			fileOpenFlag = 0;	// This is to make sure all subsequent writes are appends.
 			SaveNumericVectorToBinaryFile(fName, fileOpenFlag, v1I, 0, numValsPerTrial);
			SaveNumericVectorToBinaryFile(fName, fileOpenFlag, v0, 0, numValsPerTrial);
		} // if ( dumpRawToFile ) {

	}  // for ( long iCell = 0; iCell < N2; ++iCell )

}; // void NM::ReceptiveFieldMapSparseStreamed ( const string & fBaseName, const long & numCycles, const long & iCycle ) {

	//
	//	Save the Network to a binary file.
	//
void NM::SaveNetworkToFileSparse ( const string & fBaseName, const long & index1, const long & index2 ) {

	string indexString1 = static_cast<ostringstream*>( &(ostringstream() << index1) )->str();
	string indexString2 = static_cast<ostringstream*>( &(ostringstream() << index2) )->str();
	string fName = fBaseName + "." + indexString1 + "." + indexString2 + ".bin";

	SaveNumericVectorToBinaryFile ( fName, 1, v0, 0, numValsPerTrial );
	SaveNumericVectorToBinaryFile ( fName, 0, r0, 0, numValsPerTrial );
	SaveNumericVectorToBinaryFile ( fName, 0, v1E, 0, numValsPerTrial );
	SaveNumericVectorToBinaryFile ( fName, 0, r1E, 0, numValsPerTrial );
	SaveNumericVectorToBinaryFile ( fName, 0, v1I, 0, numValsPerTrial );
	SaveNumericVectorToBinaryFile ( fName, 0, r1I, 0, numValsPerTrial );

	fName = fBaseName + "." + indexString1 + "." + indexString2 + ".w1E0";
	WriteSparseWeightMatrixToFile ( fName, w1E0_Sparse.begin(), w1E0_Sparse.end() );

	fName = fBaseName + "." + indexString1 + "." + indexString2 + ".w1EE";
	WriteSparseWeightMatrixToFile ( fName, w1EE_Sparse.begin(), w1EE_Sparse.end()  );

	fName = fBaseName + "." + indexString1 + "." + indexString2 + ".w1IE";
	WriteSparseWeightMatrixToFile ( fName, w1IE_Sparse.begin(), w1IE_Sparse.end()  );

	fName = fBaseName + "." + indexString1 + "." + indexString2 + ".w1EI";
	WriteSparseWeightMatrixToFile ( fName, w1EI_Sparse.begin(), w1EI_Sparse.end()  );
	
}; // void NM::SaveNetworkToFileSparse ( const string &, const long &, const long & ) {

	//
	//	Save the Network to a binary file. Specialized for Focal Stim.  (Klunky)
	//
void NM::SaveNetworkToFileSparseFocalStim ( const string & fBaseName, const long & index1, const long & index2, const long & index3, const long & index4 ) {

	string indexString1 = static_cast<ostringstream*>( &(ostringstream() << index1) )->str();
	string indexString2 = static_cast<ostringstream*>( &(ostringstream() << index2) )->str();
	string indexString3 = static_cast<ostringstream*>( &(ostringstream() << index3) )->str();
	string indexString4 = static_cast<ostringstream*>( &(ostringstream() << index4) )->str();
	string fName = fBaseName + "." + indexString1 + "." + indexString2 + "." + indexString3 + "." + indexString4 + ".bin";

	SaveNumericVectorToBinaryFile ( fName, 1, v0, 0, numValsPerTrial );
	SaveNumericVectorToBinaryFile ( fName, 0, r0, 0, numValsPerTrial );
	SaveNumericVectorToBinaryFile ( fName, 0, v1E, 0, numValsPerTrial );
	SaveNumericVectorToBinaryFile ( fName, 0, r1E, 0, numValsPerTrial );
	SaveNumericVectorToBinaryFile ( fName, 0, v1I, 0, numValsPerTrial );
	SaveNumericVectorToBinaryFile ( fName, 0, r1I, 0, numValsPerTrial );

	fName = fBaseName + "." + indexString1 + "." + indexString2 + "." + indexString3 + "." + indexString4 + ".w1E0";
	WriteSparseWeightMatrixToFile ( fName, w1E0_Sparse.begin(), w1E0_Sparse.end() );

	fName = fBaseName + "." + indexString1 + "." + indexString2 + "." + indexString3 + "." + indexString4 + ".w1EE";
	WriteSparseWeightMatrixToFile ( fName, w1EE_Sparse.begin(), w1EE_Sparse.end()  );

	fName = fBaseName + "." + indexString1 + "." + indexString2 + "." + indexString3 + "." + indexString4 + ".w1IE";
	WriteSparseWeightMatrixToFile ( fName, w1IE_Sparse.begin(), w1IE_Sparse.end()  );

	fName = fBaseName + "." + indexString1 + "." + indexString2 + "." + indexString3 + "." + indexString4 + ".w1EI";
	WriteSparseWeightMatrixToFile ( fName, w1EI_Sparse.begin(), w1EI_Sparse.end()  );
	
}; // void NM::SaveNetworkToFileSparse ( const string &, const long &, const long & ) {

	//
	//	Baseline Refinement Sparse
	//
void NM::SelectiveStimSparse () {

	long iStim = 0;
	long numStim = 0;
	vector<InputPatch> inputPatches;
	vector<long> orderInputPatches;

		//	Generate the input patches.
	GenDig3SelectiveStimInputPatchListSparse ( inputPatches );

		//	Turn ON weight adaptation.
	plasticityFlag = 1;

		//	Determine a random ordering of the inputs.
	orderInputPatches.clear();
	RandomPermutation ( orderInputPatches, 0, (long) inputPatches.size() );
	numStim = (long) inputPatches.size();

		//	Apply inputs in random order.

	for ( long iStim = 0; iStim < numStim; ++iStim ) {

		for ( long iter = 0; iter < numIterPerTrial; ++iter ) {

			long iOffset = ( ( iter == 0 ) ? ( numIterPerTrialMinusOne ) : ( iter - 1 ) ) * numValsPerIter;

			SetNoise ( v0, iOffset, numValsPerIter, -noiseLevel, noiseLevel );
			AddInputPatchSparse ( iter, inputPatches [ orderInputPatches[iStim] ], (v0.begin() + iOffset) );
			Sigmoid ( v0, r0, beta, iOffset, numValsPerIter );
					
			IterateSystem5Sparse ( iter );
			//IterateSystem5SparseUpdateECellWeightsOnly ( iter );

		} // for ( long iter = 0; iter < numIterPerTrial; ++iter ) {
		
	} // for ( iStim = 0; iStim < numStim; ++iStim

		//	Turn OFF weight adaptation.
	plasticityFlag = 0;

		//	Clean up.
	CleanUpDig3RefineInputPatchList ( inputPatches );
	orderInputPatches.clear();

}; // void NM::SelectiveStimSparse ()

	//
	//	StimRespMag
	//
void NM::StimRespMag( vector<double> & x, vector<double>::iterator & xMagResp ) {
	
	double meanCorrection = (double) offsetToRFProbeInIters / durRFProbeInIters;
	vector<double>::iterator tmpMagResp = xMagResp;

			//	jCell indexes the cells on Layer C.
	for (long jCell = 0; jCell < N2; ++jCell) {

			//	Determine the average activity in the pre-stimulus time period.

		vector<double>::iterator ix = x.begin() + jCell;		
		double preTemp = 0;
		for (long k = 0; k < offsetToRFProbeInIters; ++k) {

				preTemp += *ix;
				ix += numValsPerIter;		//	Same cell, but next iteration.

		} // for ( long k = 0; k < offsetToRFProbeInIters; ++k ) {		

			  //	Determine the average activity during the stimulus time period.

		ix = x.begin() + (offsetToRFProbeInIters + 1) * numValsPerIter + jCell;
		double postTemp = 0;
		for (long k = 0; k < durRFProbeInIters; ++k) {

			postTemp += *ix;
			ix += numValsPerIter;		//	Same cell, but next iteration.

		} // for ( long k = 0; k < durRFProbeInIters; ++k ) {

			  //	Store the ratio average responses.

		*tmpMagResp = (postTemp / preTemp) * meanCorrection;
		++tmpMagResp;

	} // 		for ( long jCell = 0; jCell < N2; ++jCell ) {

}; // void NM::StimRespMag(vector<double> & x) {

	//
	//	StimRespMag
	//
void NM::StimRespMag2( vector<double> & x, vector<double>::iterator & xMagResp ) {

		//	Difference from StimRespMap is that this one takes ratio of stim response peak to pre-stim average
	
	// double meanCorrection = (double) offsetToRFProbeInIters / durRFProbeInIters;
	double preStimDur = (double) offsetToRFProbeInIters;

			//	jCell indexes the cells on Layer C.
	for (long jCell = 0; jCell < N2; ++jCell) {

			//	Determine the average activity in the pre-stimulus time period.

		vector<double>::iterator ix = x.begin() + jCell;		
		double preTemp = 0;
		for (long k = 0; k < offsetToRFProbeInIters; ++k) {

				//preTemp += *ix;
				preTemp = ( preTemp > *ix ) ? preTemp : *ix;
				ix += numValsPerIter;		//	Same cell, but next iteration.

		} // for ( long k = 0; k < offsetToRFProbeInIters; ++k ) {		

			  //	Determine the average activity during the stimulus time period.

		ix = x.begin() + (offsetToRFProbeInIters + 1) * numValsPerIter + jCell;
		double postTemp = 0;
		for (long k = 0; k < durRFProbeInIters; ++k) {

			postTemp += *ix;
			//postTemp = ( postTemp > *ix ) ? postTemp : *ix;
			ix += numValsPerIter;		//	Same cell, but next iteration.

		} // for ( long k = 0; k < durRFProbeInIters; ++k ) {

			  //	Store the ratio average responses.

		//*xMagResp = (postTemp / preTemp) * meanCorrection;
		*xMagResp = postTemp / ( preTemp / preStimDur );
		++xMagResp;

	} // 		for ( long jCell = 0; jCell < N2; ++jCell ) {

}; // void NM::StimRespMag2(vector<double> & x) {


	//
	//	SubtractVecMatMultSparse: y_bar += v_bar * M.  (M is a sparse matrix.)
	//
void NM::SubtractVecMatMultSparse ( vector<double> & y, const long & yOffset, vector<double> & v, const long & vOffset, vector<SparseWeightMatrix> & M, 
								const long & mOffset, const long & iLen ) {
	vector<double>::iterator iy = y.begin() + yOffset;
	vector<double>::iterator iyEnd = iy + iLen;
	vector<SparseWeightMatrix>::iterator im = M.begin();

	for ( ; iy != iyEnd; ++iy, ++im ) {
		
		*iy -= im->DotProduct ( (v.begin() + vOffset), mOffset );

	} // for ( ; iy != iyEnd; ++iy )
	
};  // void SubtractVecMatMultSparse ( vector<double> &, const long &, vector<double> &, const long &, vector<double> &, const long &, const long & ) {

	//
	//	SyndactylySparse
	//
void NM::SyndactylySparse () {

	long iStim = 0;
	long numStim = 0;
	vector<InputPatch> inputPatches;
	vector<long> orderInputPatches;

		//	Generate the input patches.
	//GenDig3SyndactylyInputPatchList  ( inputPatches );
	GenDig3SyndactylyInputPatchListSparse  ( inputPatches );

		//	Turn ON weight adaptation.
	plasticityFlag = 1;

		//	Determine a random ordering of the inputs.
	orderInputPatches.clear();
	RandomPermutation ( orderInputPatches, 0, (long) inputPatches.size() );
	numStim = (long) inputPatches.size();

		//	Apply inputs in random order.

	for ( long iStim = 0; iStim < numStim; ++iStim ) {

		for ( long iter = 0; iter < numIterPerTrial; ++iter ) {

			long iOffset = ( ( iter == 0 ) ? ( numIterPerTrialMinusOne ) : ( iter - 1 ) ) * numValsPerIter;

			SetNoise (  v0, iOffset, numValsPerIter, -noiseLevel, noiseLevel );
			AddInputPatchSparse ( iter, inputPatches [ orderInputPatches[iStim] ], (v0.begin() + iOffset) );
			Sigmoid ( v0, r0, beta, iOffset, numValsPerIter );
					
			IterateSystem5SparseNewNorm ( iter );

		} // for ( long iter = 0; iter < numIterPerTrial; ++iter ) {
		
	} // for ( iStim = 0; iStim < numStim; ++iStim

		//	Turn OFF weight adaptation.
	plasticityFlag = 0;

		//	Clean up.
	CleanUpDig3SyndactylyInputPatchList  ( inputPatches );
	orderInputPatches.clear();

}; // void NM::SyndactylySparse ()

	//
	//	Syndactyly Refinement Sparse 2.
	//
void NM::SyndactylySparse2 ( const string & fBaseName, const long & numCycles, const long & iCycle ) {

	// Difference from previous is feature to dump all of the
	//	time series data for each input patch.  Aids debugging and movie making.

	string indexString1 = static_cast<ostringstream*>(&(ostringstream() << numCycles))->str();
	string indexString2 = static_cast<ostringstream*>(&(ostringstream() << iCycle))->str();
	string fName = fBaseName + "." + indexString1 + "." + indexString2 + ".bin";

		// Check whether to dump raw data to file.
	long fileOpenFlag = 1;
	long dumpRawToFile = ( ( dumpRefineStepsRawToFile ) && ( ( iCycle % dumpRefineStepsRawToFileStepSize == 0 ) || ( iCycle == 1 ) ) ) ? (1) : (0);

	long iStim = 0;
	long numStim = 0;
	vector<InputPatch> inputPatches;
	vector<long> orderInputPatches;

		//	Generate the input patches.
	GenDig3SyndactylyInputPatchListSparse  ( inputPatches );

		//	Turn ON weight adaptation.
	plasticityFlag = 1;

		//	Determine a random ordering of the inputs.
	orderInputPatches.clear();
	RandomPermutation ( orderInputPatches, 0, (long) inputPatches.size() );
	numStim = (long) inputPatches.size();

		//	Apply inputs in random order.

	for ( long iStim = 0; iStim < numStim; ++iStim ) {

		for ( long iter = 0; iter < numIterPerTrial; ++iter ) {

			long iOffset = ( ( iter == 0 ) ? ( numIterPerTrialMinusOne ) : ( iter - 1 ) ) * numValsPerIter;

			SetNoise (  v0, iOffset, numValsPerIter, -noiseLevel, noiseLevel );
			AddInputPatchSparse ( iter, inputPatches [ orderInputPatches[iStim] ], (v0.begin() + iOffset) );
			Sigmoid ( v0, r0, beta, iOffset, numValsPerIter );
					
			IterateSystem5SparseNewNormRK4 ( iter );

		} // for ( long iter = 0; iter < numIterPerTrial; ++iter ) {

		if ( dumpRawToFile ) {
			SaveNumericVectorToBinaryFile(fName, fileOpenFlag, v1E, 0, numValsPerTrial);
			fileOpenFlag = 0;	// This is to make sure all subsequent writes are appends.
 			SaveNumericVectorToBinaryFile(fName, fileOpenFlag, v1I, 0, numValsPerTrial);
			SaveNumericVectorToBinaryFile(fName, fileOpenFlag, v0, 0, numValsPerTrial);
		} // if ( dumpRawToFile ) {
		
	} // for ( iStim = 0; iStim < numStim; ++iStim

		//	Turn OFF weight adaptation.
	plasticityFlag = 0;

		//	Clean up.
	CleanUpDig3RefineInputPatchList ( inputPatches );
	orderInputPatches.clear();

}; // void NM::SyndactylySparse2 ( const string & fBaseName, const long & numCycles, const long & iCycle ) {

	//
	//	SyndactylyControlSparse
	//
void NM::SyndactylyControlSparse () {

	long iStim = 0;
	long numStim = 0;
	vector<InputPatch> inputPatches;
	vector<long> orderInputPatches;

		//	Generate the input patches.
	GenDig3SyndactylyControlInputPatchList  ( inputPatches );

		//	Turn ON weight adaptation.
	plasticityFlag = 1;

		//	Determine a random ordering of the inputs.
	orderInputPatches.clear();
	RandomPermutation ( orderInputPatches, 0, (long) inputPatches.size() );
	numStim = (long) inputPatches.size();

		//	Apply inputs in random order.

	for ( long iStim = 0; iStim < numStim; ++iStim ) {

		for ( long iter = 0; iter < numIterPerTrial; ++iter ) {

			long iOffset = ( ( iter == 0 ) ? ( numIterPerTrialMinusOne ) : ( iter - 1 ) ) * numValsPerIter;

			VecCopy ( inputPatches[ orderInputPatches[iStim] ].inputPattern, iOffset, v0, iOffset, numValsPerIter );
			AddNoise (  v0, iOffset, numValsPerIter, -noiseLevel, noiseLevel );
			Sigmoid ( v0, r0, beta, iOffset, numValsPerIter );

			IterateSystem5Sparse ( iter );

		} // for ( long iter = 0; iter < numIterPerTrial; ++iter ) {
		
	} // for ( iStim = 0; iStim < numStim; ++iStim

		//	Turn OFF weight adaptation.
	plasticityFlag = 0;

		//	Clean up.
	CleanUpDig3SyndactylyInputPatchList  ( inputPatches );
	orderInputPatches.clear();

}; // void NM::SyndactylyControlSparse ()

	//
	//	Spin up the system with a few trial following initial conditions.  No plasticity during this time.
	//
void NM::SystemSpinUpSparse () {

	for ( long iTrials = 0; iTrials < numSpinUpTrials; ++iTrials ) {
			
		for ( long iter = 0; iter < numIterPerTrial; ++iter ) {
			
			IterateSystem5SparseNewNormRK4 ( iter );

		} // for ( long iter = 1; iter < numIterPerTrial; ++iter ) {

	} // for ( long iTrials = 0; iTrials < nm.numSpinUpTrials; ++iTrials ) {

}; // void NM::SystemSpinUpSparse () {

	//
	//	Uniform Refinement Sparse
	//
void NM::UniformRefinementSparse () {

	long iStim = 0;
	long numStim = 0;
	vector<InputPatch> inputPatches;
	vector<long> orderInputPatches;

		//	Generate the input patches.
	GenUniformRefineInputPatchList ( inputPatches );

		//	Turn ON weight adaptation.
	plasticityFlag = 1;

		//	Determine a random ordering of the inputs.
	orderInputPatches.clear();
	RandomPermutation ( orderInputPatches, 0, (long) inputPatches.size() );
	numStim = (long) inputPatches.size();

		//	Apply inputs in random order.

	for ( long iStim = 0; iStim < numStim; ++iStim ) {

		for ( long iter = 0; iter < numIterPerTrial; ++iter ) {

			long iOffset = ( ( iter == 0 ) ? ( numIterPerTrialMinusOne ) : ( iter - 1 ) ) * numValsPerIter;

			VecCopy ( inputPatches[ orderInputPatches[iStim] ].inputPattern, iOffset, v0, iOffset, numValsPerIter );
			AddNoise (  v0, iOffset, numValsPerIter, -noiseLevel, noiseLevel );
			Sigmoid ( v0, r0, beta, iOffset, numValsPerIter );

			IterateSystem5Sparse ( iter );

		} // for ( long iter = 0; iter < numIterPerTrial; ++iter ) {
		
	} // for ( iStim = 0; iStim < numStim; ++iStim

		//	Turn OFF weight adaptation.
	plasticityFlag = 0;

		//	Clean up.
	CleanUpDig3RefineInputPatchList ( inputPatches );
	orderInputPatches.clear();

}; // void NM::UniformRefinementSparse ()

	//
	//	WriteSpareseWeightMatrixToFile
	//		Only the weight values at the last iteration in the trial are written.
	//
void NM::WriteSparseWeightMatrixToFile ( const string & fName, vector<SparseWeightMatrix>::iterator & wPtrStart, vector<SparseWeightMatrix>::iterator & wPtrEnd ) {
	
	fstream outFile;
	outFile.open( fName.c_str(), ios::out | ios::binary );
	for ( vector<SparseWeightMatrix>::iterator wPtr = wPtrStart; wPtr != wPtrEnd; ++wPtr ) {
		wPtr->WriteToFile ( outFile );
	} // for ( vector<SparseWeightMatrix>::iterator wPtr = wPtrStart; wPtr != wPtrEnd; ++wPtr ) {
	outFile.close();

}; // void NM::WriteSparseWeightMatrixToFile ( const string & fName, vector<SparseWeightMatrix>::iterator & wPtr, vector<SparseWeightMatrix>::iterator & wPtrEnd ) { 

	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
	//
	//  END: NM class definitions that leverage sparse matrices.
	//
	//////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////
