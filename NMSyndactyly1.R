
	#
	#	NMSyndactyly1
	#

plasticityFlag = TRUE;

inStim = Dig3MapSyndact0( N, trialDurRFProbeInIters, oneSecondNumIter,
					syndactPatchSize , syndactStimDuration , syndactStimZoneID );
numStim = length ( inStim );
indexStim = sample(1:numStim);	# Want to present the stimulii in random sequence.

if ( verbose ) {
	viewFirstFew = 5;
	viewList = rep ( 0, viewFirstFew );
	for ( i in 1:viewFirstFew ) {
		whichStim = indexStim[i];
		viewList[i] = GetLin( inStim[[whichStim]]$rowPatch, inStim[[whichStim]]$columnPatch, N );
	} # for ( i in 1:viewFirstFew ) {
	viewList = sort ( viewList );
} # if ( verbose )

iStim = 1;

while ( iStim <= numStim ) {	#	Placeholder for an eventual automatic stopping criterion.

	whichStim = indexStim[iStim];
	iter = 2;
	while ( iter <= trialLengthInIters ) {

			#	Load the stimulus pattern.
		k = iter - 1;
		v0[k,] = runif(N2, min = -noiseLevel, max = noiseLevel) + inStim[[whichStim]]$inputPatch[k,];
		r0[k,] = sigmoid ( v0[k,], beta );
		
			#	Do one iteration of the system.

		source("NMDispatch.R");
		if ( wghtNormalizeFlag == 1 ) {
			w1.e.0[,,iter] = NormalizeOutputWeights1 ( w1.e.0[,,iter], n0.numPreSyn, e.wResource ); 
			w1.e.i[,,iter] = NormalizeOutputWeights1 ( w1.e.i[,,iter], n1.e.i.numPreSyn, i.wResource );
			w1.i.e[,,iter] = NormalizeOutputWeights1 ( w1.i.e[,,iter], n1.i.e.numPreSyn, e.wResource );
			w1.e.e[,,iter] = NormalizeOutputWeights1 ( w1.e.e[,,iter], n1.e.e.numPreSyn, e.wResource );
		} else if ( wghtNormalizeFlag == 2.1 ) {
			w1.e.0[,,iter] = NormalizeInputWeights1 ( w1.e.0[,,iter], e.wResource ); 
			w1.e.i[,,iter] = NormalizeInputWeights1 ( w1.e.i[,,iter], e.wResource );
			w1.i.e[,,iter] = NormalizeInputWeights1 ( w1.i.e[,,iter], i.wResource );
			w1.e.e[,,iter] = NormalizeInputWeights1 ( w1.e.e[,,iter], e.wResource );
		} else if ( wghtNormalizeFlag == 2.2 ) {
			w1.e.0[,,iter] = NormalizeInputWeights2 ( w1.e.0[,,iter], e.wResource ); 
			w1.e.i[,,iter] = NormalizeInputWeights2 ( w1.e.i[,,iter], e.wResource );
			w1.i.e[,,iter] = NormalizeInputWeights2 ( w1.i.e[,,iter], i.wResource );
			w1.e.e[,,iter] = NormalizeInputWeights2 ( w1.e.e[,,iter], e.wResource );
		} # if ( wghtNormalizeFlag ) {

			#	Increment iteration.
		iter = iter + 1;

	} # while ( iter <= trialLengthInIters )

		# 	View select time series.
		#	Show an epoch of system evolution for a few sample cells on each layer.
	if ( verbose ) {
		currentCellID = GetLin( inStim[[whichStim]]$rowPatch, inStim[[whichStim]]$columnPatch, N );
		if ( sum ( currentCellID == viewList ) ) {
			nTrials = 1; startOffsetInIter = numIter - nTrials * trialDurRFProbeInIters + 1;
			ShowTimeSeries( v0, r0, v1.e, v1.i, nTrials, trialDurRFProbeInIters, startOffsetInIter, currentCellID, knmSelect, deltaT );
		} # if ( sum ( currentCellID == viewList ) ) {
	} # if ( verbose )

		# Pop the last iteration into the first, including weights.
	v0[1,] = v0[trialLengthInIters,]; r0[1,] = r0[trialLengthInIters,];
	v1.e[1,] = v1.e[trialLengthInIters,]; r1.e[1,] = r1.e[trialLengthInIters,];
	v1.i[1,] = v1.i[trialLengthInIters,]; r1.i[1,] = r1.i[trialLengthInIters,];
	if ( plasticityFlag ) {
		w1.e.0[,,1] = w1.e.0[,,numIter]; w1.e.e[,,1] = w1.e.e[,,numIter];
		w1.e.i[,,1] = w1.e.i[,,numIter]; w1.i.e[,,1] = w1.i.e[,,numIter];
	} # if ( plasticityFlag ) {

	iStim = iStim + 1;

} # while ( iStim <= numStim ) {

#
#	Before exiting be sure to turn off plasticity.
#

plasticityFlag = FALSE;

#
#	Before exiting be sure to turn off plasticity.
#


