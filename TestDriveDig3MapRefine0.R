######################################################################################################
######################################################################################################
#
#	TestDriveDig3MapRefine0: Confirm the operation of the various input stimulus drivers.
#
######################################################################################################
######################################################################################################

testVerbose = TRUE;
if ( testVerbose ) { x11(); }

inStim = Dig3MapRefine0( N, trialLengthInIters, oneSecondNumIter, 3, noiseLevel );

inStim = Dig3MapSelStim0( N, trialDurRFProbeInIters, oneSecondNumIter,
					selectiveStimPatchSize, selectiveStimDuration, selectiveStimFactor, selectiveStimZoneID );

inStim = Dig3MapSyndact0( N, trialDurRFProbeInIters, oneSecondNumIter, syndactPatchSize , syndactStimDuration , syndactStimZoneID );

inStim = Dig3MapAmput0( N, trialDurRFProbeInIters, oneSecondNumIter, refinementAmputPatchSize, 
					refinementAmputStimDuration, kAmputDigit, kAmputZoneMax );
numStim = length(inStim);

for ( i in 1:numStim ) {

	tmp = inStim[[i]];

	if ( testVerbose ) {
		xt = tmp$inputPatch;
		txtTitle = " ";
		par ( mfrow = c ( 2, 2 ) );
		ViewSTDynamics ( xt, ((tmp$offsetPatchInIters)-1), (tmp$offsetPatchInIters), txtTitle );
		ViewSTDynamics ( xt, ((tmp$offsetPatchInIters + tmp$durationPatchInIters) - 1),
			((tmp$offsetPatchInIters + tmp$durationPatchInIters)), txtTitle );
	} # if ( testVerbose ) {

} # for ( i in 1:100 ) {

StimStats ( inStim );