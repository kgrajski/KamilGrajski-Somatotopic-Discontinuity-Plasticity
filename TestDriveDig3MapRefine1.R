######################################################################################################
######################################################################################################
#
#	TestDriveDig3MapRefine1: Generate random refinement inputs for a simulated three digit input layer.
#				No double-digit stimulus.
#
######################################################################################################
######################################################################################################

testVerbose = FALSE;
if ( testVerbose ) { x11(); }

inStim = list();

for ( i in 1:1000 ) {

	tmp = Dig3MapRefine1 ( N, trialLengthInIters, oneSecondNumIter );

	inStim[[i]] = tmp;

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