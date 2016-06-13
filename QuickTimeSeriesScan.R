#	Script to generate RF mapping inputs and see the time series response.

	for ( iCell in seq(1,N2,4) ) {	# for each cell on the input layer.

		v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
		v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
		v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;

		iter = 2;
		while ( iter <= trialDurRFProbeInIters ) {

				#	Stimulus ON/OFF.
				#		As long as the stimulus is of amplitude 1.0 in one unit, no need to call VectorNormalize.
			if ( iter >= (offsetToRFProbeInIters + 2) & iter < (timeOffRFProbeInIters + 2) ) {
				v0[iter-1,iCell] = magRFProbe + runif(1, min = -noiseLevel, max = noiseLevel);
				r0[iter-1,] = sigmoid ( v0[iter-1,], beta );
			} else if ( iter == (timeOffRFProbeInIters + 2) ) {
				v0[iter-1,iCell] = runif(1, min = -noiseLevel, max = noiseLevel);
				r0[iter-1,] = sigmoid ( v0[iter-1,], beta );
			} # if ( iter == ... )

				#	Do one iteration of the system.

			k = iter - 1;
			source("NMDispatch.R");

				#	Increment iteration.

			iter = iter + 1;

		} # while ( iter <= expLengthNumIter )

		ShowTimeSeries( v0, r0, v1.e, v1.i, nTrials, trialDurRFProbeInIters, startOffsetInIter, iCell, knmSelect, deltaT );

	} # for ( iCell = 1:N2 )
