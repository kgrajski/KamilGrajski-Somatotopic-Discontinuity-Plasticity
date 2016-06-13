	##############################################################
	##############################################################
	#
	#	RFMapExp1
	#
	#	Receptive Field Mapping
	#
	#	Given a simple two layer network.  Determine the Receptive
	#	Field for each Excitatory and Inhibitory cell in the Layer 1.
	#	For each of the input layer nodes, stimulate it and store
	#	the response of all of the cortical cells.  Post process
	#	the batch to extract the receptive field map for each
	#	cortical cell.
	#
	#	The stimulus protocol is to have a few second long
	#	'trial.'  For the first second or so just iterate the
	#	system so that a background activity level can be
	#	estimatd.  Then apply (as a ON step function) a x-second long
	#	pulse to the input layer "cell" and convert that to a firing
	#	rate.  The stimulus OFF step is exponential with time
	#	constant v0.tau (pretty fast).
	#
	#	The post processing compares the background level of
	#	activity with the activity level during and shortly after
	#	the stimulus goes off.
	#
	#	Generate map views and histograms.
	#
	##############################################################
	##############################################################

		# This will store the raw experimental data from which RF maps will be derived.
	rfMapRaw = list();
	for ( iCell in 1:N2 ) {	# for each cell on the input layer.

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
		iter = iter - 1;
		rfMapRaw[[iCell]] = list ( v0=v0[1:iter,], v1.e=v1.e[1:iter,], v1.i=v1.i[1:iter,], r0=r0[1:iter,], r1.e=r1.e[1:iter,], r1.i=r1.i[1:iter,] );

		if ( (iCell == ((trunc((N^2)/2))+1)) | (iCell == 1) ) {	
			nTrials = 1; startOffsetInIter = 1;
			if ( verbose ) {				# For the cell in the middle of the input layer, show the response of all cortical cells.
				for ( jCell in trackCells ) {
					ShowTimeSeries( v0, r0, v1.e, v1.i, nTrials, trialDurRFProbeInIters, startOffsetInIter, jCell, knmSelect, deltaT );
				} # for ( jCell in c(1, N2) ) {
			} else {
				ShowTimeSeries( v0, r0, v1.e, v1.i, nTrials, trialDurRFProbeInIters, startOffsetInIter, iCell, knmSelect, deltaT );
			} # if ( verbose )	
		} # if ( (iCell == ((trunc((N^2)/2))+1)) | (iCell == 1) | (iCell == N2) {

	} # for ( iCell = 1:N2 )

		#	Candidate for a future routine.
		#	ith row of rfMap gives the Receptive Field map for the ith cortical cell.
		#	Show the results as spatial plots.
	rfMap = GenRFMap1 ( rfMapRaw, oneSecondNumIter, offsetToRFProbeInIters, durRFProbeInIters );
	rfMap.rfAreas.e = QuantRFSize ( rfMap$r1.e.rfMap, kRFPeakToEdgeDetect );
	rfMap.rfAreas.i = QuantRFSize ( rfMap$r1.i.rfMap, kRFPeakToEdgeDetect );
	rfMap.maxResp.e = apply ( rfMap$r1.e.rfMap, 1, max );
	rfMap.maxResp.i = apply ( rfMap$r1.i.rfMap, 1, max );
	cortical.Amp.e = QuantCorticalAmp ( rfMap$r1.e.rfMap, kRFPeakToEdgeDetect );
	cortical.Amp.i = QuantCorticalAmp ( rfMap$r1.i.rfMap, kRFPeakToEdgeDetect );

	if ( verbose ) {

			#
			# Plot RF Maps
			#
			
		for ( itmp in trackCells ) {
			x11(); par(mfrow=c(2,1)); 
			txtTitle = paste("L1_E Exp", knmSelect); ViewSTDynamics ( rfMap$r1.e.rfMap, itmp, itmp, txtTitle );
			txtTitle = paste("L1_I Exp", knmSelect); ViewSTDynamics ( rfMap$r1.i.rfMap, itmp, itmp, txtTitle );
		} # for ( itmp in trackCells )

			# Plot RF size histogram(s) with fixed breaks and with default breaks for the E and I cells.
		x11(); par(mfrow=c(2,2));
		hist.rfMap.rfAreas.e = hist ( rfMap.rfAreas.e, breaks = 1:N2, plot=TRUE,
			main=paste("Receptive Field (RF) Sizes\nCortical Layer E-Cells", paste("Exp ",knmSelect), sep="\n"), ylab="Count", xlab="RF Size" );
		hist.rfMap.rfAreas.i = hist ( rfMap.rfAreas.i, breaks = 1:N2, plot=TRUE,
			main=paste("Receptive Field (RF) Sizes\nCortical Layer I-Cells", paste("Exp ",knmSelect), sep="\n"), ylab="Count", xlab="RF Size" );
		hist.rfMap.maxResp.e = hist ( rfMap.maxResp.e, plot=TRUE,
			main=paste("RF Response Magnitude\nCortical Layer E-Cells", paste("Exp ",knmSelect), sep="\n"), ylab="Count", xlab="RF Response" );
		hist.rfMap.maxResp.i = hist ( rfMap.maxResp.i, plot=TRUE,
			main=paste("RF Response Magnitude\nCortical Layer I-Cells", paste("Exp ",knmSelect), sep="\n"), ylab="Count", xlab="RF Response" );
		
		x11(); par(mfrow=c(2,2));
		hist.rfMap.rfAreas.e = hist ( rfMap.rfAreas.e, plot=TRUE,
			main=paste("Receptive Field (RF) Sizes\nCortical Layer E-Cells", paste("Exp ",knmSelect), sep="\n"), ylab="Count", xlab="RF Size" );
		hist.rfMap.rfAreas.i = hist ( rfMap.rfAreas.i, plot=TRUE,
			main=paste("Receptive Field (RF) Sizes\nCortical Layer I-Cells", paste("Exp ",knmSelect), sep="\n"), ylab="Count", xlab="RF Size" );
		hist.rfMap.maxResp.e = hist ( rfMap.maxResp.e, plot=TRUE,
			main=paste("RF Response Magnitude\nCortical Layer E-Cells", paste("Exp ",knmSelect), sep="\n"), ylab="Count", xlab="RF Response" );
		hist.rfMap.maxResp.i = hist ( rfMap.maxResp.i, plot=TRUE,
			main=paste("RF Response Magnitude\nCortical Layer I-Cells", paste("Exp ",knmSelect), sep="\n"), ylab="Count", xlab="RF Response" );
		
		x11(); par(mfrow=c(2,1));
		plot ( rfMap.maxResp.e, rfMap.rfAreas.e,
			main=paste("RF Size vs RF Response Magnitude\nCortical Layer E-Cells", paste("Exp ",knmSelect), sep="\n"),
			ylab="RF Size", xlab="RF Response Magnitude" );
		plot ( rfMap.maxResp.i, rfMap.rfAreas.i,
			main=paste("RF Size vs RF Response Magnitude\nCortical Layer I-Cells", paste("Exp ",knmSelect), sep="\n"),
			ylab="RF Size", xlab="RF Response Magnitude" );

		x11(); par(mfrow=c(2,2) );
		ShowVecAsMap ( rfMap.rfAreas.e, "rfMap.rfAreas.e" );
		ShowVecAsMap ( rfMap.rfAreas.i, "rfMap.rfAreas.i" );
		ShowVecAsMap ( rfMap.maxResp.e, "rfMap.maxResp.e" );
		ShowVecAsMap ( rfMap.maxResp.i, "rfMap.maxResp.i" );

			#	Do some cortical amplications plots
		ThreeDigitCorticalAmplification ( cortical.Amp.e, TRUE, "Magnification (E-Cells)" );
		ThreeDigitCorticalAmplification ( cortical.Amp.i, TRUE, "Magnification (I-Cells)" );

	} # if ( verbose )

