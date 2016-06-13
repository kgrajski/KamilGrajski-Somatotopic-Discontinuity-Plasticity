	#
	#	PeakStats
	#
FirstPosPeakIfAny = function ( x, iStart, iLen, numPeaks ) {

	#iStart = preStimIters;
	#iLen = stimDurIters;

	thresholdVal = 0.1;
	nups = ndowns = 2;

	pkStats = rep ( 0, 3 );

	tmp = findpeaks ( x[(iStart):(iStart+iLen-1)], npeaks=numPeaks, threshold=thresholdVal, nups=nups, ndowns=ndowns );
	if ( length(tmp) ) {
		pkStats[1] = tmp[1];
		pkStats[2] = tmp[2];
		halfHeight = 0.5 * pkStats[1];
		leftWindowLoc = max ( which( x[(iStart):(iStart+pkStats[2]-1)] <= halfHeight) );	# Left-side of peak.
		rightWindowLoc = iLen - pkStats[2] + 1;
		iTmpWindowLoc = ( which( x[(iStart+pkStats[2]):(iStart+iLen-1)] >= halfHeight) );	# Left-side of peak.
		if ( length ( iTmpWindowLoc ) ) {
			rightWindowLoc = max ( iTmpWindowLoc );
		} # if ( length ( iTmpWindowLoc ) )
		pkStats[3] = ( pkStats[2] - leftWindowLoc + 1 ) + ( rightWindowLoc );
	} # if ( length(tmp) ) 

	return ( pkStats );

} # FirstPosPeakIfAny = function ( x )

FirstNegPeakIfAny = function ( x, iStart, iLen, numPeaks ) {

	#iStart = preStimIters;
	#iLen = stimDurIters;

	x = -x;

	thresholdVal = 1.0;
	nups = ndowns = 5;

	pkStats = rep ( 0, 3 );

	tmp = findpeaks ( x[(iStart):length(x)], npeaks=numPeaks, threshold=thresholdVal, nups=nups, ndowns=ndowns );
	if ( length(tmp) ) {
		pkStats[1] = tmp[1];
		pkStats[2] = tmp[2];
		halfHeight = 0.5 * pkStats[1];
		leftWindowLoc = max ( which( x[(iStart):(iStart+pkStats[2]-1)] <= halfHeight) );	# Left-side of peak.
		rightWindowLoc = iLen - pkStats[2] + 1;
		iTmpWindowLoc = ( which( x[(iStart+pkStats[2]):(iStart+iLen-1)] >= halfHeight) );	# Left-side of peak.
		if ( length ( iTmpWindowLoc ) ) {
			rightWindowLoc = max ( iTmpWindowLoc );
		} # if ( length ( iTmpWindowLoc ) )
		pkStats[3] = ( pkStats[2] - leftWindowLoc + 1 ) + ( rightWindowLoc );
	} # if ( length(tmp) ) 
	pkStats[1] = - pkStats[1];

	return ( pkStats );

} # FirstNegPeakIfAny = function ( x )


	#
	#
	#
TSSTStatsDelta = function( base, exp ) {

	pk.e = matrix ( 0, nrow=(dim(base$pk.e)[1]), ncol=6 );
	pk.i = pk.e;

	for ( i in 1:6 ) {
		pk.e[,i] = exp$pk.e[,i] - base$pk.e[,i];
		pk.i[,i] = exp$pk.i[,i] - base$pk.i[,i];
	} # for ( i in 1:6 )
	
	return ( list ( pk.e=pk.e, pk.i=pk.i ) );
	
} # TSSTStatsDelta = function( base, exp ) {

	#
	#	TS = Time Series 
	#	ST = Spatio-temporal
	#
TSSTStatsDriver = function ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp, 
					titleTextRoot, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fName,
					xKnockBox, yKnockBox, rowSeedShowMap, colSeedShowMap, subNShowMap ) {

		#	xKnockBox, yKnockBox used to generate outline of the knockout zone.
	xKnockBox[1] = xKnockBox[1] - 0.5; xKnockBox[2] = xKnockBox[2] + 0.5;
	yKnockBox[1] = yKnockBox[1] - 0.5; yKnockBox[2] = yKnockBox[2] + 0.5;
	x.KnockBox = rbind(c(xKnockBox[1], xKnockBox[2]), c(xKnockBox[1], xKnockBox[2]),
					c(xKnockBox[1], xKnockBox[1]), c(xKnockBox[2], xKnockBox[2]));
	y.KnockBox = rbind( c ( yKnockBox[1], yKnockBox[1]), c(yKnockBox[2], yKnockBox[2]),
					c(yKnockBox[1], yKnockBox[2]), c(yKnockBox[1], yKnockBox[2]) );


	stats.base = TSSTStats ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base );
	stats.exp = TSSTStats ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.exp );
	stats.delta = TSSTStatsDelta ( stats.base, stats.exp );

	xLabText  = "Distal -> Proximal";
	yLabText = "Digit 1 -> Digit 3";
	boundaryMarks = c ( as.integer(sqrt(N2)/3)+0.5, as.integer(sqrt(N2)/3)*2+0.5 );
	skipIPlots = FALSE;

	for ( iStatType in 1:3 ) {

		if ( iStatType == 1 ) {
			tagText = "Init Pos Pk Magnitude";
		} else if ( iStatType == 2 ) {
			tagText = "Init Pos Pk Latency";
		} else if ( iStatType == 3 ) {
			tagText = "Init Pos Pk Half Width";
		} else {
			tagText = NULL;
		} # if ( iStatType == 1 )

		if ( 0 ) {
	
			if ( tiffFlag ) {
				tiff ( paste("PkAz", iStatType, fName, "tiff", sep="." ), compression="lzw", units="in", width=7.0, height=7.0, res=300 );
			} else {
				x11();
			} # if ( tiffFlag )

			par(mfrow=c(3,2) );

			titleText = paste ( "E Cell", tagText, "\n", titleTextRoot, "Stim: ", iProbeCell );	
			ShowVecAsMap2 ( stats.exp$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.exp$pk.e[,iStatType ]), max(stats.exp$pk.e[,iStatType ]) );
			abline ( h = boundaryMarks, lty=3, col=3 );
			GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

			titleText = paste ( "I Cell", tagText, "\n", titleTextRoot, "Stim: ", iProbeCell );	
			ShowVecAsMap2 ( stats.exp$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.exp$pk.i[,iStatType ]), max(stats.exp$pk.i[,iStatType ]) );
			abline ( h = boundaryMarks, lty=3, col=3 );
			GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

			titleTextRoot.base = "Baseline Refined Network";
			titleText = paste ( "E Cell", tagText, "\n", titleTextRoot.base, "Stim: ", iProbeCell );	
			ShowVecAsMap2 ( stats.base$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.base$pk.e[,iStatType ]), max(stats.base$pk.e[,iStatType ]) );
			abline ( h = boundaryMarks, lty=3, col=3 );
			GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

			titleText = paste ( "I Cell", tagText, "\n", titleTextRoot.base, "Stim: ", iProbeCell );	
			ShowVecAsMap2 ( stats.base$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.base$pk.i[,iStatType ]), max(stats.base$pk.i[,iStatType ]) );
			abline ( h = boundaryMarks, lty=3, col=3 );
			GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

			titleTextRoot.diff = "Diff Exp - Baseline";
			titleText = paste ( "E Cell", tagText, "\n", titleTextRoot.diff, "Stim: ", iProbeCell );	
			ShowVecAsMap2 ( stats.delta$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.delta$pk.e[,iStatType ]), max(stats.delta$pk.e[,iStatType ]) );
			abline ( h = boundaryMarks, lty=3, col=3 );
			GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

			titleText = paste ( "I Cell", tagText, "\n", titleTextRoot.diff, "Stim: ", iProbeCell );	
			ShowVecAsMap2 ( stats.delta$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.delta$pk.i[,iStatType ]), max(stats.delta$pk.i[,iStatType ]) );
			abline ( h = boundaryMarks, lty=3, col=3 );
			GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

			if ( tiffFlag ) {
				dev.off();
			} # if ( tiffFlag )

		} # if ( 0 )

		if ( tiffFlag ) {
			tiff ( paste("PkAz.Zoom", iStatType, fName, "tiff", sep="." ), compression="lzw", units="in", width=7.0, height=7.0, res=300 );
		} else {
			x11();
		} # if ( tiffFlag )

		par(mfrow=c(3,2) );
		
		titleText = paste ( "E Cell", tagText, "\n", titleTextRoot, "Stim: ", iProbeCell );	
		ShowVecAsMap2W ( stats.exp$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.exp$pk.e[,iStatType ]), max(stats.exp$pk.e[,iStatType ]),
						rowSeedShowMap, colSeedShowMap, subNShowMap );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleText = paste ( "I Cell", tagText, "\n", titleTextRoot, "Stim: ", iProbeCell );	
		ShowVecAsMap2W ( stats.exp$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.exp$pk.i[,iStatType ]), max(stats.exp$pk.i[,iStatType ]),
						rowSeedShowMap, colSeedShowMap, subNShowMap );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleTextRoot.base = "Baseline Refined Network";
		titleText = paste ( "E Cell", tagText, "\n", titleTextRoot.base, "Stim: ", iProbeCell );	
		ShowVecAsMap2W ( stats.base$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.base$pk.e[,iStatType ]), max(stats.base$pk.e[,iStatType ]),
						rowSeedShowMap, colSeedShowMap, subNShowMap );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleText = paste ( "I Cell", tagText, "\n", titleTextRoot.base, "Stim: ", iProbeCell );	
		ShowVecAsMap2W ( stats.base$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.base$pk.i[,iStatType ]), max(stats.base$pk.i[,iStatType ]),
						rowSeedShowMap, colSeedShowMap, subNShowMap );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleTextRoot.diff = "Diff Exp - Baseline";
		titleText = paste ( "E Cell", tagText, "\n", titleTextRoot.diff, "Stim: ", iProbeCell );	
		ShowVecAsMap2W ( stats.delta$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.delta$pk.e[,iStatType ]), max(stats.delta$pk.e[,iStatType ]),
						rowSeedShowMap, colSeedShowMap, subNShowMap );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleText = paste ( "I Cell", tagText, "\n", titleTextRoot.diff, "Stim: ", iProbeCell );	
		ShowVecAsMap2W ( stats.delta$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.delta$pk.i[,iStatType ]), max(stats.delta$pk.i[,iStatType ]),
						rowSeedShowMap, colSeedShowMap, subNShowMap );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		if ( tiffFlag ) {
			dev.off();
		} # if ( tiffFlag )

	} # 	for ( iStatType in 1:3 ) {

	for ( iStatType in 4:6 ) {

		if ( iStatType == 4 ) {
			tagText = "Neg Pk Magnitude";
			skipIPlots = TRUE;
		} else if ( iStatType == 5 ) {
			tagText = "Neg Pk Latency";
			skipIPlots = TRUE;
		} else if ( iStatType == 6 ) {
			tagText = "Neg Pk Half Width";
			skipIPlots = TRUE;
		} else {
			tagText = NULL;
		} # if ( iStatType == 1 )

		if ( 0 ) {

			if ( tiffFlag ) {
				tiff ( paste("PkAz", iStatType, fName, "tiff", sep="." ), compression="lzw", units="in", width=7.0, height=7.0, res=300 );
			} else {
				x11();
			} # if ( tiffFlag )

			par(mfcol=c(3,2) );
		
			titleText = paste ( "E Cell", tagText, "\n", titleTextRoot, "Stim: ", iProbeCell );	
			ShowVecAsMap2 ( stats.exp$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.exp$pk.e[,iStatType ]), max(stats.exp$pk.e[,iStatType ]) );
			abline ( h = boundaryMarks, lty=3, col=3 );
			GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

			titleText = paste ( "I Cell", tagText, "\n", titleTextRoot, "Stim: ", iProbeCell );	
			ShowVecAsMap2 ( stats.exp$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.exp$pk.i[,iStatType ]), max(stats.exp$pk.i[,iStatType ]) );
			abline ( h = boundaryMarks, lty=3, col=3 );
			GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

			titleTextRoot.base = "Baseline Refined Network";
			titleText = paste ( "E Cell", tagText, "\n", titleTextRoot.base, "Stim: ", iProbeCell );	
			ShowVecAsMap2 ( stats.base$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.base$pk.e[,iStatType ]), max(stats.base$pk.e[,iStatType ]) );
			abline ( h = boundaryMarks, lty=3, col=3 );
			GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

			titleTextRoot.base = "Baseline Refined Network";
			titleText = paste ( "I Cell", tagText, "\n", titleTextRoot.base, "Stim: ", iProbeCell );	
			ShowVecAsMap2 ( stats.base$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.base$pk.i[,iStatType ]), max(stats.base$pk.i[,iStatType ]) );
			abline ( h = boundaryMarks, lty=3, col=3 );
			GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

			titleTextRoot.diff = "Diff Exp - Baseline";
			titleText = paste ( "E Cell", tagText, "\n", titleTextRoot.diff, "Stim: ", iProbeCell );	
			ShowVecAsMap2 ( stats.delta$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.delta$pk.e[,iStatType ]), max(stats.delta$pk.e[,iStatType ]) );
			abline ( h = boundaryMarks, lty=3, col=3 );
			GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

			titleTextRoot.diff = "Diff Exp - Baseline";
			titleText = paste ( "I Cell", tagText, "\n", titleTextRoot.diff, "Stim: ", iProbeCell );	
			ShowVecAsMap2 ( stats.delta$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.delta$pk.i[,iStatType ]), max(stats.delta$pk.i[,iStatType ]) );
			abline ( h = boundaryMarks, lty=3, col=3 );
			GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

			if ( tiffFlag ) {
				dev.off();
			} # if ( tiffFlag )

		} # if ( 0 )

		if ( tiffFlag ) {
			tiff ( paste("PkAz.Zoom", iStatType, fName, "tiff", sep="." ), compression="lzw", units="in", width=7.0, height=7.0, res=300 );
		} else {
			x11();
		} # if ( tiffFlag )

		par(mfcol=c(3,2) );
		
		titleText = paste ( "E Cell", tagText, "\n", titleTextRoot, "Stim: ", iProbeCell );	
		ShowVecAsMap2W ( stats.exp$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.exp$pk.e[,iStatType ]), max(stats.exp$pk.e[,iStatType ]),
					rowSeedShowMap, colSeedShowMap, subNShowMap );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleText = paste ( "I Cell", tagText, "\n", titleTextRoot, "Stim: ", iProbeCell );	
		ShowVecAsMap2W ( stats.exp$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.exp$pk.i[,iStatType ]), max(stats.exp$pk.i[,iStatType ]),
					rowSeedShowMap, colSeedShowMap, subNShowMap );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleTextRoot.base = "Baseline Refined Network";
		titleText = paste ( "E Cell", tagText, "\n", titleTextRoot.base, "Stim: ", iProbeCell );	
		ShowVecAsMap2W ( stats.base$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.base$pk.e[,iStatType ]), max(stats.base$pk.e[,iStatType ]),
					rowSeedShowMap, colSeedShowMap, subNShowMap );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleTextRoot.base = "Baseline Refined Network";
		titleText = paste ( "I Cell", tagText, "\n", titleTextRoot.base, "Stim: ", iProbeCell );	
		ShowVecAsMap2W ( stats.base$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.base$pk.i[,iStatType ]), max(stats.base$pk.i[,iStatType ]),
					rowSeedShowMap, colSeedShowMap, subNShowMap );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleTextRoot.diff = "Diff Exp - Baseline";
		titleText = paste ( "E Cell", tagText, "\n", titleTextRoot.diff, "Stim: ", iProbeCell );	
		ShowVecAsMap2W ( stats.delta$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.delta$pk.e[,iStatType ]), max(stats.delta$pk.e[,iStatType ]),
					rowSeedShowMap, colSeedShowMap, subNShowMap );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleTextRoot.diff = "Diff Exp - Baseline";
		titleText = paste ( "I Cell", tagText, "\n", titleTextRoot.diff, "Stim: ", iProbeCell );	
		ShowVecAsMap2W ( stats.delta$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.delta$pk.i[,iStatType ]), max(stats.delta$pk.i[,iStatType ]),
					rowSeedShowMap, colSeedShowMap, subNShowMap );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		if ( tiffFlag ) {
			dev.off();
		} # if ( tiffFlag )

	} # 	for ( iStatType in 4:6 ) {

} # MPSTStats = function ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,

	#
	#	Compute various statistics on the raw time series
	#
TSSTStats = function ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, rData ) {

	numPeaks = 2;
	peakStats.e = matrix ( 0, nrow=N2, ncol=3*numPeaks );	# Col 1, 4: peakMagnitude; Col 2, 5: peak latency; Col 3, 6: 1/2 peak width.
	peakStats.i = matrix ( 0, nrow=N2, ncol=3*numPeaks );

	tsTmp.e = rep ( 0, numItersPerTrial);
	tsTmp.i = rep ( 0, numItersPerTrial);

	iPtrList.e = seq ( 1, N2*numItersPerTrial, N2 );
	iPtrList.i = numValsPerRFTrial + iPtrList.e;

	for ( iCell in 1:N2 ) {

				#	Using knowledge of data structure details to speed things up.
		tsTmp.e = ( rData [iPtrList.e] );
		tsTmp.i = ( rData [numValsPerRFTrial+iPtrList.e] );
		iPtrList.e = iPtrList.e + 1;

		peakStats.e[iCell,1:3] = FirstPosPeakIfAny ( tsTmp.e, preStimIters, stimDurIters, numPeaks );
		peakStats.i[iCell,1:3] = FirstPosPeakIfAny ( tsTmp.i, preStimIters, stimDurIters, numPeaks );

		peakStats.e[iCell,4:6] = FirstNegPeakIfAny ( tsTmp.e, preStimIters, stimDurIters, numPeaks );
		peakStats.i[iCell,4:6] = FirstNegPeakIfAny ( tsTmp.i, preStimIters, stimDurIters, numPeaks );

		#print ( iCell ); print ( warnings() );

	} # for ( iCell in 1:N2 )

	return ( list ( pk.e=peakStats.e, pk.i=peakStats.i ) );

} # TSSTStats = function ( N2, numValsPerRFTrial, iNumItersPerTrial, preStimIters, stimDurIters, rData ) {




