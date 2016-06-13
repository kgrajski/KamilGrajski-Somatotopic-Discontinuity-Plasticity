#
#	KnockoutFigure6.R
#
#
#	Automatically handling the control, one-side and two-side I-cell knockout variations.
#
#	KnockoutFigure6 - strictly additive to KnockoutFigure5 tailored for more
#				consolidated figures.

#	Clear the workspace.
rm(list = ls());

#	Load up some libraries and helper functions.

source ( "NMHelperFunctions.R" );

	#
	#	Additional specialized helper functions.
	#

	#
	#	Standardized plots.
	#
PlotKnockoutRFMapConsolidated = function ( ref, irefIter, exp, iexpIter, exp2, iexpIter2, exp3, iexpIter3,
						refTitleText, expTitleText, expTitleText2, expTitleText3, activationValue,
						iTrim, iFilter.E, iFilter.I, iFilter.E.2, iFilter.I.2, iFilter.E.3, iFilter.I.3,
						iRowList, iColList, x, y, tiffFlag=FALSE, iKnockOutLength=0, iExpFlag=0 ) {

		#
		#	Set up the parameters so that the "experimental zone" will be outlined on subsequent plots.
		#
	y[1] = y[1] - 0.5; y[2] = y[2] + 0.5;
	x[1] = x[1] - 0.5; x[2] = x[2];
	tmp.x = rbind(c(x[1], x[2]), c(x[1], x[2]), c(x[1], x[1]), c(x[2], x[2]));
	tmp.y = rbind(c(y[1], y[1]), c(y[2], y[2]), c(y[1], y[2]), c(y[1], y[2]));

		#	Titles and digit boundaries.
	xLabText  = "Distal -> Proximal";  yLabText = "Digit 1 -> Digit 3";
	N = as.integer ( sqrt ( dim(base$r1.i.rfMap)[1] ) );
	boundaryMarks = c ( as.integer(N/3)+0.5, as.integer(N/3)*2+0.5 );
	iRFTrackStepSize = 2;

		#
		#	RF Centroids
		#
	if ( tiffFlag ) {
		tiff ( paste("Consol.Knock", activationValue, "Centroids.Full", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )
	par ( mfrow=c(2,4) );
	ShowTopoMap1 ( ref$rfTrackData.e[[length(ref$rfTrackData.e)]], paste(paste("E-Type","Iter",irefIter,sep=" "),refTitleText, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1 ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1 ( exp2$rfTrackData.e[[length(exp2$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter2,sep=" "),expTitleText2, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1 ( exp3$rfTrackData.e[[length(exp3$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter3,sep=" "),expTitleText3, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1 ( ref$rfTrackData.i[[length(ref$rfTrackData.i )]], paste(paste("I-Type","Iter",irefIter,sep=" "),refTitleText, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1 ( exp$rfTrackData.i[[length(exp$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1 ( exp2$rfTrackData.i[[length(exp2$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter2,sep=" "),expTitleText2, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1 ( exp3$rfTrackData.i[[length(exp3$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter3,sep=" "),expTitleText3, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

	if ( tiffFlag ) {
		tiff ( paste("Consol.Knock", activationValue, "Centroids.W", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )
	par ( mfrow=c(2,4) );
	ShowTopoMap1W ( ref$rfTrackData.e[[length(ref$rfTrackData.e)]], paste(paste("E-Type","Iter",irefIter,sep=" "),refTitleText, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1W ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1W ( exp2$rfTrackData.e[[length(exp2$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter2,sep=" "),expTitleText2, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1W ( exp3$rfTrackData.e[[length(exp3$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter3,sep=" "),expTitleText3, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1W ( ref$rfTrackData.i[[length(ref$rfTrackData.i )]], paste(paste("I-Type","Iter",irefIter,sep=" "),refTitleText, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1W ( exp$rfTrackData.i[[length(exp$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1W ( exp2$rfTrackData.i[[length(exp2$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter2,sep=" "),expTitleText2, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1W ( exp3$rfTrackData.i[[length(exp3$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter3,sep=" "),expTitleText3, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

		#
		#	RF Tracks
		#
	if ( tiffFlag ) {
		tiff ( paste("Consol.Knock", activationValue, "Tracks", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )
	par ( mfrow=c(2,4) );
	ShowThreeDigitRFTrackFlexW ( ref$rfTrackData.e[[length(ref$rfTrackData.e)]], paste(paste("E-Type","Iter",irefIter,sep=" "),refTitleText, sep=""),
		TRUE, 0.5, 0, 1, iRowList, iColList, iRFTrackStepSize );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );

	ShowThreeDigitRFTrackFlexW ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""),
		TRUE, 0.5, 0, 1, iRowList, iColList, iRFTrackStepSize );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );

	ShowThreeDigitRFTrackFlexW ( exp2$rfTrackData.e[[length(exp2$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter2,sep=" "),expTitleText2, sep=""),
		TRUE, 0.5, 0, 1, iRowList, iColList, iRFTrackStepSize );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );

	ShowThreeDigitRFTrackFlexW ( exp3$rfTrackData.e[[length(exp3$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter3,sep=" "),expTitleText3, sep=""),
		TRUE, 0.5, 0, 1, iRowList, iColList, iRFTrackStepSize );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );

	ShowThreeDigitRFTrackFlexW ( ref$rfTrackData.i[[length(ref$rfTrackData.i)]], paste(paste("I-Type","Iter",irefIter,sep=" "),refTitleText, sep=""),
		TRUE, 0.5, 0, 1, iRowList, iColList, iRFTrackStepSize );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );

	ShowThreeDigitRFTrackFlexW ( exp$rfTrackData.i[[length(exp$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""),
		TRUE, 0.5, 0, 1, iRowList, iColList, iRFTrackStepSize );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );

	ShowThreeDigitRFTrackFlexW ( exp2$rfTrackData.i[[length(exp2$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter2,sep=" "),expTitleText2, sep=""),
		TRUE, 0.5, 0, 1, iRowList, iColList, iRFTrackStepSize );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );

	ShowThreeDigitRFTrackFlexW ( exp3$rfTrackData.i[[length(exp3$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter3,sep=" "),expTitleText3, sep=""),
		TRUE, 0.5, 0, 1, iRowList, iColList, iRFTrackStepSize );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

		#
		#	RF Translocation
		#
	if ( tiffFlag ) {
		tiff ( paste("Consol.Knock", activationValue, "Position_Shift", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )

	par ( mfrow=c(2,3) );
	centroidDelta.e = RFCentroidDelta ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], ref$rfTrackData.e[[length(ref$rfTrackData.e)]]);
	centroidDelta.i = RFCentroidDelta ( exp$rfTrackData.i[[length(exp$rfTrackData.i)]], ref$rfTrackData.i[[length(ref$rfTrackData.i)]]);
	centroidDelta.e[c(iTrim,iFilter.E)] = 0; centroidDelta.i[c(iTrim,iFilter.I)] = 0;

	centroidDelta.e.2 = RFCentroidDelta ( exp2$rfTrackData.e[[length(exp2$rfTrackData.e)]], ref$rfTrackData.e[[length(ref$rfTrackData.e)]]);
	centroidDelta.i.2 = RFCentroidDelta ( exp2$rfTrackData.i[[length(exp2$rfTrackData.i)]], ref$rfTrackData.i[[length(ref$rfTrackData.i)]]);
	centroidDelta.e.2[c(iTrim,iFilter.E.2)] = 0; centroidDelta.i.2[c(iTrim,iFilter.I.2)] = 0;

	centroidDelta.e.3 = RFCentroidDelta ( exp3$rfTrackData.e[[length(exp3$rfTrackData.e)]], ref$rfTrackData.e[[length(ref$rfTrackData.e)]]);
	centroidDelta.i.3 = RFCentroidDelta ( exp3$rfTrackData.i[[length(exp3$rfTrackData.i)]], ref$rfTrackData.i[[length(ref$rfTrackData.i)]]);
	centroidDelta.e.3[c(iTrim,iFilter.E.3)] = 0; centroidDelta.i.3[c(iTrim,iFilter.I.3)] = 0;
	
	ShowVecAsMap2 ( centroidDelta.e, paste(paste("E-Type RF Centroid Shift","Iter",iexpIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(centroidDelta.e), max(centroidDelta.e) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( centroidDelta.e.2, paste(paste("E-Type RF Centroid Shift","Iter",iexpIter2,sep=" "),expTitleText2, sep=""), xLabText, yLabText, min(centroidDelta.e.2), max(centroidDelta.e.2) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( centroidDelta.e.3, paste(paste("E-Type RF Centroid Shift","Iter",iexpIter3,sep=" "),expTitleText3, sep=""), xLabText, yLabText, min(centroidDelta.e.3), max(centroidDelta.e.3) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( centroidDelta.i, paste(paste("I-Type RF Centroid Shift","Iter",iexpIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(centroidDelta.i), max(centroidDelta.i) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( centroidDelta.i.2, paste(paste("I-Type RF Centroid Shift","Iter",iexpIter2,sep=" "),expTitleText2, sep=""), xLabText, yLabText, min(centroidDelta.i.2), max(centroidDelta.i.2) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( centroidDelta.i.3, paste(paste("I-Type RF Centroid Shift","Iter",iexpIter3,sep=" "),expTitleText3, sep=""), xLabText, yLabText, min(centroidDelta.i.3), max(centroidDelta.i.3) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );	

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

		#
		#	RF Translocation (with log transform for those cases (e.g., knock-in) where there is lots of movement
		#
	if ( tiffFlag ) {
		tiff ( paste("Consol.Knock", activationValue, "Position_Shift_Log10", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )

	par ( mfrow=c(2,3) );
	centroidDelta.e = log10(1+RFCentroidDelta ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], ref$rfTrackData.e[[length(ref$rfTrackData.e)]]));
	centroidDelta.i = log10(1+RFCentroidDelta ( exp$rfTrackData.i[[length(exp$rfTrackData.i)]], ref$rfTrackData.i[[length(ref$rfTrackData.i)]]));
	centroidDelta.e[c(iTrim,iFilter.E)] = 0; centroidDelta.i[c(iTrim,iFilter.I)] = 0;

	centroidDelta.e.2 = log10(1+RFCentroidDelta ( exp2$rfTrackData.e[[length(exp2$rfTrackData.e)]], ref$rfTrackData.e[[length(ref$rfTrackData.e)]]));
	centroidDelta.i.2 = log10(1+RFCentroidDelta ( exp2$rfTrackData.i[[length(exp2$rfTrackData.i)]], ref$rfTrackData.i[[length(ref$rfTrackData.i)]]));
	centroidDelta.e.2[c(iTrim,iFilter.E.2)] = 0; centroidDelta.i.2[c(iTrim,iFilter.I.2)] = 0;

	centroidDelta.e.3 = log10(1+RFCentroidDelta ( exp3$rfTrackData.e[[length(exp3$rfTrackData.e)]], ref$rfTrackData.e[[length(ref$rfTrackData.e)]]));
	centroidDelta.i.3 = log10(1+RFCentroidDelta ( exp3$rfTrackData.i[[length(exp3$rfTrackData.i)]], ref$rfTrackData.i[[length(ref$rfTrackData.i)]]));
	centroidDelta.e.3[c(iTrim,iFilter.E.3)] = 0; centroidDelta.i.3[c(iTrim,iFilter.I.3)] = 0;
	
	ShowVecAsMap2 ( centroidDelta.e, paste(paste("E-Type RF Centroid Log10 Shift","Iter",iexpIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(centroidDelta.e), max(centroidDelta.e) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( centroidDelta.e.2, paste(paste("E-Type RF Centroid Log10 Shift","Iter",iexpIter2,sep=" "),expTitleText2, sep=""), xLabText, yLabText, min(centroidDelta.e.2), max(centroidDelta.e.2) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( centroidDelta.e.3, paste(paste("E-Type RF Centroid Log10 Shift","Iter",iexpIter3,sep=" "),expTitleText3, sep=""), xLabText, yLabText, min(centroidDelta.e.3), max(centroidDelta.e.3) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( centroidDelta.i, paste(paste("I-Type RF Centroid Log10 Shift","Iter",iexpIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(centroidDelta.i), max(centroidDelta.i) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( centroidDelta.i.2, paste(paste("I-Type RF Centroid Log10 Shift","Iter",iexpIter2,sep=" "),expTitleText2, sep=""), xLabText, yLabText, min(centroidDelta.i.2), max(centroidDelta.i.2) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( centroidDelta.i.3, paste(paste("I-Type RF Centroid Log10 Shift","Iter",iexpIter3,sep=" "),expTitleText3, sep=""), xLabText, yLabText, min(centroidDelta.i.3), max(centroidDelta.i.3) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );	

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

		#
		#	RF Size
		#
	if ( tiffFlag ) {
		tiff ( paste("Consol.Knock", activationValue, "Extent", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )

	par ( mfrow=c(2,4) );
	ref.rfExtent.e = QuantRFSize ( ref$r1.e.rfMap, kRFPeakToEdgeDetect );
	ref.rfExtent.i = QuantRFSize ( ref$r1.i.rfMap, kRFPeakToEdgeDetect );

	exp.rfExtent.e = QuantRFSize ( exp$r1.e.rfMap, kRFPeakToEdgeDetect );
	exp.rfExtent.i = log10(1+QuantRFSize ( exp$r1.i.rfMap, kRFPeakToEdgeDetect ));
	exp.rfExtent.e[c(iTrim,iFilter.E)] = 0; exp.rfExtent.i[c(iTrim,iFilter.I)] = 0;

	exp2.rfExtent.e = QuantRFSize ( exp2$r1.e.rfMap, kRFPeakToEdgeDetect );
	exp2.rfExtent.i = log10(1+QuantRFSize ( exp2$r1.i.rfMap, kRFPeakToEdgeDetect ));
	exp2.rfExtent.e[c(iTrim,iFilter.E.2)] = 0; exp2.rfExtent.i[c(iTrim,iFilter.I.2)] = 0;

	exp3.rfExtent.e = QuantRFSize ( exp3$r1.e.rfMap, kRFPeakToEdgeDetect );
	exp3.rfExtent.i = log10(1+QuantRFSize ( exp3$r1.i.rfMap, kRFPeakToEdgeDetect ));
	exp3.rfExtent.e[c(iTrim,iFilter.E.3)] = 0; exp3.rfExtent.i[c(iTrim,iFilter.I.3)] = 0;

	ShowVecAsMap2 ( ref.rfExtent.e, paste(paste("E-Type RF Extent","Iter",iexpIter,sep=" "),refTitleText, sep=""), xLabText, yLabText, min(ref.rfExtent.e), max(ref.rfExtent.e) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( exp.rfExtent.e, paste(paste("E-Type RF Extent","Iter",iexpIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(exp.rfExtent.e), max(exp.rfExtent.e) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( exp2.rfExtent.e, paste(paste("E-Type RF Extent","Iter",iexpIter2,sep=" "),expTitleText2, sep=""), xLabText, yLabText, min(exp2.rfExtent.e), max(exp2.rfExtent.e) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 )
;
	ShowVecAsMap2 ( exp3.rfExtent.e, paste(paste("E-Type RF Extent","Iter",iexpIter3,sep=" "),expTitleText3, sep=""), xLabText, yLabText, min(exp3.rfExtent.e), max(exp3.rfExtent.e) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( ref.rfExtent.i, paste(paste("I-Type RF Extent","Iter",iexpIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(ref.rfExtent.i), max(ref.rfExtent.i) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( exp.rfExtent.i, paste(paste("I-Type RF Extent","Iter",iexpIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(exp.rfExtent.i), max(exp.rfExtent.i) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( exp2.rfExtent.i, paste(paste("I-Type RF Extent","Iter",iexpIter2,sep=" "),expTitleText2, sep=""), xLabText, yLabText, min(exp2.rfExtent.i), max(exp2.rfExtent.i) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( exp3.rfExtent.i, paste(paste("I-Type RF Extent","Iter",iexpIter3,sep=" "),expTitleText3, sep=""), xLabText, yLabText, min(exp3.rfExtent.i), max(exp3.rfExtent.i) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )


		#
		#	RF Size Change
		#
	if ( tiffFlag ) {
		tiff ( paste("Consol.Knock", activationValue, "ExtentDelta", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )

	par ( mfrow=c(2,3) );
	ref.rfExtent.e = QuantRFSize ( ref$r1.e.rfMap, kRFPeakToEdgeDetect );
	ref.rfExtent.i = QuantRFSize ( ref$r1.i.rfMap, kRFPeakToEdgeDetect );
	tooSmallList = which ( ref.rfExtent.e < ( mean(ref.rfExtent.e) - 3.0 * sqrt ( var ( ref.rfExtent.e ) )  ) );

	sizeDelta.e = ( QuantRFSize ( exp$r1.e.rfMap, kRFPeakToEdgeDetect ) / ref.rfExtent.e ) - 1.0;
	sizeDelta.i = ( QuantRFSize ( exp$r1.i.rfMap, kRFPeakToEdgeDetect ) / ref.rfExtent.i ) - 1.0;
	sizeDelta.e[c(iTrim,iFilter.E)] = 0; sizeDelta.i[c(iTrim,iFilter.I)] = 0;
	if ( length ( tooSmallList ) ) { sizeDelta.e[tooSmallList] = 0; }

	sizeDelta.e.2 = ( QuantRFSize ( exp2$r1.e.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( ref$r1.e.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
	sizeDelta.i.2 = ( QuantRFSize ( exp2$r1.i.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( ref$r1.i.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
	sizeDelta.e.2[c(iTrim,iFilter.E.2)] = 0; sizeDelta.i.2[c(iTrim,iFilter.I.2)] = 0;
	if ( length ( tooSmallList ) ) { sizeDelta.e.2[tooSmallList] = 0; }

	sizeDelta.e.3 = ( QuantRFSize ( exp3$r1.e.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( ref$r1.e.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
	sizeDelta.i.3 = ( QuantRFSize ( exp3$r1.i.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( ref$r1.i.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
	sizeDelta.e.3[c(iTrim,iFilter.E.3)] = 0; sizeDelta.i.3[c(iTrim,iFilter.I.3)] = 0;
	if ( length ( tooSmallList ) ) { sizeDelta.e.3[tooSmallList] = 0; }

	ShowVecAsMap2 ( sizeDelta.e, paste(paste("E-Type % RF Size Change","Iter",iexpIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(sizeDelta.e), max(sizeDelta.e) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( sizeDelta.e.2, paste(paste("E-Type % RF Size Change","Iter",iexpIter2,sep=" "),expTitleText2, sep=""), xLabText, yLabText, min(sizeDelta.e.2), max(sizeDelta.e.2) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( sizeDelta.e.3, paste(paste("E-Type % RF Size Change","Iter",iexpIter3,sep=" "),expTitleText3, sep=""), xLabText, yLabText, min(sizeDelta.e.3), max(sizeDelta.e.3) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( sizeDelta.i, paste(paste("I-Type % RF Size Change","Iter",iexpIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(sizeDelta.i), max(sizeDelta.i) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( sizeDelta.i.2, paste(paste("I-Type % RF Size Change","Iter",iexpIter2,sep=" "),expTitleText2, sep=""), xLabText, yLabText, min(sizeDelta.i.2), max(sizeDelta.i.2) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( sizeDelta.i.3, paste(paste("I-Type % RF Size Change","Iter",iexpIter3,sep=" "),expTitleText3, sep=""), xLabText, yLabText, min(sizeDelta.i.3), max(sizeDelta.i.3) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

				#
		#	RF Max Magnitude Response
		#
	if ( tiffFlag ) {
		tiff ( paste("Consol.Knock", activationValue, "MaxResp", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )
	par ( mfrow=c(2,3) );

	ref.maxresp.e = apply ( ref$r1.e.rfMap, 1, max );
	ref.maxresp.i = apply ( ref$r1.i.rfMap, 1, max );

	exp.maxresp.e = apply ( exp$r1.e.rfMap, 1, max );
	exp.maxresp.i = apply ( exp$r1.i.rfMap, 1, max );

	exp.maxresp.e.2 = apply ( exp2$r1.e.rfMap, 1, max );
	exp.maxresp.i.2 = apply ( exp2$r1.i.rfMap, 1, max );

	exp.maxresp.e.3 = apply ( exp3$r1.e.rfMap, 1, max );
	exp.maxresp.i.3 = apply ( exp3$r1.i.rfMap, 1, max );

	delta.maxresp.e = ( exp.maxresp.e / ref.maxresp.e ) - 1.0;
	delta.maxresp.i = ( exp.maxresp.i / ref.maxresp.i ) - 1.0;

	delta.maxresp.e.2 = ( exp.maxresp.e.2 / ref.maxresp.e ) - 1.0;
	delta.maxresp.i.2 = ( exp.maxresp.i.2 / ref.maxresp.i ) - 1.0;

	delta.maxresp.e.3 = ( exp.maxresp.e.3 / ref.maxresp.e ) - 1.0;
	delta.maxresp.i.3 = ( exp.maxresp.i.3 / ref.maxresp.i ) - 1.0;

	ref.maxresp.e[c(iTrim,iFilter.E)] = 0; ref.maxresp.i[c(iTrim,iFilter.I)] = 0;

	exp.maxresp.e[c(iTrim,iFilter.E)] = 0; exp.maxresp.i[c(iTrim,iFilter.I)] = 0;
	exp.maxresp.e.2[c(iTrim,iFilter.E.2)] = 0; exp.maxresp.i.2[c(iTrim,iFilter.I.2)] = 0;
	exp.maxresp.e.3[c(iTrim,iFilter.E.3)] = 0; exp.maxresp.i.3[c(iTrim,iFilter.I.3)] = 0;

	delta.maxresp.e[c(iTrim,iFilter.E)] = 0; delta.maxresp.i[c(iTrim,iFilter.I)] = 0;
	delta.maxresp.e.2[c(iTrim,iFilter.E.2)] = 0; delta.maxresp.i.2[c(iTrim,iFilter.I.2)] = 0;
	delta.maxresp.e.3[c(iTrim,iFilter.E.3)] = 0; delta.maxresp.i.3[c(iTrim,iFilter.I.3)] = 0;

	if ( 0 ) {
		ShowVecAsMap2 ( exp.maxresp.e, paste(paste("E-Type RF Mag Resp","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText,
			min(exp.maxresp.e), max(exp.maxresp.e) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

		ShowVecAsMap2 ( exp.maxresp.i, paste(paste("I-Type RF Mag Resp","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText,
			min(exp.maxresp.i), max(exp.maxresp.i) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	}

	ShowVecAsMap2 ( delta.maxresp.e, paste(paste("E-Type % RF Mag Diff","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText,
			min(delta.maxresp.e), max(delta.maxresp.e) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( delta.maxresp.e.2, paste(paste("E-Type % RF Mag Diff","Iter",irefIter,sep=" "),expTitleText2, sep=""), xLabText, yLabText,
			min(delta.maxresp.e.2), max(delta.maxresp.e.2) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( delta.maxresp.e.3, paste(paste("E-Type % RF Mag Diff","Iter",irefIter,sep=" "),expTitleText3, sep=""), xLabText, yLabText,
			min(delta.maxresp.e.3), max(delta.maxresp.e.3) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( delta.maxresp.i, paste(paste("I-Type % RF Mag Diff","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText,
			min(delta.maxresp.i), max(delta.maxresp.i) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( delta.maxresp.i.2, paste(paste("I-Type % RF Mag Diff","Iter",irefIter,sep=" "),expTitleText2, sep=""), xLabText, yLabText,
			min(delta.maxresp.i.2), max(delta.maxresp.i.2) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( delta.maxresp.i.3, paste(paste("I-Type % RF Mag Diff","Iter",irefIter,sep=" "),expTitleText3, sep=""), xLabText, yLabText,
			min(delta.maxresp.i.3), max(delta.maxresp.i.3) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

		#
		#	Intracolumnar RF Centroid Divergence
		#
	if ( tiffFlag ) {
		tiff ( paste("Consol.Knock", activationValue, "Divergence", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )
	par ( mfrow=c(2,4) );

	centroidDelta.ref = RFCentroidDelta ( ref$rfTrackData.e[[length(ref$rfTrackData.e)]], ref$rfTrackData.i[[length(ref$rfTrackData.i)]]);
	centroidDelta.exp = RFCentroidDelta ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], exp$rfTrackData.i[[length(exp$rfTrackData.i)]]);
	centroidDelta.exp.2 = RFCentroidDelta ( exp2$rfTrackData.e[[length(exp2$rfTrackData.e)]], exp2$rfTrackData.i[[length(exp2$rfTrackData.i)]]);
	centroidDelta.exp.3 = RFCentroidDelta ( exp3$rfTrackData.e[[length(exp3$rfTrackData.e)]], exp3$rfTrackData.i[[length(exp3$rfTrackData.i)]]);
	centroidDelta.ref[c(iTrim,iFilter.E,iFilter.I)] = 0;
	centroidDelta.exp[c(iTrim,iFilter.I,iFilter.E)] = 0;
	centroidDelta.exp.2[c(iTrim,iFilter.I.2,iFilter.E.2)] = 0;
	centroidDelta.exp.3[c(iTrim,iFilter.I.3,iFilter.E.3)] = 0;


	ShowVecAsMap2 ( centroidDelta.ref, paste(paste("RF Divergence","Iter",irefIter,sep=" "),refTitleText, sep=""), xLabText, yLabText, min(centroidDelta.ref), max(centroidDelta.ref) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	
	ShowVecAsMap2 ( centroidDelta.exp, paste(paste("RF Divergence","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(centroidDelta.exp), max(centroidDelta.exp) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( centroidDelta.exp.2, paste(paste("RF Divergence","Iter",irefIter,sep=" "),expTitleText2, sep=""), xLabText, yLabText, min(centroidDelta.exp.2), max(centroidDelta.exp.2) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	ShowVecAsMap2 ( centroidDelta.exp.3, paste(paste("RF Divergence","Iter",irefIter,sep=" "),expTitleText3, sep=""), xLabText, yLabText, min(centroidDelta.exp.3), max(centroidDelta.exp.3) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	centroidDelta.ref[ centroidDelta.ref > 20 ] = 0.0;
	centroidDelta.ref[ centroidDelta.ref < 1.0 ] = 0.0;
	ShowVecAsMap2 ( centroidDelta.ref, paste(paste("RF Divergence","Iter",irefIter,sep=" "),refTitleText, sep=""), xLabText, yLabText, min(centroidDelta.ref), max(centroidDelta.ref) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	
	centroidDelta.exp[ centroidDelta.exp > 20 ] = 0.0;
	centroidDelta.exp[ centroidDelta.exp < 1.0  ] = 0.0;
	ShowVecAsMap2 ( centroidDelta.exp, paste(paste("RF Divergence","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(centroidDelta.exp), max(centroidDelta.exp) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	centroidDelta.exp.2[ centroidDelta.exp.2 > 20 ] = 0.0;
	centroidDelta.exp.2[ centroidDelta.exp.2 < 1.0  ] = 0.0;
	ShowVecAsMap2 ( centroidDelta.exp.2, paste(paste("RF Divergence","Iter",irefIter,sep=" "),expTitleText2, sep=""), xLabText, yLabText, min(centroidDelta.exp.2), max(centroidDelta.exp.2) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	centroidDelta.exp.3[ centroidDelta.exp.3 > 20 ] = 0.0;
	centroidDelta.exp.3[ centroidDelta.exp.3 < 1.0  ] = 0.0;
	ShowVecAsMap2 ( centroidDelta.exp.3, paste(paste("RF Divergence","Iter",irefIter,sep=" "),expTitleText3, sep=""), xLabText, yLabText, min(centroidDelta.exp.3), max(centroidDelta.exp.3) );
	abline ( h = boundaryMarks, lty=3, col=3 ); GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

	return ( list ( centroidDelta.e = centroidDelta.e, centroidDelta.i = centroidDelta.i, sizeDelta.e = sizeDelta.e, sizeDelta.i = sizeDelta.i ) );

} # PlotKnockoutRFMapConsolidated ( ref, exp, refTitleText, expTitleText, iTrim, iFilter.E, iFilter.I, iColList ) {


RFFlyOver = function ( ref, exp, N, iRowParm, iColParm, rfExpZone.x, rfExpZone.y ) {
	
	iRowParm = as.integer ( iRowParm );
	iColParm = as.integer ( iColParm );

	numColors = 128;
	xLabText  = "Distal -> Proximal";  yLabText = "Digit 1 -> Digit 3";
	graphicalPanelMarks = ( N + 0.5 );
	digitBorderMarks = c((N/3)+0.5, (2*N/3)+0.5);
	digitBorderMarks = c ( digitBorderMarks, N+digitBorderMarks );

	zmin = Inf; zmax = -Inf;

	#negFlag = matrix ( 0, nrow=100, ncol=3 ); negFlagCount = 1;

	for ( iCol in seq ( iColParm[1], iColParm[2], iColParm[3] ) ) {

		for ( iRow in seq ( iRowParm[1], iRowParm[2], iRowParm[3] ) ) {

			iCell = GetLin ( iRow, iCol, N );
			z = c ( ref$r1.i.rfMap[iCell, ], ref$r1.e.rfMap[iCell, ], exp$r1.i.rfMap[iCell, ], exp$r1.e.rfMap[iCell, ] );
			tmin = min ( z ); tmax = max ( z );
			#if ( min(z) < 0 ) {
				#negFlag[negFlagCount,1] = iRow;
				#negFlag[negFlagCount,2] = iCol;
				#negFlag[negFlagCount,3] = min(z);
				#negFlagCount = negFlagCount + 1;
			#} # if ( min(z) < 0 )
			if ( tmin < zmin ) { zmin = tmin; }
			if ( tmax > zmax ) { zmax = tmax; }

		} # for ( iRow in seq ( iRowParm[1], iRowParm[2], iRowParm[3] ) ) {

	} # for ( iCol in seq ( iColParm[1], iColParm[2], iColParm[3] ) ) {
	#zlim = log10 (1 + c ( zmin, zmax ) );
	zlim = ( c ( zmin, zmax ) );


	x11();
	saveHTML ( {

		for ( iCol in seq ( iColParm[1], iColParm[2], iColParm[3] ) ) {

			for ( iRow in seq ( iRowParm[1], iRowParm[2], iRowParm[3] ) ) {

				iCell = GetLin ( iRow, iCol, N );
				#z = log10 ( 1 + InterleaveForDisplay ( ref$r1.i.rfMap[iCell, ], ref$r1.e.rfMap[iCell, ],
				z = (InterleaveForDisplay ( ref$r1.i.rfMap[iCell, ], ref$r1.e.rfMap[iCell, ],
					exp$r1.i.rfMap[iCell, ], exp$r1.e.rfMap[iCell, ] ) );
				titleText = paste ( "Cortical Column # ", iCell, sep="" );
				image.plot( c(1:(2*N)), c(1:(2*N)), matrix ( (z), nrow=2*N, ncol=2*N, byrow=FALSE ),
					col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
					zlim = zlim, main=titleText, xlab=xLabText, ylab=yLabText );			
				abline(h=digitBorderMarks, col=3, lty=3, lwd=0.5 );
				abline(h=graphicalPanelMarks, v=graphicalPanelMarks, col=1, lty=1, lwd=2 );
				GenOutline4X( N, rfExpZone.x, rfExpZone.y, TRUE );

			} # for ( iRow in seq ( iRowParm[1], iRowParm[2], iRowParm[3] ) ) {

		} # for ( iCol in seq ( iColParm[1], iColParm[2], iColParm[3] ) ) {

	}, interval = 1 );  # saveGIF

} # RFFlyOver = function ( ref, exp, ... ) {

	#
	#
RFFlyOverA = function ( ref, exp, N, iRowParm, iColParm, rfExpZone.x, rfExpZone.y, kRFPeakToEdgeDetect ) {
	
	iRowParm = as.integer ( iRowParm );
	iColParm = as.integer ( iColParm );

	numColors = 128;
	xLabText  = "Distal -> Proximal";  yLabText = "Digit 1 -> Digit 3";
	graphicalPanelMarks = ( N + 0.5 );
	digitBorderMarks = c((N/3)+0.5, (2*N/3)+0.5);
	digitBorderMarks = c ( digitBorderMarks, N+digitBorderMarks );

	x11();
	saveHTML ( {
		for ( iRow in seq ( iRowParm[1], iRowParm[2], iRowParm[3] ) ) {

			for ( iCol in seq ( iColParm[1], iColParm[2], iColParm[3] ) ) {

				iCell = GetLin ( iRow, iCol, N );

				tmp1 = QuantRFSizeA ( ref$r1.i.rfMap, kRFPeakToEdgeDetect );
				t.ref.r1.i = ref$r1.i.rfMap[iCell, ];
				t.ref.r1.i[tmp1] = 0.0;

				tmp1 = QuantRFSizeA ( ref$r1.e.rfMap, kRFPeakToEdgeDetect );
				t.ref.r1.e = ref$r1.e.rfMap[iCell, ];
				t.ref.r1.e[tmp1] = 0.0;

				tmp1 = QuantRFSizeA ( ref$r1.i.rfMap, kRFPeakToEdgeDetect );
				t.exp.r1.i = exp$r1.i.rfMap[iCell, ];
				t.exp.r1.i[tmp1] = 0.0;

				tmp1 = QuantRFSizeA ( ref$r1.e.rfMap, kRFPeakToEdgeDetect );
				t.exp.r1.e = exp$r1.e.rfMap[iCell, ];
				t.exp.r1.e[tmp1] = 0.0;

				z = InterleaveForDisplay ( t.ref.r1.i, t.ref.r1.e, t.exp.r1.i, t.exp.r1.e );

				zlim = c ( min(z), max(z) );
				titleText = paste ( "Cortical Column # ", iCell, sep="" );

				image.plot( c(1:(2*N)), c(1:(2*N)), matrix ( z, nrow=2*N, ncol=2*N, byrow=FALSE ),
					col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
					zlim = zlim, main=titleText, xlab=xLabText, ylab=yLabText );
			
				abline(h=digitBorderMarks, col=3, lty=3, lwd=0.5 );
				abline(h=graphicalPanelMarks, v=graphicalPanelMarks, col=1, lty=1, lwd=2 );
				GenOutline4X( N, rfExpZone.x, rfExpZone.y, TRUE );

			} # for ( iCol in seq ( iColStart, iColEnd, iColStep ) {

		} # for ( iRow in seq ( iRowStart, iRowEnd, iRowStep ) {
	}, interval = 1 );  # saveGIF

} # A = function ( ref, exp, ... ) {

#
#	MAIN
#

		#
		#	Global constants.
		#
N = 75;
N2 = N * N;
kRFPeakToEdgeDetect = 0.5;
iTrim = EdgeTrimCellList ( seq ( 1, N2 ), N, 2 );	# Trim the outer-most edges.
boundaryMarks = c ( as.integer(N/3)+0.5, as.integer(N/3)*2+0.5 );
iRFTrackStepSize = 3;

makeMovie = FALSE;
tiffFlag = FALSE;

		#
		#	Detailed parameters describing the knockouts.
		#
activationValueList = c ( "0.2" );
iKnockLengthList = c ( 8 );

activationValue = "0.2";
iKnockLength = 8;

		#
		#	START: Generate Figures
		#
	iBase = iKnockLength;

		#	Set the directory where the experimental data is sitting.
	tmpText = paste ( "F:/NMLab/Working/T.75.7.Lesion.x", activationValue, ".", sep="" );
	fDir = paste ( tmpText, iKnockLength, "/", sep="" );

		#
		#
		#
	iKnockOffset = ceiling( ( N / 2 ) ) - 4;
	controlTrack = c ( ceiling( ( 2 * N / 3) ), N + 1 - ceiling( N / 6 ) );
	rfExpZone.x = c ( iKnockOffset, iKnockOffset + iKnockLength - 1 );
	rowParmsRFFlyOver = c ( iKnockOffset - 3, iKnockOffset + iKnockLength + 3, 1 );
	iRowList = seq ( rfExpZone.x[1]-1, rfExpZone.x[2]+1, 1 );
	#iRowList = c ( seq ( 24, 32, 2 ), seq ( 42, 50, 2 ) );

		#
		#	Get the Baseline Refinement RF Map (Use the Placebo output from Knock exps).
		#
	fRoot = "BorderKnockout_Control_Placebo.RFMap";
	fileRootName = paste(fDir, fRoot, sep="//");
	base = GetRFMapData2A ( fileRootName, iBase, 0, 0, 0, kRFPeakToEdgeDetect, N2 );

		#
		#	1.  Generate Consolidated Knockout Figure: CONTROL E, I and E&I
		#
	rfExpZone.y = rfExpZone.y.ctl = c ( ceiling(N/6), ceiling(N/6)+1 );
	iColList = iColList.knockout = c ( controlTrack, rfExpZone.y[1]-1, rfExpZone.y[2]+1 );
	iColList.knockin = c ( controlTrack );
	itmp = (rfExpZone.y[1]-1) * N + iKnockOffset;

	fRoot = "BorderKnockout_Control_I.RFMap"; fileRootName = paste(fDir, fRoot, sep="\\");
	iStart = 1; iEnd = 1; iStepSize = 1; iExpFlag = 1;
	test.ctl.I = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
	test.ctl.I.iFilter.E = NULL;
	test.ctl.I.iFilter.I = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1 ) );

	fRoot = "BorderKnockout_Control_E.RFMap"; fileRootName = paste(fDir, fRoot, sep="\\");
	iStart = 2; iEnd = 2; iStepSize = 1; iExpFlag = 2;
	test.ctl.E = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
	test.ctl.E.iFilter.I = NULL;
	test.ctl.E.iFilter.E = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );

	fRoot = "BorderKnockout_Control_EI.RFMap"; fileRootName = paste(fDir, fRoot, sep="\\");
	iStart = 3; iEnd = 3; iStepSize = 1; iExpFlag = 3;
	test.ctl.EI = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
	test.ctl.EI.iFilter.E = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
	test.ctl.EI.iFilter.I = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );

	tmp = PlotKnockoutRFMapConsolidated ( base, 15, test.ctl.E, 15, test.ctl.I, 15, test.ctl.EI, 15,
			"\nBaseline Network", "\nWithin-Rep E Only", "\nWithin-Rep I Only", "\nWithin-Rep E and I ",
			activationValue, iTrim, test.ctl.E.iFilter.E, test.ctl.E.iFilter.I, test.ctl.I.iFilter.E, test.ctl.I.iFilter.I,
			test.ctl.EI.iFilter.E, test.ctl.EI.iFilter.I, iRowList, iColList, rfExpZone.x, rfExpZone.y, tiffFlag, iKnockLength, iExpFlag );

	if ( makeMovie ) {
		rowParmsRFFlyOver = c ( rfExpZone.x[1]-10, rfExpZone.x[2]+10, 1 );
		colParmsRFFlyOver = c ( as.integer(N/6)-10, as.integer(N/6)+10, 1 );
		RFFlyOver ( base, test.ctl.E, N, rowParmsRFFlyOver, colParmsRFFlyOver, rfExpZone.x, rfExpZone.y );
	} # if ( makeMovie )

		#
		#	2.  Generate Consolidated Knockout Figure: One-Side E, I and E&I
		#
	#activationValue = "-2.0"; iKnockLength = 4;
	rfExpZone.y = rfExpZone.y.oneside = c ( round((N/3),0)-1, round((N/3),0)-0 );
	iColList = c ( controlTrack, ( ceiling ( N / 3 ) - 2 ),  ( ceiling ( N / 3 ) + 1 ) );
	itmp = (rfExpZone.y[1]-1) * N + iKnockOffset;

	fRoot = "BorderKnockout_OneSide_I.RFMap"; fileRootName = paste(fDir, fRoot, sep="\\");
	iStart = 4; iEnd = 4; iStepSize = 1; iExpFlag = 4;
	test.oneside.I = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
	test.oneside.I.iFilter.E = NULL;
	test.oneside.I.iFilter.I = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1 ) );

	fRoot = "BorderKnockout_OneSide_E.RFMap";	fileRootName = paste(fDir, fRoot, sep="\\");
	iStart = 5; iEnd = 5; iStepSize = 1; iExpFlag = 5;
	test.oneside.E = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
	test.oneside.E.iFilter.I = NULL;
	test.oneside.E.iFilter.E = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );

	fRoot = "BorderKnockout_OneSide_EI.RFMap"; fileRootName = paste(fDir, fRoot, sep="\\");
	iStart = 6; iEnd = 6; iStepSize = 1; iExpFlag = 6;
	test.oneside.EI = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
	test.oneside.EI.iFilter.E = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
	test.oneside.EI.iFilter.I = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
	
	tmp = PlotKnockoutRFMapConsolidated ( base, 15, test.oneside.E, 15, test.oneside.I, 15, test.oneside.EI, 15,
			"\nBaseline Network", "\nBoundary E Only", "\nBoundary I Only", "\nBoundary E and I ",
			activationValue, iTrim, test.oneside.E.iFilter.E, test.oneside.E.iFilter.I, test.oneside.I.iFilter.E, test.oneside.I.iFilter.I,
			test.oneside.EI.iFilter.E, test.oneside.EI.iFilter.I, iRowList, iColList, rfExpZone.x, rfExpZone.y, tiffFlag, iKnockLength, iExpFlag );

		#
		#	END: Generate Figures
		#


