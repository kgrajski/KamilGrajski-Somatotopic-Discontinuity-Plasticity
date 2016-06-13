#
#	KnockoutFigure4.R
#
#
#	Automatically handling the control, one-side and two-side I-cell knockout variations.
#

#	Clear the workspace.
rm(list = ls());

#	Load up some libraries and helper functions.

source ( "NMHelperFunctions.R" );

	#
	#	Standardized plots.
	#
PlotKnockoutRFMap = function ( ref, irefIter, exp, iexpIter, refTitleText, expTitleText,
						iTrim, iFilter.E, iFilter.I, iRowList, iColList, x, y, tiffFlag=FALSE, iKnockOutLength=0, iExpFlag=0 ) {

		#	Set up the parameters so that the "experimental zone" will be outlined on subsequent plots.

	y[1] = y[1] - 0.5; y[2] = y[2] + 0.5;
	x[1] = x[1] - 0.5; x[2] = x[2] + 0.5;
	tmp.x = rbind(c(x[1], x[2]), c(x[1], x[2]), c(x[1], x[1]), c(x[2], x[2]));
	tmp.y = rbind(c(y[1], y[1]), c(y[2], y[2]), c(y[1], y[2]), c(y[1], y[2]));

	xLabText  = "Distal -> Proximal";  yLabText = "Digit 1 -> Digit 3";
	N = as.integer ( sqrt ( dim(base$r1.i.rfMap)[1] ) );
	boundaryMarks = c ( as.integer(N/3)+0.5, as.integer(N/3)*2+0.5 );
	iRFTrackStepSize = 2;

		#	RF Centroids
	if ( tiffFlag ) {
		tiff ( paste("Knock", "Centroids", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )
	par ( mfcol=c(2,2) );
	ShowTopoMap1 ( ref$rfTrackData.e[[length(ref$rfTrackData.e)]], paste(paste("E-Type","Iter",irefIter,sep=" "),refTitleText, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1 ( ref$rfTrackData.i[[length(ref$rfTrackData.i )]], paste(paste("I-Type","Iter",irefIter,sep=" "),refTitleText, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1 ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowTopoMap1 ( exp$rfTrackData.i[[length(exp$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""), FALSE, 0.5, 0 );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

		#	RF Tracks
	if ( tiffFlag ) {
		tiff ( paste("Knock", "Tracks", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )
	par ( mfcol=c(2,2) );
	ShowThreeDigitRFTrackFlex ( ref$rfTrackData.e[[length(ref$rfTrackData.e)]], paste(paste("E-Type","Iter",irefIter,sep=" "),refTitleText, sep=""),
		TRUE, 0.5, 0, 1, iRowList, iColList, iRFTrackStepSize );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowThreeDigitRFTrackFlex ( ref$rfTrackData.i[[length(ref$rfTrackData.i)]], paste(paste("I-Type","Iter",irefIter,sep=" "),refTitleText, sep=""),
		TRUE, 0.5, 0, 1, iRowList, iColList, iRFTrackStepSize );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowThreeDigitRFTrackFlex ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""),
		TRUE, 0.5, 0, 1, iRowList, iColList, iRFTrackStepSize );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	ShowThreeDigitRFTrackFlex ( exp$rfTrackData.i[[length(exp$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""),
		TRUE, 0.5, 0, 1, iRowList, iColList, iRFTrackStepSize );
	GenOutline1X ( tmp.x, tmp.y, 2, 1, 2 );
	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

		#	RF Translocation
	if ( tiffFlag ) {
		tiff ( paste("Knock", "Position_Size", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )
	par ( mfcol=c(2,2) );
	centroidDelta.e = RFCentroidDelta ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], ref$rfTrackData.e[[length(ref$rfTrackData.e)]]);
	centroidDelta.i = RFCentroidDelta ( exp$rfTrackData.i[[length(exp$rfTrackData.i)]], ref$rfTrackData.i[[length(ref$rfTrackData.i)]]);
	centroidDelta.e[c(iTrim,iFilter.E)] = 0;
	centroidDelta.i[c(iTrim,iFilter.I)] = 0;
	#zmin = min ( centroidDelta.e, centroidDelta.i ); zmax = max ( centroidDelta.e, centroidDelta.i );
	ShowVecAsMap2 ( centroidDelta.e, paste(paste("E-Type RF Centroid Shift","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(centroidDelta.e), max(centroidDelta.e) );
	abline ( h = boundaryMarks, lty=3, col=3 );
	GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( centroidDelta.i, paste(paste("I-Type RF Centroid Shift","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(centroidDelta.i), max(centroidDelta.i) );
	abline ( h = boundaryMarks, lty=3, col=3 );
	GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

		#	RF Size Change
	#ref.rfSize = QuantRFSize ( exp$r1.e.rfMap, kRFPeakToEdgeDetect ); exp.rfSize = QuantRFSize ( ref$r1.e.rfMap, kRFPeakToEdgeDetect );
	#ref.rfSize.mean = mean ( ref.rfSize ); ref.rfSize.sd = sqrt ( var ( ref.rfSize ) );
	#iSizeFilter = seq(1,N2)[( ref.rfSize < ( ref.rfSize.mean - 2.5*ref.rfSize.sd ) ) | ( ref.rfSize > ( ref.rfSize.mean + 2.5*ref.rfSize.sd ) )];

	sizeDelta.e = ( QuantRFSize ( exp$r1.e.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( ref$r1.e.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
	sizeDelta.i = ( QuantRFSize ( exp$r1.i.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( ref$r1.i.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
	sizeDelta.e[c(iTrim,iFilter.E)] = 0;
	sizeDelta.i[c(iTrim,iFilter.I)] = 0;
	ShowVecAsMap2 ( sizeDelta.e, paste(paste("E-Type % RF Size Change","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(sizeDelta.e), max(sizeDelta.e) );
	abline ( h = boundaryMarks, lty=3, col=3 );
	GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	ShowVecAsMap2 ( sizeDelta.i, paste(paste("I-Type % RF Size Change","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(sizeDelta.i), max(sizeDelta.i) );
	abline ( h = boundaryMarks, lty=3, col=3 );
	GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

		#	Intracolumnar RF Centroid Divergence
	if ( tiffFlag ) {
		tiff ( paste("Knock", "Divergence", iKnockOutLength, iExpFlag, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )
	par ( mfcol=c(2,2) );
	centroidDelta.ref = RFCentroidDelta ( ref$rfTrackData.e[[length(ref$rfTrackData.e)]], ref$rfTrackData.i[[length(ref$rfTrackData.i)]]);
	centroidDelta.exp = RFCentroidDelta ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], exp$rfTrackData.i[[length(exp$rfTrackData.i)]]);
	centroidDelta.ref[c(iTrim,iFilter.E,iFilter.I)] = 0;
	centroidDelta.exp[c(iTrim,iFilter.I,iFilter.E)] = 0;

	ShowVecAsMap2 ( centroidDelta.ref, paste(paste("Baseline RF Divergence","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(centroidDelta.ref), max(centroidDelta.ref) );
	abline ( h = boundaryMarks, lty=3, col=3 );
	GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	
	ShowVecAsMap2 ( centroidDelta.exp, paste(paste("Knockout RF Divergence","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(centroidDelta.exp), max(centroidDelta.exp) );
	abline ( h = boundaryMarks, lty=3, col=3 );
	GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );

	centroidDelta.ref[ centroidDelta.ref > 20 ] = 0.0;
	ShowVecAsMap2 ( centroidDelta.ref, paste(paste("Baseline RF Divergence","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(centroidDelta.ref), max(centroidDelta.ref) );
	abline ( h = boundaryMarks, lty=3, col=3 );
	GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	
	centroidDelta.exp[ centroidDelta.exp > 20 ] = 0.0;
	ShowVecAsMap2 ( centroidDelta.exp, paste(paste("Knockout RF Divergence","Iter",irefIter,sep=" "),expTitleText, sep=""), xLabText, yLabText, min(centroidDelta.exp), max(centroidDelta.exp) );
	abline ( h = boundaryMarks, lty=3, col=3 );
	GenOutline1X ( tmp.x, tmp.y, "white", 1, 0.5 );
	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

	return ( list ( centroidDelta.e = centroidDelta.e, centroidDelta.i = centroidDelta.i, sizeDelta.e = sizeDelta.e, sizeDelta.i = sizeDelta.i ) );

} # PlotKnockoutRFMap ( ref, exp, refTitleText, expTitleText, iTrim, iFilter.E, iFilter.I, iColList ) {

GenOutline1X = function ( tmp.x, tmp.y, kCol, kLty, kLwd ) {
	
	lines ( tmp.x, tmp.y, col=kCol, lty=kLty, lwd=kLwd );

} # GenOutline1X = function ( tmp.x, tmp.y, kCol, kLty, kLwd ) {

GenOutline4X = function ( N, x, y, doPlot ) {
	
	y[1] = y[1] - 0.5;
	y[2] = y[2] + 0.5;
	x[1] = x[1] - 0.5;
	x[2] = x[2] + 0.5;

	tmp.x = rbind(c(x[1], x[2]), c(x[1], x[2]), c(x[1], x[1]), c(x[2], x[2]));
	tmp.y = rbind(c(y[1], y[1]), c(y[2], y[2]), c(y[1], y[2]), c(y[1], y[2]));

	if ( doPlot ) {
		GenOutline1X ( tmp.x, tmp.y, "white", 1, .5 );
		GenOutline1X ( tmp.x, N+tmp.y, "white", 1, 0.5 );
		GenOutline1X ( N+tmp.x, tmp.y, "white", 1, 0.5 );
		GenOutline1X ( N+tmp.x, N+tmp.y, "white", 1, 0.5 );
	} # if ( doPlot )

} # GenOutline4X = function ( N, x, y, doPlot ) {

RFFlyOver = function ( ref, exp, N, iRowParm, iColParm, rfExpZone.x, rfExpZone.y ) {
	
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
				z = InterleaveForDisplay ( ref$r1.i.rfMap[iCell,], ref$r1.e.rfMap[iCell, ], exp$r1.i.rfMap[iCell, ], exp$r1.e.rfMap[iCell, ] );
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

} # RFFlyOver = function ( ref, exp, ... ) {

#
#	MAIN
#

		#
		#	Global constants.
		#
N = 45;
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
iBase = iKnockLength = 4;

		#
		#
		#
iKnockOffset = round( ( N / 2 ), 0 ) - 4;
controlTrack = c ( round( ( 2 * N / 3), 0 ), N + 1 - round( N / 6, 0 ) );
rfExpZone.x = c ( iKnockOffset + 1, iKnockOffset + iKnockLength );
rowParmsRFFlyOver = c ( iKnockOffset - 3, iKnockOffset + iKnockLength + 3, 1 );
iRowList = seq ( rfExpZone.x[1]-2, rfExpZone.x[2]+2, 2 );

		#	Set the directory where the experimental data is sitting.
fDir = "D:/NMLab/S.45.7.Knockout.2XRFProbeMag/";
fDir = "D:/NMLab/S.45.7.Knockout/";

		#
		#	Get the Baseline Refinement RF Map.
		#
fRoot = "Base.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
base = GetRFMapData2A ( fileRootName, 15, 15, 15, 15, kRFPeakToEdgeDetect, N2 );

		##############################################################
		##############################################################
		#
		#	0.  Get the Placebo
		#
		##############################################################
		##############################################################

fRoot = "BorderKnockout_Control_Placebo.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 0; iEnd = 0; iStepSize = 1; iExpFlag = 0;
placebo = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );

			#	Set up and do the basic RF Map plots
iColList = c ( controlTrack, as.integer ( N / 6 - 1 ), as.integer ( N / 6 + 2 ) );
rfExpZone.y = c ( as.integer(N/6), as.integer(N/6)+1 );
rfstats.placebo = PlotKnockoutRFMap ( base, 15, placebo, 15, "\nBaseline Refinement", "\nControl Placebo", iTrim,
				NULL, NULL, iRowList, iColList, rfExpZone.x, rfExpZone.y, tiffFlag, iKnockLength, iExpFlag );

			#	Set up and do the RF "fly over"
if ( makeMovie ) {
	colParmsRFFlyOver = c ( as.integer(N/6)-3, as.integer(N/6)+4, 1 );
	RFFlyOver ( base, placebo, N, rowParmsRFFlyOver, colParmsRFFlyOver, rfExpZone.x, rfExpZone.y );
} # if ( makeMovie )

summary ( rfstats.placebo$centroidDelta.e )
summary ( rfstats.placebo$centroidDelta.i )
summary ( rfstats.placebo$sizeDelta.e )
summary ( rfstats.placebo$sizeDelta.i )

		##############################################################
		##############################################################
		#
		#	1.  Get the Knockout RF Map - CONTROL I Only.
		#
		##############################################################
		##############################################################

fRoot = "BorderKnockout_Control_I.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 1; iEnd = 1; iStepSize = 1; iExpFlag = 1;
test.ctl.I = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );

			#	Set up and do the basic RF Map plots
rfExpZone.y = c ( round((N/6),0), round((N/6),0)+1 );			# 	Partially describes the zone that was manipulated.
iColList = c ( controlTrack, as.integer ( N / 6 - 1 ), as.integer ( N / 6 + 2 ) );
itmp = as.integer ( N/6 ) * N + 1 + iKnockOffset;
iFilter.E = NULL;
iFilter.I = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
rfstats.test.ctl.I = PlotKnockoutRFMap ( base, 15, test.ctl.I, 15, "\nBaseline Refinement", "\nControl I Only", iTrim,
				iFilter.E, iFilter.I, iRowList, iColList, rfExpZone.x, rfExpZone.y, tiffFlag, iKnockLength, iExpFlag );

			#	Set up and do the RF "fly over"
if ( makeMovie ) {
	colParmsRFFlyOver = c ( as.integer(N/6)-2, as.integer(N/3) + 1, 1 );	#	Partially describes the fly-over zone.
	RFFlyOver ( base, test.ctl.I, N, rowParmsRFFlyOver, colParmsRFFlyOver, rfExpZone.x, rfExpZone.y );
} # if ( makeMovie )

summary ( rfstats.test.ctl.I$centroidDelta.e )
summary ( rfstats.test.ctl.I$centroidDelta.i )
summary ( rfstats.test.ctl.I$sizeDelta.e )
summary ( rfstats.test.ctl.I$sizeDelta.i )

		##############################################################
		##############################################################
		#
		#	2.  Get the Knockout RF Map - CONTROL E Only.
		#		
		##############################################################
		##############################################################

fRoot = "BorderKnockout_Control_E.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 2; iEnd = 2; iStepSize = 1; iExpFlag = 2;
test.ctl.E = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );

			#	Set up and do the basic RF Map plots
rfExpZone.y = c ( round((N/6),0), round((N/6),0)+1 );			# 	Partially describes the zone that was manipulated.
iColList = c ( controlTrack, as.integer ( N / 6 - 1 ), as.integer ( N / 6 + 2 ) );
iFilter.I = NULL;
itmp = as.integer ( N/6 ) * N + 1 + iKnockOffset;
iFilter.E = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
rfstats.test.ctl.E = PlotKnockoutRFMap ( base, 15, test.ctl.E, 15, "\nBaseline Refinement", "\nControl E Only", iTrim,
				iFilter.E, iFilter.I, iRowList, iColList, rfExpZone.x, rfExpZone.y, tiffFlag, iKnockLength, iExpFlag );

			#	Set up and do the RF "fly over"
if ( makeMovie ) {
	colParmsRFFlyOver = c ( as.integer(N/6)-2, as.integer(N/3) + 1, 1 );	#	Partially describes the fly-over zone.
	RFFlyOver ( base, test.ctl.E, N, rowParmsRFFlyOver, colParmsRFFlyOver, rfExpZone.x, rfExpZone.y );
} # if ( makeMovie )

		##############################################################
		##############################################################
		#
		#	3.  Get the Knockout RF Map - CONTROL E and I.
		#
		##############################################################
		##############################################################

fRoot = "BorderKnockout_Control_EI.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 3; iEnd = 3; iStepSize = 1; iExpFlag = 3;
test.ctl.EI = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );

			#	Set up and do the basic RF Map plots
rfExpZone.y = c ( round((N/6),0), round((N/6),0)+1 );			# 	Partially describes the zone that was manipulated.
iColList = c ( controlTrack, as.integer ( N / 6 - 1 ), as.integer ( N / 6 + 2 ) );
itmp = as.integer ( N/6 ) * N + 1 + iKnockOffset;
iFilter.E = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
iFilter.I = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
rfstats.test.ctl.EI = PlotKnockoutRFMap ( base, 15, test.ctl.EI, 15, "\nBaseline Refinement", "\nControl E and I", iTrim,
				iFilter.E, iFilter.I, iRowList, iColList, rfExpZone.x, rfExpZone.y, tiffFlag, iKnockLength, iExpFlag );
			
			#	Set up and do the RF "fly over"
if ( makeMovie ) {
	colParmsRFFlyOver = c ( as.integer(N/6)-2, as.integer(N/3) + 1, 1 );	#	Partially describes the fly-over zone.
	RFFlyOver ( base, test.ctl.EI, N, rowParmsRFFlyOver, colParmsRFFlyOver, rfExpZone.x, rfExpZone.y );
} # if ( makeMovie )
		##############################################################
		##############################################################
		#
		#	4.  Get the Knockout RF Map - ONE-SIDE I.
		#
		##############################################################
		##############################################################

fRoot = "BorderKnockout_OneSide_I.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 4; iEnd = 4; iStepSize = 1; iExpFlag = 4;
test.oneside.I = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );

			#	Set up and do the basic RF Map plots
rfExpZone.y = c ( round((N/3),0)-1, round((N/3),0)-0 );					# 	Partially describes the zone that was manipulated.
iColList = c ( controlTrack, 16, 17 );
iColList = c ( controlTrack, as.integer ( N / 3 - 1 ), as.integer ( N / 3 + 0 ) );
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 3 ) );
iFilter.E = NULL;
itmp = ( as.integer ( N/3 ) - 2) * N + 1 + iKnockOffset;
iFilter.I = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
rfstats.test.oneside.I = PlotKnockoutRFMap ( base, 15, test.oneside.I, 15, "\nBaseline Refinement", "\nBoundary One Side I Only", iTrim,
					iFilter.E, iFilter.I, iRowList, iColList, rfExpZone.x, rfExpZone.y, tiffFlag, iKnockLength, iExpFlag );
			
			#	Set up and do the RF "fly over"
if ( makeMovie ) {
	colParmsRFFlyOver = c ( as.integer(N/3) - (N/6), as.integer(N/3) + (N/6), 1 );	#	Partially describes the fly-over zone.
	RFFlyOver ( base, test.oneside.I, N, rowParmsRFFlyOver, colParmsRFFlyOver, rfExpZone.x, rfExpZone.y );
} # if ( makeMovie )

		##############################################################
		##############################################################
		#
		#	5.  Get the Knockout RF Map - ONE-SIDE E.
		#
		##############################################################
		##############################################################

fRoot = "BorderKnockout_OneSide_E.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 5; iEnd = 5; iStepSize = 1; iExpFlag = 5;
test.oneside.E = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
		
			#	Set up and do the basic RF Map plots
rfExpZone.y = c ( round((N/3),0)-1, round((N/3),0)-0 );					# 	Partially describes the zone that was manipulated.
iColList = c ( controlTrack, 16, 17 ); iColList = c ( controlTrack, 14, 15 );
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 3 ) );
iFilter.I = NULL;
itmp = ( as.integer ( N/3 ) - 2) * N + 1 + iKnockOffset;
iFilter.E = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
rfstats.test.oneside.E = PlotKnockoutRFMap ( base, 15, test.oneside.E, 15, "\nBaseline Refinement", "\nBoundary One Side E Only", iTrim,
					iFilter.E, iFilter.I, iRowList, iColList, rfExpZone.x, rfExpZone.y, tiffFlag, iKnockLength, iExpFlag );
			
			#	Additional longitudinal tracks.
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 2 ) );

			#	Set up and do the RF "fly over"
 if ( makeMovie )	{		
	colParmsRFFlyOver = c ( as.integer(N/3) - (N/6), as.integer(N/3) + (N/6), 1 );	#	Partially describes the fly-over zone.
	RFFlyOver ( base, test.oneside.E, N, rowParmsRFFlyOver, colParmsRFFlyOver, rfExpZone.x, rfExpZone.y );
} # if ( makeMovie )

		##############################################################
		##############################################################
		#
		#	6.  Get the Knockout RF Map - ONE-SIDE E and I.
		#
		##############################################################
		##############################################################

fRoot = "BorderKnockout_OneSide_EI.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 6; iEnd = 6; iStepSize = 1; iExpFlag = 6;
test.oneside.EI = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );

			#	Set up and do the basic RF Map plots
rfExpZone.y = c ( round((N/3),0)-1, round((N/3),0)-0 );					# 	Partially describes the zone that was manipulated.
iColList = c ( controlTrack, 16, 17 ); iColList = c ( controlTrack, 14, 15 );
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 3 ) );
itmp = ( as.integer ( N/3 ) - 2) * N + 1 + iKnockOffset;
iFilter.I = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
iFilter.E = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
rfstats.test.oneside.EI = PlotKnockoutRFMap ( base, 15, test.oneside.EI, 15, "\nBaseline Refinement", "\nBoundary One Side E and I", iTrim,
					iFilter.E, iFilter.I, iRowList, iColList, rfExpZone.x, rfExpZone.y, tiffFlag, iKnockLength, iExpFlag );
			
			#	Set up and do the RF "fly over"
if ( makeMovie ) {
	colParmsRFFlyOver = c ( as.integer(N/3) - (N/6), as.integer(N/3) + (N/6), 1 );	#	Partially describes the fly-over zone.
	RFFlyOver ( base, test.oneside.EI, N, rowParmsRFFlyOver, colParmsRFFlyOver, rfExpZone.x, rfExpZone.y );
} # if ( makeMovie )

		##############################################################
		##############################################################
		#
		#	7.  Get the Knockout RF Map - TWO-SIDE I.
		#
		##############################################################
		##############################################################

fRoot = "BorderKnockout_TwoSide_I.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 7; iEnd = 7; iStepSize = 1; iExpFlag = 7;
test.twoside.I = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );

			#	Set up and do the basic RF Map plots
rfExpZone.y = c ( round((N/3),0)-1, round((N/3),0)+2 );					# 	Partially describes the zone that was manipulated.
iColList = c ( controlTrack, 16, 17 ); iColList = c ( controlTrack, 14, 15 );
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 3 ) );
itmp = ( as.integer ( N/3 ) - 2) * N + 1 + iKnockOffset;
iFilter.E = NULL;
iFilter.I = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1 ),
			seq ( itmp + 2*N, itmp + 2*N + iKnockLength - 1 ), seq ( itmp + 3*N, itmp + 3*N + iKnockLength - 1 ) );
rfstats.test.twoside.I = PlotKnockoutRFMap ( base, 15, test.twoside.I, 15, "\nBaseline Refinement", "\nBoundary Two Side I Only", iTrim,
					iFilter.E, iFilter.I, iRowList, iColList, rfExpZone.x, rfExpZone.y, tiffFlag, iKnockLength, iExpFlag );
			
			#	Set up and do the RF "fly over"
if ( makeMovie ) {
	colParmsRFFlyOver = c ( as.integer(N/3) - (N/6), as.integer(N/3) + (N/6), 1 );	#	Partially describes the fly-over zone.
	RFFlyOver ( base, test.twoside.I, N, rowParmsRFFlyOver, colParmsRFFlyOver, rfExpZone.x, rfExpZone.y );
} # if ( makeMovie )
		
		##############################################################
		##############################################################
		#
		#	8.  Get the Knockout RF Map - TWO-SIDE E.
		#
		##############################################################
		##############################################################

fRoot = "BorderKnockout_TwoSide_E.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 8; iEnd = 8; iStepSize = 1; iExpFlag = 8;
test.twoside.E = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
			
			#	Set up and do the basic RF Map plots
rfExpZone.y = c ( round((N/3),0)-1, round((N/3),0)+2 );					# 	Partially describes the zone that was manipulated.
iColList = c ( controlTrack, 16, 17 ); iColList = c ( controlTrack, 14, 15 );
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 3 ) );
itmp = ( as.integer ( N/3 ) - 2) * N + 1 + iKnockOffset;
iFilter.I = NULL;
iFilter.E = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1 ),
			seq ( itmp + 2*N, itmp + 2*N + iKnockLength - 1 ), seq ( itmp + 3*N, itmp + 3*N + iKnockLength - 1 ) );
rfstats.test.twoside.I = PlotKnockoutRFMap ( base, 15, test.twoside.E, 15, "\nBaseline Refinement", "\nBoundary Two Side E Only", iTrim,
					iFilter.E, iFilter.I, iRowList, iColList, rfExpZone.x, rfExpZone.y, tiffFlag, iKnockLength, iExpFlag );
			
			#	Set up and do the RF "fly over"
if ( makeMovie ) {
	colParmsRFFlyOver = c ( as.integer(N/3) - (N/6), as.integer(N/3) + (N/6), 1 );	#	Partially describes the fly-over zone.
	RFFlyOver ( base, test.twoside.E, N, rowParmsRFFlyOver, colParmsRFFlyOver, rfExpZone.x, rfExpZone.y );
} # if ( makeMovie )
		
		##############################################################
		##############################################################
		#
		#	9.  Get the Knockout RF Map - TWO-SIDE E and I.
		#
		##############################################################
		##############################################################

fRoot = "BorderKnockout_TwoSide_EI.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 9; iEnd = 9; iStepSize = 1; iExpFlag = 9;
test.twoside.EI = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
test.twoside.EI = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, 0.1, N2 );

			#	Set up and do the basic RF Map plots
rfExpZone.y = c ( round((N/3),0)-1, round((N/3),0)+2 );					# 	Partially describes the zone that was manipulated.
iColList = c ( controlTrack, 16, 17 ); iColList = c ( controlTrack, 14, 15 );
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 3 ) );
itmp = ( as.integer ( N/3 ) - 2) * N + 1 + iKnockOffset;
iFilter.E = iFilter.I = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1 ),
			seq ( itmp + 2*N, itmp + 2*N + iKnockLength - 1 ), seq ( itmp + 3*N, itmp + 3*N + iKnockLength - 1 ) );
rfstats.test.twoside.I = PlotKnockoutRFMap ( base, 15, test.twoside.EI, 15, "\nBaseline Refinement", "\nBoundary Two Side E and I", iTrim,
					iFilter.E, iFilter.I, iRowList, iColList, rfExpZone.x, rfExpZone.y, tiffFlag, iKnockLength, iExpFlag );
			
			#	Set up and do the RF "fly over"
if ( makeMovie ) {
	colParmsRFFlyOver = c ( as.integer(N/3) - (N/6), as.integer(N/3) + (N/6), 1 );	#	Partially describes the fly-over zone.
	RFFlyOver ( base, test.twoside.EI, N, rowParmsRFFlyOver, colParmsRFFlyOver, rfExpZone.x, rfExpZone.y );
} # if ( makeMovie )	
		

