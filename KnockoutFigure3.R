#
#	KnockoutFigure3.R
#

#
#	Automatically handling the control, one-side and two-side I-cell knockout variations.
#

#	Clear the workspace.
rm(list = ls());

	#
	#	Standardized plots.
	#
PlotKnockoutRFMap = function ( ref, irefIter, exp, iexpIter, refTitleText, expTitleText, iTrim, iFilter.E, iFilter.I, iColList ) {

	N = as.integer ( sqrt ( dim(base$r1.i.rfMap)[1] ) );
	boundaryMarks = c ( as.integer(N/3)+0.5, as.integer(N/3)*2+0.5 );

		#	RF Centroids
	x11(); par ( mfcol=c(2,2) );
	ShowTopoMap1 ( ref$rfTrackData.e[[length(ref$rfTrackData.e)]], paste(paste("E-Type","Iter",irefIter,sep=" "),refTitleText, sep=""), FALSE, 0.5, 0 );
	ShowTopoMap1 ( ref$rfTrackData.i[[length(ref$rfTrackData.i )]], paste(paste("I-Type","Iter",irefIter,sep=" "),refTitleText, sep=""), FALSE, 0.5, 0 );
	ShowTopoMap1 ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""), FALSE, 0.5, 0 );
	ShowTopoMap1 ( exp$rfTrackData.i[[length(exp$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""), FALSE, 0.5, 0 );

		#	RF Tracks
	x11(); par(mfcol=c(2,2));
	ShowThreeDigitRFTrackFlex ( ref$rfTrackData.e[[length(ref$rfTrackData.e)]], paste(paste("E-Type","Iter",irefIter,sep=" "),refTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	ShowThreeDigitRFTrackFlex ( ref$rfTrackData.i[[length(ref$rfTrackData.i)]], paste(paste("I-Type","Iter",irefIter,sep=" "),refTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	ShowThreeDigitRFTrackFlex ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], paste(paste("E-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	ShowThreeDigitRFTrackFlex ( exp$rfTrackData.i[[length(exp$rfTrackData.i)]], paste(paste("I-Type","Iter",iexpIter,sep=" "),expTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	
		#	RF Translocation
	x11(); par(mfrow=c(2,2));
	centroidDelta.e = RFCentroidDelta ( exp$rfTrackData.e[[length(exp$rfTrackData.e)]], ref$rfTrackData.e[[length(ref$rfTrackData.e)]]);
	centroidDelta.i = RFCentroidDelta ( exp$rfTrackData.i[[length(exp$rfTrackData.i)]], ref$rfTrackData.i[[length(ref$rfTrackData.i)]]);
	centroidDelta.e[c(iTrim,iFilter.E)] = 0;
	centroidDelta.i[c(iTrim,iFilter.I)] = 0;
	#zmin = min ( centroidDelta.e, centroidDelta.i ); zmax = max ( centroidDelta.e, centroidDelta.i );
	ShowVecAsMap2 ( centroidDelta.e, paste(paste("E-Type RF Centroid Shift","Iter",irefIter,sep=" "),expTitleText, sep=""), "xLabText", "yLabText", min(centroidDelta.e), max(centroidDelta.e) );
	abline ( h = boundaryMarks, lty=3, col=3 );
	ShowVecAsMap2 ( centroidDelta.i, paste(paste("I-Type RF Centroid Shift","Iter",irefIter,sep=" "),expTitleText, sep=""), "xLabText", "yLabText", min(centroidDelta.i), max(centroidDelta.i) );
	abline ( h = boundaryMarks, lty=3, col=3 );

		#	RF Size Change
	sizeDelta.e = ( QuantRFSize ( exp$r1.e.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( ref$r1.e.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
	sizeDelta.i = ( QuantRFSize ( exp$r1.i.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( ref$r1.i.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
	sizeDelta.e[c(iTrim,iFilter.E)] = 0;
	sizeDelta.i[c(iTrim,iFilter.I)] = 0;
	#zmin = min ( sizeDelta.e, sizeDelta.i ); zmax = max ( sizeDelta.e, sizeDelta.i );
	ShowVecAsMap2 ( sizeDelta.e, paste(paste("E-Type % RF Size Change","Iter",irefIter,sep=" "),expTitleText, sep=""), "xLabText", "yLabText", min(sizeDelta.e), max(sizeDelta.e) );
	abline ( h = boundaryMarks, lty=3, col=3 );
	ShowVecAsMap2 ( sizeDelta.i, paste(paste("I-Type % RF Size Change","Iter",irefIter,sep=" "),expTitleText, sep=""), "xLabText", "yLabText", min(sizeDelta.i), max(sizeDelta.i) );
	abline ( h = boundaryMarks, lty=3, col=3 );

	return ( list ( centroidDelta.e = centroidDelta.e, centroidDelta.i = centroidDelta.i, sizeDelta.e = sizeDelta.e, sizeDelta.i = sizeDelta.i ) );

} # PlotKnockoutRFMap ( ref, exp, refTitleText, expTitleText, iTrim, iFilter.E, iFilter.I, iColList ) {


#	Load up some libraries and helper functions.
require(ggplot2);
require(reshape2);
require(ellipse);
require(stats4);
source ( "NMHelperFunctions.R" );

#
#	MAIN
#

		#
		#	Global constants.
		#
N = 45;
N2 = N * N;
kRFPeakToEdgeDetect = 0.5;
iTrim = c ( 1:N, seq ( N2, N2 - N + 1, -1 ), seq ( N, N2, N ), seq ( 1, N2, N ) );	# Trim the outer-most edges.

		#
		#	Detailed parameters describing the knockouts.
		#
iKnockOffset = as.integer ( N / 2 ) - 4;
controlTrack = as.integer ( 2 * N / 3);
iBase = iKnockLength = 12;

		#	Set the directory where the experimental data is sitting.
fDir = "D:/NMLab/S.45.7.Knockout - Copy/";

		#
		#	Get the Baseline Refinement RF Map.
		#
fRoot = "Base.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
base = GetRFMapData2A ( fileRootName, 15, 15, 15, 15, kRFPeakToEdgeDetect, N2 );

		#
		#	0.  Get the Placebo
		#
fRoot = "BorderKnockout_Control_Placebo.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 0; iEnd = 0; iStepSize = 1;
placebo = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iColList = c ( controlTrack, as.integer ( N / 6 - 1 ), as.integer ( N / 6 + 2 ) );
rfstats.placebo = PlotKnockoutRFMap ( base, 15, placebo, 15, "\nBaseline Refinement", "\nControl Placebo", iTrim, NULL, NULL, iColList );

		#
		#	1.  Get the Knockout RF Map - CONTROL I Only.
		#
fRoot = "BorderKnockout_Control_I.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 1; iEnd = 1; iStepSize = 1;
test.ctl.I = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iColList = c ( controlTrack, 10, as.integer ( N / 6 ), as.integer ( N / 6 ) );
iColList = c ( controlTrack, 10 );
itmp = as.integer ( N/6 ) * N + 1 + iKnockOffset;
iFilter.E = NULL;
iFilter.I = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
rfstats.test.ctl.I = PlotKnockoutRFMap ( base, 15, test.ctl.I, 15, "\nBaseline Refinement", "\nControl I Only", iTrim, iFilter.E, iFilter.I, iColList );
summary ( rfstats.test.ctl.I$centroidDelta.e )
summary ( rfstats.test.ctl.I$centroidDelta.i )
summary ( rfstats.test.ctl.I$sizeDelta.e )
summary ( rfstats.test.ctl.I$sizeDelta.i )

		#
		#	2.  Get the Knockout RF Map - CONTROL E Only.
		#
fRoot = "BorderKnockout_Control_E.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 2; iEnd = 2; iStepSize = 1;
test.ctl.E = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iColList = c ( controlTrack, as.integer ( N / 6 - 1 ), as.integer ( N / 6 + 2 ) );
iFilter.I = NULL;
itmp = as.integer ( N/6 ) * N + 1 + iKnockOffset;
iFilter.E = c ( seq ( itmp, itmp + iKnockLength - 1 ), seq ( itmp + N, itmp + N + iKnockLength - 1  ) );
rfstats.test.ctl.E = PlotKnockoutRFMap ( base, 15, test.ctl.E, 15, "\nBaseline Refinement", "\nControl E Only", iTrim, iFilter.E, iFilter.I, iColList );
summary ( rfstats.test.ctl.E$centroidDelta.e )
summary ( rfstats.test.ctl.E$centroidDelta.i )
summary ( rfstats.test.ctl.E$sizeDelta.e )
summary ( rfstats.test.ctl.E$sizeDelta.i )
		#
		#	3.  Get the Knockout RF Map - CONTROL E and I.
		#
fRoot = "BorderKnockout_Control_EI.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 3; iEnd = 3; iStepSize = 1;
test.ctl.EI = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iColList = c ( controlTrack, as.integer ( N / 6 - 1 ), as.integer ( N / 6 + 2 ) );
rfstats.test.ctl.EI = PlotKnockoutRFMap ( base, 15, test.ctl.EI, 15, "\nBaseline Refinement", "\nControl E and I", iTrim, iFilter.E, iFilter.I, iColList );

		#
		#	4.  Get the Knockout RF Map - ONE-SIDE I.
		#
fRoot = "BorderKnockout_OneSide_I.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 4; iEnd = 4; iStepSize = 1;
test.oneside.I = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 2 ) );
rfstats.test.oneside.I = PlotKnockoutRFMap ( base, 15, test.oneside.I, 15, "\nBaseline Refinement", "\nBoundary One Side I Only", iTrim, iFilter.E, iFilter.I, iColList );

		#
		#	5.  Get the Knockout RF Map - ONE-SIDE E.
		#
fRoot = "BorderKnockout_OneSide_E.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 5; iEnd = 5; iStepSize = 1;
test.oneside.E = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 2 ) );
rfstats.test.oneside.E = PlotKnockoutRFMap ( base, 15, test.oneside.E, 15, "\nBaseline Refinement", "\nBoundary One Side E Only", iTrim, iFilter.E, iFilter.I, iColList );

		#
		#	6.  Get the Knockout RF Map - ONE-SIDE E and I.
		#
fRoot = "BorderKnockout_OneSide_EI.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 6; iEnd = 6; iStepSize = 1;
test.oneside.EI = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 2 ) );
rfstats.test.oneside.EI = PlotKnockoutRFMap ( base, 15, test.oneside.EI, 15, "\nBaseline Refinement", "\nBoundary One Side E and I", iTrim, iFilter.E, iFilter.I, iColList );

		#
		#	7.  Get the Knockout RF Map - TWO-SIDE I.
		#
fRoot = "BorderKnockout_TwoSide_I.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 7; iEnd = 7; iStepSize = 1;
test.twoside.I = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 2 ) );
rfstats.test.twoside.I = PlotKnockoutRFMap ( base, 15, test.twoside.I, 15, "\nBaseline Refinement", "\nBoundary Two Sides I Only", iTrim, iFilter.E, iFilter.I, iColList );

		#
		#	8.  Get the Knockout RF Map - TWO-SIDE E.
		#
fRoot = "BorderKnockout_TwoSide_I.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 7; iEnd = 7; iStepSize = 1;
test.twoside.E = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 2 ) );
rfstats.test.twoside.E = PlotKnockoutRFMap ( base, 15, test.twoside.E, 15, "\nBaseline Refinement", "\nBoundary Two Sides E Only", iTrim, iFilter.E, iFilter.I, iColList );

		#
		#	9.  Get the Knockout RF Map - TWO-SIDE E and I.
		#
fRoot = "BorderKnockout_TwoSide_EandI.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 9; iEnd = 9; iStepSize = 1;
test.twoside.EI = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iColList = c ( controlTrack, as.integer ( N / 3 - 2 ), as.integer ( N / 3 + 2 ) );
rfstats.test.twoside.EI = PlotKnockoutRFMap ( base, 15, test.twoside.EI, 15, "\nBaseline Refinement", "\nBoundary Two Sides E and I", iTrim, iFilter.E, iFilter.I, iColList );


	#
	#	2.C:	Spatial map RF Centroid shifts.
	#
dTmp.e = RFCentroidDelta ( test.ctl.eandi$rfTrackData.e[[iTop.test.ctl]],  base$rfTrackData.e[[iTop.base]]);
dTmp.i = RFCentroidDelta ( test.ctl.eandi$rfTrackData.i[[iTop.test.ctl]],  base$rfTrackData.i[[iTop.base]]);
iFilter = ( as.integer ( N / 6 ) + 1 ) * N + iKnockOffset + 1; iFilter = seq ( iFilter, iFilter + iKnockLength - 1 );
dTmp.e[c(iTrim,iFilter)] = 0; dTmp.i[c(iTrim,iFilter)] = 0;
x11(); par(mfrow=c(2,2));
ShowVecAsMap2 ( dTmp.e, "titleText", "xLabText", "yLabText", min(dTmp.e), max(dTmp.e) );
ShowVecAsMap2 ( dTmp.i, "titleText", "xLabText", "yLabText", min(dTmp.i), max(dTmp.i) );

dTmp.e = ( QuantRFSize ( test.ctl.eandi$r1.e.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( base$r1.e.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
dTmp.i = ( QuantRFSize ( test.ctl.eandi$r1.i.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( base$r1.i.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
dTmp.e[c(iTrim,iFilter)] = 0; dTmp.i[c(iTrim,iFilter)] = 0;
ShowVecAsMap2 ( dTmp.e, "titleText", "xLabText", "yLabText", min(dTmp.e), max(dTmp.e) );
ShowVecAsMap2 ( dTmp.i, "titleText", "xLabText", "yLabText", min(dTmp.i), max(dTmp.i) );

	#
	#	3.A:	Baseline v Knockout One-Side
	#
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = "\nOne-Sided I-Type Knockout ";
iWhich = iTop.test.oneside; ShowTopoMap1 ( test.oneside$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.test.oneside; ShowTopoMap1 ( test.oneside$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

	#
	#	3.B:	Baseline v Knockout One-Side. Longitudinal Tracks.
	#
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 3 ) + 1 );
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );

subTitleText = "\nOne-Sided I-Type Knockout ";
iWhich = iTop.test.oneside; ShowThreeDigitRFTrackFlex ( test.oneside$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.oneside; ShowThreeDigitRFTrackFlex ( test.oneside$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 3 ) - 1 );
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );

subTitleText = "\nOne-Sided I-Type Knockout ";
iWhich = iTop.test.oneside; ShowThreeDigitRFTrackFlex ( test.oneside$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.oneside; ShowThreeDigitRFTrackFlex ( test.oneside$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );

controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 3 ) - 2 );
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );

subTitleText = "\nOne-Sided I-Type Knockout ";
iWhich = iTop.test.oneside; ShowThreeDigitRFTrackFlex ( test.oneside$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.oneside; ShowThreeDigitRFTrackFlex ( test.oneside$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	
	#
	#	3.C:	Spatial map RF Centroid shifts.
	#
dTmp.e = RFCentroidDelta ( test.oneside$rfTrackData.e[[iTop.test.ctl]],  base$rfTrackData.e[[iTop.base]]);
dTmp.i = RFCentroidDelta ( test.oneside$rfTrackData.i[[iTop.test.ctl]],  base$rfTrackData.i[[iTop.base]]);
iTrim = c ( 1:N, seq ( N2, N2 - N + 1, -1 ), seq ( N, N2, N ), seq ( 1, N2, N ) );	# Trim the edges.
iFilter = ( as.integer ( N / 3 ) - 1 ) * N + iKnockOffset + 1; iFilter = seq ( iFilter, iFilter + iKnockLength - 1 );
dTmp.e[c(iTrim)] = 0; dTmp.i[c(iTrim,iFilter)] = 0;
x11(); par(mfrow=c(2,2));
ShowVecAsMap2 ( dTmp.e, "titleText", "xLabText", "yLabText", min(dTmp.e), max(dTmp.e) );
ShowVecAsMap2 ( dTmp.i, "titleText", "xLabText", "yLabText", min(dTmp.i), max(dTmp.i) );

dTmp.e = ( QuantRFSize ( test.oneside$r1.e.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( base$r1.e.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
dTmp.i = ( QuantRFSize ( test.oneside$r1.i.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( base$r1.i.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
dTmp.e[c(iTrim)] = 0; dTmp.i[c(iTrim,iFilter)] = 0;
ShowVecAsMap1 ( dTmp.e, "titleText", "xLabText", "yLabText" );
ShowVecAsMap1 ( dTmp.i, "titleText", "xLabText", "yLabText" );

x11(); par(mfrow=c(2,2));
iCell = iFilter[1] + N; ShowVecAsMap1 ( base$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] + N; ShowVecAsMap1 ( base$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] + N; ShowVecAsMap1 ( test.oneside$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] + N; ShowVecAsMap1 ( test.oneside$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );

x11(); par(mfrow=c(2,2));
iCell = iFilter[1] - 1; ShowVecAsMap1 ( base$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] - 1; ShowVecAsMap1 ( base$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] - 1; ShowVecAsMap1 ( test.oneside$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] - 1; ShowVecAsMap1 ( test.oneside$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );



	#
	#	4.A:	Baseline v Knockout Two-Side I only
	#
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = "\nTwo-Sided I-Type Knockout ";
iWhich = iTop.test.twoside; ShowTopoMap1 ( test.twoside$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.test.twoside; ShowTopoMap1 ( test.twoside$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

	#
	#	4.B:	Baseline v Knockout Two-Side. Longitudinal Tracks.
	#
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 3 ) ); 	# Knocked out track
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
subTitleText = "\nTwo-Sided I-Type Knockout ";
iWhich = iTop.test.twoside; ShowThreeDigitRFTrackFlex ( test.twoside$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.twoside; ShowThreeDigitRFTrackFlex ( test.twoside$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 3 ) + 1 );	# Knocked out track
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
subTitleText = "\nTwo-Sided I-Type Knockout ";
iWhich = iTop.test.twoside; ShowThreeDigitRFTrackFlex ( test.twoside$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.twoside; ShowThreeDigitRFTrackFlex ( test.twoside$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 3 ) + 2 );
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
subTitleText = "\nTwo-Sided I-Type Knockout ";
iWhich = iTop.test.twoside; ShowThreeDigitRFTrackFlex ( test.twoside$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.twoside; ShowThreeDigitRFTrackFlex ( test.twoside$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 3 ) - 1 );
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
subTitleText = "\nTwo-Sided I-Type Knockout ";
iWhich = iTop.test.twoside; ShowThreeDigitRFTrackFlex ( test.twoside$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.twoside; ShowThreeDigitRFTrackFlex ( test.twoside$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	
	#
	#	4.C:	Spatial map RF Centroid shifts.
	#
dTmp.e = RFCentroidDelta ( test.twoside$rfTrackData.e[[iTop.test.ctl]],  base$rfTrackData.e[[iTop.base]]);
dTmp.i = RFCentroidDelta ( test.twoside$rfTrackData.i[[iTop.test.ctl]],  base$rfTrackData.i[[iTop.base]]);
iTrim = c ( 1:N, seq ( N2, N2 - N + 1, -1 ), seq ( N, N2, N ), seq ( 1, N2, N ) );	# Trim the edges.
iFilter0 = ( as.integer ( N / 3 ) - 1 ) * N + iKnockOffset + 1; iFilter0 = seq ( iFilter0, iFilter0 + iKnockLength - 1 );
iFilter1 = ( as.integer ( N / 3 ) ) * N + iKnockOffset + 1; iFilter1 = seq ( iFilter1, iFilter1 + iKnockLength - 1 );
dTmp.e[c(iTrim)] = 0; dTmp.i[c(iTrim,iFilter0,iFilter1)] = 0;
x11(); par(mfrow=c(2,2));
ShowVecAsMap2 ( dTmp.e, "titleText", "xLabText", "yLabText", min(dTmp.e), max(dTmp.e) );
ShowVecAsMap2 ( dTmp.i, "titleText", "xLabText", "yLabText", min(dTmp.i), max(dTmp.i) );

dTmp.e = ( QuantRFSize ( test.twoside$r1.e.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( base$r1.e.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
dTmp.i = ( QuantRFSize ( test.twoside$r1.i.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( base$r1.i.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
dTmp.e[c(iTrim)] = 0; dTmp.i[c(iTrim,iFilter0,iFilter1)] = 0;
ShowVecAsMap2 ( dTmp.e, "titleText", "xLabText", "yLabText", min(dTmp.e), max(dTmp.e) );
ShowVecAsMap2 ( dTmp.i, "titleText", "xLabText", "yLabText", min(dTmp.i), max(dTmp.i) );

x11(); par(mfrow=c(2,2));
iCell = iFilter0[ (iKnockLength/2) ]; ShowVecAsMap1 ( base$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter0[ (iKnockLength/2) ]; ShowVecAsMap1 ( base$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter0[ (iKnockLength/2) ]; ShowVecAsMap1 ( test.twoside$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter0[ (iKnockLength/2) ]; ShowVecAsMap1 ( test.twoside$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );

x11(); par(mfrow=c(2,2));
iCell = iFilter0[ (iKnockLength/2) ] - N; ShowVecAsMap1 ( base$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter0[ (iKnockLength/2) ] - N; ShowVecAsMap1 ( base$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter0[ (iKnockLength/2) ] - N; ShowVecAsMap1 ( test.twoside$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter0[ (iKnockLength/2) ] - N; ShowVecAsMap1 ( test.twoside$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );

x11(); par(mfrow=c(2,2));
iCell = iFilter0[ (iKnockLength/2) ] - 2*N; ShowVecAsMap1 ( base$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter0[ (iKnockLength/2) ] - 2*N; ShowVecAsMap1 ( base$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter0[ (iKnockLength/2) ] - 2*N; ShowVecAsMap1 ( test.twoside$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter0[ (iKnockLength/2) ] - 2*N; ShowVecAsMap1 ( test.twoside$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );

x11(); par(mfrow=c(2,2));
iCell = iFilter0[ (iKnockLength/2) ] + 2*N; ShowVecAsMap1 ( base$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter0[ (iKnockLength/2) ] + 2*N; ShowVecAsMap1 ( base$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter0[ (iKnockLength/2) ] + 2*N; ShowVecAsMap1 ( test.twoside$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter0[ (iKnockLength/2) ] + 2*N; ShowVecAsMap1 ( test.twoside$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );

	#
	#	5.A:	Baseline v Knockout Two-Side E and I
	#
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = "\nTwo-Sided E and I-Type Knockout ";
iWhich = iTop.test.twoside.eandi; ShowTopoMap1 ( test.twoside.eandi$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.test.twoside.eandi; ShowTopoMap1 ( test.twoside.eandi$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

	#
	#	5.B:	Baseline v Knockout Two-Side. Longitudinal Tracks.
	#
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 3 ) ); 	# Knocked out track
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
subTitleText = "\nTwo-Sided I-Type Knockout ";
iWhich = iTop.test.twoside.eandi; ShowThreeDigitRFTrackFlex ( test.twoside.eandi$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.twoside.eandi; ShowThreeDigitRFTrackFlex ( test.twoside.eandi$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 3 ) + 1 );	# Knocked out track
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
subTitleText = "\nTwo-Sided E and I-Type Knockout ";
iWhich = iTop.test.twoside.eandi; ShowThreeDigitRFTrackFlex ( test.twoside.eandi$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.twoside.eandi; ShowThreeDigitRFTrackFlex ( test.twoside.eandi$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 3 ) + 2 );
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
subTitleText = "\nTwo-Sided E and I-Type Knockout ";
iWhich = iTop.test.twoside.eandi; ShowThreeDigitRFTrackFlex ( test.twoside.eandi$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.twoside.eandi; ShowThreeDigitRFTrackFlex ( test.twoside.eandi$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 3 ) - 1 );
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
subTitleText = "\nTwo-Sided E and I-Type Knockout ";
iWhich = iTop.test.twoside.eandi; ShowThreeDigitRFTrackFlex ( test.twoside.eandi$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.twoside.eandi; ShowThreeDigitRFTrackFlex ( test.twoside.eandi$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	
	#
	#	5.C:	Spatial map RF Centroid shifts.
	#
dTmp.e = RFCentroidDelta ( test.twoside.eandi$rfTrackData.e[[iTop.test.ctl]],  base$rfTrackData.e[[iTop.base]]);
dTmp.i = RFCentroidDelta ( test.twoside.eandi$rfTrackData.i[[iTop.test.ctl]],  base$rfTrackData.i[[iTop.base]]);
iTrim = c ( 1:N, seq ( N2, N2 - N + 1, -1 ), seq ( N, N2, N ), seq ( 1, N2, N ) );	# Trim the edges.
iFilter0 = ( as.integer ( N / 3 ) - 1 ) * N + iKnockOffset + 1; iFilter0 = seq ( iFilter0, iFilter0 + iKnockLength - 1 );
iFilter1 = ( as.integer ( N / 3 ) ) * N + iKnockOffset + 1; iFilter1 = seq ( iFilter1, iFilter1 + iKnockLength - 1 );
dTmp.e[c(iTrim,iFilter0,iFilter1)] = 0; dTmp.i[c(iTrim,iFilter0,iFilter1)] = 0;
x11(); par(mfrow=c(2,2));
ShowVecAsMap2 ( dTmp.e, "titleText", "xLabText", "yLabText", min(dTmp.e), max(dTmp.e) );
ShowVecAsMap2 ( dTmp.i, "titleText", "xLabText", "yLabText", min(dTmp.i), max(dTmp.i) );

dTmp.e = ( QuantRFSize ( test.twoside.eandi$r1.e.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( base$r1.e.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
dTmp.i = ( QuantRFSize ( test.twoside.eandi$r1.i.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( base$r1.i.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
dTmp.e[c(iTrim,iFilter0,iFilter1)] = 0; dTmp.i[c(iTrim,iFilter0,iFilter1)] = 0;
ShowVecAsMap2 ( dTmp.e, "titleText", "xLabText", "yLabText", min(dTmp.e), max(dTmp.e) );
ShowVecAsMap2 ( dTmp.i, "titleText", "xLabText", "yLabText", min(dTmp.i), max(dTmp.i) );

ShowVecAsMap2 ( (abs(dTmp.e) >= pRFAreaChangeThreshold), "titleText", "xLabText", "yLabText", 0, 1 );
ShowVecAsMap2 ( (abs(dTmp.i) >= pRFAreaChangeThreshold), "titleText", "xLabText", "yLabText", 0, 1 );

x11(); par(mfrow=c(2,2));
iCell = iFilter[1] - N; ShowVecAsMap1 ( base$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] - N; ShowVecAsMap1 ( base$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] - N; ShowVecAsMap1 ( test.oneside$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] - N; ShowVecAsMap1 ( test.oneside$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );

x11(); par(mfrow=c(2,2));
iCell = iFilter[1] + 2 * N; ShowVecAsMap1 ( base$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] + 2 * N; ShowVecAsMap1 ( base$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] + 2 * N; ShowVecAsMap1 ( test.oneside$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] + 2 * N; ShowVecAsMap1 ( test.oneside$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );



iFilter = ( as.integer ( N / 6 ) + 1 ) * N + iKnockOffset + 1; iFilter = seq ( iFilter, iFilter + iKnockLength - 1 );



		
