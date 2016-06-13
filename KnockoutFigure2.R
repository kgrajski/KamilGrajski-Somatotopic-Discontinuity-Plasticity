#
#	KnockoutFigure2.R
#

#
#	Automatically handling the control, one-side and two-side I-cell knockout variations.
#

#	Standardized plots.

#	There may be further consolidation into summary figures for publication.

#	Clear the workspace.
rm(list = ls());

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
		#	Some additional variables needed if/when NMExpZone isn't the source.
kRFPeakToEdgeDetect = 0.50;
N = 45;
N2 = N * N;
kRFPeakToEdgeDetect = 0.5;
pRFAreaChangeThreshold = 0.10;

		#	Set the directory where the experimental data is sitting.
fDir = "C:/Users/KAG/Documents/NMLab/Simulations/S.2.12/"
fDir = "D:/NMLab/S.45.7.Knockout/";

		#
		#	Get the Baseline Refinement RF Map.
		#
fRoot = "Base.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iBase = 15; iStart = 15; iEnd = 15; iStepSize = 10;
base = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iTop.base = length ( base$rfTrackData.e );

iKnockOffset = as.integer ( N / 2 ) - 4;
iBase = iKnockLength = 12;

		#
		#	Get the Placebo
		#
fRoot = "BorderKnockout_Control_Placebo.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 0; iEnd = 0; iStepSize = 1;
placebo = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iTop.placebo = length ( placebo $rfTrackData.e );

		#
		#	Get the Knockout RF Map - CONTROL I Only.
		#
fRoot = "BorderKnockout_Control_I.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 1; iEnd = 1; iStepSize = 1;
test.ctl = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iTop.test.ctl = length ( test.ctl$rfTrackData.e );

		#
		#	Get the Knockout RF Map - CONTROL E and I.
		#
fRoot = "BorderKnockout_Control_EandI.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 2; iEnd = 2; iStepSize = 1;
test.ctl.eandi = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iTop.test.ctl.eandi = length ( test.ctl.eandi$rfTrackData.e );

		#
		#	Get the Knockout RF Map - ONE-SIDE.
		#
fRoot = "BorderKnockout_OneSide_I.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 3; iEnd = 3; iStepSize = 3;
test.oneside = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iTop.test.oneside = length ( test.oneside$rfTrackData.e );

		#
		#	Get the Knockout RF Map - TWO-SIDE.
		#
fRoot = "BorderKnockout_TwoSide_I.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 0; iEnd = 0; iStepSize = 1;
test.twoside = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iTop.test.twoside = length ( test.twoside$rfTrackData.e );

		#
		#	Get the Knockout RF Map - TWO-SIDE.
		#
fRoot = "BorderKnockout_TwoSide_EandI.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iStart = 0; iEnd = 0; iStepSize = 1;
test.twoside.eandi = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iTop.test.twoside.eandi  = length ( test.twoside.eandi $rfTrackData.e );

	#
	#	0: Placebo Run.
	#

	#
	#	0.A:	Baseline v Placebo
	#
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = "\nControl I-Type Knockout ";
iWhich = iTop.placebo; ShowTopoMap1 ( placebo$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.placebo; ShowTopoMap1 ( placebo$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
	
	#
	#	0.B:	Baseline v Placebo. Longitudinal Tracks.
	#
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 6 + 1 ) );
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );

subTitleText = "\nControl I-Type Knockout ";
iWhich = iTop.placebo; ShowThreeDigitRFTrackFlex ( placebo$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.placebo; ShowThreeDigitRFTrackFlex ( placebo$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
		
	#
	#	0.C:	Baseline v Placebo Spatial map RF Centroid and RF size shifts.
	#
dTmp.e = RFCentroidDelta ( placebo$rfTrackData.e[[iTop.test.ctl]],  base$rfTrackData.e[[iTop.base]]);
dTmp.i = RFCentroidDelta ( placebo$rfTrackData.i[[iTop.test.ctl]],  base$rfTrackData.i[[iTop.base]]);
iTrim = c ( 1:N, seq ( N2, N2 - N + 1, -1 ), seq ( N, N2, N ), seq ( 1, N2, N ) );	# Trim the edges.
iFilter = ( as.integer ( N / 6 ) + 1 ) * N + iKnockOffset + 1; iFilter = seq ( iFilter, iFilter + iKnockLength - 1 );
dTmp.e[c(iTrim,iFilter)] = 0; dTmp.i[c(iTrim,iFilter)] = 0;
x11(); par(mfrow=c(2,2));
ShowVecAsMap1 ( dTmp.e, "titleText", "xLabText", "yLabText" );
ShowVecAsMap1 ( dTmp.i, "titleText", "xLabText", "yLabText" );

dTmp.e = ( QuantRFSize ( placebo$r1.e.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( base$r1.e.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
dTmp.i = ( QuantRFSize ( placebo$r1.i.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( base$r1.i.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
dTmp.e[c(iTrim,iFilter)] = 0; dTmp.i[c(iTrim,iFilter)] = 0;
ShowVecAsMap1 ( dTmp.e, "titleText", "xLabText", "yLabText" );
ShowVecAsMap1 ( dTmp.i, "titleText", "xLabText", "yLabText" );

x11(); par(mfrow=c(2,2));
iCell = iFilter[1] - 1; ShowVecAsMap1 ( base$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] - 1; ShowVecAsMap1 ( base$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] - 1; ShowVecAsMap1 ( placebo$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[1] - 1; ShowVecAsMap1 ( placebo$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );

	#
	#	1.A:	Baseline v Knockout Control I Only
	#
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = "\nControl I-Type Knockout ";
iWhich = iTop.test.ctl; ShowTopoMap1 ( test.ctl$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.test.ctl; ShowTopoMap1 ( test.ctl$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
	
	#
	#	1.B:	Baseline v Knockout Control I Only. Longitudinal Tracks.
	#
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 6 + 1 ) );
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );

subTitleText = "\nControl I-Type Knockout ";
iWhich = iTop.test.ctl; ShowThreeDigitRFTrackFlex ( test.ctl$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.ctl; ShowThreeDigitRFTrackFlex ( test.ctl$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
		
	#
	#	1.C:	Spatial map RF Centroid and RF size shifts.
	#
dTmp.e = RFCentroidDelta ( test.ctl$rfTrackData.e[[iTop.test.ctl]],  base$rfTrackData.e[[iTop.base]]);
dTmp.i = RFCentroidDelta ( test.ctl$rfTrackData.i[[iTop.test.ctl]],  base$rfTrackData.i[[iTop.base]]);
iTrim = c ( 1:N, seq ( N2, N2 - N + 1, -1 ), seq ( N, N2, N ), seq ( 1, N2, N ) );	# Trim the edges.
iFilter0 = ( as.integer ( N / 6 ) + 1 ) * N + iKnockOffset + 1; iFilter0 = seq ( iFilter0, iFilter0 + iKnockLength - 1 );
iFilter1 = ( as.integer ( N / 6 ) ) * N + iKnockOffset + 1; iFilter1 = seq ( iFilter1, iFilter1 + iKnockLength - 1 );
dTmp.e[c(iTrim)] = 0; dTmp.i[c(iTrim,iFilter0,iFilter1)] = 0;
x11(); par(mfrow=c(2,2));
ShowVecAsMap2 ( dTmp.e, "titleText", "xLabText", "yLabText", min(c(dTmp.e, dTmp.i)), max(c(dTmp.e, dTmp.i)) );
ShowVecAsMap2 ( dTmp.i, "titleText", "xLabText", "yLabText", min(c(dTmp.e, dTmp.i)), max(c(dTmp.e, dTmp.i)) );

dTmp.e = ( QuantRFSize ( test.ctl$r1.e.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( base$r1.e.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
dTmp.i = ( QuantRFSize ( test.ctl$r1.i.rfMap, kRFPeakToEdgeDetect ) / QuantRFSize ( base$r1.i.rfMap, kRFPeakToEdgeDetect ) ) - 1.0;
dTmp.e[c(iTrim)] = 0; dTmp.i[c(iTrim,iFilter0,iFilter1)] = 0;
ShowVecAsMap1 ( dTmp.e, "titleText", "xLabText", "yLabText" );
ShowVecAsMap1 ( dTmp.i, "titleText", "xLabText", "yLabText" );

x11(); par(mfrow=c(2,2));
iCell = iFilter[iKnockLength/2] - N; ShowVecAsMap1 ( base$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[iKnockLength/2] - N; ShowVecAsMap1 ( base$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[iKnockLength/2] - N; ShowVecAsMap1 ( test.ctl$r1.e.rfMap[iCell,], "titleText", "xLabText", "yLabText" );
iCell = iFilter[iKnockLength/2] - N; ShowVecAsMap1 ( test.ctl$r1.i.rfMap[iCell,], "titleText", "xLabText", "yLabText" );

	#
	#	2.A:	Baseline v Knockout Control E and I Type
	#
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = "\nControl E and I-Type Knockout ";
iWhich = iTop.test.ctl.eandi; ShowTopoMap1 ( test.ctl.eandi$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.test.ctl.eandi; ShowTopoMap1 ( test.ctl.eandi$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
	
	#
	#	2.B:	Baseline v Knockout Control E and I Type. Longitudinal Tracks.
	#
controlTrack = as.integer ( 2 * N / 3);
iColList = c ( controlTrack, as.integer ( N / 6 + 2 ) );
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );

subTitleText = "\nControl E and I-Type Knockout ";
iWhich = iTop.test.ctl.eandi; ShowThreeDigitRFTrackFlex ( test.ctl.eandi$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
iWhich = iTop.test.ctl.eandi; ShowThreeDigitRFTrackFlex ( test.ctl.eandi$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList );
	
	#
	#	2.C:	Spatial map RF Centroid shifts.
	#
dTmp.e = RFCentroidDelta ( test.ctl.eandi$rfTrackData.e[[iTop.test.ctl]],  base$rfTrackData.e[[iTop.base]]);
dTmp.i = RFCentroidDelta ( test.ctl.eandi$rfTrackData.i[[iTop.test.ctl]],  base$rfTrackData.i[[iTop.base]]);
iTrim = c ( 1:N, seq ( N2, N2 - N + 1, -1 ), seq ( N, N2, N ), seq ( 1, N2, N ) );	# Trim the edges.
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







		
