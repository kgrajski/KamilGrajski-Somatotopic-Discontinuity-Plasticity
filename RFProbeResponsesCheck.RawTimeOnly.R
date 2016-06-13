#
#	RFProbeResponsesCheck.R
#

#
#	Use this in conjunction with RF Probe Raw data files from CPP
#

#
#	Quick visual check on time series generated during RF mapping
#

#
#	Eventually need to develop some automated routines to make sure
#	there isn't funny stuff happening.
#

#{ rm(list = ls()); }

#
#	Helper Functions
#

source ( "NMHelperFunctions.R" );
library (animation);

#
#	Three-panel time series plot: r1E, r1I, r0.
#

QuickCheckTimeSeries = function ( r1E, r1I, rZ, iCell, deltaT, plotFlag, ylim ) {

	tVal = 1:(length(r1E)) * deltaT;

	if ( plotFlag ) {
		x11();
		par(mfrow=c(3,1));
	} # if ( plotFlag )		

	plot( tVal, r1E, type="p", col=1, ylim=ylim, ylab="Firing Rate", xlab="Time (sec)");
	title(main=paste("Cortical Layer E-Cell Firing Rate",paste("Cell ",iCell),sep='\n'));

	plot( tVal, r1I, type="p", col=1, ylim=ylim, ylab="Firing Rate", xlab="Time (sec)");
	title(main=paste("Cortical Layer I-Cell Firing Rate",paste("Cell ",iCell),sep='\n'));

	plot( tVal, rZ, type="p", col=1, ylim=ylim, ylab="Firing Rate", xlab="Time (sec)");
	title(main=paste("Input Layer Firing Rate",paste("Cell ",iCell),sep='\n'));

} # QuickCheckTimeSeries = function ( r1E, r1I, r0, cellID, deltaT, ... ) {

QuickCheckTimeSeries2 = function ( tRange, r1E, r1I, rZ, iCell, deltaT, plotFlag, titleText, ylabText ) {

	tVal = tRange * deltaT;

	if ( plotFlag ) {
		x11();
		par(mfrow=c(3,1));
	} # if ( plotFlag )		

	plot( tVal, r1E, type="p", col=1, ylab=ylabText, xlab="Time (sec)");
	title ( main = paste("Cortical Layer E-Cell", titleText, paste("Cell ",iCell),sep='\n'));

	plot( tVal, r1I, type="p", col=1, ylab=ylabText, xlab="Time (sec)");
	title(main=paste("Cortical Layer I-Cell", titleText, paste("Cell ",iCell),sep='\n'));

	plot( tVal, rZ, type="p", col=1, ylab=ylabText, xlab="Time (sec)");
	title(main=paste("Input Layer Cell", titleText, paste("Cell ",iCell),sep='\n'));

} # QuickCheckTimeSeries2 = function ( r1E, r1I, r0, cellID, deltaT, ... ) {

GoodToShow = function ( N, iCandidate ) {

	candidateList = TrimEdgesFromCellList ( seq ( 1, N2, 1 ), N, 3 );
	show = FALSE;
	if ( sum ( candidateList == iCandidate ) ) { show = TRUE; }

	#if ( iCandidate <= N ) { show = FALSE; }
	#if ( iCandidate > (N^2 - N) ) { show = FALSE; }
	#if ( iCandidate %% N == 0 ) { show = FALSE; }
	#if ( iCandidate %% N == 1 ) { show = FALSE; }

	return ( show );

} # GoodToShow = function ( N, iCandidate ) {

#
#	Main
#

source("NMHelperFunctions.R");

	#
	#	Set some constants
	#
N = 45;
deltaT = 0.001;
trialDuration = 0.500;

N2 = N * N;
numItersPerTrial = 500;
numValsPerRFTrial = numItersPerTrial * N * N;

	#
	#	Read the RF Probe Raw Data
	#

#fName = "D:/NMLab/Keep - Syndactyly Exps/T.7/RF Raw Output For Completed Baseline v1.0/Base.RFMap.Raw.15.1.bin";
fName = "F:/NMLab/Working/S.45.7.Lesion.4/BorderKnockout_Control_Placebo.Raw.4.0.bin";

	#
	#	Spatial plots
	#
digitBoundaryMarks = c ( as.integer(N/3)+0.5, as.integer(N/3)*2+0.5 );
boundaryMarkSet = c ( digitBoundaryMarks, N+digitBoundaryMarks, 2*N+digitBoundaryMarks );
minZ = 0; maxZ = 1;

iShowTheseRows = seq ( 15, 30 );
iShowTheseCells = c ( 1, 13*N+iShowTheseRows, 14*N+iShowTheseRows, 15*N+iShowTheseRows, 16*N+iShowTheseRows );
iShowTheseCells = c ( 1, 13*N+iShowTheseRows, 16*N+iShowTheseRows );


iTargetCell = 652;

finfo = file.info ( fName ); toread = file ( fName, "rb" );
for ( iProbeCell in 1:max(iShowTheseCells) ) {
	alldata = readBin ( toread, "double", size=8, n=(3*numValsPerRFTrial), endian="little" );

	if ( sum(iShowTheseCells == iProbeCell) ) {
		startOffset.e = 1; startOffset.i = numValsPerRFTrial + 1;
		startOffset.0 = 2 * numValsPerRFTrial + 1;
		v1e = alldata [ (startOffset.e):(startOffset.e + numValsPerRFTrial - 1) ];
		v1i = alldata [ (startOffset.i):(startOffset.i + numValsPerRFTrial - 1) ];
		v0 = alldata [ (startOffset.0):(startOffset.0 + numValsPerRFTrial - 1) ];
		r1e = sigmoid ( v1e, 4 ); r1i = sigmoid ( v1i, 4 ); r0 = sigmoid ( v0, 4 );

		x11(); par(mfcol=c(3,2));

		iCellWhich = seq ( iTargetCell, numValsPerRFTrial, N2 );

		ylim = c ( min( v1e[iCellWhich], v1i[iCellWhich], v0[iCellWhich] ), max( v1e[iCellWhich], v1i[iCellWhich], v0[iCellWhich] ) );
		titleText = ylabText = "Avg. Membrane Potential";
		QuickCheckTimeSeries2 ( seq(1,numItersPerTrial), v1e[iCellWhich], v1i[iCellWhich], v0[iCellWhich], iTargetCell, deltaT, FALSE, titleText, ylabText );
		ylim = c ( 0, 1 );	
		titleText = ylabText = "Avg. Spiking Rate";
		QuickCheckTimeSeries2 ( seq(1,numItersPerTrial), r1e[iCellWhich], r1i[iCellWhich], r0[iCellWhich], iTargetCell, deltaT, FALSE, titleText, ylabText );

		Sys.sleep(1.0);
	} # 	if ( sum(iShowTheseCells == iProbeCell) ) {

} # for ( iProbeCell in 1:max(iShowTheseCells) ) {
close ( toread );

