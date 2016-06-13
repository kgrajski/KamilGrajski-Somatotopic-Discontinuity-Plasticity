#
#	RFProbeResponsesCheck.Time.R
#

#
#	Use this in conjunction wxlaith RF Probe Raw data files from CPP
#

#
#	Check time series generated during RF mapping in knockout vs baseline
#

{ rm(list = ls()); }

	#
	#	Helper Functions
	#

source ( "NMHelperFunctions.R" );
library (animation);

ExtractRFProbeRawData = function ( fName, N2, numValsPerRFTrial, iCellList ) {

	v1e = matrix ( 0, nrow=numValsPerRFTrial, ncol=length(iCellList) );
	v1i = matrix ( 0, nrow=numValsPerRFTrial, ncol=length(iCellList) );
	v0 = matrix ( 0, nrow=numValsPerRFTrial, ncol=length(iCellList) );
	iCount = 1;

	finfo = file.info ( fName ); toread = file ( fName, "rb" );
	for ( iProbeCell in 1:max(iCellList) ) {
		alldata = readBin ( toread, "double", size=8, n=(3*numValsPerRFTrial), endian="little" );
		if ( sum( iCellList == iProbeCell ) > 0 ) {
			startOffset.e = 1; startOffset.i = numValsPerRFTrial + 1;
			startOffset.0 = 2 * numValsPerRFTrial + 1;
			v1e[,iCount] = alldata [ (startOffset.e):(startOffset.e + numValsPerRFTrial - 1) ];
			v1i[,iCount] = alldata [ (startOffset.i):(startOffset.i + numValsPerRFTrial - 1) ];
			v0[,iCount] = alldata [ (startOffset.0):(startOffset.0 + numValsPerRFTrial - 1) ];
		} # if ( sum( iCellList == iProbeCell ) > 0 ) {
	} # for ( iProbeCell in 1:max(iCellList) ) {
	close ( toread );

	return ( list ( v1e=v1e, v1i=v1i, v0=v0 ) );

} # ExtractRFProbeRawData = function ( fName, N2, numValsPerRFTrial ) {

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
#	Extract the vt values for the given iCell from the CPP generated RF Raw data file.
#	Compute the rt values.
#
TSResponse = function ( N2, numValsPerRFTrial, iCell, xt, eScale ) {

	iCellList = seq ( iCell, numValsPerRFTrial, N2 );
	startOffset.e = 1; startOffset.i = numValsPerRFTrial + 1;
	startOffset.0 = 2 * numValsPerRFTrial + 1;

	vt.e = ( xt [ (startOffset.e):(startOffset.e + numValsPerRFTrial - 1) ] ) / eScale;
	vt.e = vt.e[iCellList];

	vt.i = xt [ (startOffset.i):(startOffset.i + numValsPerRFTrial - 1) ];
	vt.i = vt.i[iCellList];

	vt.0 = xt [ (startOffset.0):(startOffset.0 + numValsPerRFTrial - 1) ];
	vt.0 = vt.0[iCellList];


	rt.e = sigmoid ( vt.e, 4 );
	rt.i = sigmoid ( vt.i, 4 );
	rt.0 = sigmoid ( vt.0, 4 );

	vt.min = min ( vt.e, vt.i, vt.0 );
	vt.max = max ( vt.e, vt.i, vt.0 );
	rt.min = min ( rt.e, rt.i, rt.0 );
	rt.max = max ( rt.e, rt.i, rt.0 );

	return ( list ( vt.e=vt.e, vt.i=vt.i, vt.0=vt.0, vt.min=vt.min, vt.max=vt.max,
			rt.e=rt.e, rt.i=rt.i, rt.0=rt.0, rt.min=rt.min, rt.max=rt.max ) );

} # TSResponse = function ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale ) {


#
#	Prepare a baseline vs experimental dataset for the given iCell
#
ExpRefTSResponse = function ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale ) {

	xt.base = TSResponse ( N2, numValsPerRFTrial, iCell, alldata.base, eScale );
	xt.exp = TSResponse ( N2, numValsPerRFTrial, iCell, alldata.exp, eScale );
	
	v1e.delta = xt.exp$vt.e - xt.base$vt.e;
	v1i.delta = xt.exp$vt.i - xt.base$vt.i;
	v0.delta = xt.exp$vt.0 - xt.base$vt.0;
	
	vt.min = min ( xt.base$vt.min, xt.exp$vt.min ); vt.max = max ( xt.base$vt.max, xt.exp$vt.max );
	delta.min = min ( v1e.delta, v1i.delta, v0.delta ); delta.max = max ( v1e.delta, v1i.delta, v0.delta );

	return ( list ( v1e.base=xt.base$vt.e, v1i.base=xt.base$vt.i, v0.base=xt.base$vt.0,
				v1e.exp=xt.exp$vt.e, v1i.exp=xt.exp$vt.i, v0.exp=xt.exp$vt.0,
				v1e.delta=v1e.delta, v1i.delta=v1i.delta, v0.delta=v0.delta,
				vt.min=vt.min, vt.max=vt.max, delta.min=delta.min, delta.max=delta.max ) );

} # TSResponse = function ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale ) {

#
#
#
KnockoutTimeSeriesPlot = function ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleTextRoot, subTitleTextRoot, xlabTextRoot, ylabTextRoot ) {

	x11(); par ( mfcol=c(3,4) );

	iRecordCellsList = c ( iProbeCell, iProbeCell + N, iProbeCell + 2*N, iProbeCell + 3*N );

	for ( iCell in iRecordCellsList ) {

		tmp = ExpRefTSResponse ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale );

		titleText = paste ( titleTextRoot, "\nStim: ", iProbeCell, "Rec: ", iCell );
		xlabText = xlabTextRoot;
		ylabText = ylabTextRoot;

		ylim = c(tmp$vt.min, tmp$vt.max);
		ylim = c(-2, 1);
		plot ( tmp$v1e.exp, type="l", ylim=ylim, col=3, lty=1,
			main=titleText, xlab=xlabText, ylab=ylabText );
		lines ( tmp$v1i.exp, col=2, lty=1 );
		lines ( tmp$v0.exp, col=1, lty=1 );

		titleText = paste ( "Baseline Refined", "\nStimLoc: ", iProbeCell, "Record: ", iCell );
		plot ( tmp$v1e.base, type="l", ylim=ylim, col=3, lty=1,
			main=titleText, xlab=xlabText, ylab=ylabText );
		lines ( tmp$v1i.base, col=2, lty=1 );
		lines ( tmp$v0.base, col=1, lty=1 );

		titleText = paste ( "Diff. Exp - Control", "\nStimLoc: ", iProbeCell, "Record: ", iCell );
		plot ( tmp$v1e.delta, type="l", ylim=c(tmp$delta.min, tmp$delta.max), col=3, lty=1,
			main=titleText, xlab=xlabText, ylab=ylabText );
		lines ( tmp$v1i.delta, col=2, lty=1 );
		lines ( tmp$v0.delta, col=1, lty=1 );

	} # for ( iCell in iRecordCellsList ) {

} # KnockoutTimeSeriesPlot = function ( N, numValsPerRFTrial, iProbeCell...

#
#
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

#
#
#
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

#
#
#
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
numItersPerTrial = 350;
numValsPerRFTrial = numItersPerTrial * N * N;
titleText = "Two-Sided I-Only Knockout";
subTitleTextRoot = "TBD";
xlabTextRoot = "Msec.";
ylabTextRoot = "v(t)";
eScale = 1.0;

	#
	#	Read the RF Probe Raw Data
	#
fDir = "D:/NMLab/S.45.7.Knockout/WithSelectedRawRF/";
fDir = "D:/NMLab/S.45.7.Knockout/WithSelectedRawRF/";
fName.exp = "BorderKnockout_Control_I.RFMap.Raw.4.1.bin";
#fName.base = "BorderKnockout_Control_Placebo.Raw.4.0.bin";
fName.base = "Base.RFMap.Raw.15.1.bin";

fName.base = paste ( fDir, fName.base, sep="" );
fName.exp = paste ( fDir, fName.exp, sep="" );

	#
	#	Time Series: Exp vs Base plots
	#
finfo.base = file.info ( fName.base );
finfo.exp = file.info ( fName.exp );

toread.base = file ( fName.base, "rb" );
toread.exp = file ( fName.exp, "rb" );

tmp = seq ( 685-N, 710-N, 1 );
iProbeCellList = c ( tmp, tmp-45, tmp-(3*45), tmp+45, tmp+(2*45), tmp+(3*45) );
iProbeCellList = tmp;
maxProbeCell = max ( iProbeCellList );

for ( iProbeCell in 1:maxProbeCell ) {

	alldata.base = readBin ( toread.base, "double", size=8, n=(3*numValsPerRFTrial), endian="little" );
	alldata.exp = readBin ( toread.exp, "double", size=8, n=(3*numValsPerRFTrial), endian="little" );

	if ( sum(iProbeCellList == iProbeCell) ) {

		KnockoutTimeSeriesPlot ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot );
		
	} # if ( runif(1) <= pCut )

} # for ( iProbeCell in ...

close ( toread.base );
close ( toread.exp );

