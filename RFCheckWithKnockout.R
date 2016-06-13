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

fDir = "D:/NMLab/Working/S.45.7B.Knockout.4/";

fName = "BorderKnockout_OneSide_E.RFMap.Raw.4.5.bin";

fName = paste ( fDir, fName, sep="" );
yLabText="Bot: Skin Layer.  Mid: Cortical I.  Top: Cortical E."; xLabText="Distal -> Proximal";
numColors = 128;

	#
	#	Spatial plots
	#
layerMarks = c(N+0.5, 2*N+0.5);
digitBoundaryMarks = c ( as.integer(N/3)+0.5, as.integer(N/3)*2+0.5 );
boundaryMarkSet = c ( digitBoundaryMarks, N+digitBoundaryMarks, 2*N+digitBoundaryMarks );
minZ = 0; maxZ = 1;
pCut = 0.02;
nCount = 0;

iShowTheseRows = seq ( 15, 30 );
iShowTheseCells = c ( 13*N+iShowTheseRows, 14*N+iShowTheseRows, 15*N+iShowTheseRows, 16*N+iShowTheseRows );
iShowTheseCells = c ( 13*N+iShowTheseRows, 14*N+iShowTheseRows );
iShowTheseCells = c ( 640, 651, 652, 653 );

saveHTML ( { 
	finfo = file.info ( fName ); toread = file ( fName, "rb" );
	for ( iProbeCell in 1:max(iShowTheseCells) ) {
		alldata = readBin ( toread, "double", size=8, n=(3*numValsPerRFTrial), endian="little" );
		if ( sum( iShowTheseCells == iProbeCell ) > 0 ) {
			nCount = nCount + 1;
			startOffset.e = 1; startOffset.i = numValsPerRFTrial + 1;
			startOffset.0 = 2 * numValsPerRFTrial + 1;
			v1e = alldata [ (startOffset.e):(startOffset.e + numValsPerRFTrial - 1) ];
			v1i = alldata [ (startOffset.i):(startOffset.i + numValsPerRFTrial - 1) ];
			v0 = alldata [ (startOffset.0):(startOffset.0 + numValsPerRFTrial - 1) ];
			r1e = sigmoid ( v1e, 4 ); r1i = sigmoid ( v1i, 4 ); r0 = sigmoid ( v0, 4 );
			for ( iTime in seq(90,150) ) {

				titleText = paste ( "Response at Each Layer to RF Probe # ", iProbeCell, "\nTime=", iTime*deltaT, sep="" );
				
				w = c ( r0[(iTime - 1) * N2 + seq ( 1:N2 )], r1i[(iTime - 1) * N2 + seq ( 1:N2 )],
					r1e[(iTime - 1) * N2 + seq ( 1:N2 )] );

				image.plot ( c(1:N), c(1:(3*N)), matrix ( log10(w), nrow=N, ncol=3*N, byrow=FALSE ),
					col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
					zlim = c ( log10(minZ+1e-6), log10(maxZ) ), main=titleText, xlab=xLabText, ylab=yLabText );

				abline(h=layerMarks, col=1, lty=1 );
				abline(h=c(boundaryMarkSet), col=5, lty=2, lwd=0.25 );

				Sys.sleep(1.0);

			} # for ( iTime in seq(1,numItersPerTrial) ) {
		} # if ( runif(1) <= pCut )
	} # for ( iCell in 1:N^2 )
	close ( toread );
}, interval = deltaT );  # saveGIF

	#
	#	Time series plots
	#

for ( iCell in seq( 665, 66 ) ) {
	x11(); par(mfcol=c(3,2));

	tRange = 1:500;
	iCellWhich = seq ( iCell, numValsPerRFTrial, N2 )[tRange];
	ylim = c ( min( v1e[iCellWhich], v1i[iCellWhich], v0[iCellWhich] ), max( v1e[iCellWhich], v1i[iCellWhich], v0[iCellWhich] ) );
	
	#QuickCheckTimeSeries ( v1e[iCellWhich], v1i[iCellWhich], v0[iCellWhich], iCell, deltaT, FALSE, ylim );
	titleText = ylabText = "Avg. Membrane Potential";
	QuickCheckTimeSeries2 ( tRange, v1e[iCellWhich], v1i[iCellWhich], v0[iCellWhich], iCell, deltaT, FALSE, titleText, ylabText );

	ylim = c ( 0, 1 );
	
	#QuickCheckTimeSeries ( r1e[iCellWhich], r1i[iCellWhich], r0[iCellWhich], iCell, deltaT, FALSE, ylim );
	titleText = ylabText = "Avg. Spiking Rate";
	QuickCheckTimeSeries2 ( tRange, r1e[iCellWhich], r1i[iCellWhich], r0[iCellWhich], iCell, deltaT, FALSE, titleText, ylabText );

	Sys.sleep(1.0);

} # for ( iCell in 1:N^2 )

x11(); par(mfrow=c(3,1));
plot ( r1e[iCellWhich], type="l" );
plot ( r1i[iCellWhich], type="l" );
plot ( r0[iCellWhich], type="l" );

saveHTML ( { 
	finfo = file.info ( fName ); toread = file ( fName, "rb" );
	for ( iProbeCell in 1:649 ) {
		alldata = readBin ( toread, "double", size=8, n=(3*numValsPerRFTrial), endian="little" );
		#if ( GoodToShow ( N, iProbeCell ) && (runif(1) <= pCut) ) {
		if ( iProbeCell == 649 ) {
			nCount = nCount + 1;
			startOffset.e = 1; startOffset.i = numValsPerRFTrial + 1;
			startOffset.0 = 2 * numValsPerRFTrial + 1;
			v1e = alldata [ (startOffset.e):(startOffset.e + numValsPerRFTrial - 1) ];
			v1i = alldata [ (startOffset.i):(startOffset.i + numValsPerRFTrial - 1) ];
			v0 = alldata [ (startOffset.0):(startOffset.0 + numValsPerRFTrial - 1) ];
			minZ = min ( v0, v1i, v1e ); maxZ = max ( v0, v1i, v1e );
			#r1e = sigmoid ( v1e, 4 ); r1i = sigmoid ( v1i, 4 ); r0 = sigmoid ( v0, 4 );
			for ( iTime in seq(1,numItersPerTrial) ) {
				titleText = paste ( "Response at Each Layer to RF Probe # ", iProbeCell, "\nTime=", iTime*deltaT, sep="" );
				w = c ( v0[(iTime - 1) * N2 + seq ( 1:N2 )], v1i[(iTime - 1) * N2 + seq ( 1:N2 )],
					v1e[(iTime - 1) * N2 + seq ( 1:N2 )] );
				image.plot ( c(1:N), c(1:(3*N)), matrix ( w, nrow=N, ncol=3*N, byrow=FALSE ),
					col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
				zlim = c ( minZ, maxZ ), main=titleText, xlab=xLabText, ylab=yLabText );
				abline(h=layerMarks, col=1, lty=1 );
				abline(h=c(), col=1, lty=1 );
			} # for ( iTime in seq(1,numItersPerTrial) ) {
		} # if ( runif(1) <= pCut )
	} # for ( iCell in 1:N^2 )
	close ( toread );
}, interval = deltaT );  # saveGIF

	#
	#	Grab the weights.
	#

fRoot = "SysInit.25.0.w1E0"; init.w1.e0 = e0 = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = "SysInit.25.0.w1EE"; init.w1.ee = ee = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = "SysInit.25.0.w1EI"; init.w1.ei = ei = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = "SysInit.25.0.w1IE"; init.w1.ie = ie = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
tmp = GetInputWeights ( e0, 1 );



	#	Play area....
			x11(); par(mfrow=c(2,2));

			for ( iTime in seq(251,350) ) {

				titleText = paste ( "Response at Each Layer to RF Probe # ", iProbeCell, "\nTime=", iTime*deltaT, sep="" );
				
				w = c ( r0[(iTime - 1) * N2 + seq ( 1:N2 )], r1i[(iTime - 1) * N2 + seq ( 1:N2 )],
					r1e[(iTime - 1) * N2 + seq ( 1:N2 )] );

				par(mfrow=c(2,2));
				image.plot ( c(1:N), c(1:(N)), matrix ( (r1e[(iTime - 1) * N2 + seq ( 1:N2 )]), nrow=N, ncol=N, byrow=FALSE ),
					col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
					main=titleText, xlab=xLabText, ylab=yLabText );

				image.plot ( c(1:N), c(1:(N)), matrix ( (r1i[(iTime - 1) * N2 + seq ( 1:N2 )]), nrow=N, ncol=N, byrow=FALSE ),
					col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
					main=titleText, xlab=xLabText, ylab=yLabText );

				image.plot ( c(1:N), c(1:(N)), matrix ( (r0[(iTime - 1) * N2 + seq ( 1:N2 )]), nrow=N, ncol=N, byrow=FALSE ),
					col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
					main=titleText, xlab=xLabText, ylab=yLabText );


				abline(h=layerMarks, col=1, lty=1 );
				abline(h=c(boundaryMarkSet), col=5, lty=2, lwd=0.25 );

				Sys.sleep(1.0);

			} # for ( iTime in seq(1,numItersPerTrial) ) {


