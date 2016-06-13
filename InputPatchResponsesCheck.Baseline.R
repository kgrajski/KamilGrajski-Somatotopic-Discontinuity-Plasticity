#
#	InputPatchResponsesCheck.R
#

#
#	Use this in conjunction with RF Probe Raw data files from CPP
#

#
#	Quick visual check on time series generated during "training"
#

#
#	Eventually need to develop some automated routines to make sure
#	there isn't funny stuff happening.
#

rm(list = ls());
library(animation);

#
#	Helper Functions
#

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


#
#	Main
#

source("NMHelperFunctions.R");

	#	Set some constants
N = 33;
N2 = N * N;
deltaT = 0.005;
trialDuration = 0.350;
numItersPerTrial = as.integer ( trialDuration * 1 / deltaT );
numValsPerTrial = numItersPerTrial * N * N;

	#	Expected file size (for sanity checking)
	#		For 3-digit baseline:
digitWidth = N / 3;
patchSize = 7;
numPatches.baseline = 3 * ( digitWidth - patchSize ) * ( N - patchSize );
expectedLen.baseline = 3 * numValsPerTrial * numPatches.baseline;

		#
		#	Spatial plots
		#
startIter = 1;
endIter = numItersPerTrial;
layerMarks = c(N+0.5, 2*N+0.5);
startOffset.e = 1;
startOffset.i = numValsPerTrial + 1;
startOffset.0 = 2 * numValsPerTrial + 1;
yLabText="Bot: Skin Layer.  Mid: Cortical I.  Top: Cortical E."; xLabText="Distal -> Proximal";
zLim = c(0, 1);
numColors = 128;

	#
	#	Read some file information.  Typical sized files will be
	#	too big to read in.  Need to read in just what will be shown.
	#
fDir = "E:/NMLab/Simulations/Syndactyly/S.33.7.5.100/";
fDir = "E:/NMLab/Simulations/Syndactyly_5.1/S.33.7.25.100//";
fDir = "E:/NMLab/Simulations/Syndactyly_0.0025_0.00025/S.33.1.100.100/"
#fDir = "E:/NMLab/Simulations/NM5.2 Test Area/"
fDir = "E:/NMLab/Simulations/NM5.3_Test_Area/M.050.1.0//";
fDir = "D:/NMLab/Simulations/Syndactyly/N.100.1.0//";
fDir = "E:/NMLab/Simulations/NM5.5/S.1/"

fName = "Base.RawPatch.20.20.bin";

fName = paste ( fDir, fName, sep="" );
finfo = file.info ( fName );
toread = file ( fName, "rb" );
numStimToShow = 5;
saveHTML ( { 
	for ( iStim in 1:numStimToShow ) {
		alldata = readBin ( toread, "double", size=8, n=(3*numValsPerTrial), endian="little" );
		v1e = ( alldata [ (startOffset.e):(startOffset.e + numValsPerTrial - 1) ]);
		v1i = ( alldata [ (startOffset.i):(startOffset.i + numValsPerTrial - 1) ] );
		v0 = ( alldata [ (startOffset.0):(startOffset.0 + numValsPerTrial - 1) ] );

		r1e = sigmoid ( v1e, 4 );
		r1i = sigmoid ( v1i, 4 );
		r0 = sigmoid ( v0, 4 );

		r0 = r0 / max(r0);
		for ( iTime in seq(1,numItersPerTrial) ) {
			titleText = paste ( "Response at Each Layer to Input Patch # ", iStim, "\nTime=", iTime*deltaT, sep="" );
			w = c ( r0[(iTime - 1) * N2 + seq ( 1:N2 )], r1i[(iTime - 1) * N2 + seq ( 1:N2 )],
				r1e[(iTime - 1) * N2 + seq ( 1:N2 )] );
			image ( c(1:N), c(1:(3*N)), matrix ( w, nrow=N, ncol=3*N, byrow=FALSE ),
				col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
				zlim = zLim, main=titleText, xlab=xLabText, ylab=yLabText );
			abline(h=layerMarks, col=1, lty=1 );
			abline(h=c(), col=1, lty=1 );
		} # for ( iCell in 1:N^2 )
	} # for ( iProbeCell = 1:N2 )
}, interval = 0.1 );  # saveGIF
close ( toread );

	#
	#	Time series plots
	#

x11(); par(mfcol=c(3,2))
for ( iCell in seq(818,818) ) {
	iCellWhich = seq ( iCell, numValsPerTrial, N2 );
	zLim = c ( min ( v1e[iCellWhich], v1i[iCellWhich], v0[iCellWhich] ), max ( v1e[iCellWhich], v1i[iCellWhich], v0[iCellWhich] ) );
	QuickCheckTimeSeries ( v1e[iCellWhich], v1i[iCellWhich], v0[iCellWhich], iCell, deltaT, FALSE, zLim );

	zLim = c ( min ( r1e[iCellWhich], r1i[iCellWhich], r0[iCellWhich] ), max ( r1e[iCellWhich], r1i[iCellWhich], r0[iCellWhich] ) );
	QuickCheckTimeSeries ( r1e[iCellWhich], r1i[iCellWhich], r0[iCellWhich], iCell, deltaT, FALSE, zLim );
	Sys.sleep(0.5);

} # for ( iCell in 1:N^2 )

x11(); par(mfrow=c(3,1));
plot ( r1e[iCellWhich] );
plot ( r1i[iCellWhich] );
plot ( r0[iCellWhich] );

x11(); par(mfrow=c(3,1));
plot ( v1e[iCellWhich] );
plot ( v1i[iCellWhich] );
plot ( v0[iCellWhich] );


	