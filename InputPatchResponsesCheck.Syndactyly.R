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

	plot( tVal, r1E, type="l", col=1, ylim=ylim, ylab="Firing Rate", xlab="Time (sec)");
	title(main=paste("Cortical Layer E-Cell Firing Rate",paste("Cell ",iCell),sep='\n'));

	plot( tVal, r1I, type="l", col=1, ylim=ylim, ylab="Firing Rate", xlab="Time (sec)");
	title(main=paste("Cortical Layer I-Cell Firing Rate",paste("Cell ",iCell),sep='\n'));

	plot( tVal, rZ, type="l", col=1, ylim=ylim, ylab="Firing Rate", xlab="Time (sec)");
	title(main=paste("Input Layer Firing Rate",paste("Cell ",iCell),sep='\n'));

} # QuickCheckTimeSeries = function ( r1E, r1I, r0, cellID, deltaT, ... ) {


#
#	Main
#

source("NMHelperFunctions.R");

	#	Set some constants
N = 30;
N2 = N * N;
deltaT = 0.01;
trialDuration = 1.0;
numItersPerTrial = as.integer ( trialDuration * 1 / deltaT );
numValsPerTrial = numItersPerTrial * N * N;

	#	Expected file size (for sanity checking)
	#		For 3-digit baseline:
digitWidth = N / 3;
patchSize = 4;
numPatches.baseline = 3 * ( digitWidth - patchSize ) * ( N - patchSize );
numStim = expectedLen.baseline = 3 * numValsPerTrial * numPatches.baseline;

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
fDir = "D:/NMLab/S.30.4.Control/";
fDir = "D:/NMLab/S.30.4.Test/";
fName = "Base.RawPatch.25.1.bin";
fName = paste ( fDir, fName, sep="" );
finfo = file.info ( fName );
toread = file ( fName, "rb" );
numStimToShow = 20;
saveHTML ( { 
	for ( iStim in 1:numStimToShow ) {
		alldata = readBin ( toread, "double", size=8, n=(3*numValsPerTrial), endian="little" );
		r1e = alldata [ (startOffset.e):(startOffset.e + numValsPerTrial - 1) ];
		r1i = alldata [ (startOffset.i):(startOffset.i + numValsPerTrial - 1) ];
		r0 = alldata [ (startOffset.0):(startOffset.0 + numValsPerTrial - 1) ];	
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
}, interval = 0.1, htmlfile="PatchFig3.Test" );  # saveGIF
close ( toread );

	#
	#	Time series plots
	#

for ( iCell in 1:N^2 ) {

	iCellWhich = seq ( iCell, numValsPerTrial, N2 );
	QuickCheckTimeSeries ( r1e[iCellWhich], r1i[iCellWhich], r0[iCellWhich], iCell, deltaT, FALSE, ylim );
	Sys.sleep ( 1.0 );

} # for ( iCell in 1:N^2 )

	#
	#	No longer used, but may come in handy one day.
	#

#titleText = paste ( "RF Probe # ", iProbeCell, " E\nTime=", iTime*deltaT, sep="" );
#ShowVecAsMap2 ( r1e[(iTime - 1) * N2 + seq ( 1:N2 )], titleText, xLabText, yLabText, minZ, maxZ );
#titleText = paste ( "RF Probe # ", iProbeCell, " I\nTime=", iTime*deltaT, sep="" );
#ShowVecAsMap2 ( r1i[(iTime - 1) * N2 + seq ( 1:N2 )], titleText, xLabText, yLabText, minZ, maxZ );
#titleText = paste ( "RF Probe # ", iProbeCell, " S\nTime=", iTime*deltaT, sep="" );
#ShowVecAsMap2 ( r0[(iTime - 1) * N2 + seq ( 1:N2 )], titleText, xLabText, yLabText, minZ, maxZ );