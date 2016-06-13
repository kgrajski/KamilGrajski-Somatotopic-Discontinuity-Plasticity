#
#	WeightEvolutionMovie.R
#

#
#	Use this in conjunction with RF Probe Raw data files from CPP
#

#
#	Quick visual check on weight time series
#

#
#	Eventually need to develop some automated routines to make sure
#	there isn't funny stuff happening.
#

rm(list = ls());

#
#	Helper Functions
#

source ( "NMHelperFunctions.R" );

InterleaveForDisplay = function ( w1, w2, w3, w4 ) {

	N = as.integer ( sqrt ( length ( w1 ) ) );
	wLen = 4 * N * N;
	w = rep ( 0, wLen );

	iLower = 1;
	iUpper = as.integer ( wLen/2 + 1 );

	for ( i in 1:N ) {

		w[iLower:(iLower+N-1)] = w1[((i-1)*N+1):(i*N)];
		iLower = iLower + N;

		w[iLower:(iLower+N-1)] = w3[((i-1)*N+1):(i*N)];
		iLower = iLower + N;

		w[iUpper:(iUpper+N-1)] = w2[((i-1)*N+1):(i*N)];
		iUpper = iUpper + N;

		w[iUpper:(iUpper+N-1)] = w4[((i-1)*N+1):(i*N)];
		iUpper = iUpper + N;

	} # for ( i in (2 * N) ) {

	return ( w );
	
} # InterleaveForDisplay = function ( w1, w2, w3, w4 ) {

#
#	Main
#


	#
	#	Set some constants
	#
N = 33;
deltaT = 0.005;
trialDuration = 0.25;

N2 = N * N;
numItersPerTrial = as.integer ( trialDuration * 1 / deltaT );
numValsPerRFTrial = numItersPerTrial * N * N;

	#
	#	Read the RF Probe Raw Data
	#

fDir = "D:/NMLab/S.21.4.Test/";
fDir = "E:/NMLab/Simulations/Syndactyly_5.1/S.33.7.5.75//";
fName = "Base.RFMap.Raw.25.25.bin";

fName = paste ( fDir, fName, sep="" );

iBase = 25;
iStart = 0;
iEnd = 25;
iStepSize = iCycle = 5;

	#
	#	Spatial plots
	#

yLabText="TBD"; xLabText="TBD";

numColors = 64;
rSeq = as.integer ( 1055 );

x11();

#saveHTML ( {

	for ( iProbeCell in rSeq ) {

		for ( iCycle in seq ( iStart, iEnd, iStepSize ) ) {

			titleText = paste ( "Cortical Column # ", iProbeCell, "\nCycle=", iCycle, sep="" );

			fRoot = paste("Base", iBase, iCycle, "w1E0", sep="."); E0 = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
			fRoot = paste("Base", iBase, iCycle, "w1EE", sep="."); EE = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
			fRoot = paste("Base", iBase, iCycle, "w1IE", sep="."); IE = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
			fRoot = paste("Base", iBase, iCycle, "w1EI", sep="."); EI = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );

			w.e0 = GetInputWeights ( E0, iProbeCell ); w.ee = GetInputWeights ( EE, iProbeCell );
			w.ie = GetInputWeights ( IE, iProbeCell ); w.ei = GetInputWeights ( EI, iProbeCell );

			#t1 = t2 = t3 = t4 = rep ( 0, N2 ); iTestLoc = (N2 / 2);
			#t1[iTestLoc] = 1.0; t2[iTestLoc] = 2.0; t3[iTestLoc] = 3.0; t4[iTestLoc] = 4.0;
			#w = InterleaveForDisplay ( t1, t2, t3, t4 );
			w = InterleaveForDisplay ( w.e0, w.ee, w.ie, w.ei );

			zlim = c ( min(w), max(w) );

			#image ( c(1:(N)), c(1:(N)), matrix ( t1, nrow=N, ncol=N, byrow=FALSE ),
			#col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
			#zlim = zlim, main=titleText, xlab=xLabText, ylab=yLabText );

			image ( c(1:(2*N)), c(1:(2*N)), matrix ( w, nrow=2*N, ncol=2*N, byrow=FALSE ),
					col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
					zlim = zlim, main=titleText, xlab=xLabText, ylab=yLabText );
			
			abline(h=(N+0.5), col=1, lty=1 );
			abline(v=(N+0.5), col=1, lty=1 );

			Sys.sleep ( 1.0 );
		
		} # for ( iCycle in seq ( iStart, iEnd, iStepSize ) {

	} # for ( iProbeCell in rSeq )

#}, interval = deltaT, htmlfile = "RFchk002.html" );  # saveGIF







