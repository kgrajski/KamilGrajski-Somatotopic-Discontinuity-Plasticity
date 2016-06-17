#
#	SyndactylyPaperFigureGenerator.FinalPubs.Fig2.R
#

#
#	This script is an ugly cut & paste from a bunch of other scripts that make more A:B:C comparison plots.
#

#	Clear the workspace.
rm(list = ls());

#
source ( "NMHelperFunctions.R" );

#
#	Local help functions.
#

######################################################################################################
######################################################################################################
#
#	Publication Quality Figure Showing Membrane Potential and Firing Rate 
#
######################################################################################################
######################################################################################################

RFProbeTrialTimeSeriesPlot = function ( N2, numValsPerRFTrial, iProbeCell, iStart, iEnd,
								alldata.base, alldata.exp, alldata.exp2, alldata.exp3, eScale,
								titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag,
								deltaTInMsec ) {
	#alldata.exp = alldata.base; alldata.exp2 = alldata.base; alldata.exp3 = alldata.base;
	titleText = "RF Probe Trial\nAvg. Membrane Potential";
	xlabText = xlabTextRoot;
	ylabText = ylabTextRoot;

	ts.exp = ExpRefTSResponse ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale );

	if ( tiffFlag ) { tiff(filename="Fig2.tiff",compression="lzw",units="in",width=7.0,height=7.0,res=600); } else { x11(); }
	par(mfrow=c(2,2));
	titleTextToPrint = paste ( titleText );
	xVal = iStart:iEnd;
	yMin = min ( ts.exp$v1e.base[xVal], ts.exp$v1i.base, ts.exp$v0.base );
	yMax = max ( ts.exp$v1e.base[xVal], ts.exp$v1i.base, ts.exp$v0.base ) + max(ts.exp$v0.base[xVal]);
	ylim = c ( yMin, yMax );
	plot ( xVal*deltaTInMsec, ts.exp$v1e.base[xVal], type="l", cex=0.25, ylim=ylim, col=1,
		main=titleTextToPrint, xlab=xlabText, ylab="" );
	lines ( xVal*deltaTInMsec, ts.exp$v1i.base[xVal], col=2, cex=0.25 );
	axis(1); mtext(side=2, "Ve, Vi", cex=0.75, line=2);
	abline ( h=0, lty=4, col=5 );
	par(new=TRUE);
	plot ( xVal*deltaTInMsec, ts.exp$v0.base[xVal], type="l", cex=0.1, col=3, lty=1, xaxt="n", yaxt="n", xlab="", ylab="" );
	#abline ( h=0, lty=4, col=5 );
	axis(4); mtext(side=4, "Vs", cex=0.75, line=0);
	legend ( 102, 0.45, c("Ve", "Vi", "Vs"), lty=c(1,1,1), lwd=c(1,1,1), col=c(1,2,3), cex=0.60);

	titleText = "RF Probe Trial\nAvg. Spiking Rate";
	xlabText = xlabTextRoot;
	ylabText = ylabTextRoot;
	titleTextToPrint = paste ( titleText );
	xVal = iStart:iEnd;
	yMin = min ( ts.exp$v1e.base[xVal], ts.exp$v1i.base, ts.exp$v0.base );
	yMax = max ( ts.exp$v1e.base[xVal], ts.exp$v1i.base, ts.exp$v0.base ) + max(ts.exp$v0.base[xVal]);
	ylim = c ( yMin, yMax );
	ylim = sigmoid ( ylim, 4.0);
	plot ( xVal*deltaTInMsec, sigmoid(ts.exp$v1e.base[xVal],4), type="l", cex=0.25, ylim=ylim, col=1,
		main=titleTextToPrint, xlab=xlabText, ylab="" );
	lines ( xVal*deltaTInMsec, sigmoid(ts.exp$v1i.base[xVal],4), col=2, cex=0.25 );
	lines ( xVal*deltaTInMsec, sigmoid(ts.exp$v0.base[xVal],4), col=3, cex=0.25 );
	axis(1); mtext(side=2, "Re, Ri, Rs", cex=0.75, line=2);

	legend ( 125, 0.45, c("Re", "Ri", "Rs"), lty=c(1,1,1), lwd=c(1,1,1), col=c(1,2,3), cex=0.60);

	if ( tiffFlag ) { dev.off(); }

} # RFProbeTrialTimeSeriesPlot = function ( N2, numValsPerRFTrial, iProbeCell, ...



	#
	#	Additional helper functions.
	#

	#
	#	Set some constants
	#
N = 45;
N2 = N * N;
deltaT = 0.001;
trialDuration = 0.350;
numItersPerTrial = 350;
numValsPerRFTrial = numItersPerTrial * N * N;
tiffFlag = FALSE;
initialPass = FALSE;
iShortPlot = FALSE;

	#
	#	Identify the pool of recording sites of potential interest.
	#
iShowTheseRows = seq ( 18, 32 );
iCellList = c ( 333 );

	#
	#	Identify the speciufic RF probes of interest.
	#
iSetProbeList = c ( 333 );

	#
	#	Set up titles
	#
eScale = 1.0;

	#
	#	Directory storing the Raw RF Probe Files.
	#

fDir = "D:/NMLab/Keep - Syndactyly Exps/T.7/RF Raw Output For Completed Baseline v1.0/";
fName.base = paste ( "Base.RFMap.Raw.15.1.bin" );

	#
	#	Set up final final names and read the data.
	#
fName.base = paste ( fDir, fName.base, sep="" );

	#
	#	Get the subset of data of interest and save it.
	#	The files can be quite large so have to split doing the control and exps modes.
	#
base = ExtractRFProbeRawData ( fName.base, N2, numValsPerRFTrial, iCellList );

		#
		#	Time series difference plots.
		#
iProbe = 1;	
iProbeCell = iSetProbeList[iProbe];
indexToRFProbeRawData = which ( iCellList == iProbeCell );
alldata.base = c ( base$v1e[,indexToRFProbeRawData], base$v1i[,indexToRFProbeRawData], base$v0[,indexToRFProbeRawData] );

titleText = "RF Probe Trial";
subTitleText = NULL;
xlabTextRoot = "Msec.";
ylabTextRoot = "Vs(t), Vi(t)";

deltaTInMsec = 1;
iStart = 90;
iEnd = 180;
eScale = 1.0;
tiffFlag = TRUE;
RFProbeTrialTimeSeriesPlot ( N2, numValsPerRFTrial, iProbeCell, iStart, iEnd,
						alldata.base, alldata.base, alldata.base, alldata.base, eScale,
						titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag,
						deltaTInMsec );



