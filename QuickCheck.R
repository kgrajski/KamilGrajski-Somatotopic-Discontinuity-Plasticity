######################################################################################################
######################################################################################################
#
#	Long-term state topo map evolution and various quantitative measure displays.
#
#	Run once (which is slow) to gather up all of the topo map data into a single list.
#
#	That way can rerun or do various other experiments faster.
#
#######################################################################################################
#######################################################################################################

#
#
#	Optimized for SYNDACTYLY.
#
#	There are 5 stages:
#
#	Base
#	SyndactCtl
#	SyndactReleaseCtl
#	SyndactExp
#	SyndactReleaseExp
#

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
kRFPeakToEdgeDetect = 0.5;

		#	Set the directory where the experimental data is sitting.
N = 75;

fDir = "D:/NMLab/Working/T.75.7.Rand.3/";

		#	Set some additional values.
N2 = N * N;

		#
		#	Get the Baseline Refinement data.
		#
fRoot = "Base.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");

iBase = 15;
iStart = 0;
iEnd = 5;
iStepSize = 5;
base = GetRFMapData2 ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iTop.base = length ( base$rfTrackData.e );

fRoot = "CheckInputStim.Baseline";
fileRootName = paste(fDir, fRoot, sep="\\");
#base.stimCount = GetStimCountData1 ( fileRootName, N2 );
base.stimCount = GetStimCountData ( fileRootName, N2 );

		#	Figure 0a.  Stimcount plots.
#if ( tiffFlag ) { tiff("Fig3.tif",compression="lzw",units="in",width=7.0,height=7.0,res=300); } else { x11(); }
x11(); par(mfrow=c(2,2));
titleText = paste ( "Stimulation Count\nBaseline Cycles\nMin(Blue)=", min(base.stimCount), " Max(Red)=", max(base.stimCount),sep="" );
xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";
ShowVecAsMap1 ( base.stimCount, titleText, xlab, ylab );
titleText = paste ( "Stimulation Count\nSyndactyly Cycles\nMin(Blue)=", min(syndact.stimCount), " Max(Red)=", max(syndact.stimCount),sep="" );
ShowVecAsMap1 ( syndact.stimCount, titleText, xlab, ylab );
#if ( tiffFlag ) { dev.off(); }


		#	Figure 1a.  Plot of E-Cell Centroids
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );


		#	Figure 1b.  Plot of I-Cell Centroids

subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

x11();
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
x11();
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );







