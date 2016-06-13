4######################################################################################################
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
#numIterPerTrial = 150;

g0 = 3;

		#	Set the directory where the experimental data is sitting.
E = 65;
N = 15;
P = "5.EOnly";
fDir = paste("E:\\NMLab\\E",E,".",N,".",P,sep="");
fDir = "../";
fDir = "D:/NMLab/S.21.4.Control/"
fDir = "D:/NMLab/S.21.4.Test.1/"
fDir = "D:/NMLab/S.27.4.Test.2/"
fDir = "D:/NMLab/S.15.3.Control.Reverse/"

		#	Set some additional values.
N2 = N * N;

		#
		#	Get the Baseline Refinement data.
		#
fRoot = "Base.RFMap";
fRoot = "SyndactExp.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");

iBase = 25;
iStart = 1;
iEnd = 18;
base = GetRFMapData ( fileRootName, iBase, iStart, iEnd, kRFPeakToEdgeDetect, N2 );

fRoot = "CheckInputStim.Baseline";
fRoot = "CheckInputStim.Syndact";
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
iWhich = iEnd; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );


		#	Figure 1b.  Plot of I-Cell Centroids

subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iEnd; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );










