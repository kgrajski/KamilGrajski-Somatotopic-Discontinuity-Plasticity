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
#if ( !alreadyCleared ) { rm(list = ls()); }
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
#numIterPerTrial = 250;

g0 = 3;
N = 45;
N2 = N * N;

		#	Set the directory where the experimental data is sitting.

fDir = "E:/NMLab/Simulations/T.14/";

		#
		#	Get the Baseline Refinement data.
		#
fRoot = "Base.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");

iBase = 5;
iSequence = c ( 0, iBase );
base = GetRFMapData3 ( fileRootName, iBase, iSequence, kRFPeakToEdgeDetect, N2 );
iTop.base = length(base$rfTrackData.e);

fRoot = "CheckInputStim.Baseline";
fileRootName = paste(fDir, fRoot, sep="\\");
#base.stimCount = GetStimCountData1 ( fileRootName, N2 );
base.stimCount = GetStimCountData ( fileRootName, N2 );

		#
		#	Get the Syndactyly Experiment data
		#
fRoot = "SyndactExp.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iBase = 5;
iSequence = c ( 1, 5 );
syndact = GetRFMapData3 ( fileRootName, iBase, iSequence, kRFPeakToEdgeDetect, N2 );
iTop.syndact = length(syndact$rfTrackData.e);

fRoot = "CheckInputStim.Syndact";
fileRootName = paste(fDir, fRoot, sep="\\");
#syndact.stimCount = GetStimCountData1 ( fileRootName, N2 );
syndact.stimCount = GetStimCountData ( fileRootName, N2 );

		#
		#	Get the Syndactyly Release data.
		#
fRoot = "SyndactReleaseExp.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iBase = 5;
iSequence = c ( 1, iBase );
syndactRel = GetRFMapData3 ( fileRootName, iBase, iSequence, kRFPeakToEdgeDetect, N2 );
iTop.syndactRel = length(syndactRel$rfTrackData.e);


		#	Figure 0a.  Stimcount plots.
#if ( tiffFlag ) { tiff("Fig3.tif",compression="lzw",units="in",width=7.0,height=7.0,res=300); } else { x11(); }
x11(); par(mfrow=c(2,2));
titleText = paste ( "Stimulation Count\nBaseline Cycles\nMin(Blue)=", min(base.stimCount), " Max(Red)=", max(base.stimCount),sep="" );
xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";
ShowVecAsMap1 ( base.stimCount, titleText, xlab, ylab );
x11(); par(mfrow=c(2,2));
titleText = paste ( "Stimulation Count\nSyndactyly Cycles\nMin(Blue)=", min(syndact.stimCount), " Max(Red)=", max(syndact.stimCount),sep="" );
ShowVecAsMap1 ( syndact.stimCount, titleText, xlab, ylab );
#if ( tiffFlag ) { dev.off(); }


		#	Figure 1a.  Plot of E-Cell Centroids

x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = "\nSyndactyly";
iWhich = 1; ShowTopoMap1 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.syndact; ShowTopoMap1 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = "\nSyndactyly Release";
iWhich = 1; ShowTopoMap1 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.syndactRel; ShowTopoMap1 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

		#	Figure 1b.  Plot of I-Cell Centroids
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = "\nSyndactyly";
iWhich = 1; ShowTopoMap1 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.syndact; ShowTopoMap1 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = "\nSyndactyly Release";
iWhich = 1; ShowTopoMap1 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.syndactRel; ShowTopoMap1 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

		#	Figure 2a.	Plot of E-Cell Tracks RF Extents
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowThreeDigitRFTrack1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.base; ShowThreeDigitRFTrack1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = "\nSyndactyly";
iWhich = 1; ShowThreeDigitRFTrack1 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndact; ShowThreeDigitRFTrack1 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = "\nSyndactyly Release";
iWhich = 1; ShowThreeDigitRFTrack1 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack1 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

		#	Figure 2b.	Plot of I-Cell Tracks RF Extents
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowThreeDigitRFTrack1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.base; ShowThreeDigitRFTrack1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = "\nSyndactyly";
iWhich = 1; ShowThreeDigitRFTrack1 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndact; ShowThreeDigitRFTrack1 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = "\nSyndactyly Release";
iWhich = 1; ShowThreeDigitRFTrack1 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack1 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

		#	Figure 3a.	Plot of E-Cell Tracks RF Extents II
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowThreeDigitRFTrack2 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.base; ShowThreeDigitRFTrack2 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = "\nSyndactyly";
iWhich = 1; ShowThreeDigitRFTrack2 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndact; ShowThreeDigitRFTrack2 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = "\nSyndactyly Release";
iWhich = 1; ShowThreeDigitRFTrack2 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack2 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

		#	Figure 3b.	Plot of I-Cell Tracks RF Extents II
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowThreeDigitRFTrack2 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.base; ShowThreeDigitRFTrack2 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = "\nSyndactyly";
iWhich = 1; ShowThreeDigitRFTrack2 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndact; ShowThreeDigitRFTrack2 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = "\nSyndactyly Release";
iWhich = 1; ShowThreeDigitRFTrack2 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack2 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

		#	Figure 4a.	Plot of E-Cell Tracks RF Extents III (Dist, Mid, Prox - across digits)
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowThreeDigitRFTrack3 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );
iWhich = iTop.base; ShowThreeDigitRFTrack3 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );

subTitleText = "\nSyndactyly";
iWhich = 1; ShowThreeDigitRFTrack3 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );
iWhich = iTop.syndact; ShowThreeDigitRFTrack3 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );

subTitleText = "\nSyndactyly Release";
iWhich = 1; ShowThreeDigitRFTrack3 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack3 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );

		#	Figure 4b.	Plot of I-Cell Tracks RF Extents II
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowThreeDigitRFTrack3 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );
iWhich = iTop.base; ShowThreeDigitRFTrack3 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );

subTitleText = "\nSyndactyly";
iWhich = 1; ShowThreeDigitRFTrack3 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );
iWhich = iTop.syndact; ShowThreeDigitRFTrack3 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );

subTitleText = "\nSyndactyly Release";
iWhich = 1; ShowThreeDigitRFTrack3 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack3 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );

		#	Figure 5.  Identify: RF Divergence Spatial Plot

iWhich = 1; div.base.init = IntraColumnarRFDivergence ( base$rfTrackData.e[[iWhich]], base$rfTrackData.i[[iWhich]] );
iWhich = iTop.base; div.base.final = IntraColumnarRFDivergence ( base$rfTrackData.e[[iWhich]], base$rfTrackData.i[[iWhich]] );

iWhich = 1; div.syndact.init = IntraColumnarRFDivergence ( syndact$rfTrackData.e[[iWhich]], syndact$rfTrackData.i[[iWhich]] );
iWhich = iTop.syndact; div.syndact.final = IntraColumnarRFDivergence ( syndact$rfTrackData.e[[iWhich]], syndact$rfTrackData.i[[iWhich]] );

iWhich = 1; div.syndactRel.init = IntraColumnarRFDivergence ( syndactRel$rfTrackData.e[[iWhich]], syndactRel$rfTrackData.i[[iWhich]] );
iWhich = iTop.syndactRel; div.syndactRel.final = IntraColumnarRFDivergence ( syndactRel$rfTrackData.e[[iWhich]], syndactRel$rfTrackData.i[[iWhich]] );

minDiv = min ( c ( div.base.init, div.base.final, div.syndact.init, div.syndact.final, div.syndactRel.init, div.syndactRel.final ) );
maxDiv = max ( c ( div.base.init, div.base.final, div.syndact.init, div.syndact.final, div.syndactRel.init, div.syndactRel.final ) );

x11(); par ( mfrow = c ( 3,2 ) );
xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

iWhich = 1; titleText = paste ( "Intracolumnar Receptive Field Divergence\nBaseline Refinement\nIter ", (iWhich-1), sep="" );
ShowVecAsMap2 ( div.base.init, titleText, xlab, ylab, minDiv, maxDiv );

iWhich = iEnd; titleText = paste ( "Intracolumnar Receptive Field Divergence\nBaseline Refinement\nIter ", (iWhich-1), sep="" );
ShowVecAsMap2 ( div.base.final, titleText, xlab, ylab, minDiv, maxDiv );

iWhich = 1; titleText = paste ( "Intracolumnar Receptive Field Divergence\nSyndactyly\nIter ", (iWhich), sep="" );
ShowVecAsMap2 ( div.syndact.init, titleText, xlab, ylab, minDiv, maxDiv );

iWhich = iEnd; titleText = paste ( "Intracolumnar Receptive Field Divergence\nSyndactyly\nIter ", (iWhich), sep="" );
ShowVecAsMap2 ( div.syndact.final, titleText, xlab, ylab, minDiv, maxDiv );

iWhich = 1; titleText = paste ( "Intracolumnar Receptive Field Divergence\nSyndactyly Release\nIter ", (iWhich), sep="" );
ShowVecAsMap2 ( div.syndactRel.init, titleText, xlab, ylab, minDiv, maxDiv ); 

iWhich = iEnd; titleText = paste ( "Intracolumnar Receptive Field Divergence\nSyndactyly Release\nIter ", (iWhich), sep="" );
ShowVecAsMap2 ( div.syndactRel.final, titleText, xlab, ylab, minDiv, maxDiv );

		#	Identify: Optional Histogram Plots of RF Divergence
x11(); par ( mfrow = c ( 3,2 ) );
hist ( div.base.init, freq = FALSE, breaks = seq ( 0, maxDiv+0.5, 0.5 ) );
hist ( div.base.final, freq = FALSE, breaks = seq ( 0, maxDiv+0.5, 0.5 ) );
hist ( div.syndact.init, freq = FALSE, breaks = seq ( 0, maxDiv+0.5, 0.5 ) );
hist ( div.syndact.final, freq = FALSE, breaks = seq ( 0, maxDiv+0.5, 0.5 ) );
hist ( div.syndactRel.init, freq = FALSE, breaks = seq ( 0, maxDiv+0.5, 0.5 ) );
hist ( div.syndactRel.final, freq = FALSE, breaks = seq ( 0, maxDiv+0.5, 0.5 ) );

		#	Characterize: Select a few cells for a closer look.
x11(); par(mfrow=c(3,2)); ShowVecAsMap1 ( (div.base.final > sqrt(2.0)), titleText, xlab, ylab );
div.base.final.cellList = which ( div.base.final > sqrt(2.0) );

		#
digitWidth = (N / 3);
div.base.final.cellList = c ( (digitWidth * N) + 6 - N, (digitWidth * N) + 6, (digitWidth * 2 * N) + 6 + 3 * N);
symbolChoice = c ( 1, 0, 2 );

		#	Topo Quant figures

x11(); par(mfrow=c(2,1));

base.q.e = TopoQuantSeries ( N, g0, base$rfTrackData.e );
base.q.i = TopoQuantSeries ( N, g0, base$rfTrackData.i );

syndact.q.e = TopoQuantSeries ( N, g0, syndact$rfTrackData.e );
syndact.q.i = TopoQuantSeries ( N, g0, syndact$rfTrackData.i );

syndactRel.q.e = TopoQuantSeries ( N, g0, syndactRel$rfTrackData.e );
syndactRel.q.i = TopoQuantSeries ( N, g0, syndactRel$rfTrackData.i );

x1.e = c ( base.q.e$q1, syndact.q.e$q1, syndactRel.q.e$q1 );
x1.i = c ( base.q.i$q1, syndact.q.i$q1, syndactRel.q.i$q1 );

ylim = c ( min ( x1.e, x1.i ), max ( x1.e, x1.i ) );
plot ( x1.e, type="l", col=1, ylim=ylim, xlab="Iteration", ylab="C Metric",
	main="(F,G) Topography Measure vs Refinement Iterations\nE Cell Layer (Black); I Cell Layer (Red)" );
lines ( x1.i, col=2 );
abline ( v=c(26,51), lty=4, col=4 );

x2.e = c ( base.q.e$q2, syndact.q.e$q2, syndactRel.q.e$q2 );
x2.i = c ( base.q.i$q2, syndact.q.i$q2, syndactRel.q.i$q2 );

ylim = c ( min ( x2.e, x2.i ), max ( x2.e, x2.i ) );
plot ( x2.e, type="l", col=1, ylim=ylim, xlab="Iteration", ylab="C Metric",
	main="(F,G) Topography Measure vs Refinement Iterations\nE Cell Layer (Black); I Cell Layer (Red)" );
lines ( x2.i, col=2 );
abline ( v=c(26,51), lty=4, col=4 );

#
#	Get Weights
#

fRoot.first = paste ( "Base", iBase, "0", sep="." );
fRoot.last = paste ( "Base", iBase, iBase, sep="." );
init.w1.e0 = GetSparseWeightMatrix ( paste(fDir, paste(fRoot.first,"w1E0",sep="."), sep="\\") );
init.w1.ee = GetSparseWeightMatrix ( paste(fDir, paste(fRoot.first,"w1EE",sep="."), sep="\\") );
init.w1.ei = GetSparseWeightMatrix ( paste(fDir, paste(fRoot.first,"w1EI",sep="."), sep="\\") );
init.w1.ie = GetSparseWeightMatrix ( paste(fDir, paste(fRoot.first,"w1IE",sep="."), sep="\\") );

refine.w1.e0 = GetSparseWeightMatrix ( paste(fDir, paste(fRoot.last,"w1E0",sep="."), sep="\\") );
refine.w1.ee = GetSparseWeightMatrix ( paste(fDir, paste(fRoot.last,"w1EE",sep="."), sep="\\") );
refine.w1.ei = GetSparseWeightMatrix ( paste(fDir, paste(fRoot.last,"w1EI",sep="."), sep="\\") );
refine.w1.ie = GetSparseWeightMatrix ( paste(fDir, paste(fRoot.last,"w1IE",sep="."), sep="\\") );

fRoot.last = paste ( "SyndactExp", iBase, iBase, sep="." );
syndact.w1.e0 = GetSparseWeightMatrix ( paste(fDir, paste(fRoot.last,"w1E0",sep="."), sep="\\")  );
syndact.w1.ee = GetSparseWeightMatrix ( paste(fDir, paste(fRoot.last,"w1EE",sep="."), sep="\\")  );
syndact.w1.ei = GetSparseWeightMatrix ( paste(fDir, paste(fRoot.last,"w1EI",sep="."), sep="\\") );
syndact.w1.ie = GetSparseWeightMatrix ( paste(fDir, paste(fRoot.last,"w1IE",sep="."), sep="\\") );


