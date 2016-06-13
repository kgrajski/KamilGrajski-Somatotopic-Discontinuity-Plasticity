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
#	Optimized for SELECTIVE STIMULATION.
#
#	There are 5 stages:

#	Clear the workspace.
if ( !alreadyCleared ) { rm(list = ls()); }

#	Load up some libraries and helper functions.
require(ggplot2);
require(reshape2);
require(ellipse);
require(stats4);
require(fields);
source ( "NMHelperFunctions.R" );

#
#	MAIN
#

		#
		#	Some additional variables needed if/when NMExpZone isn't the source.
kRFPeakToEdgeDetect = 0.5;
numIterPerTrial = 150;
g0 = 3;

		#
		#	Get the Baseline Refinement data.
		#
fRoot = "Base.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");

iBase = 25; iStart = 0; iEnd = 25;
base = GetRFMapData ( fileRootName, iBase, iStart, iEnd, kRFPeakToEdgeDetect, N2 );

fRoot = "CheckInputStim.Baseline";
fileRootName = paste(fDir, fRoot, sep="\\");
base.stimCount = GetStimCountData ( fileRootName, N2 );

		#
		#	Get the SelStim Experiment data.
		#
fRoot = "SelStimExp.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iBase = 25; iStart = 1; iEnd = 25;
#selStim = GetRFMapData ( fileRootName, iBase, iStart, iEnd, kRFPeakToEdgeDetect, N2 );
selStim = GetRFMapData1 ( fileRootName, iBase, iStart, iEnd, expSectorID, expFactor, kRFPeakToEdgeDetect, N2 );


fRoot = "CheckInputStim.SelStim";
fRoot2 = paste ( fRoot, expSectorID, expFactor, sep="." );
fileRootName = paste(fDir, fRoot, sep="\\");
fileRootName = paste(fDir, fRoot2, sep="\\");
selStim.stimCount = GetStimCountData ( fileRootName, N2 );

		#
		#	Get the SelStim Control data.
		#
fRoot = "SelStimCtl.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iBase = 25; iStart = 1; iEnd = 25;
selStim.ctl = GetRFMapData ( fileRootName, iBase, iStart, iEnd, kRFPeakToEdgeDetect, N2 );

		#	Figure 0a.  Stimcount plots.
x11(); par(mfrow=c(2,2));
#titleText = paste ( "Stimulation Count\nBaseline Cycles\nMin(Blue)=", min(base.stimCount), " Max(Red)=", max(base.stimCount),sep="" );
titleText = paste ( "Stimulation Count\nBaseline Cycles", sep="" );
xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";
ShowVecAsMapContour ( base.stimCount, titleText, xlab, ylab );
abline(h=c(11.5,22.5),lty=3,col=605);
#titleText = paste ( "Stimulation Count\nFocal Stimulation Cycles\nMin(Blue)=", min(selStim.stimCount), " Max(Red)=", max(selStim.stimCount),sep="" );
titleText = paste ( "Stimulation Count\nFocal Stimulation Cycles", sep="" );
ShowVecAsMapContour ( selStim.stimCount, titleText, xlab, ylab );
abline(h=c(11.5,22.5),lty=3,col=605);

iFirst = 1; iLast = 25; iLast.base = 26
iTmp = 25;

		#	Figure 1.  CTL: Plot of E-Cell & I-Cell Centroids
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement I";
iWhich = iLast.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iTmp,sep=" "), subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iLast.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iTmp,sep=" "), subTitleText, sep=""), FALSE, 0.5, 0 );
subTitleText = "\nBaseline Refinement II";
iWhich = iFirst; ShowTopoMap1 ( selStim.ctl$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "), subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iFirst; ShowTopoMap1 ( selStim.ctl$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "), subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iLast; ShowTopoMap1 ( selStim.ctl$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "), subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iLast; ShowTopoMap1 ( selStim.ctl$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "), subTitleText, sep=""), FALSE, 0.5, 0 );

		#	Figure 2.  Exp: Plot of E-Cell & I-Cell Centroids
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = iLast.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iTmp,sep=" "), subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iLast.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iTmp,sep=" "), subTitleText, sep=""), FALSE, 0.5, 0 );
subTitleText = "\nFocal Stimulation";
iWhich = iFirst; ShowTopoMap1 ( selStim$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "), subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iFirst; ShowTopoMap1 ( selStim$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "), subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iLast; ShowTopoMap1 ( selStim$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "), subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iLast; ShowTopoMap1 ( selStim$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "), subTitleText, sep=""), FALSE, 0.5, 0 );

		#	Figure 2A: CTL: Spatial plot of quantitative measure of RF translocation.
minMoveDist = sqrt(2.0)-0.1;
minMoveDist = 0.9
x11(); par( mfrow=c(2,2) );
mainText.e = "RF Centroid Translocation\nFocal Stim vs Focal Stim Ctl\nE-Type ";
mainText.i = "RF Centroid Translocation\nFocal Stim vs Focal Stim Ctl\nI-Type ";
xlabText="Distal -> Proximal"; ylabText="Digit 1 -> Digit 3";
RFTransLocPlot ( selStim$rfTrackData.e[[iLast]]$rfCentroid, selStim.ctl$rfTrackData.e[[iLast]]$rfCentroid, mainText.e, xlabText, ylabText, minMoveDist );
RFTransLocPlot ( selStim$rfTrackData.i[[iLast]]$rfCentroid, selStim.ctl$rfTrackData.i[[iLast]]$rfCentroid, mainText.i, xlabText, ylabText, minMoveDist );

		#	Figure 3.	CTL: Plot of E-Cell Tracks RF Extents
#x11(); par(mfrow=c(3,2));
#subTitleText = "\nBaseline Refinement I";
#iWhich = iLast.base; ShowThreeDigitRFTrack3 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iTmp,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );
#iWhich = iLast.base; ShowThreeDigitRFTrack3 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iTmp,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );
#subTitleText = "\nBaseline Refinement II";
#iWhich = iFirst; ShowThreeDigitRFTrack3 ( selStim.ctl$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );
#iWhich = iFirst; ShowThreeDigitRFTrack3 ( selStim.ctl$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );
#iWhich = iLast; ShowThreeDigitRFTrack3 ( selStim.ctl$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );
#iWhich = iLast; ShowThreeDigitRFTrack3 ( selStim.ctl$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0 );

		#	Figure 4.	EXP: Plot of E-Cell Tracks RF Extents
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = iLast.base; ShowThreeDigitRFTrack3 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iTmp,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 2 );
iWhich = iLast.base; ShowThreeDigitRFTrack3 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iTmp,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 2 );
subTitleText = "\nFocal Stimulation";
iWhich = iFirst; ShowThreeDigitRFTrack3 ( selStim$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 2 );
iWhich = iFirst; ShowThreeDigitRFTrack3 ( selStim$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 2 );
iWhich = iLast; ShowThreeDigitRFTrack3 ( selStim$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 2 );
iWhich = iLast; ShowThreeDigitRFTrack3 ( selStim$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 2 );

if ( 0 ) {
		#	Figure 5.	CTL: Plot of E-Cell Tracks RF Extents
x11(); par(mfrow=c(3,2));
iColList = c ( (N/3)+5, 2*(N/3)+5 );
			#	iColList: used to highlight specific longitudinal tracks; which are of interest depends
			#	on the experimental Sector ID
iColList = c ( 19 );
subTitleText = "\nBaseline Refinement I";
iWhich = iLast.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iTmp,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList, 2 );
iWhich = iLast.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iTmp,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList, 2 );
subTitleText = "\nBaseline Refinement II";
iWhich = iFirst; ShowThreeDigitRFTrackFlex ( selStim.ctl$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList, 2 );
iWhich = iFirst; ShowThreeDigitRFTrackFlex ( selStim.ctl$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList, 2 );
iWhich = iLast; ShowThreeDigitRFTrackFlex ( selStim.ctl$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList, 2 );
iWhich = iLast; ShowThreeDigitRFTrackFlex ( selStim.ctl$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList, 2 );

		#	Figure 6.	EXP: Plot of E-Cell Tracks RF Extents
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = iLast.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iTmp,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList, 2 );
iWhich = iLast.base; ShowThreeDigitRFTrackFlex ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iTmp,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList, 2 );
subTitleText = "\nFocal Stimulation";
iWhich = iFirst; ShowThreeDigitRFTrackFlex ( selStim$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList, 2 );
iWhich = iFirst; ShowThreeDigitRFTrackFlex ( selStim$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList, 2 );
iWhich = iLast; ShowThreeDigitRFTrackFlex ( selStim$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList, 2 );
iWhich = iLast; ShowThreeDigitRFTrackFlex ( selStim$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1, iColList, 2 );
		
	#	Figure 7.  Identify: RF Divergence Spatial Plot

iWhich = 1; div.base.init = IntraColumnarRFDivergence ( base$rfTrackData.e[[iWhich]], base$rfTrackData.i[[iWhich]] );
iWhich = 26; div.base.final = IntraColumnarRFDivergence ( base$rfTrackData.e[[iWhich]], base$rfTrackData.i[[iWhich]] );

iWhich = 1; div.selStim.init = IntraColumnarRFDivergence ( selStim$rfTrackData.e[[iWhich]], selStim$rfTrackData.i[[iWhich]] );
iWhich = 25; div.selStim.final = IntraColumnarRFDivergence ( selStim$rfTrackData.e[[iWhich]], selStim$rfTrackData.i[[iWhich]] );

iWhich = 1; div.selStim.ctl.init = IntraColumnarRFDivergence ( selStim.ctl$rfTrackData.e[[iWhich]], selStim.ctl$rfTrackData.i[[iWhich]] );
iWhich = 25; div.selStim.ctl.final = IntraColumnarRFDivergence ( selStim.ctl$rfTrackData.e[[iWhich]], selStim.ctl$rfTrackData.i[[iWhich]] );

minDiv = min ( c ( div.base.init, div.base.final, div.selStim.init, div.selStim.final, div.selStim.ctl.init, div.selStim.ctl.final ) );
maxDiv = max ( c ( div.base.init, div.base.final, div.selStim.init, div.selStim.final, div.selStim.ctl.init, div.selStim.ctl.final ) );

x11(); par ( mfrow = c ( 3,2 ) );
xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

iWhich = 1; titleText = paste ( "Intracolumnar Receptive Field Divergence\nBaseline Refinement\nIter ", (iWhich-1), sep="" );
ShowVecAsMap2 ( div.base.init, titleText, xlab, ylab, minDiv, maxDiv );

iWhich = 26; titleText = paste ( "Intracolumnar Receptive Field Divergence\nBaseline Refinement\nIter ", (iWhich-1), sep="" );
ShowVecAsMap2 ( div.base.final, titleText, xlab, ylab, minDiv, maxDiv );

iWhich = 1; titleText = paste ( "Intracolumnar Receptive Field Divergence\nFocal Stimulation\nIter ", (iWhich), sep="" );
ShowVecAsMap2 ( div.selStim.init, titleText, xlab, ylab, minDiv, maxDiv );

iWhich = 25; titleText = paste ( "Intracolumnar Receptive Field Divergence\nFocal Stimulation\nIter ", (iWhich), sep="" );
ShowVecAsMap2 ( div.selStim.final, titleText, xlab, ylab, minDiv, maxDiv );

iWhich = 1; titleText = paste ( "Intracolumnar Receptive Field Divergence\nFocal Stimulation Control\nIter ", (iWhich), sep="" );
ShowVecAsMap2 ( div.selStim.ctl.init, titleText, xlab, ylab, minDiv, maxDiv ); 

iWhich = 25; titleText = paste ( "Intracolumnar Receptive Field Divergence\nFocal Stimulation Control\nIter ", (iWhich), sep="" );
ShowVecAsMap2 ( div.selStim.ctl.final, titleText, xlab, ylab, minDiv, maxDiv );

		#	Identify: Optional Histogram Plots of RF Divergence
x11(); par ( mfrow = c ( 3,2 ) );
hist ( div.base.init, freq = FALSE, breaks = seq ( 0, maxDiv+0.5, 0.5 ) );
hist ( div.base.final, freq = FALSE, breaks = seq ( 0, maxDiv+0.5, 0.5 ) );
hist ( div.selStim.init, freq = FALSE, breaks = seq ( 0, maxDiv+0.5, 0.5 ) );
hist ( div.selStim.final, freq = FALSE, breaks = seq ( 0, maxDiv+0.5, 0.5 ) );
hist ( div.selStim.ctl.init, freq = FALSE, breaks = seq ( 0, maxDiv+0.5, 0.5 ) );
hist ( div.selStim.ctl.final, freq = FALSE, breaks = seq ( 0, maxDiv+0.5, 0.5 ) );

		#	Characterize: Select a few cells for a closer look.
		#	A few different options.
x11(); par(mfrow=c(3,2)); ShowVecAsMap1 ( (div.selStim.final > sqrt(2.0)), titleText, xlab, ylab );
div.base.final.cellList = which ( div.selStim.final > sqrt(2.0) );
div.base.final.cellList = deepLookList = seq ( 1, N2, 1 );
symbolChoice = c ( 1, 0, 2 );

		#	Topo Quant figures

base.q.e = TopoQuantSeries ( N, g0, base$rfTrackData.e );
base.q.i = TopoQuantSeries ( N, g0, base$rfTrackData.i );

selStim.q.e = TopoQuantSeries ( N, g0, selStim$rfTrackData.e );
selStim.q.i = TopoQuantSeries ( N, g0, selStim$rfTrackData.i );

selStim.ctl.q.e = TopoQuantSeries ( N, g0, selStim.ctl$rfTrackData.e );
selStim.ctl.q.i = TopoQuantSeries ( N, g0, selStim.ctl$rfTrackData.i );

x1.e = c ( base.q.e$q1, selStim.q.e$q1, selStim.ctl.q.e$q1 );
x1.i = c ( base.q.i$q1, selStim.q.i$q1, selStim.ctl.q.i$q1 );

x11(); par(mfrow=c(2,1));
ylim = c ( min ( x1.e, x1.i ), max ( x1.e, x1.i ) );
plot ( x1.e, type="l", col=1, ylim=ylim, xlab="Iteration", ylab="C Metric",
	main="(F,G) Topography Measure vs Refinement Iterations\nE Cell Layer (Black); I Cell Layer (Red)" );
lines ( x1.i, col=2 );
abline ( v=c(26,51), lty=4, col=4 );
abline ( h=0.16, lty=4, col=4 );

x2.e = c ( base.q.e$q2, selStim.q.e$q2, selStim.ctl.q.e$q2 );
x2.i = c ( base.q.i$q2, selStim.q.i$q2, selStim.ctl.q.i$q2 );

ylim = c ( min ( x2.e, x2.i ), max ( x2.e, x2.i ) );
plot ( x2.e, type="l", col=1, ylim=ylim, xlab="Iteration", ylab="C Metric",
	main="(F,G) Topography Measure vs Refinement Iterations\nE Cell Layer (Black); I Cell Layer (Red)" );
lines ( x2.i, col=2 );
abline ( v=c(26,51), lty=4, col=4 );

} # if ( 0 ) Some of these figures are more useful than others.














