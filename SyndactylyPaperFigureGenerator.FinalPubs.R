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
g0 = 3;
N = 45;
N2 = N * N;

		#	Set the directory where the experimental data is sitting.

fDir = "E:/NMLab/Working/T.45.7.Rand.1/";
fDir = "D:/NMLab/Keep - Syndactyly Exps/T.7/";

		
		#	Get the Baseline Refinement data.
		#
fRoot = "Base.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");

iBase = 15;
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
iBase = 15;
iSequence = c ( 1, iBase );
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
iBase = 15;
iSequence = c ( 1, iBase );
syndactRel = GetRFMapData3 ( fileRootName, iBase, iSequence, kRFPeakToEdgeDetect, N2 );
iTop.syndactRel = length(syndactRel$rfTrackData.e);


		#	Figure 0a.  Stimcount plots.
#if ( tiffFlag ) { tiff("Fig3.tif",compression="lzw",units="in",width=7.0,height=7.0,res=300); } else { x11(); }
x11(); par(mfrow=c(2,2));
#titleText = paste ( "Stimulation Count\nBaseline Cycles\nMin(Blue)=", min(base.stimCount), " Max(Red)=", max(base.stimCount),sep="" );
titleText = paste ( "Stimulation Count\nBaseline Cycles",sep="" );
xlab="Dist.->Prox."; ylab="D1 -> D3";
ShowVecAsMap ( base.stimCount, titleText);
#x11(); par(mfrow=c(2,2));
titleText = paste ( "Stimulation Count\nSyndactyly Cycles",sep="" );
ShowVecAsMap1 ( syndact.stimCount, titleText, xlab, ylab );
#if ( tiffFlag ) { dev.off(); }


		#	Figure 1a.  Plot of E-Cell Centroids

x11(); par(mfrow=c(2,2), pin=c(2,2));
#x11(); layout( matrix(seq(1,6),3,2,byrow=TRUE), widths=c(lcm(4), lcm(4), lcm(4), lcm(4), lcm(4), lcm(4)), heights=c(lcm(4), lcm(4), lcm(4), lcm(4), lcm(4), lcm(4)), respect=FALSE );

subTitleText = " Baseline";
iWhich = 1; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",0,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = " Syndactyly";
#iWhich = 1; ShowTopoMap1 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.syndact; ShowTopoMap1 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

#x11(); par(mfrow=c(2,2));
subTitleText = " Release";
#iWhich = 1; ShowTopoMap1 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.syndactRel; ShowTopoMap1 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

		#	Figure 1b.  Plot of I-Cell Centroids
x11(); par(mfrow=c(2,2), pin=c(2,2));
subTitleText = " Baseline";
#iWhich = 1; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",0,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = " Syndactyly";
iWhich = 1; ShowTopoMap1 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.syndact; ShowTopoMap1 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

#x11(); par(mfrow=c(2,2));
subTitleText = " Release";
#iWhich = 1; ShowTopoMap1 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.syndactRel; ShowTopoMap1 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = "\nBaseline Refinement";
x11(); iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

		#	Figure 2a.	Plot of E-Cell Tracks RF Extents
x11(); par(mfrow=c(2,2), pin=c(2,2));
subTitleText = " Baseline";
iWhich = 1; ShowThreeDigitRFTrack1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",0,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
iWhich = iTop.base; ShowThreeDigitRFTrack1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );

subTitleText = " Syndactyly";
#iWhich = 1; ShowThreeDigitRFTrack1 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndact; ShowThreeDigitRFTrack1 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );

subTitleText = " Release";
#iWhich = 1; ShowThreeDigitRFTrack1 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack1 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );

		#	Figure 2b.	Plot of I-Cell Tracks RF Extents
x11(); par(mfrow=c(2,2), pin=c(2,2));
subTitleText = " Baseline";
iWhich = 1; ShowThreeDigitRFTrack1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",0,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.base; ShowThreeDigitRFTrack1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = " Syndactyly";
#iWhich = 1; ShowThreeDigitRFTrack1 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndact; ShowThreeDigitRFTrack1 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = " Release";
#iWhich = 1; ShowThreeDigitRFTrack1 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack1 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

		#	Figure 3a.	Plot of E-Cell Tracks RF Extents II
x11(); par(mfrow=c(2,2));
subTitleText = " Baseline";
iWhich = 1; ShowThreeDigitRFTrack2 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.base; ShowThreeDigitRFTrack2 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = " Syndactyly";
#iWhich = 1; ShowThreeDigitRFTrack2 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndact; ShowThreeDigitRFTrack2 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = " Release";
#iWhich = 1; ShowThreeDigitRFTrack2 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack2 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

		#	Figure 3b.	Plot of I-Cell Tracks RF Extents II
x11(); par(mfrow=c(2,2), pin=c(2,2));
subTitleText = " Baseline";
iWhich = 1; ShowThreeDigitRFTrack3A ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
iWhich = iTop.base; ShowThreeDigitRFTrack3A ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );

subTitleText = " Syndactyly";
#iWhich = 1; ShowThreeDigitRFTrack3A ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
iWhich = iTop.syndact; ShowThreeDigitRFTrack3A ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );

subTitleText = " Release";
#iWhich = 1; ShowThreeDigitRFTrack3A ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack3A ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );

		#	Figure 4a.	Plot of E-Cell Tracks RF Extents III (Dist, Mid, Prox - across digits)
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowThreeDigitRFTrack3 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
iWhich = iTop.base; ShowThreeDigitRFTrack3 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 2 );

subTitleText = "\nSyndactyly";
iWhich = 1; ShowThreeDigitRFTrack3 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
iWhich = iTop.syndact; ShowThreeDigitRFTrack3 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 2 );

subTitleText = "\nSyndactyly Release";
iWhich = 1; ShowThreeDigitRFTrack3 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iWhich,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack3 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );

		#	Figure 4b.	Plot of I-Cell Tracks RF Extents II
x11(); par(mfrow=c(3,2));
subTitleText = "\nBaseline Refinement";
iWhich = 1; ShowThreeDigitRFTrack3 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",0,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
iWhich = iTop.base; ShowThreeDigitRFTrack3 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 2 );

subTitleText = "\nSyndactyly";
iWhich = 1; ShowThreeDigitRFTrack3 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",0,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
iWhich = iTop.syndact; ShowThreeDigitRFTrack3 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 2 );

subTitleText = "\nSyndactyly Release";
iWhich = 1; ShowThreeDigitRFTrack3 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",0,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack3 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );

		#
		#	Figure 5.  Identify: RF Divergence Spatial Plot
		#

iWhich = 1; div.base.init = IntraColumnarRFDivergence ( base$rfTrackData.e[[iWhich]], base$rfTrackData.i[[iWhich]] );
iWhich = iTop.base; div.base.final = IntraColumnarRFDivergence ( base$rfTrackData.e[[iWhich]], base$rfTrackData.i[[iWhich]] );

iWhich = 1; div.syndact.init = IntraColumnarRFDivergence ( syndact$rfTrackData.e[[iWhich]], syndact$rfTrackData.i[[iWhich]] );
iWhich = iTop.syndact; div.syndact.final = IntraColumnarRFDivergence ( syndact$rfTrackData.e[[iWhich]], syndact$rfTrackData.i[[iWhich]] );

iWhich = 1; div.syndactRel.init = IntraColumnarRFDivergence ( syndactRel$rfTrackData.e[[iWhich]], syndactRel$rfTrackData.i[[iWhich]] );
iWhich = iTop.syndactRel; div.syndactRel.final = IntraColumnarRFDivergence ( syndactRel$rfTrackData.e[[iWhich]], syndactRel$rfTrackData.i[[iWhich]] );

minDiv = min ( c ( div.base.init, div.base.final, div.syndact.init, div.syndact.final, div.syndactRel.init, div.syndactRel.final ) );
maxDiv = max ( c ( div.base.init, div.base.final, div.syndact.init, div.syndact.final, div.syndactRel.init, div.syndactRel.final ) );

minDiv = 0; maxDiv = 5;

x11();
#par ( mfrow = c ( 2,2 ) );
par(mfrow=c(2,2),pin=c(2.5,2.5),oma=c( 0,3,0,0 ),mar=c(4,2,3,1)+0.1);
xlab="Distal -> Proximal"; ylab="D1 -> D3";
xlab=""; ylab=""
boundaryMarks = c((N/3)+0.5, (2*N/3)+0.5);
titleText = paste ( "RF Divergence\nBaseline Iter ", (0), sep="" );
ShowVecAsMap2 ( div.base.init, titleText, xlab, ylab, minDiv, maxDiv );
abline(h=boundaryMarks, lty=4, col=3);
mtext(side=1, "Distal -> Proximal", cex=0.75, line=2);
mtext(side=2, "D1 -> D3", cex=0.75, line=2);

titleText = paste ( "RF Divergence\nBaseline Iter ", (iBase), sep="" );
ShowVecAsMap2ImageOnly ( div.base.final, titleText, xlab, ylab, minDiv, maxDiv );
abline(h=boundaryMarks, lty=4, col=3);

titleText = paste ( "RF Divergence\nSyndactyly Iter ", (iBase), sep="" );
ShowVecAsMap2ImageOnly ( div.syndact.final, titleText, xlab, ylab, minDiv, maxDiv );
abline(h=boundaryMarks, lty=4, col=3);

titleText = paste ( "RF Divergence\nRelease Iter ", (iBase), sep="" );
ShowVecAsMap2ImageOnly ( div.syndactRel.final, titleText, xlab, ylab, minDiv, maxDiv );
abline(h=boundaryMarks, lty=4, col=3);

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

#
#	Deep dive on RFs and Weights
#

init.r1.e.rfMap = base$r1.e.rfMap[[1]];
init.r1.i.rfMap = base$r1.i.rfMap[[1]];
final.r1.e.rfMap = base$r1.e.rfMap[[iTop.base]];
final.r1.i.rfMap = base$r1.i.rfMap[[iTop.base]];
syndact.r1.e.rfMap = syndact$r1.e.rfMap[[iTop.syndact]];
syndact.r1.i.rfMap = syndact$r1.i.rfMap[[iTop.syndact]];

	#
	#	RF Property Scatterplots as a health check.
	#

init.rfMap.rfAreas.e = QuantRFSize ( init.r1.e.rfMap, kRFPeakToEdgeDetect );
init.rfMap.rfAreas.i = QuantRFSize ( init.r1.i.rfMap, kRFPeakToEdgeDetect );
final.rfMap.rfAreas.e = QuantRFSize ( final.r1.e.rfMap, kRFPeakToEdgeDetect );
final.rfMap.rfAreas.i = QuantRFSize ( final.r1.i.rfMap, kRFPeakToEdgeDetect );
syndact.rfMap.rfAreas.e = QuantRFSize ( syndact.r1.e.rfMap, kRFPeakToEdgeDetect );
syndact.rfMap.rfAreas.i = QuantRFSize ( syndact.r1.i.rfMap, kRFPeakToEdgeDetect );

init.rfMap.maxResp.e = apply ( init.r1.e.rfMap, 1, max );
init.rfMap.maxResp.i = apply ( init.r1.i.rfMap, 1, max );
final.rfMap.maxResp.e = apply ( final.r1.e.rfMap, 1, max );
final.rfMap.maxResp.i = apply ( final.r1.i.rfMap, 1, max );
syndact.rfMap.maxResp.e = apply ( syndact.r1.e.rfMap, 1, max );
syndact.rfMap.maxResp.i = apply ( syndact.r1.i.rfMap, 1, max );

x11(); par(mfrow=c(2,2));
initText = "Initial Random"; finalText = "Baseline Refined";
xlim.e = c ( min(init.rfMap.maxResp.e, final.rfMap.maxResp.e), max(init.rfMap.maxResp.e, final.rfMap.maxResp.e) );
ylim.e = c ( min(init.rfMap.rfAreas.e, final.rfMap.rfAreas.e), max(init.rfMap.rfAreas.e, final.rfMap.rfAreas.e) );
plot ( init.rfMap.maxResp.e, init.rfMap.rfAreas.e, xlim=xlim.e, ylim=ylim.e, pch=".", col=1, cex=2,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer E-Cells\n",initText,"(.); ",finalText,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.e, final.rfMap.rfAreas.e, pch="*", col=2, cex=1);

xlim.i = c ( min(init.rfMap.maxResp.i, final.rfMap.maxResp.i), max(init.rfMap.maxResp.i, final.rfMap.maxResp.i) );
ylim.i = c ( min(init.rfMap.rfAreas.i, final.rfMap.rfAreas.i), max(init.rfMap.rfAreas.i, final.rfMap.rfAreas.i) );
plot ( init.rfMap.maxResp.i, init.rfMap.rfAreas.i, xlim=xlim.i, ylim=ylim.i, pch=".", col=1, cex=2,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer I-Cells\n",initText,"(.); ",finalText,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.i, final.rfMap.rfAreas.i, pch="*", col=2, cex=1);

	#
	#	Deep Dive RF Plots
	#

init.r1.e.rfMap = base$r1.e.rfMap[[1]];
init.r1.i.rfMap = base$r1.i.rfMap[[1]];
final.r1.e.rfMap = base$r1.e.rfMap[[iTop.base]];
final.r1.i.rfMap = base$r1.i.rfMap[[iTop.base]];
syndact.r1.e.rfMap = syndact$r1.e.rfMap[[iTop.syndact]];
syndact.r1.i.rfMap = syndact$r1.i.rfMap[[iTop.syndact]];

#div.base.final.cellList = which ( div.base.final > 2*sqrt(2.0) );
div.base.final.cellList = c ( 423, 648 );
#div.base.final.cellList = seq ( 13*45+1, 17*45, 3 );
boundaryMarks = c((N/3)+0.5, (2*N/3)+0.5);

		#
		#	Initial Random -> Baseline
		#
tiffFlag = FALSE;
for ( iCell in div.base.final.cellList ) {

	if ( tiffFlag ) {
		tiff ( paste("DeepRF", iCell, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )
	par(mfrow=c(3,2),pin=c(1.5, 1.5),oma=c( 0,3,0,3 ),mar=c(0,2,3,1)+0.1);
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";
	xlab=""; ylab="";

	iWhich = 1; titleText = paste ( " E Cell RF Extent\nBaseline Iter ", (iWhich-1), " Col ", iCell, sep="" );
	rfdat = init.r1.e.rfMap[iCell,];
	ShowVecAsMap1 ( rfdat, titleText, xlab, ylab, FALSE, "n", "s" );
	abline(h=boundaryMarks, lty=4, col=3);
	mtext(side=1, "Distal -> Proximal", cex=0.75, line=0);
	mtext(side=2, "D1 -> D3", cex=0.75, line=2);
	
	rfdat = init.r1.i.rfMap[iCell,];	
	iWhich = 1;	titleText = paste ( " I Cell RF Extent\nBaseline Iter ", (iWhich-1), " Col ", iCell, sep="" );
	ShowVecAsMap1 ( rfdat, titleText, xlab, ylab, FALSE, "n", "n" );
	abline ( h = boundaryMarks, lty=3, col=3 );

	rfdat = final.r1.e.rfMap[iCell,];
	iWhich = iBase; titleText = paste ( "Baseline Iter ", iWhich, " Col ", iCell, sep="" );
	ShowVecAsMap1 ( rfdat, titleText, xlab, ylab, FALSE, "n", "n" );
	abline ( h = boundaryMarks, lty=3, col=3 );

	rfdat = final.r1.i.rfMap[iCell,];
	iWhich = iBase; titleText = paste ( "Baseline Iter ", iWhich, " Col ", iCell, sep="" );
	ShowVecAsMap1 ( rfdat, titleText, xlab, ylab, FALSE, "n", "n" );
	abline ( h = boundaryMarks, lty=3, col=3 );

	rfdat = syndact.r1.e.rfMap[iCell,];
	iWhich = iBase; titleText = paste ( "Syndact Iter ", iWhich, " Col ", iCell, sep="" );
	ShowVecAsMap1 ( rfdat, titleText, xlab, ylab, FALSE, "n", "n" );
	abline ( h = boundaryMarks[2], lty=3, col=3 );

	rfdat = syndact.r1.i.rfMap[iCell,];
	iWhich = iBase; titleText = paste ( "Syndact Iter ", iWhich, " Col ", iCell, sep="" );
	ShowVecAsMap1 ( rfdat, titleText, xlab, ylab, FALSE, "n", "n" );
	abline ( h = boundaryMarks[2], lty=3, col=3 );

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

} # for ( iCell in div.base.final.cellList ) {

		#
		#	Baseline -> Syndactyly
		#
for ( iCell in div.base.final.cellList ) {

	x11(); par(mfcol=c(2,2));
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

	iWhich = iBase;
	titleText = paste ( "Column # ", iCell, " E Cell RF Extent\nBaseline Refinement Iter ", iWhich, sep="" );
	ShowVecAsMap1 ( final.r1.e.rfMap[iCell,], titleText, xlab, ylab );
	abline ( h = boundaryMarks, lty=3, col=3 );

	iWhich = iBase;
	titleText = paste ( "Column # ", iCell, " I Cell RF Extent\nBaseline Refinement Iter ", iWhich, sep="" );
	ShowVecAsMap1 ( final.r1.i.rfMap[iCell,], titleText, xlab, ylab );
	abline ( h = boundaryMarks, lty=3, col=3 );

	iWhich = iBase;
	titleText = paste ( "Column # ", iCell, " E Cell RF Extent\nSyndactyly Iter ", iWhich, sep="" );
	ShowVecAsMap1 ( syndact.r1.e.rfMap[iCell,], titleText, xlab, ylab );
	abline ( h = boundaryMarks, lty=3, col=3 );

	iWhich = iBase;
	titleText = paste ( "Column # ", iCell, " I Cell RF Extent\nSyndactyly Iter ", iWhich, sep="" );
	ShowVecAsMap1 ( syndact.r1.i.rfMap[iCell,], titleText, xlab, ylab );
	abline ( h = boundaryMarks, lty=3, col=3 );

} # for ( iCell in div.base.final.cellList ) {

	#
	# 	Weights
	#
digitBorderline = c ( 0, 4.5 );
iBorder = 1;

for ( iCell in div.base.final.cellList ) {

	if ( tiffFlag ) {
		tiff ( paste("W", iCell, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )
	par(mfrow=c(3,4),pin=c(0.8,0.8),oma=c( 0,0,0,3 ),mar=c(4,2,3,1)+0.1);
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";
	xlab=""; ylab="";
	iHorizontal = FALSE;

		#	Initial Conditions

	iWhich = 1;	titleText = paste (" E<-S Initial\nColumn ", iCell, sep="" );
	vTmp = GetInputWeights ( init.w1.e0, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1NoLog ( vTmp, titleText, xlab, ylab, iHorizontal, "s", "s" );
	#mtext(titleText,side=3,line=0,font=2,cex=1.0);
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 3, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 1;
	titleText = paste (" E<-E Initial " );
	vTmp = GetInputWeights ( init.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1NoLog ( vTmp, titleText, "", "", iHorizontal );

	iWhich = 1;
	titleText = paste (" E<-I Initial " );
	vTmp = GetInputWeights ( init.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1NoLog ( vTmp, titleText, "", "", iHorizontal );

	iWhich = 1;
	titleText = paste (" I<-E Initial " );
	vTmp = GetInputWeights ( init.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1NoLog ( vTmp, titleText, "", "", iHorizontal );

		#	Refined Map

	iWhich = 16;
	titleText = titleText = paste (" E<-S Refined " );
	vTmp = GetInputWeights ( refine.w1.e0, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1NoLog ( vTmp, titleText, "", "", iHorizontal );
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 3, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 16;
	titleText = paste (" E<-E Refined " );
	vTmp = GetInputWeights ( refine.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1NoLog ( vTmp, titleText, "", "", iHorizontal );

	iWhich = 16;
	titleText = paste (" E<-I Refined " );
	vTmp = GetInputWeights ( refine.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1NoLog ( vTmp, titleText, "", "", iHorizontal );

	iWhich = 16;
	titleText = paste (" I<-E Refined " );
	vTmp = GetInputWeights ( refine.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1NoLog ( vTmp, titleText, "", "", iHorizontal );

		#	Syndactyly

	iWhich = 15;
	titleText = paste (" E<-S Synd. " );
	vTmp = GetInputWeights ( syndact.w1.e0, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1NoLog ( vTmp, titleText, "", "", iHorizontal );
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 3, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 15;
	titleText = paste (" E<-E Synd. " );
	vTmp = GetInputWeights ( syndact.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1NoLog ( vTmp, titleText, "", "", iHorizontal );

	iWhich = 15;
	titleText = paste (" E<-I Synd. " );
	vTmp = GetInputWeights ( syndact.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1NoLog ( vTmp, titleText, "", "", iHorizontal );

	iWhich = 15;
	titleText = paste (" I<-E Synd. " );
	vTmp = GetInputWeights ( syndact.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1NoLog ( vTmp, titleText, "", "", iHorizontal );

	iBorder = iBorder + 1;

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

} # for ( iCell in div.base.final.cellList ) {






