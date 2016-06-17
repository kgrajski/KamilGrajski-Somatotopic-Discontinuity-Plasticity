#
#	This script started out as one thing and has wound up as the
#	generator of the Figures for the draft publication and some
#	additional sanity check plots that could be used in a talk
#	as back up material, but not in a paper.
#

#
#	Kamil A. Grajski
#	Copyright 2016
#

#
#
#	Optimized for SYNDACTYLY paper.
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
		#	tiffFlag is ready to publish
		#
tiffFlag = TRUE;

		#
		#	Constants
		#
kRFPeakToEdgeDetect = 0.5;
g0 = 3;
N = 45;
N2 = N * N;

		#
		#	Set the directory where the experimental data is sitting.
		#	
fDir = "E:/NMLab/Working/T.45.7.Rand.1/";
fDir = "D:/NMLab/Keep - Syndactyly Exps/T.7/";

		#
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

		#
		#	Figure 1.0.  Stimcount plots.
		#
if ( tiffFlag ) { tiff(filename="Fig1.tiff",compression="lzw",units="in",width=7.0,height=7.0,res=600); } else { x11(); }
par(mfrow=c(2,2));
xlab="Dist.->Prox."; ylab="D1 -> D3";
titleText = paste ( "Stimulation Count\nBaseline Cycles",sep="" );
ShowVecAsMap1G ( base.stimCount, 32, titleText, xlab, ylab);
titleText = paste ( "Stimulation Count\nSyndactyly Cycles",sep="" );
ShowVecAsMap1G ( syndact.stimCount, 32, titleText, xlab, ylab, FALSE, "n", "n" );
if ( tiffFlag ) { dev.off(); }

		#
		#	Figure 2.0.  Different script (time series).
		#

		#
		#	Figure 3.0.  Plot of E-Cell Centroids
		#
if ( tiffFlag ) { tiff(filename="Fig3.tiff",compression="lzw",units="in",width=7.0,height=7.0,res=600); } else { x11(); }
par(mfrow=c(2,2), pin=c(2,2));
subTitleText = " Baseline";
iWhich = 1; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",0,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
subTitleText = " Syndactyly";
iWhich = iTop.syndact; ShowTopoMap1 ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
subTitleText = " Release";
iWhich = iTop.syndactRel; ShowTopoMap1 ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
if ( tiffFlag ) { dev.off(); }

		#
		#	Figure 4.0.  Plot of I-Cell Centroids
		#
if ( tiffFlag ) { tiff(filename="Fig4.tiff",compression="lzw",units="in",width=7.0,height=7.0,res=600); } else { x11(); }
par(mfrow=c(2,2), pin=c(2,2));
subTitleText = " Baseline";
iWhich = 1; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",0,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
subTitleText = " Syndactyly";
iWhich = iTop.syndact; ShowTopoMap1 ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
subTitleText = " Release";
iWhich = iTop.syndactRel; ShowTopoMap1 ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
if ( tiffFlag ) { dev.off(); }

		#
		#	Figure 5.0.	 Plot of E-Cell Tracks RF Extents
		#
if ( tiffFlag ) { tiff(filename="Fig5.tiff",compression="lzw",units="in",width=7.0,height=7.0,res=600); } else { x11(); }
par(mfrow=c(2,2), pin=c(2,2));
subTitleText = " Baseline";
iWhich = 1; ShowThreeDigitRFTrack1E ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",0,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
iWhich = iTop.base; ShowThreeDigitRFTrack1E ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
subTitleText = " Syndactyly";
iWhich = iTop.syndact; ShowThreeDigitRFTrack1E ( syndact$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
subTitleText = " Release";
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack1E ( syndactRel$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
if ( tiffFlag ) { dev.off(); }

		#
		#	Figure 6.0.	 Plot of I-Cell Tracks RF Extents
		#
if ( tiffFlag ) { tiff("Fig6.tiff",compression="lzw",units="in",width=7.0,height=7.0,res=600); } else { x11(); }
par(mfrow=c(2,2), pin=c(2,2));
subTitleText = " Baseline";
iWhich = 1; ShowThreeDigitRFTrack1E ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",0,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.base; ShowThreeDigitRFTrack1E ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
subTitleText = " Syndactyly";
iWhich = iTop.syndact; ShowThreeDigitRFTrack1E ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
subTitleText = " Release";
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack1E ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
if ( tiffFlag ) { dev.off(); }

		#
		#	Figure 7.0.	 Plot of I-Cell Tracks RF Extents II
		#
if ( tiffFlag ) { tiff(filename="Fig7.tiff",compression="lzw",units="in",width=7.0,height=7.0,res=600); } else { x11(); }
par(mfrow=c(2,2), pin=c(2,2));
subTitleText = " Baseline";
iWhich = 1; ShowThreeDigitRFTrack3AE ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iWhich-1,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
iWhich = iTop.base; ShowThreeDigitRFTrack3AE ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
subTitleText = " Syndactyly";
iWhich = iTop.syndact; ShowThreeDigitRFTrack3AE ( syndact$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
subTitleText = " Release";
iWhich = iTop.syndactRel; ShowThreeDigitRFTrack3AE ( syndactRel$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iBase,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 3 );
if ( tiffFlag ) { dev.off(); }

		#
		#	Figure 8.0.  Identify: RF Divergence Spatial Plot
		#
if ( tiffFlag ) { tiff(filename="Fig8.tiff",compression="lzw",units="in",width=7.0,height=7.0,res=600); } else { x11(); }
iWhich = 1; div.base.init = IntraColumnarRFDivergence ( base$rfTrackData.e[[iWhich]], base$rfTrackData.i[[iWhich]] );
iWhich = iTop.base; div.base.final = IntraColumnarRFDivergence ( base$rfTrackData.e[[iWhich]], base$rfTrackData.i[[iWhich]] );

iWhich = 1; div.syndact.init = IntraColumnarRFDivergence ( syndact$rfTrackData.e[[iWhich]], syndact$rfTrackData.i[[iWhich]] );
iWhich = iTop.syndact; div.syndact.final = IntraColumnarRFDivergence ( syndact$rfTrackData.e[[iWhich]], syndact$rfTrackData.i[[iWhich]] );

iWhich = 1; div.syndactRel.init = IntraColumnarRFDivergence ( syndactRel$rfTrackData.e[[iWhich]], syndactRel$rfTrackData.i[[iWhich]] );
iWhich = iTop.syndactRel; div.syndactRel.final = IntraColumnarRFDivergence ( syndactRel$rfTrackData.e[[iWhich]], syndactRel$rfTrackData.i[[iWhich]] );

minDiv = min ( c ( div.base.init, div.base.final, div.syndact.init, div.syndact.final, div.syndactRel.init, div.syndactRel.final ) );
maxDiv = max ( c ( div.base.init, div.base.final, div.syndact.init, div.syndact.final, div.syndactRel.init, div.syndactRel.final ) );

minDiv = 0; maxDiv = 5;
par(mfrow=c(2,2),pin=c(2.5,2.5),oma=c( 0,3,0,0 ),mar=c(4,2,3,1)+0.1);
xlab="Distal -> Proximal"; ylab="D1 -> D3";
xlab=""; ylab=""
boundaryMarks = c((N/3)+0.5, (2*N/3)+0.5);
titleText = paste ( "RF Divergence\nBaseline Iter ", (0), sep="" );
ShowVecAsMap2G ( div.base.init, 32, titleText, xlab, ylab, minDiv, maxDiv );
abline(h=boundaryMarks, lty=4, col=3);
mtext(side=1, "Distal -> Proximal", cex=0.75, line=0);
mtext(side=2, "D1 -> D3", cex=0.75, line=0);

titleText = paste ( "RF Divergence\nBaseline Iter ", (iBase), sep="" );
ShowVecAsMap2GImageOnly ( div.base.final, 32, titleText, xlab, ylab, minDiv, maxDiv );
abline(h=boundaryMarks, lty=4, col=3);

titleText = paste ( "RF Divergence\nSyndactyly Iter ", (iBase), sep="" );
ShowVecAsMap2GImageOnly ( div.syndact.final, 32, titleText, xlab, ylab, minDiv, maxDiv );
#abline(h=boundaryMarks, lty=4, col=3);

titleText = paste ( "RF Divergence\nRelease Iter ", (iBase), sep="" );
ShowVecAsMap2GImageOnly ( div.syndactRel.final, 32, titleText, xlab, ylab, minDiv, maxDiv );
abline(h=boundaryMarks, lty=4, col=3);
if ( tiffFlag ) { dev.off(); }

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

	#
	#	Initial Random -> Baseline
	#
tiffFlag = TRUE;
useLog = TRUE;
numColors = 32;
boundaryMarks = c((N/3)+0.5, (2*N/3)+0.5);
div.base.final.cellList = c ( 333, 648 );
for ( iCell in div.base.final.cellList ) {

	if ( tiffFlag ) {
		tiff ( filename=paste("RF",iCell,"tiff",sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )
	par(mfrow=c(3,2),pin=c(1.5, 1.5),oma=c( 0,3,0,3 ),mar=c(0,2,3,1)+0.1);
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";
	xlab=""; ylab="";

	minE = min ( init.r1.e.rfMap[iCell,], final.r1.e.rfMap[iCell,], syndact.r1.e.rfMap[iCell,] );
	maxE = max ( init.r1.e.rfMap[iCell,], final.r1.e.rfMap[iCell,], syndact.r1.e.rfMap[iCell,] );
	minI = min ( init.r1.i.rfMap[iCell,], final.r1.i.rfMap[iCell,], syndact.r1.i.rfMap[iCell,] );
	maxI = max ( init.r1.i.rfMap[iCell,], final.r1.i.rfMap[iCell,], syndact.r1.i.rfMap[iCell,] );

	iWhich = 1; titleText = paste ( " E Cell RF Extent\nBaseline Iter ", (iWhich-1), " Col ", iCell, sep="" );
	rfdat = init.r1.e.rfMap[iCell,];
	ShowVecAsMap2G ( rfdat, numColors, titleText, xlab, ylab, minE, maxE, FALSE, "n", "s", useLog );
	abline(h=boundaryMarks, lty=3, col="white", lwd=0.5);
	mtext(side=1, "Distal -> Proximal", cex=0.75, line=0);
	mtext(side=2, "D1 -> D3", cex=0.75, line=2);
	
	rfdat = init.r1.i.rfMap[iCell,];	
	iWhich = 1;	titleText = paste ( " I Cell RF Extent\nBaseline Iter ", (iWhich-1), " Col ", iCell, sep="" );
	ShowVecAsMap2G ( rfdat, numColors, titleText, xlab, ylab, minI, maxI, FALSE, "n", "n", useLog );
	abline ( h = boundaryMarks, lty=3, col=3 );

	rfdat = final.r1.e.rfMap[iCell,];
	iWhich = iBase; titleText = paste ( "Baseline Iter ", iWhich, " Col ", iCell, sep="" );
	ShowVecAsMap2G ( rfdat, numColors, titleText, xlab, ylab, minE, maxE, FALSE, "n", "n", useLog );
	abline ( h = boundaryMarks, lty=3, col=3 );

	rfdat = final.r1.i.rfMap[iCell,];
	iWhich = iBase; titleText = paste ( "Baseline Iter ", iWhich, " Col ", iCell, sep="" );
	ShowVecAsMap2G ( rfdat, numColors, titleText, xlab, ylab, minI, maxI, FALSE, "n", "n", useLog );
	abline ( h = boundaryMarks, lty=3, col=3 );

	rfdat = syndact.r1.e.rfMap[iCell,];
	iWhich = iBase; titleText = paste ( "Syndact Iter ", iWhich, " Col ", iCell, sep="" );
	ShowVecAsMap2G ( rfdat, numColors, titleText, xlab, ylab, minE, maxE, FALSE, "n", "n", useLog );
	abline ( h = boundaryMarks[2], lty=3, col=3 );

	rfdat = syndact.r1.i.rfMap[iCell,];
	iWhich = iBase; titleText = paste ( "Syndact Iter ", iWhich, " Col ", iCell, sep="" );
	ShowVecAsMap2G ( rfdat, numColors, titleText, xlab, ylab, minI, maxI, FALSE, "n", "n", useLog );
	abline ( h = boundaryMarks[2], lty=3, col=3 );

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

} # for ( iCell in div.base.final.cellList ) {

	#
	# 	Weights
	#
digitBorderline = c ( 0, 4.5 );
iBorder = 1;
useLog = FALSE;
numColors = 32;
for ( iCell in div.base.final.cellList ) {

	if ( tiffFlag ) {
		tiff ( filename=paste("W", iCell, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
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
	ShowVecAsMap1G ( vTmp, numColors, titleText, xlab, ylab, iHorizontal, "s", "s" );
	#mtext(titleText,side=3,line=0,font=2,cex=1.0);
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 3, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 1;
	titleText = paste (" E<-E Initial " );
	vTmp = GetInputWeights ( init.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1G ( vTmp, numColors, titleText, "", "", iHorizontal, "n", "n" );

	iWhich = 1;
	titleText = paste (" E<-I Initial " );
	vTmp = GetInputWeights ( init.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1G ( vTmp, numColors, titleText, "", "", iHorizontal, "n", "n" );

	iWhich = 1;
	titleText = paste (" I<-E Initial " );
	vTmp = GetInputWeights ( init.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1G ( vTmp, numColors, titleText, "", "", iHorizontal, "n", "n" );

		#	Refined Map

	iWhich = 16;
	titleText = titleText = paste (" E<-S Refined " );
	vTmp = GetInputWeights ( refine.w1.e0, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1G ( vTmp, numColors, titleText, "", "", iHorizontal, "n", "n"  );
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 3, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 16;
	titleText = paste (" E<-E Refined " );
	vTmp = GetInputWeights ( refine.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1G( vTmp, numColors, titleText, "", "", iHorizontal, "n", "n"  );

	iWhich = 16;
	titleText = paste (" E<-I Refined " );
	vTmp = GetInputWeights ( refine.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1G ( vTmp, numColors, titleText, "", "", iHorizontal, "n", "n"  );

	iWhich = 16;
	titleText = paste (" I<-E Refined " );
	vTmp = GetInputWeights ( refine.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1G ( vTmp, numColors, titleText, "", "", iHorizontal, "n", "n"  );

		#	Syndactyly

	iWhich = 15;
	titleText = paste (" E<-S Synd. " );
	vTmp = GetInputWeights ( syndact.w1.e0, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1G ( vTmp, numColors, titleText, "", "", iHorizontal, "n", "n" );
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 3, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 15;
	titleText = paste (" E<-E Synd. " );
	vTmp = GetInputWeights ( syndact.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1G ( vTmp, numColors, titleText, "", "", iHorizontal, "n", "n" );

	iWhich = 15;
	titleText = paste (" E<-I Synd. " );
	vTmp = GetInputWeights ( syndact.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1G ( vTmp, numColors, titleText, "", "", iHorizontal, "n", "n" );

	iWhich = 15;
	titleText = paste (" I<-E Synd. " );
	vTmp = GetInputWeights ( syndact.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1G ( vTmp, numColors, titleText, "", "", iHorizontal, "n", "n" );

	iBorder = iBorder + 1;

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

} # for ( iCell in div.base.final.cellList ) {






