
######################################################################################################
######################################################################################################
#
#	CorticalMagnificationAnalysis.R
#
#	Automate the analysis of cortical magnification factor and receptive field extent analysis.
#	Visualization and tests of statistical significance.
#
######################################################################################################
######################################################################################################

#	Clear the workspace.
if ( !alreadyCleared ) { rm(list = ls()); }

source ( "NMHelperFunctions.R" );
source ( "Dig3SectorLocs1.R" );
source ( "D3LongAxisLocs1.R" );
source ( "D3CorticalMagnificationPlots.R" );
source ( "D3RFDataBoxPlots.R" );
source ( "InvCortMagRuleCheck.R" );
source ( "MapSectorNumberToLabel.R" );
source ( "D3MannWhitneyTestFrontEnd.R" );
source ( "LabelledCDF.R" );
source ( "MapMagToCellList.R" );

		#
		#	Some additional variables needed if/when NMExpZone isn't the source.
N = 45;
N2 = N * N;
kRFPeakToEdgeDetect = 0.5;

		#
		#

E = 65;
P = "4.x64Check";
P = 4;
P = 5;

E = 65;
N = 45;
P = "5.EOnly";
#fDir = paste("E:\\NMLab\\E",E,".",N,".",P,sep="");
#fDir = "../";
#fDir =" E:/NMLab/Simulations/T.7/";

	######################################################
	#
	#	Get the network with random initial conditions.
	#
	######################################################
fName = paste(fDir, "\\", "Base.RFMap.15.0.bin", sep="" );
initText = "Refined";

	finfo = file.info ( fName );
	toread = file ( fName, "rb" );
	alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
	close ( toread );	
	itmp = length ( alldata ) / 2;
	N2 = as.integer ( sqrt ( itmp ) );
	init.r1.e.rfMap = r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
	init.r1.i.rfMap = r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

init.cortical.Amp.e = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
init.cortical.Amp.i = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );
init.rfMap.rfAreas.e = QuantRFSize ( r1.e.rfMap, kRFPeakToEdgeDetect );
init.rfMap.rfAreas.i = QuantRFSize ( r1.i.rfMap, kRFPeakToEdgeDetect );
init.rfMap.maxResp.e = apply ( r1.e.rfMap, 1, max );
init.rfMap.maxResp.i = apply ( r1.i.rfMap, 1, max );


	######################################################
	#
	#	Get the refined network.
	#
	######################################################
fName = paste(fDir, "\\", "Base.RFMap.15.15.bin", sep="" );
finalText = "Refined";

	finfo = file.info ( fName );
	toread = file ( fName, "rb" );
	alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
	close ( toread );	
	itmp = length ( alldata ) / 2;
	N2 = as.integer ( sqrt ( itmp ) );
	final.r1.e.rfMap = r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
	final.r1.i.rfMap = r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

final.cortical.Amp.e = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
final.cortical.Amp.i = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );
final.rfMap.rfAreas.e = QuantRFSize ( r1.e.rfMap, kRFPeakToEdgeDetect );
final.rfMap.rfAreas.i = QuantRFSize ( r1.i.rfMap, kRFPeakToEdgeDetect );
final.rfMap.maxResp.e = apply ( r1.e.rfMap, 1, max );
final.rfMap.maxResp.i = apply ( r1.i.rfMap, 1, max );

	######################################################
	#
	#	Get the syndactyly network.
	#
	######################################################
fName = paste(fDir, "\\", "SyndactExp.RFMap.15.15.bin", sep="" );
finalText2 = "Syndactyly";

	finfo = file.info ( fName );
	toread = file ( fName, "rb" );
	alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
	close ( toread );	
	itmp = length ( alldata ) / 2;
	N2 = as.integer ( sqrt ( itmp ) );
	final2.r1.e.rfMap = r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
	final2.r1.i.rfMap = r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

final2.cortical.Amp.e = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
final2.cortical.Amp.i = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );
final2.rfMap.rfAreas.e = QuantRFSize ( r1.e.rfMap, kRFPeakToEdgeDetect );
final2.rfMap.rfAreas.i = QuantRFSize ( r1.i.rfMap, kRFPeakToEdgeDetect );
final2.rfMap.maxResp.e = apply ( r1.e.rfMap, 1, max );
final2.rfMap.maxResp.i = apply ( r1.i.rfMap, 1, max );

	######################################################
	#
	#	Get the syndactyly control network.
	#
	######################################################
fName = paste(fDir, "\\", "SyndactCtl.RFMap.20.20.bin", sep="" );
finalText2Ctl = "SyndactylyCtl";

	finfo = file.info ( fName );
	toread = file ( fName, "rb" );
	alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
	close ( toread );	
	itmp = length ( alldata ) / 2;
	N2 = as.integer ( sqrt ( itmp ) );
	final2Ctl.r1.e.rfMap = r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
	final2Ctl.r1.i.rfMap = r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

final2Ctl.cortical.Amp.e = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
final2Ctl.cortical.Amp.i = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );
final2Ctl.rfMap.rfAreas.e = QuantRFSize ( r1.e.rfMap, kRFPeakToEdgeDetect );
final2Ctl.rfMap.rfAreas.i = QuantRFSize ( r1.i.rfMap, kRFPeakToEdgeDetect );
final2Ctl.rfMap.maxResp.e = apply ( r1.e.rfMap, 1, max );
final2Ctl.rfMap.maxResp.i = apply ( r1.i.rfMap, 1, max );
	

	######################################################
	#
	#	Get the syndactyly release network.
	#
	######################################################
fName = paste(fDir, "\\", "SyndactReleaseExp.RFMap.20.20.bin", sep="" );
finalText3 = "SyndactReleaseExp";

	finfo = file.info ( fName );
	toread = file ( fName, "rb" );
	alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
	close ( toread );	
	itmp = length ( alldata ) / 2;
	N2 = as.integer ( sqrt ( itmp ) );
	final3.r1.e.rfMap = r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
	final3.r1.i.rfMap = r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

final3.cortical.Amp.e = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
final3.cortical.Amp.i = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );
final3.rfMap.rfAreas.e = QuantRFSize ( r1.e.rfMap, kRFPeakToEdgeDetect );
final3.rfMap.rfAreas.i = QuantRFSize ( r1.i.rfMap, kRFPeakToEdgeDetect );
final3.rfMap.maxResp.e = apply ( r1.e.rfMap, 1, max );
final3.rfMap.maxResp.i = apply ( r1.i.rfMap, 1, max );

	######################################################
	#
	#	Get the syndactyly release control network.
	#
	######################################################

fName = paste(fDir, "\\", "SyndactReleaseCtl.RFMap.20.20.bin", sep="" );
finalText3Ctl = "SyndactReleaseCtl";

	finfo = file.info ( fName );
	toread = file ( fName, "rb" );
	alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
	close ( toread );	
	itmp = length ( alldata ) / 2;
	N2 = as.integer ( sqrt ( itmp ) );
	final3Ctl.r1.e.rfMap = r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
	final3Ctl.r1.i.rfMap = r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

final3Ctl.cortical.Amp.e = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
final3Ctl.cortical.Amp.i = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );
final3Ctl.rfMap.rfAreas.e = QuantRFSize ( r1.e.rfMap, kRFPeakToEdgeDetect );
final3Ctl.rfMap.rfAreas.i = QuantRFSize ( r1.i.rfMap, kRFPeakToEdgeDetect );
final3Ctl.rfMap.maxResp.e = apply ( r1.e.rfMap, 1, max );
final3Ctl.rfMap.maxResp.i = apply ( r1.i.rfMap, 1, max );


	######################################################
	#
	#	Data Visualizations.
	# 	Scatterplot: layer C E-Cell Area vs Resp.
	#
	######################################################
x11();
xlim.e = c ( min(init.rfMap.maxResp.e,final.rfMap.maxResp.e), max(init.rfMap.maxResp.e,final.rfMap.maxResp.e) );
ylim.e = c ( min(init.rfMap.rfAreas.e,final.rfMap.rfAreas.e), max(init.rfMap.rfAreas.e,final.rfMap.rfAreas.e) );
plot ( init.rfMap.maxResp.e, init.rfMap.rfAreas.e, xlim=xlim.e, ylim=ylim.e, pch=".", col=1, cex=4,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer E-Cells\n",initText,"(.); ",finalText,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.e, final.rfMap.rfAreas.e, pch="*", col=2, cex=2 );
delta.area.e = delta.maxResp.e = slope.e = rep ( 0, N2 );
for ( i in 1:N2 ) {
	#lines ( c( init.rfMap.maxResp.e[i], final.rfMap.maxResp.e[i] ), c( init.rfMap.rfAreas.e[i], final.rfMap.rfAreas.e[i] ) )
	delta.area.e[i] = ( final.rfMap.rfAreas.e[i] - init.rfMap.rfAreas.e[i] );
	delta.maxResp.e[i] = ( final.rfMap.maxResp.e[i] - init.rfMap.maxResp.e[i] );
	slope.e[i] = ( final.rfMap.rfAreas.e[i] - init.rfMap.rfAreas.e[i] ) / ( final.rfMap.maxResp.e[i] - init.rfMap.maxResp.e[i] );
} # for ( i in 1:N2 ) {

	##############################################################
	#
	#	Scatterplot: spatial mapping of scatterplot quadrants.
	#
	#		E Cells
	#
	##############################################################
x11();
par(mfrow=c(2,2));

ofInterest = which ( (final.rfMap.rfAreas.e > mean(final.rfMap.rfAreas.e)) & (final.rfMap.maxResp.e < mean(final.rfMap.maxResp.e) ) );
tmp = rep ( 0,N2); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Above Av RF Extent\nBelow Average Resp\nE Cells" );

ofInterest = which ( (final.rfMap.rfAreas.e > mean(final.rfMap.rfAreas.e)) & (final.rfMap.maxResp.e > mean(final.rfMap.maxResp.e) ) );
tmp = rep ( 0,N2); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Above Av RF Extent\nAbove Average Resp\nE Cells" );

ofInterest = which ( (final.rfMap.rfAreas.e < mean(final.rfMap.rfAreas.e)) & (final.rfMap.maxResp.e < mean(final.rfMap.maxResp.e) ) );
tmp = rep ( 0,N2); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Below Av RF Extent\nBelow Average Resp\nE Cells" );

ofInterest = which ( (final.rfMap.rfAreas.e < mean(final.rfMap.rfAreas.e)) & (final.rfMap.maxResp.e > mean(final.rfMap.maxResp.e) ) );
tmp = rep ( 0,N2); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Below Av RF Extent\nAbove Average Resp\nE Cells" );

	######################################################
	#
	# 	Scatterplot: Layer C I-Cell Area vs Resp.
	#
	######################################################
x11();
xlim.i = c ( min(init.rfMap.maxResp.i,final.rfMap.maxResp.i), max(init.rfMap.maxResp.i,final.rfMap.maxResp.i) );
ylim.i = c ( min(init.rfMap.rfAreas.i,final.rfMap.rfAreas.i), max(init.rfMap.rfAreas.i,final.rfMap.rfAreas.i) );
plot ( init.rfMap.maxResp.i, init.rfMap.rfAreas.i, xlim=xlim.i, ylim=ylim.i, pch=".", col=1, cex=4,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer I-Cells\n",initText,"(.); ",finalText,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.i, final.rfMap.rfAreas.i, pch="*", col=2, cex=2 );

#points ( final.rfMap.maxResp.i[div.base.final.cellList], final.rfMap.rfAreas.i[div.base.final.cellList],
	#pch="o", col=3, cex=3 );

delta.area.i = delta.maxResp.i = slope.i = rep ( 0, N2 );
for ( i in 1:N2 ) {
	#lines ( c( init.rfMap.maxResp.i[i], final.rfMap.maxResp.i[i] ), c( init.rfMap.rfAreas.i[i], final.rfMap.rfAreas.i[i] ) )
	delta.area.i[i] = ( final.rfMap.rfAreas.i[i] - init.rfMap.rfAreas.i[i] );
	delta.maxResp.i[i] = ( final.rfMap.maxResp.i[i] - init.rfMap.maxResp.i[i] );
	slope.i[i] = ( final.rfMap.rfAreas.i[i] - init.rfMap.rfAreas.i[i] ) / ( final.rfMap.maxResp.i[i] - init.rfMap.maxResp.i[i] );
} # for ( i in 1:N2 ) {


	##############################################################
	#
	#	Scatterplot: spatial mapping of scatterplot quadrants.
	#
	#		I Cells
	#
	##############################################################
x11();
par(mfrow=c(2,2));

ofInterest = which ( (final.rfMap.rfAreas.i > mean(final.rfMap.rfAreas.i)) & (final.rfMap.maxResp.i < mean(final.rfMap.maxResp.i) ) );
tmp = rep ( 0,N2); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Above Av RF Extent\nBelow Average Resp\nI Cells" );

ofInterest = which ( (final.rfMap.rfAreas.i > mean(final.rfMap.rfAreas.i)) & (final.rfMap.maxResp.i > mean(final.rfMap.maxResp.i) ) );
tmp = rep ( 0,N2); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Above Av RF Extent\nAbove Average Resp\nI Cells" );

ofInterest = which ( (final.rfMap.rfAreas.i < mean(final.rfMap.rfAreas.i)) & (final.rfMap.maxResp.i < mean(final.rfMap.maxResp.i) ) );
tmp = rep ( 0,N2); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Below Av RF Extent\nBelow Average Resp\nI Cells" );

ofInterest = which ( (final.rfMap.rfAreas.i < mean(final.rfMap.rfAreas.i)) & (final.rfMap.maxResp.i > mean(final.rfMap.maxResp.i) ) );
tmp = rep ( 0,N2); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Below Av RF Extent\nAbove Average Resp\nI Cells" );




	######################################################################
	#
 	# 	Scatterplot: Combined final state scatter plot E-cell, I-cell
	#
	######################################################################
x11();
xlim = c ( min ( xlim.i, xlim.e ), max ( xlim.i, xlim.e ) );
ylim = c ( min ( ylim.i, ylim.e ), max ( ylim.i, ylim.e ) );
plot ( final.rfMap.maxResp.e, final.rfMap.rfAreas.e, xlim=xlim, ylim=ylim, pch="*", col=1,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer E-Cells\nBaseline(Black *); Syndactyly( Red *)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.i, final.rfMap.rfAreas.i, pch="*", col=2 );

	############################################################################################
	#
	#	START
	#
	#	(INVERSE) CORTICAL MAGNIFICATION ANALYSIS
	#	This means partitioning the data by relative frequency of stimulation rather than
	#	sector or longitudonal axis.  It is the cleanest test.
	#
	############################################################################################
	
		#
		#	Extract from the stimulus pattern the relative frequencies of stimulation.
		#

fRoot = "CheckInputStim.Baseline";
fileRootName = paste(fDir, fRoot, sep="\\");
#mapStimCount = base.stimCount = GetStimCountData1 ( fileRootName, N2 );
mapStimCount = base.stimCount = GetStimCountData ( fileRootName, N2 );

fRoot = "CheckInputStim.Syndact";
fileRootName = paste(fDir, fRoot, sep="\\");
#mapStimCount.syndact = GetStimCountData1 ( fileRootName, N2 );
mapStimCount.syndact = GetStimCountData ( fileRootName, N2 );

		#
		#	Draw some pictures to visualize the data.
		#
xIDSubset = c ( 10, 20, 25 );

		#
		#	CDF: E Cell CORTICAL MAGNIFICATION & RF EXTENT:
		#		Initial (init) vs Baseline (final)
		#		Baseline (final) vs Syndactyly (final2)
		#
init.CMag.e = apply ( init.cortical.Amp.e, 2, sum ); final.CMag.e = apply ( final.cortical.Amp.e, 2, sum );
final2.CMag.e = apply ( final2.cortical.Amp.e, 2, sum ); final3.CMag.e = apply ( final3.cortical.Amp.e, 2, sum );

x11(); par(mfrow=c(3,2));

xLabelText = "Cortical Magnification"; yLabelText = "CDF"; zTypeText = "E Cell"
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
subText = "Init.(Dash); Base(Solid); 10x(B); 20x(R); 25x(G)";
CorticalMagCDF ( mapStimCount, init.CMag.e, final.CMag.e, titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

xLabelText = "RF Extent"; yLabelText = "CDF";
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
InvCorticalMagCDF ( mapStimCount, init.cortical.Amp.e, init.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e,
				titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

xLabelText = "Cortical Magnification"; yLabelText = "CDF"; zTypeText = "E Cell"
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
subText = "Base.(Dash); S-Rel.(Solid); 10x(B); 20x(R); 25x(G)";
CorticalMagCDF ( mapStimCount, final.CMag.e, final3.CMag.e, titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

xLabelText = "RF Extent"; yLabelText = "CDF";
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
InvCorticalMagCDF ( mapStimCount, final.cortical.Amp.e, final.rfMap.rfAreas.e, final3.cortical.Amp.e, final3.rfMap.rfAreas.e,
				titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );


		#
		#	CDF: CDF: I Cell CORTICAL MAGNIFICATION & RF EXTENT: Initial (init) vs Baseline (final)
		#
init.CMag.i = apply ( init.cortical.Amp.i, 2, sum ); final.CMag.i = apply ( final.cortical.Amp.i, 2, sum );
final2.CMag.i = apply ( final2.cortical.Amp.i, 2, sum ); final3.CMag.i = apply ( final3.cortical.Amp.i, 2, sum );

x11(); par(mfrow=c(2,2));

xLabelText = "Cortical Magnification"; yLabelText = "CDF"; zTypeText = "I Cell"
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
CorticalMagCDF ( mapStimCount, init.CMag.i, final.CMag.i, titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );


xLabelText = "RF Extent"; yLabelText = "CDF";
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
InvCorticalMagCDF ( mapStimCount, init.cortical.Amp.i, init.rfMap.rfAreas.i, final.cortical.Amp.i, final.rfMap.rfAreas.i,
				titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

 xLabelText = "Cortical Magnification"; yLabelText = "CDF"; zTypeText = "I Cell"
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
CorticalMagCDF ( mapStimCount, init.CMag.i, final.CMag.i, titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );


xLabelText = "RF Extent"; yLabelText = "CDF";
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
InvCorticalMagCDF ( mapStimCount, init.cortical.Amp.i, init.rfMap.rfAreas.i, final.cortical.Amp.i, final.rfMap.rfAreas.i,
				titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

		#
		#	WMM: CORTICAL MAGNIFICATION
		#
pValue = 1e-20;

init.CMag.MWTest.e = MWTestForCorticalMagnification ( mapStimCount, init.CMag.e, init.CMag.e, pValue );
init.CMag.MWTest.e$pValues[init.CMag.MWTest.e$acceptAltHyp,]

final.CMag.MWTest.e = MWTestForCorticalMagnification ( mapStimCount, final.CMag.e, final.CMag.e, pValue  );
final.CMag.MWTest.e$pValues[final.CMag.MWTest.e$acceptAltHyp,]

init.final.CMag.MWTest.e = MWTestForCorticalMagnification ( mapStimCount, init.CMag.e, final.CMag.e, pValue  );
init.final.CMag.MWTest.e$pValues[init.final.CMag.MWTest.e$acceptAltHyp,]

		#
		#	WMM: RECEPTIVE FIELD EXTENT.
		#
init.rfAreas.MWTest.e = MWTestForInvCorticalMagnification ( mapStimCount, init.cortical.Amp.e, init.rfMap.rfAreas.e, init.cortical.Amp.e, init.rfMap.rfAreas.e, pValue );
init.rfAreas.MWTest.e$pValues[init.rfAreas.MWTest.e$acceptAltHyp,]

final.rfAreas.MWTest.e = MWTestForInvCorticalMagnification ( mapStimCount, final.cortical.Amp.e, final.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e, pValue );
final.rfAreas.MWTest.e$pValues[final.rfAreas.MWTest.e$acceptAltHyp,]

init.final.rfAreas.MWTest.e = MWTestForInvCorticalMagnification ( mapStimCount, init.cortical.Amp.e, init.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e, pValue );
init.final.rfAreas.MWTest.e$pValues[init.final.rfAreas.MWTest.e$acceptAltHyp,]

	############################################################################################
	#
	#	END
	#	(INVERSE) CORTICAL MAGNIFICATION ANALYSIS
	#
	############################################################################################

	######################################################
	#
	#	CDF: Sector Pair-wise plot to compare Init, Final
	#		Receptive Field Extent
	#
	######################################################

x11(); par(mfrow=c(3,3));
mainText1 = "CDF of RF Extent Sector ";
mainText2 = paste("\n",initText," vs ", finalText, sep="");
mainText3 = "\nC Layer E Cells";
for ( iDigit in seq ( 3, 1, -1 ) ) {
	for ( iSector in seq ( 1, 3, 1 ) ) {
		iWhich = iDigit + (iSector - 1)*3;
		tmpLabel = MapSectorNumberToLabel ( iSector );
		D3PairSectorsRFDataCDF ( init.cortical.Amp.e, init.rfMap.rfAreas.e, iWhich, final.cortical.Amp.e, final.rfMap.rfAreas.e, iWhich,
							paste(mainText1, paste(iDigit, tmpLabel, sep=""), mainText2, mainText3, sep="" ),
							paste("Init.(Dash); Final(Solid)"), paste("Receptive Field Extent"), paste("CDF"), FALSE );
	} # for ( iSector in seq ( 1, 3, 1 ) ) {
} # for ( iDigit in seq ( 3, 1, -1 ) ) {

x11(); par(mfrow=c(3,3));
mainText1 = "CDF of RF Extent Sector ";
mainText2 = paste("\n",initText," vs ", finalText, sep="");
mainText3 = "\nC Layer I Cells";
for ( iDigit in seq ( 3, 1, -1 ) ) {
	for ( iSector in seq ( 1, 3, 1 ) ) {
		iWhich = iDigit + (iSector - 1)*3;
		tmpLabel = MapSectorNumberToLabel ( iSector );
		D3PairSectorsRFDataCDF ( init.cortical.Amp.i, init.rfMap.rfAreas.i, iWhich, final.cortical.Amp.i, final.rfMap.rfAreas.i, iWhich,
							paste(mainText1, paste(iDigit, tmpLabel, sep=""), mainText2, mainText3, sep="" ),
							paste("Init.(Dash); Final(Solid)"), paste("Receptive Field Extent"), paste("CDF"), FALSE );
	} # for ( iSector in seq ( 1, 3, 1 ) ) {
} # for ( iDigit in seq ( 3, 1, -1 ) ) {

	#######################################################
	#
	#	MWW: Sector pair-wise test of statistical diffs.
	#
	#######################################################

pLevel = 1e-5;
mww.tmp = D3AllSectorsMannWhitneyTest ( init.cortical.Amp.e, init.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e, pLevel );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]

mww.tmp = D3AllSectorsMannWhitneyTest ( final.cortical.Amp.e, final.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e, pLevel );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]

mww.tmp = D3AllSectorsMannWhitneyTest ( init.cortical.Amp.e, init.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e, pLevel );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]

	#####################################################################
	#
	#	CDF: Longitudonal Tracks Pair-wise plot to compare Init, Final
	#		Receptive Field Extent.
	#
	#####################################################################

x11(); par(mfcol=c(3,2));
mainText1 = "CDF of RF Extent Long. Axis ";
mainText2 = paste("\n",initText," vs ", finalText, sep="");
mainText3 = "\nC Layer E Cells";
subText = "Init.(Dash); Final(Solid); D1(B); D2(R); D3(G)";
for ( iLongAxis in seq ( as.integer(N/3), 5, -1 ) ) {
	D3AllLongAxisRFDataCDF ( iLongAxis, init.cortical.Amp.e, init.rfMap.rfAreas.e,
						final.cortical.Amp.e, final.rfMap.rfAreas.e,
						paste(mainText1, iLongAxis, mainText2, mainText3, sep="" ),
						subText, paste("Receptive Field Extent"), paste("CDF"), FALSE );
	#abline(v=25,lty=4,col=4);
	#abline(v=30,lty=4,col=5);
} # for ( iDigit in seq ( 3, 1, -1 ) ) {

x11(); par(mfcol=c(3,2));
mainText1 = "CDF of RF Extent Long. Axis ";
mainText2 = paste("\n",initText," vs ", finalText, sep="");
mainText3 = "\nC Layer I Cells";
subText = "Init.(Dash); Final(Solid); D1(B); D2(R); D3(G)";
subText = "Init.(Dash( iLongAxis in seq ( as.integer(N/3), 5, -1 ) )";
for ( iLongAxis in seq ( as.integer(N/3), 5, -1 ) ) {

	D3AllLongAxisRFDataCDF ( iLongAxis, init.cortical.Amp.i, init.rfMap.rfAreas.i,
						final.cortical.Amp.i, final.rfMap.rfAreas.i,
						paste(mainText1, iLongAxis, mainText2, mainText3, sep="" ),
						subText, paste("Receptive Field Extent"), paste("CDF"), FALSE );
	#abline(v=25,lty=4,col=4);
	#abline(v=30,lty=4,col=5);
} # for ( iDigit in seq ( 3, 1, -1 ) ) {

	#######################################################
	#
	#	START: Longitudonal Axis Cortical Magnification Plots
	#
	#######################################################

iTrim = 2;

init.sTable.e = GenD3LongAxisSummaryTable ( init.cortical.Amp.e, init.rfMap.rfAreas.e, mapStimCount, iTrim );
final.sTable.e = GenD3LongAxisSummaryTable ( final.cortical.Amp.e, final.rfMap.rfAreas.e, mapStimCount, iTrim );
final2.sTable.e = GenD3LongAxisSummaryTable ( final2.cortical.Amp.e, final2.rfMap.rfAreas.e, mapStimCount.syndact, iTrim );
final3.sTable.e = GenD3LongAxisSummaryTable ( final3.cortical.Amp.e, final3.rfMap.rfAreas.e, mapStimCount, iTrim );
all.sTable.e = round( cbind ( init.sTable.e, final.sTable.e, final2.sTable.e, final3.sTable.e ), 1 );
all.sTable.exp.e = round( cbind ( final.sTable.e, final2.sTable.e, final3.sTable.e ), 1 );
write.table ( all.sTable.exp.e, file="all.sTable.exp.syndactyly.e.txt", col.names=FALSE, row.names=FALSE, sep="\t" );

init.sTable.i = GenD3LongAxisSummaryTable ( init.cortical.Amp.i, init.rfMap.rfAreas.i, mapStimCount, iTrim );
final.sTable.i = GenD3LongAxisSummaryTable ( final.cortical.Amp.i, final.rfMap.rfAreas.i, mapStimCount, iTrim );
final2.sTable.i = GenD3LongAxisSummaryTable ( final2.cortical.Amp.i, final2.rfMap.rfAreas.i, mapStimCount.syndact, iTrim );
final3.sTable.i = GenD3LongAxisSummaryTable ( final3.cortical.Amp.i, final3.rfMap.rfAreas.i, mapStimCount, iTrim );
all.sTable.i = round( cbind ( init.sTable.i, final.sTable.i, final2.sTable.i, final3.sTable.i ), 1 );
all.sTable.exp.i = round( cbind ( final.sTable.i, final2.sTable.i, final3.sTable.i ), 1 );
write.table ( all.sTable.exp.i, file="all.sTable.exp.syndactyly.i.txt", col.names=FALSE, row.names=FALSE, sep="\t" );

final2Ctl.sTable.e = GenD3LongAxisSummaryTable ( final2Ctl.cortical.Amp.e, final2Ctl.rfMap.rfAreas.e, mapStimCount.syndact, iTrim );
final3Ctl.sTable.e = GenD3LongAxisSummaryTable ( final3Ctl.cortical.Amp.e, final3Ctl.rfMap.rfAreas.e, mapStimCount, iTrim );
allCtl.sTable.e = round( cbind ( init.sTable.e, final.sTable.e, final2Ctl.sTable.e, final3Ctl.sTable.e ), 1 );
allCtl.sTable.exp.e = round( cbind ( final.sTable.e, final2Ctl.sTable.e, final3Ctl.sTable.e ), 1 );
write.table ( allCtl.sTable.exp.e, file="allCtl.sTable.exp.syndactyly.e.txt", col.names=FALSE, row.names=FALSE, sep="\t" );

final2Ctl.sTable.i = GenD3LongAxisSummaryTable ( final2Ctl.cortical.Amp.i, final2Ctl.rfMap.rfAreas.i, mapStimCount.syndact, iTrim );
final3Ctl.sTable.i = GenD3LongAxisSummaryTable ( final3Ctl.cortical.Amp.i, final3Ctl.rfMap.rfAreas.i, mapStimCount, iTrim );
allCtl.sTable.i = round( cbind ( init.sTable.i, final.sTable.i, final2Ctl.sTable.i, final3Ctl.sTable.i ), 1 );
allCtl.sTable.exp.i = round( cbind ( final.sTable.i, final2Ctl.sTable.i, final3Ctl.sTable.i ), 1 );
write.table ( allCtl.sTable.exp.i, file="allCtl.sTable.exp.syndactyly..txt", col.names=FALSE, row.names=FALSE, sep="\t" );




for ( iDigit in seq ( 1 ) ) {

	if ( iDigit == 1 ) { iAxisList = seq ( 1, 10 ); }
	if ( iDigit == 2 ) { iAxisList = c ( 1, 2, 5  ); }

	for ( iAxis in iAxisList ) {

		x11(); par ( mfcol = c ( 3,2 ) );
		xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

		iWhich = 25; titleText = paste ( "E Cortical Magnification\nBaseline Refinement\nIter ", iWhich, ", Digit ", iDigit, ", Axis ", iAxis, sep="" );
		ShowVecAsMap1 ( GetLongAxisCorticalMag ( final.cortical.Amp.e, N, iDigit, iAxis, iTrim ), titleText, xlab, ylab );

		iWhich = 25; titleText = paste ( "E Cortical Magnification\nSyndactyly Refinement\nIter ", iWhich, ", Digit ", iDigit, ", Axis ", iAxis, sep="" );
		ShowVecAsMap1 ( GetLongAxisCorticalMag ( final2.cortical.Amp.e, N, iDigit, iAxis, iTrim ), titleText, xlab, ylab );

		iWhich = 25; titleText = paste ( "E Cortical Magnification\nSyndactyly Release Refinement\nIter ", iWhich, ", Digit ", iDigit, ", Axis ", iAxis, sep="" );
		ShowVecAsMap1 ( GetLongAxisCorticalMag ( final3.cortical.Amp.e, N, iDigit, iAxis, iTrim ), titleText, xlab, ylab );

		iWhich = 25; titleText = paste ( "I Cortical Magnification\nBaseline Refinement\nIter ", iWhich, ", Digit ", iDigit, ", Axis ", iAxis, sep="" );
		ShowVecAsMap1 ( GetLongAxisCorticalMag ( final.cortical.Amp.i, N, iDigit, iAxis, iTrim ), titleText, xlab, ylab );

		iWhich = 25; titleText = paste ( "I Cortical Magnification\nSyndactyly Refinement\nIter ", iWhich, ", Digit ", iDigit, ", Axis ", iAxis, sep="" );
		ShowVecAsMap1 ( GetLongAxisCorticalMag ( final2.cortical.Amp.i, N, iDigit, iAxis, iTrim ), titleText, xlab, ylab );

		iWhich = 25; titleText = paste ( "I Cortical Magnification\nSyndactyly Release Refinement\nIter ", iWhich, ", Digit ", iDigit, ", Axis ", iAxis, sep="" );
		ShowVecAsMap1 ( GetLongAxisCorticalMag ( final3.cortical.Amp.i, N, iDigit, iAxis, iTrim ), titleText, xlab, ylab );


	} # for ( iLongAxis in c ( 9, 10, 11, 12 ) ) {
} # for ( iDigit in seq ( 1, 2 ) )); Final(Solid); D1(B); D2(R); D3(G)";

	#######################################################
	#
	#	END: Longitudonal Axis Cortical Magnification Plots
	#
	#######################################################

	#####################################################################
	#
	#	MWW: Longitudonal Tracks Pair-wise plot to compare Init, Final
	#		Receptive Field Extent.
	#
	#####################################################################

pLevel = 0.0001;
mww.tmp = D3AllLongAxesMannWhitneyTest ( init.cortical.Amp.e, init.rfMap.rfAreas.e, init.cortical.Amp.e, init.rfMap.rfAreas.e, pLevel );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]

mww.tmp = D3AllLongAxesMannWhitneyTest ( final.cortical.Amp.e, final.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e, pLevel );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]

mww.tmp = D3AllLongAxesMannWhitneyTest ( init.cortical.Amp.e, init.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e, pLevel );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]



		#
		#	The following does CDF controlling for position only.  Call: D3AllLongAxesMannWhitneyTest2 .
		#
pLevel = 0.0001; trimEnds = 2;

		#	Experimental digit 1 & 2: E: baseline & syndactyly

mww.tmp = D3AllLongAxesMannWhitneyTest2 ( final.cortical.Amp.e, final.rfMap.rfAreas.e, final2.cortical.Amp.e, final2.rfMap.rfAreas.e, pLevel, trimEnds );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]
sameCol = ( mww.tmp$pValues[,1] == mww.tmp$pValues[,2] );
base.syndact.Mags.e = mww.tmp$pValues[sameCol,c(1,2,4,5)]

		#	Experimental digit 1 & 2: I: baseline & syndactyly
mww.tmp = D3AllLongAxesMannWhitneyTest2 ( final.cortical.Amp.i, final.rfMap.rfAreas.i, final2.cortical.Amp.i, final2.rfMap.rfAreas.i, pLevel, trimEnds );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]
sameCol = mww.tmp$pValues[,1] == mww.tmp$pValues[,2];
base.syndact.Mags.i = mww.tmp$pValues[sameCol,c(1,2,4,5)]

		#	Experimental digit 1 & 2: E: baseline & syndactyly release
mww.tmp = D3AllLongAxesMannWhitneyTest2 ( final.cortical.Amp.e, final.rfMap.rfAreas.e, final3.cortical.Amp.e, final3.rfMap.rfAreas.e, pLevel, trimEnds );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]
sameCol = mww.tmp$pValues[,1] == mww.tmp$pValues[,2];
base.syndactRel.Mags.e = mww.tmp$pValues[sameCol,c(1,2,4,5)]

		#	Experimental digit 1 & 2: I: baseline & syndactyly release
mww.tmp = D3AllLongAxesMannWhitneyTest2 ( final.cortical.Amp.i, final.rfMap.rfAreas.i, final3.cortical.Amp.i, final3.rfMap.rfAreas.i, pLevel, trimEnds );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]
sameCol = mww.tmp$pValues[,1] == mww.tmp$pValues[,2];
base.syndactRel.Mags.i = mww.tmp$pValues[sameCol,c(1,2,4,5)]



		#	Experimental digit 1 & 2: E: baseline & syndactyly release ctl
mww.tmp = D3AllLongAxesMannWhitneyTest2 ( final.cortical.Amp.e, final.rfMap.rfAreas.e, final3Ctl.cortical.Amp.e, final3Ctl.rfMap.rfAreas.e, pLevel, trimEnds );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]
sameCol = mww.tmp$pValues[,1] == mww.tmp$pValues[,2];
base.syndactRel.Mags.e = mww.tmp$pValues[sameCol,c(1,2,4,5)]

		#	Experimental digit 1 & 2: I: baseline & syndactyly release ctl
mww.tmp = D3AllLongAxesMannWhitneyTest2 ( final.cortical.Amp.i, final.rfMap.rfAreas.i, final3Ctl.cortical.Amp.i, final3Ctl.rfMap.rfAreas.i, pLevel, trimEnds );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]
sameCol = mww.tmp$pValues[,1] == mww.tmp$pValues[,2];
base.syndactRel.Mags.i = mww.tmp$pValues[sameCol,c(1,2,4,5)]

mags.e = cbind ( seq(1,30), base.syndact.Mags.e[,3], base.syndact.Mags.e[,4], base.syndactRel.Mags.e[,4] );
mags.i = cbind ( seq(1,30), base.syndact.Mags.i[,3], base.syndact.Mags.i[,4], base.syndactRel.Mags.i[,4] );

		#
		#	The following does CDF controlling for position and stimulation count.  Call: D3AllLongAxesMannWhitneyTest1 
		#

pLevel = 0.0001;
trimEnds = 4;
final.cutOff = final2.cutOff = final3.cutOff = 25;


		#	Init vs Init: E Cells
	mww.tmp = NULL; mww.tmp = D3AllLongAxesMannWhitneyTest1 ( 	mapStimCount, final.cutOff, init.cortical.Amp.e, init.rfMap.rfAreas.e,
										mapStimCount, final.cutOff, init.cortical.Amp.e, init.rfMap.rfAreas.e, pLevel, trimEnds );
	mww.tmp$pValues[mww.tmp$acceptAltHyp,]

		#	Init vs Init: I Cells
	mww.tmp = NULL; mww.tmp = D3AllLongAxesMannWhitneyTest1 ( 	mapStimCount, final.cutOff, init.cortical.Amp.i, init.rfMap.rfAreas.i,
										mapStimCount, final.cutOff, init.cortical.Amp.i, init.rfMap.rfAreas.i, pLevel, trimEnds );
	mww.tmp$pValues[mww.tmp$acceptAltHyp,]

		#	Base vs Base: E Cells
	mww.tmp = NULL; mww.tmp = D3AllLongAxesMannWhitneyTest1 ( 	mapStimCount, final.cutOff, final.cortical.Amp.e, final.rfMap.rfAreas.e,
										mapStimCount, final.cutOff, final.cortical.Amp.e, final.rfMap.rfAreas.e, pLevel, trimEnds );
	mww.tmp$pValues[mww.tmp$acceptAltHyp,]

		#	Base vs Base: I Cells
	mww.tmp = NULL; mww.tmp = D3AllLongAxesMannWhitneyTest1 ( 	mapStimCount, final.cutOff, final.cortical.Amp.i, final.rfMap.rfAreas.i,
										mapStimCount, final.cutOff, final.cortical.Amp.i, final.rfMap.rfAreas.i, pLevel, trimEnds );
	mww.tmp$pValues[mww.tmp$acceptAltHyp,]


		#	Control Digit: Base vs Syndact: E Cells
	mww.tmp = NULL; mww.tmp = D3AllLongAxesMannWhitneyTest1 ( 	mapStimCount, final.cutOff, final.cortical.Amp.e, final.rfMap.rfAreas.e,
										mapStimCount.syndact, final2.cutOff, final2.cortical.Amp.e, final2.rfMap.rfAreas.e, pLevel, trimEnds );
	mww.tmp$pValues[mww.tmp$acceptAltHyp,]

		#	Control Digit: Base vs Syndact: I Cells
	mww.tmp = NULL; mww.tmp = D3AllLongAxesMannWhitneyTest1 ( 	mapStimCount, final.cutOff, final.cortical.Amp.i, final.rfMap.rfAreas.i,
										mapStimCount.syndact, final2.cutOff, final2.cortical.Amp.i, final2.rfMap.rfAreas.i, pLevel, trimEnds );
	mww.tmp$pValues[mww.tmp$acceptAltHyp,]


		#	Control Digit: Base vs Syndact Rel: E Cells
	mww.tmp = NULL; mww.tmp = D3AllLongAxesMannWhitneyTest1 ( 	mapStimCount, final.cutOff, final.cortical.Amp.e, final.rfMap.rfAreas.e,
										mapStimCount, final3.cutOff, final3.cortical.Amp.e, final3.rfMap.rfAreas.e, pLevel, trimEnds );
	mww.tmp$pValues[mww.tmp$acceptAltHyp,]

		#	Control Digit: Base vs Syndact Rel: I Cells
	mww.tmp = NULL; mww.tmp = D3AllLongAxesMannWhitneyTest1 ( 	mapStimCount, final.cutOff, final.cortical.Amp.i, final.rfMap.rfAreas.i,
										mapStimCount, final3.cutOff, final3.cortical.Amp.i, final3.rfMap.rfAreas.i, pLevel, trimEnds );
	mww.tmp$pValues[mww.tmp$acceptAltHyp,]



		#	Control Digit: Base vs Syndact Rel Ctl: E Cells
	mww.tmp = NULL; mww.tmp = D3AllLongAxesMannWhitneyTest1 ( 	mapStimCount, final.cutOff, final.cortical.Amp.e, final.rfMap.rfAreas.e,
										mapStimCount, final3.cutOff, final3Ctl.cortical.Amp.e, final3Ctl.rfMap.rfAreas.e, pLevel, trimEnds );
	mww.tmp$pValues[mww.tmp$acceptAltHyp,]

		#	Control Digit: Base vs Syndact Rel Ctl: I Cells
	mww.tmp = NULL; mww.tmp = D3AllLongAxesMannWhitneyTest1 ( 	mapStimCount, final.cutOff, final.cortical.Amp.i, final.rfMap.rfAreas.i,
										mapStimCount, final3.cutOff, final3Ctl.cortical.Amp.i, final3Ctl.rfMap.rfAreas.i, pLevel, trimEnds );
	mww.tmp$pValues[mww.tmp$acceptAltHyp,]




