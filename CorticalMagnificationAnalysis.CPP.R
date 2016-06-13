
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

rm(list = ls());

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
N = 30;
N2 = N * N;
kRFPeakToEdgeDetect = 0.5;
trialDurRFProbeInIters = 150;
oneSecondNumIter = 100;
refinementPatchSize = 3;
		#
		#

E = 65;
P = 4;
fDir = paste("E:\\NMLab\\E",E,".",N,".",P,sep="");

	######################################################
	#
	#	Get the network with random initial conditions.
	#
	######################################################
#fName = paste(fDir, "\\", "Base.RFMap.25.0.bin", sep="" );
fName = paste(fDir, "\\", "Base.RFMap.25.25.bin", sep="" );

initText = "Refined";

	finfo = file.info ( fName );
	toread = file ( fName, "rb" );
	alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
	close ( toread );	
	itmp = length ( alldata ) / 2;
	N2 = as.integer ( sqrt ( itmp ) );
	r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
	r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

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
fName = paste(fDir, "\\", "SyndactExp.RFMap.25.25.bin", sep="" );
#fName = paste(fDir, "\\", "Base.RFMap.25.25.bin", sep="" );
finalText = "Syndactyly";

	finfo = file.info ( fName );
	toread = file ( fName, "rb" );
	alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
	close ( toread );	
	itmp = length ( alldata ) / 2;
	N2 = as.integer ( sqrt ( itmp ) );
	r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
	r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

final.cortical.Amp.e = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
final.cortical.Amp.i = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );
final.rfMap.rfAreas.e = QuantRFSize ( r1.e.rfMap, kRFPeakToEdgeDetect );
final.rfMap.rfAreas.i = QuantRFSize ( r1.i.rfMap, kRFPeakToEdgeDetect );
final.rfMap.maxResp.e = apply ( r1.e.rfMap, 1, max );
final.rfMap.maxResp.i = apply ( r1.i.rfMap, 1, max );
	
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
plot ( init.rfMap.maxResp.i, init.rfMap.rfAreas.i, xlim=xlim.i, ylim=ylim.i, pch=".", col=3, cex=3,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer I-Cells\n",initText,"(.); ",finalText,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.i, final.rfMap.rfAreas.i, pch="*", col=4 );
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
#xStim = Dig3MapRefine0( N, trialDurRFProbeInIters, oneSecondNumIter, refinementPatchSize, 0.25 );
#mapStimCount = xStim[[length(xStim)]]$stimCount;

fRoot = "CheckInputStim.Baseline";
fileRootName = paste(fDir, fRoot, sep="\\");
mapStimCount = base.stimCount = GetStimCountData ( fileRootName, N2 );

		#
		#	Draw some pictures to visualize the data.
		#
xIDSubset = c ( 4, 8, 16 );
xIDSubset = c ( 5, 10, 25 );
#xIDSubset = c ( 1, 3, 5 );

		#
		#	CDF: CORTICAL MAGNIFICATION
		#
init.CMag.e = apply ( init.cortical.Amp.e, 2, sum );
final.CMag.e = apply ( final.cortical.Amp.e, 2, sum );

x11(); par(mfrow=c(2,2));
xLabelText = "Cortical Magnification";
yLabelText = "CDF";
zTypeText = "E Cell"
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
subText = "Init.(Dash); Final(Solid); 5x(B); 10x(R); 20x(G)";
#subText = "Init.(Dash); Final(Solid); 2x(B); 4x(R); 9x(G)";
CorticalMagCDF ( mapStimCount, init.CMag.e, final.CMag.e, titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

		#
		#	CDF: RECEPTIVE FIELD EXTENT
		#
xLabelText = "RF Extent";
yLabelText = "CDF";
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
subText = "Init.(Dash); Final(Solid); 5x(B); 10x(R); 20x(G)";
#subText = "Init.(Dash); Final(Solid); 2x(B); 4x(R); 9x(G)";
InvCorticalMagCDF ( mapStimCount, init.cortical.Amp.e, init.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e,
				titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

		#
		#	CDF: CORTICAL MAGNIFICATION
		#
init.CMag.i = apply ( init.cortical.Amp.i, 2, sum );
final.CMag.i = apply ( final.cortical.Amp.i, 2, sum );


xLabelText = "Cortical Magnification";
yLabelText = "CDF";
zTypeText = "I Cell"
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
subText = "Init.(Dash); Final(Solid); 5x(B); 10x(R); 20x(G)";
#subText = "Init.(Dash); Final(Solid); 2x(B); 4x(R); 9x(G)";
CorticalMagCDF ( mapStimCount, init.CMag.i, final.CMag.i, titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

		#
		#	CDF: RECEPTIVE FIELD EXTENT
		#
xLabelText = "RF Extent";
yLabelText = "CDF";
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
subText = "Init.(Dash); Final(Solid); 5x(B); 10x(R); 20x(G)";
#subText = "Init.(Dash); Final(Solid); 2x(B); 4x(R); 9x(G)";
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
for ( iLongAxis in seq ( as.integer(N/3), 5, -1 ) ) {
	D3AllLongAxisRFDataCDF ( iLongAxis, init.cortical.Amp.i, init.rfMap.rfAreas.i,
						final.cortical.Amp.i, final.rfMap.rfAreas.i,
						paste(mainText1, iLongAxis, mainText2, mainText3, sep="" ),
						subText, paste("Receptive Field Extent"), paste("CDF"), FALSE );
	#abline(v=25,lty=4,col=4);
	#abline(v=30,lty=4,col=5);
} # for ( iDigit in seq ( 3, 1, -1 ) ) {

	#####################################################################
	#
	#	MWW: Longitudonal Tracks Pair-wise plot to compare Init, Final
	#		Receptive Field Extent.
	#
	#####################################################################
pLevel = 1e-5;
mww.tmp = D3AllLongAxesMannWhitneyTest ( init.cortical.Amp.e, init.rfMap.rfAreas.e, init.cortical.Amp.e, init.rfMap.rfAreas.e, pLevel );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]

mww.tmp = D3AllLongAxesMannWhitneyTest ( final.cortical.Amp.e, final.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e, pLevel );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]

mww.tmp = D3AllLongAxesMannWhitneyTest ( init.cortical.Amp.e, init.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e, pLevel );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]

