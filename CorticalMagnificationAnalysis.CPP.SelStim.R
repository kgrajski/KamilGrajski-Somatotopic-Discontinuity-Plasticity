
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

if ( !alreadyCleared ) { rm(list = ls()); }

		#
		#	Some additional variables needed if/when NMExpZone isn't the source.
kRFPeakToEdgeDetect = 0.5;
trialDurRFProbeInIters = 150;
oneSecondNumIter = 100;
		#
		#

	###############################################################
	#
	#	Get the baseline refined network.  Call it baseline I.
	#
	###############################################################

fName = paste(fDir, "\\", "Base.RFMap.25.25.bin", sep="" );
initText = "Baseline";

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


	###############################################################
	#
	#	Get the experimental network.
	#
	###############################################################

fName = paste(fDir, "\\", "SelStimExp.RFMap.25.25.bin", sep="" );
fName = paste(fDir, "\\", "SelStimExp.RFMap.25.25.",expSectorID,".",expFactor,".bin", sep="" );
finalText = "FocalStim";

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


	###############################################################
	#
	#	Get the control network.  Call it Baseline II
	#
	###############################################################

fName = paste(fDir, "\\", "SelStimCtl.RFMap.25.25.bin", sep="" );
finalText2 = "FocalStim.Ctl";

	finfo = file.info ( fName );
	toread = file ( fName, "rb" );
	alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
	close ( toread );	
	itmp = length ( alldata ) / 2;
	N2 = as.integer ( sqrt ( itmp ) );
	final.ctl.r1.e.rfMap = r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
	final.ctl.r1.i.rfMap = r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

final.ctl.cortical.Amp.e = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
final.ctl.cortical.Amp.i = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );
final.ctl.rfMap.rfAreas.e = QuantRFSize ( r1.e.rfMap, kRFPeakToEdgeDetect );
final.ctl.rfMap.rfAreas.i = QuantRFSize ( r1.i.rfMap, kRFPeakToEdgeDetect );
final.ctl.rfMap.maxResp.e = apply ( r1.e.rfMap, 1, max );
final.ctl.rfMap.maxResp.i = apply ( r1.i.rfMap, 1, max );

	
	###############################################################
	#
	#	Data Visualizations.
	#
	###############################################################

		#
		#	The following scatterplot compares the last stage of baseline I
		#	vs last stage of focal stimulation for cortical E cells.
		#

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

		#
		#	The following scatterplot compares the last stage of baseline II
		#	vs last stage of focal stimulation for cortical E cells.
		#

x11();
xlim.e = c ( min(final.ctl.rfMap.maxResp.e,final.rfMap.maxResp.e), max(final.ctl.rfMap.maxResp.e,final.rfMap.maxResp.e) );
ylim.e = c ( min(final.ctl.rfMap.rfAreas.e,final.rfMap.rfAreas.e), max(final.ctl.rfMap.rfAreas.e,final.rfMap.rfAreas.e) );
plot ( final.ctl.rfMap.maxResp.e, final.ctl.rfMap.rfAreas.e, xlim=xlim.e, ylim=ylim.e, pch=".", col=1, cex=4,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer E-Cells\n",initText,"(.); ",finalText,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.e, final.rfMap.rfAreas.e, pch="*", col=2, cex=2 );
delta.area.e = delta.maxResp.e = slope.e = rep ( 0, N2 );
for ( i in 1:N2 ) {
	#lines ( c( final.ctl.rfMap.maxResp.e[i], final.rfMap.maxResp.e[i] ), c( final.ctl.rfMap.rfAreas.e[i], final.rfMap.rfAreas.e[i] ) )
	delta.area.e[i] = ( final.rfMap.rfAreas.e[i] - final.ctl.rfMap.rfAreas.e[i] );
	delta.maxResp.e[i] = ( final.rfMap.maxResp.e[i] - final.ctl.rfMap.maxResp.e[i] );
	slope.e[i] = ( final.rfMap.rfAreas.e[i] - final.ctl.rfMap.rfAreas.e[i] ) / ( final.rfMap.maxResp.e[i] - final.ctl.rfMap.maxResp.e[i] );
} # for ( i in 1:N2 ) {

		#
		#	The following figure displays the spatial distribution of the four
		#	quadrants of the scatterplot of focal stimulation cortical E cell results.
		#

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

		#
		#	The following scatterplot compares the last stage of baseline I
		#	vs last stage of focal stimulation for cortical I cells.
		#

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

		#
		#	The following scatterplot compares the last stage of baseline II
		#	vs last stage of focal stimulation for cortical I cells.
		#

x11();
xlim.i = c ( min(final.ctl.rfMap.maxResp.i,final.rfMap.maxResp.i), max(final.ctl.rfMap.maxResp.i,final.rfMap.maxResp.i) );
ylim.i = c ( min(final.ctl.rfMap.rfAreas.i,final.rfMap.rfAreas.i), max(final.ctl.rfMap.rfAreas.i,final.rfMap.rfAreas.i) );
plot ( final.ctl.rfMap.maxResp.i, final.ctl.rfMap.rfAreas.i, xlim=xlim.i, ylim=ylim.i, pch=".", col=3, cex=3,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer I-Cells\n",initText,"(.); ",finalText,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.i, final.rfMap.rfAreas.i, pch="*", col=4 );
delta.area.i = delta.maxResp.i = slope.i = rep ( 0, N2 );
for ( i in 1:N2 ) {
	#lines ( c( final.ctl.rfMap.maxResp.i[i], final.rfMap.maxResp.i[i] ), c( final.ctl.rfMap.rfAreas.i[i], final.rfMap.rfAreas.i[i] ) )
	delta.area.i[i] = ( final.rfMap.rfAreas.i[i] - final.ctl.rfMap.rfAreas.i[i] );
	delta.maxResp.i[i] = ( final.rfMap.maxResp.i[i] - final.ctl.rfMap.maxResp.i[i] );
	slope.i[i] = ( final.rfMap.rfAreas.i[i] - final.ctl.rfMap.rfAreas.i[i] ) / ( final.rfMap.maxResp.i[i] - final.ctl.rfMap.maxResp.i[i] );
} # for ( i in 1:N2 ) {

		#
		#	The following figure displays the spatial distribution of the four
		#	quadrants of the scatterplot of focal stimulation cortical I cell results.
		#

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


	###############################################################
	#
	#	(INVERSE) CORTICAL MAGNIFICATION ANALYSIS
	#
	###############################################################

		#
		#	Data visualization: by input layer node analysis.
		#	Display scatterplots and spatial maps of cortical
		#	magnification (CM) and receptive field sizes.
		#

iTrimRings = 1;
delta.CM.min = 0.1;
delta.MeanRFArea.min = 0.1;
mainTextIn = "Focal Stim. vs Focal Stim. Control";
expMain = CMPointByPoint ( final.ctl.cortical.Amp.e, final.cortical.Amp.e, final.ctl.cortical.Amp.i, final.cortical.Amp.i,
			final.ctl.rfMap.rfAreas.e, final.rfMap.rfAreas.e, final.ctl.rfMap.rfAreas.i, final.rfMap.rfAreas.i,
			iTrimRings, delta.CM.min, delta.MeanRFArea.min, mainTextIn );

delta.CM.min = 0.1;
delta.MeanRFArea.min = 0.1;
mainTextIn = "Focal Stim. Control vs Baseline Final";
ctlMain = CMPointByPoint ( init.cortical.Amp.e, final.ctl.cortical.Amp.e, init.cortical.Amp.i, final.ctl.cortical.Amp.i,
			init.rfMap.rfAreas.e, final.ctl.rfMap.rfAreas.e, init.rfMap.rfAreas.i, final.ctl.rfMap.rfAreas.i,
			iTrimRings, delta.CM.min, delta.MeanRFArea.min, mainTextIn );

delta.CM.min = 0.1;
delta.MeanRFArea.min = 0.1;
mainTextIn = "Focal Stim. vs Baseline Final";
ctlSecondary = CMPointByPoint ( init.cortical.Amp.e, final.cortical.Amp.e, init.cortical.Amp.i, final.cortical.Amp.i,
			init.rfMap.rfAreas.e, final.rfMap.rfAreas.e, init.rfMap.rfAreas.i, final.rfMap.rfAreas.i,
			iTrimRings, delta.CM.min, delta.MeanRFArea.min, mainTextIn );










	
		#
		#	Extract relative frequencies of stimulation.
		#

fRoot = "CheckInputStim.Baseline";
fileRootName = paste(fDir, fRoot, sep="\\");
mapStimCount = base.stimCount = init.stimCount = GetStimCountData ( fileRootName, N2 );

fRoot = "CheckInputStim.SelStim";
fileRootName = paste(fDir, fRoot, sep="\\");
#final.stimCount = GetStimCountData ( fileRootName, N2 );
final.stimCount = selStim.stimCount;



		#
		#	Set some partition values.
		#
xIDSubset = c ( 9, 18, 27 );
#xIDSubset = c ( 5, 10, 25 );
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


	######################################################
	#
	#	CDF: Sector Pair-wise plot to compare Init, Final
	#		Receptive Field Extent
	#
	#	This is a special version to focus just on the mid phalange and
	#	compare final vs final.ctl.
	#
	######################################################
x11(); par(mfcol=c(3,2));
mainText1 = "CDF of RF Extent Sector ";
mainText2 = paste("\n",initText," vs ", finalText, sep="");
mainText3 = "\nC Layer E Cells";
for ( iDigit in seq ( 3, 1, -1 ) ) {
	for ( iSector in seq ( 2, 2, 1 ) ) {
		iWhich = iDigit + (iSector - 1)*3;
		tmpLabel = MapSectorNumberToLabel ( iSector );
		D3PairSectorsRFDataCDF ( final.ctl.cortical.Amp.e, final.ctl.rfMap.rfAreas.e, iWhich, final.cortical.Amp.e, final.rfMap.rfAreas.e, iWhich,
							paste(mainText1, paste(iDigit, tmpLabel, sep=""), mainText2, mainText3, sep="" ),
							paste("Init.(Dash); Final(Solid)"), paste("Receptive Field Extent"), paste("CDF"), FALSE );
	} # for ( iSector in seq ( 1, 3, 1 ) ) {
} # for ( iDigit in seq ( 3, 1, -1 ) ) {

mainText1 = "CDF of RF Extent Sector ";
mainText2 = paste("\n",initText," vs ", finalText, sep="");
mainText3 = "\nC Layer I Cells";
for ( iDigit in seq ( 3, 1, -1 ) ) {
	for ( iSector in seq ( 2, 2, 1 ) ) {
		iWhich = iDigit + (iSector - 1)*3;
		tmpLabel = MapSectorNumberToLabel ( iSector );
		D3PairSectorsRFDataCDF ( final.ctl.cortical.Amp.i, final.ctl.rfMap.rfAreas.i, iWhich, final.cortical.Amp.i, final.rfMap.rfAreas.i, iWhich,
							paste(mainText1, paste(iDigit, tmpLabel, sep=""), mainText2, mainText3, sep="" ),
							paste("Init.(Dash); Final(Solid)"), paste("Receptive Field Extent"), paste("CDF"), FALSE );
	} # for ( iSector in seq ( 1, 3, 1 ) ) {
} # for ( iDigit in seq ( 3, 1, -1 ) ) {



	#######################################################
	#
	#	MWW: Sector pair-wise test of statistical diffs.
	#
	#######################################################
pLevel = 1e-1;
mww.tmp = D3AllSectorsMannWhitneyTest ( init.cortical.Amp.e, init.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e, pLevel );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]

mww.tmp = D3AllSectorsMannWhitneyTest ( final.cortical.Amp.e, final.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e, pLevel );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]

mww.tmp = D3AllSectorsMannWhitneyTest ( final.ctl.cortical.Amp.e, final.ctl.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e, pLevel );
mww.tmp$pValues[mww.tmp$acceptAltHyp,]

	#####################################################################
	#
	#	CDF: Longitudonal Tracks Pair-wise plot to compare Init, Final
	#		Receptive Field Extent.
	#
	#####################################################################
x11(); par(mfcol=c(3,3));
mainText1 = "CDF of RF Extent Long. Axis ";
mainText2 = paste("\n",initText," vs ", finalText, sep="");
mainText3 = "\nC Layer E Cells";
subText = "Init.(Dash); Final(Solid); D1(B); D2(R); D3(G)";
for ( iLongAxis in seq ( as.integer(N/3), 2, -1 ) ) {
	D3AllLongAxisRFDataCDF ( iLongAxis, init.cortical.Amp.e, init.rfMap.rfAreas.e,
						final.cortical.Amp.e, final.rfMap.rfAreas.e,
						paste(mainText1, iLongAxis, mainText2, mainText3, sep="" ),
						subText, paste("Receptive Field Extent"), paste("CDF"), FALSE );
	#abline(v=25,lty=4,col=4);
	#abline(v=30,lty=4,col=5);
} # for ( iDigit in seq ( 3, 1, -1 ) ) {

x11(); par(mfcol=c(3,3));
mainText1 = "CDF of RF Extent Long. Axis ";
mainText2 = paste("\n",initText," vs ", finalText, sep="");
mainText3 = "\nC Layer I Cells";
subText = "Init.(Dash); Final(Solid); D1(B); D2(R); D3(G)";
for ( iLongAxis in seq ( as.integer(N/3), 2, -1 ) ) {
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

	#####################################################################
	#####################################################################
	#
	#	RF Overlap Analysis.	
	#
	#####################################################################
	#####################################################################

init.ovlap.long.e = RFOverlap ( init.cortical.Amp.e, TRUE );
init.ovlap.cross.e = RFOverlap ( init.cortical.Amp.e, FALSE );
final.ovlap.long.e = RFOverlap ( final.cortical.Amp.e, TRUE );
final.ovlap.cross.e = RFOverlap ( final.cortical.Amp.e, FALSE );

init.ovlap.long.i = RFOverlap ( init.cortical.Amp.i, TRUE );
init.ovlap.cross.i = RFOverlap ( init.cortical.Amp.i, FALSE );
final.ovlap.long.i = RFOverlap ( final.cortical.Amp.i, TRUE );
final.ovlap.cross.i = RFOverlap ( final.cortical.Amp.i, FALSE );

x11(); par(mfrow=c(2,2));
mainTxt="Pairwise RF Extent % Overlap"; xlabTxt = "ith Pairwise Comparison"; ylabTxt = "Proportion Overlap";
PlotRFOverlap ( init.ovlap.long.e, paste(mainTxt,"\nBaseline Refined\nE Cells; Long. Axis",sep=""), xlabTxt, ylabTxt );
PlotRFOverlap ( init.ovlap.long.i, paste(mainTxt,"\nBaseline Refined\nI Cells; Long. Axis",sep=""), xlabTxt, ylabTxt);
PlotRFOverlap ( final.ovlap.long.e, paste(mainTxt,"\nFocal Stimulation\nE Cells; Long. Axis",sep=""), xlabTxt, ylabTxt );
PlotRFOverlap ( final.ovlap.long.i, paste(mainTxt,"\nFocal Stimulation\nI Cells; Long. Axis",sep=""), xlabTxt, ylabTxt );

x11(); par(mfrow=c(2,2));
mainTxt="Pairwise RF Extent % Overlap"; xlabTxt = "ith Pairwise Comparison"; ylabTxt = "Proportion Overlap";
PlotRFOverlap ( init.ovlap.cross.e, paste(mainTxt,"\nBaseline Refined\nE Cells; Trans. Axis",sep=""), xlabTxt, ylabTxt );
PlotRFOverlap ( init.ovlap.cross.i, paste(mainTxt,"\nBaseline Refined\nI Cells; Trans. Axis",sep=""), xlabTxt, ylabTxt);
PlotRFOverlap ( final.ovlap.cross.e, paste(mainTxt,"\nFocal Stimulation\nE Cells; Trans. Axis",sep=""), xlabTxt, ylabTxt );
PlotRFOverlap ( final.ovlap.cross.i, paste(mainTxt,"\nFocal Stimulation\nI Cells; Trans. Axis",sep=""), xlabTxt, ylabTxt );



