
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
source ( "LabelledCDF.UseDep.R" );
source ( "MapMagToCellList.R" );

	######################################################
	#
	#	Get the network with random initial conditions.
	#
	######################################################
load(file="Run.25.25");
#(file="Run.100.0");


		# Earlier versions of RFMapExp1 didn't do this.
rfMap.maxResp.e = apply ( rfMap$r1.e.rfMap, 1, max );
rfMap.maxResp.i = apply ( rfMap$r1.i.rfMap, 1, max );

init.cortical.Amp.e = cortical.Amp.e;
init.cortical.Amp.i = cortical.Amp.i;
init.rfMap.maxResp.e = rfMap.maxResp.e;
init.rfMap.rfAreas.e = rfMap.rfAreas.e;
init.rfMap.maxResp.i = rfMap.maxResp.i;
init.rfMap.rfAreas.i = rfMap.rfAreas.i;

	######################################################
	#
	#	Get the refined network.
	#
	######################################################
load(file="SelStim.7.10.10");
#load(file="Run.100.50");
#load(file="Run.100.100");

		# 	Earlier versions of RFMapExp1 didn't do this.
rfMap.maxResp.e = apply ( rfMap$r1.e.rfMap, 1, max );
rfMap.maxResp.i = apply ( rfMap$r1.i.rfMap, 1, max );

final.cortical.Amp.e = cortical.Amp.e;
final.cortical.Amp.i = cortical.Amp.i;
final.rfMap.maxResp.e = rfMap.maxResp.e;
final.rfMap.rfAreas.e = rfMap.rfAreas.e;
final.rfMap.maxResp.i = rfMap.maxResp.i;
final.rfMap.rfAreas.i = rfMap.rfAreas.i;
	
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
	main=paste("RF Extent vs RF Peak Response\nCortical Layer E-Cells\nInitial(.); Final(*)"),
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
	##############################################################
x11();
par(mfrow=c(2,2));

ofInterest = which ( (final.rfMap.rfAreas.e > mean(final.rfMap.rfAreas.e)) & (final.rfMap.maxResp.e < mean(final.rfMap.maxResp.e) ) );
tmp = rep ( 0,225); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Above Av RF Extent; Below Average Resp" );

ofInterest = which ( (final.rfMap.rfAreas.e > mean(final.rfMap.rfAreas.e)) & (final.rfMap.maxResp.e > mean(final.rfMap.maxResp.e) ) );
tmp = rep ( 0,225); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Above Av RF Extent; Above Average Resp" );

ofInterest = which ( (final.rfMap.rfAreas.e < mean(final.rfMap.rfAreas.e)) & (final.rfMap.maxResp.e < mean(final.rfMap.maxResp.e) ) );
tmp = rep ( 0,225); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Below Av RF Extent; Below Average Resp" );

ofInterest = which ( (final.rfMap.rfAreas.e < mean(final.rfMap.rfAreas.e)) & (final.rfMap.maxResp.e > mean(final.rfMap.maxResp.e) ) );
tmp = rep ( 0,225); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Below Av RF Extent; Above Average Resp" );

	######################################################
	#
	# 	Scatterplot: Layer C I-Cell Area vs Resp.
	#
	######################################################
x11();
xlim.i = c ( min(init.rfMap.maxResp.i,final.rfMap.maxResp.i), max(init.rfMap.maxResp.i,final.rfMap.maxResp.i) );
ylim.i = c ( min(init.rfMap.rfAreas.i,final.rfMap.rfAreas.i), max(init.rfMap.rfAreas.i,final.rfMap.rfAreas.i) );
plot ( init.rfMap.maxResp.i, init.rfMap.rfAreas.i, xlim=xlim.i, ylim=ylim.i, pch=".", col=3, cex=3,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer I-Cells\nInitial(.); Final(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.i, final.rfMap.rfAreas.i, pch="*", col=4 );
delta.area.i = delta.maxResp.i = slope.i = rep ( 0, N2 );
for ( i in 1:N2 ) {
	#lines ( c( init.rfMap.maxResp.i[i], final.rfMap.maxResp.i[i] ), c( init.rfMap.rfAreas.i[i], final.rfMap.rfAreas.i[i] ) )
	delta.area.i[i] = ( final.rfMap.rfAreas.i[i] - init.rfMap.rfAreas.i[i] );
	delta.maxResp.i[i] = ( final.rfMap.maxResp.i[i] - init.rfMap.maxResp.i[i] );
	slope.i[i] = ( final.rfMap.rfAreas.i[i] - init.rfMap.rfAreas.i[i] ) / ( final.rfMap.maxResp.i[i] - init.rfMap.maxResp.i[i] );
} # for ( i in 1:N2 ) {

	######################################################################
	#
 	# 	Scatterplot: Combined final state scatter plot E-cell, I-cell
	#
	######################################################################
x11();
xlim = c ( min ( xlim.i, xlim.e ), max ( xlim.i, xlim.e ) );
ylim = c ( min ( ylim.i, ylim.e ), max ( ylim.i, ylim.e ) );
plot ( final.rfMap.maxResp.e, final.rfMap.rfAreas.e, xlim=xlim, ylim=ylim, pch="*", col=1,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer E-Cells\nInitial(Black *); Final( Red *)"),
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
source ( "NMHelperFunctions.R" );
#xStim = Dig3MapRefine0( N, trialDurRFProbeInIters, oneSecondNumIter, refinementPatchSize, 0.25 );
xStim = Dig3MapSelStim0( N, trialDurRFProbeInIters, oneSecondNumIter, selectiveStimPatchSize, selectiveStimDuration, selectiveStimFactor, selectiveStimZoneID );
mapStimCount = xStim[[length(xStim)]]$stimCount;

		#
		#	Draw some pictures to visualize the data.
		#
xIDSubset = c ( 2, 4, 8 );
#xIDSubset = c ( 1, 3, 5 );

		#
		#	CDF: E Cell CORTICAL MAGNIFICATION
		#
init.CMag.e = apply ( init.cortical.Amp.e, 2, sum );
final.CMag.e = apply ( final.cortical.Amp.e, 2, sum );

x11(); par(mfrow=c(2,2));
xLabelText = "Cortical Magnification";
yLabelText = "CDF";
zTypeText = "E Cell (Final)"
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
subText = "Init.(Dash); Final(Solid); 2x(B); 4x(R); >8x(G)";
CorticalMagCDF.UseDep ( mapStimCount, init.CMag.e, final.CMag.e, titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

		#
		#	CDF: E Cell RECEPTIVE FIELD EXTENT
		#
xLabelText = "RF Extent";
yLabelText = "CDF";
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
subText = "Init.(Dash); Final(Solid); 2x(B); 4x(R); >8x(G)";
InvCorticalMagCDF.UseDep ( mapStimCount, init.cortical.Amp.e, init.rfMap.rfAreas.e, final.cortical.Amp.e, final.rfMap.rfAreas.e,
				titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

		#
		#	CDF: I Cell CORTICAL MAGNIFICATION
		#
init.CMag.i = apply ( init.cortical.Amp.i, 2, sum );
final.CMag.i = apply ( final.cortical.Amp.i, 2, sum );

x11(); par(mfrow=c(2,2));
xLabelText = "Cortical Magnification";
yLabelText = "CDF";
zTypeText = "I Cell (Final)"
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
subText = "Init.(Dash); Final(Solid); 2x(B); 4x(R); >8x(G)";
CorticalMagCDF.UseDep ( mapStimCount, init.CMag.i, final.CMag.i, titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

		#
		#	CDF: I Cell RECEPTIVE FIELD EXTENT
		#
xLabelText = "RF Extent";
yLabelText = "CDF";
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
subText = "Init.(Dash); Final(Solid); 2x(B); 4x(R); >8x(G)";
InvCorticalMagCDF.UseDep ( mapStimCount, init.cortical.Amp.i, init.rfMap.rfAreas.i, final.cortical.Amp.e, final.rfMap.rfAreas.e,
				titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

		#
		#	WMM: CORTICAL MAGNIFICATION
		#
pValue = 1e-5;

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
mainText2 = "\nInitial Conds. vs Refined";
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
mainText2 = "\nInitial Conds. vs Refined";
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
mainText2 = "\nInitial Conds. vs Refined";
mainText3 = "\nC Layer E Cells";
subText = "Init.(Dash); Final(Solid); D1(B); D2(R); D3(G)";
for ( iLongAxis in seq ( (N/3), 1, -1 ) ) {
	D3AllLongAxisRFDataCDF ( iLongAxis, init.cortical.Amp.e, init.rfMap.rfAreas.e,
						final.cortical.Amp.e, final.rfMap.rfAreas.e,
						paste(mainText1, iLongAxis, mainText2, mainText3, sep="" ),
						subText, paste("Receptive Field Extent"), paste("CDF"), FALSE );
	abline(v=25,lty=4,col=4);
	abline(v=30,lty=4,col=5);
} # for ( iDigit in seq ( 3, 1, -1 ) ) {

x11(); par(mfcol=c(3,2));
mainText1 = "CDF of RF Extent Long. Axis ";
mainText2 = "\nInitial Conds. vs Refined";
mainText3 = "\nC Layer I Cells";
subText = "Init.(Dash); Final(Solid); D1(B); D2(R); D3(G)";
for ( iLongAxis in seq ( (N/3), 1, -1 ) ) {
	D3AllLongAxisRFDataCDF ( iLongAxis, init.cortical.Amp.i, init.rfMap.rfAreas.i,
						final.cortical.Amp.i, final.rfMap.rfAreas.i,
						paste(mainText1, iLongAxis, mainText2, mainText3, sep="" ),
						subText, paste("Receptive Field Extent"), paste("CDF"), FALSE );
	abline(v=25,lty=4,col=4);
	abline(v=30,lty=4,col=5);
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

