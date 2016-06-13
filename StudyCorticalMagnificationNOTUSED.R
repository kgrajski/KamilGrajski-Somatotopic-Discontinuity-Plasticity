
######################################################################################################
######################################################################################################
#
#	StudyCorticalMagnification.R
#
#	Do some scatterplots: RF extent vs response magnitude.
#	Do some plots to look at cortical magnification, inverse magnification.#
#
######################################################################################################
######################################################################################################

rm(list = ls());

library(ggplot2);
library(reshape2);

source ( "NMHelperFunctions.R" );
source ( "Dig3SectorLocs1.R" );
source ( "D3LongAxisLocs1.R" );
source ( "D3CorticalMagnificationPlots.R" );
source ( "D3RFDataBoxPlots.R" );
source ( "InvCortMagRuleCheck.R" );
source ( "MapSectorNumberToLabel.R" );
source ( "D3MannWhitneyTestFrontEnd.R" );
source ( "LabelledCDF.R" );

	#
	#	Get the network with random initial conditions.
	#
load(file="Run.25.0");

		# Earlier versions of RFMapExp1 didn't do this.
rfMap.maxResp.e = apply ( rfMap$r1.e.rfMap, 1, max );
rfMap.maxResp.i = apply ( rfMap$r1.i.rfMap, 1, max );

init.cortical.Amp.e = cortical.Amp.e;
init.cortical.Amp.i = cortical.Amp.i;
init.rfMap.maxResp.e = rfMap.maxResp.e;
init.rfMap.rfAreas.e = rfMap.rfAreas.e;
init.rfMap.maxResp.i = rfMap.maxResp.i;
init.rfMap.rfAreas.i = rfMap.rfAreas.i;

	#
	#	Get the refined network.
	#
load(file="Run.25.25");

		# 	Earlier versions of RFMapExp1 didn't do this.
rfMap.maxResp.e = apply ( rfMap$r1.e.rfMap, 1, max );
rfMap.maxResp.i = apply ( rfMap$r1.i.rfMap, 1, max );

final.cortical.Amp.e = cortical.Amp.e;
final.cortical.Amp.i = cortical.Amp.i;
final.rfMap.maxResp.e = rfMap.maxResp.e;
final.rfMap.rfAreas.e = rfMap.rfAreas.e;
final.rfMap.maxResp.i = rfMap.maxResp.i;
final.rfMap.rfAreas.i = rfMap.rfAreas.i;

	#
	# 	Scatterplot of Layer C E-Cell Area vs Resp.
	#
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

	#
	#	Draw plots showing the spatial distributions of sectors of the scatterplot.
	#
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

	#
	#	Determine which receptive field extents got bigger, smaller, stayed the same.
	#	Generate summaries of how RF extent changed.
	#	Generate summaries of how the response magnitude depending upon whether the
	#	RF extent got bigger, smaller, or stayed the same.
	#
area.lt0.e = ( delta.area.e<0 );
area.eq0.e = ( delta.area.e==0 );
area.gt0.e = ( delta.area.e>0 );
area.count = c ( sum(area.lt0.e), sum(area.eq0.e), sum(area.gt0.e) );
area.count.pct = 100 * area.count/N2;

summary ( -delta.area.e[area.lt0.e]/init.rfMap.rfAreas.e[area.lt0.e] )
summary ( delta.area.e[area.eq0.e] )
summary ( delta.area.e[area.gt0.e]/init.rfMap.rfAreas.e[area.gt0.e] )
summary ( delta.maxResp.e[area.lt0.e]/init.rfMap.maxResp.e[area.lt0.e] )
summary ( delta.maxResp.e[area.gt0.e]/init.rfMap.maxResp.e[area.gt0.e] )
summary ( delta.maxResp.e[area.eq0.e]/init.rfMap.maxResp.e[area.eq0.e] )

		#
		# 	Scatterplot of Layer C I-Cell Area vs Resp.
		#
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

		#
 		# 	Combined final state scatter plot E-cell, I-cell
		#
x11();
xlim = c ( min ( xlim.i, xlim.e ), max ( xlim.i, xlim.e ) );
ylim = c ( min ( ylim.i, ylim.e ), max ( ylim.i, ylim.e ) );
plot ( final.rfMap.maxResp.e, final.rfMap.rfAreas.e, xlim=xlim, ylim=ylim, pch="*", col=1,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer E-Cells\nInitial(Black *); Final( Red *)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.i, final.rfMap.rfAreas.i, pch="*", col=2 );


		#
		#	By Sector Cortical Magnification Plots
		#
D3CorticalMagnificationPlots ( final.rfMap.rfAreas.e, "Cortical Layer E Cells" );
D3CorticalMagnificationPlots ( final.rfMap.rfAreas.i, "Cortical Layer I Cells" );

D3PairedCorticalMagnificationPlots ( init.rfMap.rfAreas.e, final.rfMap.rfAreas.e, "Init C E", "Final C E" );
D3PairedCorticalMagnificationPlots ( init.rfMap.rfAreas.i, final.rfMap.rfAreas.i, "Init C I", "Final C I" );
D3PairedCorticalMagnificationPlots ( final.rfMap.rfAreas.i, final.rfMap.rfAreas.e, "Final C I", "Final C E" );




			#	Boxplots: Look at all of the sectors on a hand.
x11(); par(mfrow=c(2,1));
D3AllSectorsRFDataDistr ( init.rfMap.rfAreas.e, "RF Extent", "Init C E", FALSE, ylimMax=50 );
abline ( h = mean ( init.rfMap.rfAreas.e ), lty=4, col=1 );
abline ( h = max ( init.rfMap.rfAreas.e ), lty=4, col=5 );
D3AllSectorsRFDataDistr ( final.rfMap.rfAreas.e, "RF Extent", "Final C E", FALSE, ylimMax=50 );
abline ( h = mean ( init.rfMap.rfAreas.e ), lty=4, col=1 );
abline ( h = max ( init.rfMap.rfAreas.e ), lty=4, col=5 );

x11(); par(mfrow=c(2,1));
D3AllSectorsRFDataDistr ( init.rfMap.maxResp.e, "RF Max Resp", "Init C E", FALSE );
D3AllSectorsRFDataDistr ( final.rfMap.maxResp.e, "RF Max Resp", "Final C E", FALSE );

x11(); par(mfrow=c(2,1));
D3AllSectorsRFDataDistr ( init.rfMap.rfAreas.i, "RF Extent", "Init C I", FALSE );
D3AllSectorsRFDataDistr ( final.rfMap.rfAreas.i, "RF Extent", "Final C I", FALSE );

			#	Boxplots: Sector Pair-wise compare Init, Final, CE.
x11(); par(mfrow=c(3,3));
for ( iDigit in seq ( 3, 1, -1 ) ) {
	for ( iSector in seq ( 1, 3, 1 ) ) {
		iWhich = iDigit + (iSector - 1)*3;
		xLabel = MapSectorNumberToLabel ( iSector );
		D3PairSectorsRFDataDistr ( init.rfMap.rfAreas.e, iWhich, final.rfMap.rfAreas.e, iWhich, "RF Extent", 
			paste("IC CE"," D",iDigit,xLabel,sep=""), paste("Fin CE"," D",iDigit,xLabel,sep=""), FALSE );
	} # for ( iSector in seq ( 1, 3, 1 ) ) {
} # for ( iDigit in seq ( 3, 1, -1 ) ) {

			#	Boxplots: Sector Pair-wise compare Init, Final, CI.
x11(); par(mfrow=c(3,3));
for ( iDigit in seq ( 3, 1, -1 ) ) {
	for ( iSector in seq ( 1, 3, 1 ) ) {
		iWhich = iDigit + (iSector - 1)*3;
		xLabel = MapSectorNumberToLabel ( iSector );
		D3PairSectorsRFDataDistr ( init.rfMap.rfAreas.i, iWhich, final.rfMap.rfAreas.i, iWhich, "RF Extent", 
			paste("IC CI"," D",iDigit,xLabel,sep=""), paste("Fin CI"," D",iDigit,xLabel,sep=""), FALSE );
	} # for ( iSector in seq ( 1, 3, 1 ) ) {
} # for ( iDigit in seq ( 3, 1, -1 ) ) {


			#	Boxplots: Sector Pair-wise compare Final E, Final I.
x11(); par(mfrow=c(3,3));
for ( iDigit in seq ( 3, 1, -1 ) ) {
	for ( iSector in seq ( 1, 3, 1 ) ) {
		iWhich = iDigit + (iSector - 1)*3;
		xLabel = MapSectorNumberToLabel ( iSector );
		D3PairSectorsRFDataDistr ( final.rfMap.rfAreas.e, iWhich, final.rfMap.rfAreas.i, iWhich, "RF Extent",
			paste("FIN CI"," D",iDigit,xLabel,sep=""), paste("FIN CE"," D",iDigit,xLabel,sep=""), FALSE );
	} # for ( iSector in seq ( 1, 3, 1 ) ) {
} # for ( iDigit in seq ( 3, 1, -1 ) ) {

		#
		#	Longitudonal Axis: Some examples.		#

			#	Boxplots: Compare indiviual longitudonal tracks, such two different tracsk as on a single digit.
D3PairLongAxisRFDataDistr( final.rfMap.rfAreas.e, 1, 1, final.rfMap.rfAreas.e, 1, 3, "Extent", "Final C E", "Final C E" );

			#	Boxplots: Compare indiviual longitudonal tracks, such same track different conditions.
D3PairLongAxisRFDataDistr( init.rfMap.rfAreas.e, 2, 3, final.rfMap.rfAreas.e, 2, 3, "Extent", "Init C E", "Final C E" );

			#	Boxplots: Look at all of the tracks on a single digit.
x11(); par(mfrow=c(3,1));
D3SingleDigitLongAxisRFDataDistr( init.rfMap.rfAreas.e, 3, "RF Extent", "Init C E", FALSE  );
D3SingleDigitLongAxisRFDataDistr( init.rfMap.rfAreas.e, 2, "RF Extent", "Init C E", FALSE  );
D3SingleDigitLongAxisRFDataDistr( init.rfMap.rfAreas.e, 1, "RF Extent", "Init C E", FALSE  );

x11(); par(mfrow=c(3,1));
D3SingleDigitLongAxisRFDataDistr( final.rfMap.rfAreas.e, 3, "RF Extent", "Final C E", FALSE  );
D3SingleDigitLongAxisRFDataDistr( final.rfMap.rfAreas.e, 2, "RF Extent", "Final C E", FALSE  );
D3SingleDigitLongAxisRFDataDistr( final.rfMap.rfAreas.e, 1, "RF Extent", "Final C E", FALSE  );

x11(); par(mfrow=c(3,1));
D3SingleDigitLongAxisRFDataDistr( init.rfMap.rfAreas.i, 3, "RF Extent", "Init C I", FALSE  );
D3SingleDigitLongAxisRFDataDistr( init.rfMap.rfAreas.i, 2, "RF Extent", "Init C I", FALSE  );
D3SingleDigitLongAxisRFDataDistr( init.rfMap.rfAreas.i, 1, "RF Extent", "Init C I", FALSE  );

x11(); par(mfrow=c(3,1));
D3SingleDigitLongAxisRFDataDistr( final.rfMap.rfAreas.i, 3, "RF Extent", "Final C I", FALSE  );
D3SingleDigitLongAxisRFDataDistr( final.rfMap.rfAreas.i, 2, "RF Extent", "Final C I", FALSE  );
D3SingleDigitLongAxisRFDataDistr( final.rfMap.rfAreas.i, 1, "RF Extent", "Final C I", FALSE  );

x11(); par(mfrow=c(2,2));
D3SingleDigitLongAxisRFDataDistr( init.rfMap.maxResp.e, 3, "RF Resp", "Init C E", FALSE  );
D3SingleDigitLongAxisRFDataDistr( final.rfMap.maxResp.e, 3, "RF Resp", "Final C E", FALSE  );
D3SingleDigitLongAxisRFDataDistr( init.rfMap.rfAreas.e, 3, "RF Extent", "Init C E", FALSE  );
D3SingleDigitLongAxisRFDataDistr( final.rfMap.rfAreas.e, 3, "RF Extent", "Final C E", FALSE  );

x11(); par(mfrow=c(2,2));
D3SingleDigitLongAxisRFDataDistr( init.rfMap.maxResp.e, 2, "RF Resp", "Init C E", FALSE  );
D3SingleDigitLongAxisRFDataDistr( final.rfMap.maxResp.e, 2, "RF Resp", "Final C E", FALSE  );
D3SingleDigitLongAxisRFDataDistr( init.rfMap.rfAreas.e, 2, "RF Extent", "Init C E", FALSE  );
D3SingleDigitLongAxisRFDataDistr( final.rfMap.rfAreas.e, 2, "RF Extent", "Final C E", FALSE  );


	############################################################################################
	############################################################################################
	#
	#	Examine the data for evidence for the (inverse) cortical magnification rule.
	#
	############################################################################################
	############################################################################################

		#	Extract from the stimulus pattern the relative frequencies of stimulation.
source ( "NMHelperFunctions.R" );
xStim = Dig3MapRefine0( N, trialDurRFProbeInIters, oneSecondNumIter, refinementPatchSize, 0.25 );
mapStimCount = xStim[[length(xStim)]]$stimCount;

		#	Draw some pictures to visualize the data.
xIDSubset = c ( 2, 4, 8 );

			#	Generate a plot of a low, mid and high relative frequency of stimulation
			#	CORTICAL MAGNIFICATION & RECEPTIVE FIELD EXTENT
init.CMag.e = apply ( init.cortical.Amp.e, 2, sum );
final.CMag.e = apply ( final.cortical.Amp.e, 2, sum );

x11(); par(mfrow=c(2,1));
xLabelText = "Cortical Magnification";
yLabelText = "CDF";
zTypeText = "E Cell (Final)"
titleText = paste("CDF of ", xLabelText, "\nBy S Layer Rel. Freq. Stim.\n", zTypeText, sep="" );
subText = "";
LabelledCDF ( mapStimCount, init.CMag.e, final.CMag.e, titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

			#	RECEPTIVE FIELD EXTENT
xLabelText = "RF Extent";
yLabelText = "CDF";
titleText = paste("CDF of ", xLabelText, "\nBy Stim. Count\nE Cell (Final)", sep="" );
subText = "Initial (Dashed); Final (Solid); 2x (Black); 4x (Red); 8x (Green)";
LabelledCDF ( mapStimCount, init.rfMap.rfAreas.e, final.rfMap.rfAreas.e, titleText, subText, xLabelText, yLabelText, xIDSubset, FALSE );

		#	The following could be eventually upgraded to a routine.
		#	Do Wilcox-Mann-Whitney test on cortical magnification partitioned by relative frequency of stimulation.
pValue = 0.001;

init.CMag.MWTest.e = D3RelFreqStimMannWhitneyTest ( mapStimCount, init.CMag.e, init.CMag.e, pValue );
init.CMag.MWTest.e$pValues[init.CMag.MWTest.e$acceptAltHyp,]

final.CMag.MWTest.e = D3RelFreqStimMannWhitneyTest ( mapStimCount, init.CMag.e, init.CMag.e, pValue  );
final.CMag.MWTest.e$pValues[final.CMag.MWTest.e$acceptAltHyp,]

init.final.CMag.MWTest.e = D3RelFreqStimMannWhitneyTest ( mapStimCount, init.CMag.e, init.CMag.e, pValue  );
init.final.CMag.MWTest.e$pValues[init.final.CMag.MWTest.e$acceptAltHyp,]

		#	Do Wilcox-Mann-Whitney test on RF extent partitioned by relative frequency of stimulation.
init.rfAreas.MWTest.e = D3RelFreqStimMannWhitneyTest ( mapStimCount, init.rfMap.rfAreas.e, init.rfMap.rfAreas.e, pValue );
init.rfAreas.MWTest.e$pValues[init.rfAreas.MWTest.e$acceptAltHyp,]

final.rfAreas.MWTest.e = D3RelFreqStimMannWhitneyTest ( mapStimCount, final.rfMap.rfAreas.e, final.rfMap.rfAreas.e, pValue  );
final.rfAreas.MWTest.e$pValues[final.rfAreas.MWTest.e$acceptAltHyp,]

init.final.rfAreas.MWTest.e = D3RelFreqStimMannWhitneyTest ( mapStimCount, init.rfMap.rfAreas.e, final.rfMap.rfAreas.e, pValue  );
init.final.rfAreas.MWTest.e$pValues[init.final.rfAreas.MWTest.e$acceptAltHyp,]

		#	Compute the cortical magnifications sorted by rel stim freq(u, sd).
		#	Compute the RF extent sorted by rel stim freq(u, sd).
		#	Do some Wilcox tests on RF extens partitioned by rel stim freq.
init.invCMag.Chk.e = InvCortMagRuleCheck ( mapStimCount, init.cortical.Amp.e, 0 );
final.invCMag.Chk.e = InvCortMagRuleCheck ( mapStimCount, final.cortical.Amp.e, 0 );
init.invCMag.Chk.i = InvCortMagRuleCheck ( mapStimCount, init.cortical.Amp.i, 0 );
final.invCMag.Chk.i = InvCortMagRuleCheck ( mapStimCount, final.cortical.Amp.i, 0 );
ofInterest = which ( (final.rfMap.rfAreas.e < mean(final.rfMap.rfAreas.e)) & (final.rfMap.maxResp.e > mean(final.rfMap.maxResp.e) ) );
tmp = rep ( 0,225); tmp[ofInterest] = 1;
ShowVecAsMap ( tmp, "Below Av RF Extent; Above Average Resp" );





