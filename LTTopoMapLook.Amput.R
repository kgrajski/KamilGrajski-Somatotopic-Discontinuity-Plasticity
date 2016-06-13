######################################################################################################
######################################################################################################
#
#	Long-term state topo map evolution and various quantitative measure displays.
#
#	Run once (which is slow) to gather up all of the topo map data into a single list.
#
#	That way can rerun or do various other experiments faster.
#
#
#	Modified to ease analysis of the Amputation Experiments.
#
#######################################################################################################
######################################################################################################

#rm(list = ls());
source ( "NMHelperFunctions.R" );

require(ggplot2);
require(reshape2);
require(ellipse);
require(stats4);

iBase = 25;
iStart = 0;
iEnd = 25;

	#	Allocate arrays to compute and capture receptive field data from each iteration.
evol.rfTrackData.e = list();
evol.rfTrackData.i = list();

evol.cortical.Amp.e = list()
evol.cortical.Amp.i = list()

evol.quant.topoMap.e = rep ( 0, iEnd - iStart + 1 );
evol.quant.topoMap.i = rep ( 0, iEnd - iStart + 1 );

	#	Get the baseline topo map result.
iCount = 1;
load ( file = paste ( "Run.25.25" ) );
source ( "NMHelperFunctions.R" );
evol.rfTrackData.e[[iCount]] = ExtractRFPositionData ( rfMap$r1.e.rfMap, kRFPeakToEdgeDetect );
evol.rfTrackData.i[[iCount]] = ExtractRFPositionData ( rfMap$r1.i.rfMap, kRFPeakToEdgeDetect );
evol.cortical.Amp.e[[iCount]] = QuantCorticalAmp ( rfMap$r1.e.rfMap, kRFPeakToEdgeDetect );
evol.cortical.Amp.i[[iCount]] = QuantCorticalAmp ( rfMap$r1.i.rfMap, kRFPeakToEdgeDetect );
iCount = iCount + 1;
	
	#	Now get the set of post-amputation refinement results.
for ( iRefinements in iStart:iEnd ) {
		
	load ( file = paste ("Amput", kAmputDigit, kAmputZoneMax, numAmputRefinements, iRefinements, sep="." ));

	evol.rfTrackData.e[[iCount]] = ExtractRFPositionData ( rfMap$r1.e.rfMap, kRFPeakToEdgeDetect );
	evol.rfTrackData.i[[iCount]] = ExtractRFPositionData ( rfMap$r1.i.rfMap, kRFPeakToEdgeDetect );

	evol.cortical.Amp.e[[iCount]] = QuantCorticalAmp ( rfMap$r1.e.rfMap, kRFPeakToEdgeDetect );
	evol.cortical.Amp.i[[iCount]] = QuantCorticalAmp ( rfMap$r1.i.rfMap, kRFPeakToEdgeDetect );

	iCount = iCount + 1;	
	
} # for ( iRefinement in iStart:iEnd ) {
numCount = iCount - 1;

	#	Store the aggregated data so that it needs to be done only once.
	#	Several of the additional quantitative metrics are derived from the 
save ( file="Amput.RData", list=c("evol.rfTrackData.e", "evol.rfTrackData.i", "evol.cortical.Amp.e", "evol.cortical.Amp.i" ) );
load ( file="Amput.RData" );

	#	Do a succession of topo map displays.  There are two sets.
	#	First set: display first and last iteration topo map and topo map "tracks".
numCount = length ( evol.rfTrackData.e );
for ( iCount in c(1,2,numCount) ) {

	ShowTopoMap ( evol.rfTrackData.e[[iCount]], paste("E-Type","Iter",iCount-1,sep=" "), TRUE, 0.5 );
	ShowThreeDigitRFTrack0 ( evol.rfTrackData.e[[iCount]], paste("E-Type","Iter",iCount-1,sep=" "), TRUE, 0.5 );

	ShowTopoMap ( evol.rfTrackData.i[[iCount]], paste("I-Type","Iter",iCount-1,sep=" "), TRUE, 0.5 );
	ShowThreeDigitRFTrack0 ( evol.rfTrackData.i[[iCount]], paste("I-Type","Iter",iCount-1,sep=" "), TRUE, 0.5 );

} # for ( iRefinement in iStart:iEnd ) {

	#	Second set: for each of several E- and I- cells show the progression of RF size over iterations.
x11(); par(mfcol=c(3,2));
iStart = 1; iEnd = numCount ; iStep = 1;
for ( tmpCell in c ( ceiling(N/2), ceiling(N2/2), N2-ceiling(N/2) ) ) {
	ShowTopoMapOverTimeForOneCell ( evol.rfTrackData.e, tmpCell, iStart, iEnd, iStep,
		paste("E-Type","Iter",tmpCell,sep=" "), TRUE, 0.5, TRUE );
} # for ( tmpCell in c ( 1, ceiling(N2/2), N2 ) ) {

for ( tmpCell in c ( ceiling(N/2), ceiling(N2/2), N2-ceiling(N/2) ) ) {
	ShowTopoMapOverTimeForOneCell ( evol.rfTrackData.i, tmpCell, iStart, iEnd, iStep,
		paste("I-Type","Cell",tmpCell,sep=" "), TRUE, 0.5, TRUE );
} # for ( tmpCell in c ( 1, ceiling(N2/2), N2 ) ) {

	#	Now do some topo map quantitative measurements and display results for E and I cells.
	#	First method: fancy style from the literature
x11(); par(mfrow=c(2,1));
for ( iCount in 1:numCount) {

	evol.quant.topoMap.e[iCount] = TopoQuant1 ( N, g0, evol.rfTrackData.e[[iCount]] );
	evol.quant.topoMap.i[iCount] = TopoQuant1 ( N, g0, evol.rfTrackData.i[[iCount]] );

} # for ( iRefinement in iStart:iEnd ) {

ylim = c ( min(c(evol.quant.topoMap.e, evol.quant.topoMap.i)), max(c(evol.quant.topoMap.e, evol.quant.topoMap.i)) );
plot ( evol.quant.topoMap.e, type="l", col=1, ylim=ylim, xlab="Iteration", ylab="C Metric",
	main="(F,G) Topography Measure vs Refinement Iterations\nE Cell Layer (Black); I Cell Layer (Red)" );
lines ( evol.quant.topoMap.i, col=2 );

	#	Second method: R^2 type style
for ( iCount in 1:numCount) {

	evol.quant.topoMap.e[iCount] = TopoQuant2 ( N, evol.rfTrackData.e[[iCount]] );
	evol.quant.topoMap.i[iCount] = TopoQuant2 ( N, evol.rfTrackData.i[[iCount]] );

} # for ( iRefinement in iStart:iEnd ) {
ylim = c ( min(c(evol.quant.topoMap.e, evol.quant.topoMap.i)), max(c(evol.quant.topoMap.e, evol.quant.topoMap.i)) );
plot ( evol.quant.topoMap.e, type="l", col=1, ylim=ylim, xlab="Iteration", ylab="C Metric",
	main="R^2-Type Topography Measure vs Refinement Iterations\nE Cell Layer (Black); I Cell Layer (Red)" );
lines ( evol.quant.topoMap.i, col=2 );

	#	Now plot the evolution of cortical amplification.
cAmp.e = cAmp.i = matrix ( 0, nrow=numCount, ncol=9 );
for ( iCount in 1:numCount ) {
	cAmp.e[iCount,] = ThreeDigitCorticalAmplification ( evol.cortical.Amp.e[[iCount]], FALSE, "" );
	cAmp.i[iCount,] = ThreeDigitCorticalAmplification ( evol.cortical.Amp.i[[iCount]], FALSE, "" );
} # for ( iCount in 1:numCount )

x11(); par(mfrow=c(2,1));
ylim = c ( 0, max(cAmp.e) );
plot( cAmp.e[,1], type="l", lty=1, col=1, ylim=ylim, main="Evolution Input Layer Amplification in E-Cells",
xlab="Iterations", ylab="Area", sub="Digits: 3(Solid), 2(Dash), 1(Dot); Sector: Distal(Blk), Mid(Red), Prox(Green)" );
lines ( cAmp.e[,2], type="l", lty=1, col=2 ); lines ( cAmp.e[,3], type="l", lty=1, col=3 );
lines( cAmp.e[,4], type="l", lty=2, col=1 ); lines ( cAmp.e[,5], type="l", lty=2, col=2 ); lines ( cAmp.e[,6], type="l", lty=2, col=3 ); 
lines( cAmp.e[,7], type="l", lty=3, col=1 ); lines ( cAmp.e[,8], type="l", lty=3, col=2 ); lines ( cAmp.e[,9], type="l", lty=3, col=3 );
abline(h=49, lty=3, col=4);


ylim = c ( 0, max(cAmp.i) );
plot( cAmp.i[,1], type="l", lty=1, col=1, ylim=ylim, main="Evolution Input Layer Amplification in E-Cells",
xlab="Iterations", ylab="Area",sub="Digits: 3(Solid), 2(Dash), 1(Dot); Sector: Distal(Blk), Mid(Red), Prox(Green)" );
lines ( cAmp.i[,2], type="l", lty=1, col=2 ); lines ( cAmp.i[,3], type="l", lty=1, col=3 );
lines( cAmp.i[,4], type="l", lty=2, col=1 ); lines ( cAmp.i[,5], type="l", lty=2, col=2 ); lines ( cAmp.i[,6], type="l", lty=2, col=3 ); 
lines( cAmp.i[,7], type="l", lty=3, col=1 ); lines ( cAmp.i[,8], type="l", lty=3, col=2 ); lines ( cAmp.i[,9], type="l", lty=3, col=3 ); 
abline(h=49, lty=3, col=4);











