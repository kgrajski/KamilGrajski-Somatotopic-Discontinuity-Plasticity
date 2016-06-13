######################################################################################################
######################################################################################################
#
#	TopoQuant1: Return a scalar representing neighborhood preservation
#			between the "skin" input layer and an upper layer.
#
#			This routine implements the case where F and G are both
#			Euclidean distances.  F will be the pairwise distances
#			representing the coordinates of each input layer cell.
#			G will be the pairwise distances over the positions of
#			the receptive field centers.
#
#			All the distances are 2D.
#	
#	Input:
#		N - this implies the size and positions to use for computing input layer
#		g0 - a standard deviation figure to use for defining neighborhood
#		rfTopo - this is the list returned by ExtractRFPositionData which includes
#				rfCentroid and rfCovar.
#
#		This routine uses only the rfCentroid (Euclidean distance).
#		
#
#	Output:
#		C - a scalar
#
######################################################################################################
######################################################################################################

TopoQuant1= function ( N, g0, rfTopo ) {

		#	Generate cell positions for the N x N input layer.

	N2 = N^2;
	inData = matrix ( 0, nrow=N2, ncol=2 );
	iTmp = 1;
	for ( iCell in 1:N2 ) {

		inData[iTmp, 1] = GetRow ( iCell, N );
		inData[iTmp, 2] = GetCol ( iCell, N );
		iTmp = iTmp + 1;

	} # for ( icol = 1:N )
	tmpF = dist ( inData, method="euclidean" );
	tmpF[(tmpF > g0 * sqrt(2.0))] = 0.

		#	Generate the pairwise distances for the input positions of RF centroids.
	tmpG = dist ( rfTopo$rfCentroid, method="euclidean" );
	tmpG[(tmpG > g0 * sqrt(2.0))] = 0.

	qVal = sum ( tmpF * tmpG ) / (N2 * (N2 - 1) );
	#qVal = sum ( tmpF * tmpG );

	return ( qVal );

} # TopoQuant1= function ( x ) { {




######################################################################################################
######################################################################################################
#
#	Long-term state topo map Quantitative evolution display.
#
#	Run once (which is slow) to gather up all of the topo map data into a single list.
#	That way can rerun or do various other experiments faster.
#
#######################################################################################################
######################################################################################################

rm(list = ls());
source ( "NMHelperFunctions.R" );

require(ggplot2);
require(reshape2);
require(ellipse);
require(stats4);

iBase = 25;
iStart = 0;
iEnd = 25;
#iStep = 5;

evol.rfTrackData.e = list();
evol.rfTrackData.i = list();

evol.quant.topoMap.e = rep ( 0, iEnd - iStart + 1 );
evol.quant.topoMap.i = rep ( 0, iEnd - iStart + 1 );

iCount = 1;
for ( iRefinement in iStart:iEnd ) {

	fileName = paste ( "Run.", iBase, ".", iRefinement, sep="" );
	load ( file = fileName );

	rfTopoMap.e = ExtractRFPositionData ( rfMap$r1.e.rfMap, kRFPeakToEdgeDetect );
	rfTopoMap.i = ExtractRFPositionData ( rfMap$r1.i.rfMap, kRFPeakToEdgeDetect );

	evol.rfTrackData.e[[iCount]] = rfTopoMap.e;
	evol.rfTrackData.i[[iCount]] = rfTopoMap.i;

	iCount = iCount + 1;	
	
} # for ( iRefinement in iStart:iEnd ) 
numCount = iCount - 1;


for ( iCount in 1:numCount) {

	evol.quant.topoMap.e[iCount] = TopoQuant1 ( N, g0, evol.rfTrackData.e[[iCount]] );
	evol.quant.topoMap.i[iCount] = TopoQuant1 ( N, g0, evol.rfTrackData.i[[iCount]] );

} # for ( iRefinement in iStart:iEnd ) {


x11(); par(mfrow=c(2,1));
ylim = c ( min(c(evol.quant.topoMap.e, evol.quant.topoMap.i)), max(c(evol.quant.topoMap.e, evol.quant.topoMap.i)) );
plot ( evol.quant.topoMap.e, type="l", col=1, ylim=ylim, xlab="Iteration", ylab="C Metric",
	main="Topography Measure vs Refinement Iterations\nE Cell Layer (Black); I Cell Layer (Red)" );
lines ( evol.quant.topoMap.i, col=2 );






