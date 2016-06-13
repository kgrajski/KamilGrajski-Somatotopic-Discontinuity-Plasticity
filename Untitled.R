######################################################################################################
######################################################################################################
#
#	Plot Receptive Field Overlap Results
#
#
#	Input:
#		ovLap - a matrix where each column is of length N-1 corresponding to the successive
#				pairwise proportion overlap of receptive fields
#
######################################################################################################
######################################################################################################

PlotRFOverlap = function ( ovLap, mainTxt, xlabTxt, ylabTxt ) {

	numPlots = dim ( ovLap )[2];

	ylimExtra = 10;
	ylim = c ( min ( ovLap ) - ylimExtra, max ( ovLap, 100 ) );
	plot ( ovLap[,1], ylim=ylim, type="l", main=mainTxt, xlab=xlabTxt, ylab=ylabTxt );

	for ( i in 2:numPlots ) {

		lines ( ovLap[,i] );

	} # for ( i in 2:numPlots )	

} # PlotRFOverlap = function ( ovLap, mainTxt, xlabTxt, ylabTxt ) {


