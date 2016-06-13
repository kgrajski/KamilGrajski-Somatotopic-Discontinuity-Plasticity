#
#	Long-term state variable evolution display.
#

######################################################################################################
######################################################################################################
#
#	PlotAllWeights: Plotting of the time evolution of weights.
#
######################################################################################################
######################################################################################################

PlotAllWeights = function ( tw, twWhichCell, titleText ) {

	nRow = dim(tw)[1];
	nCol = dim(tw)[2];
	tmpX = 1:nCol;

	whichRow = which( apply ( tw, 1, sum ) != 0 );
	ylim = c ( min ( tw[whichRow,] ), max ( tw[whichRow,] ) );
	plot ( tmpX, tw[whichRow[1],], ylim=ylim, type="l", col=1,
		main=paste(titleText, "Evolution Into Cell", twWhichCell),
		xlab="Iteration", ylab="Weight" );
	for ( i in 2:length(whichRow) ) {
		lines ( tmpX, tw[whichRow[i],] );
	} # for ( i in 2:nRow )

} # PlotAllWeights = function ( tw, twWhichCell ) {


######################################################################################################
######################################################################################################
#
#	WeightEvolAll: Drive plotting of the time evolution of weights.  The data may have to be
#				gathered across a bunch of files first.
#
######################################################################################################
######################################################################################################

WeightEvolAll = function ( N2, numIter, iwBase, iwStart, iwEnd, twWhichCell ) {

	numFiles = iwEnd - iwStart + 1;
	t.w1.e.0 = matrix ( 0, nrow = N2, ncol = ( numIter * numFiles ) );
	t.w1.e.e = matrix ( 0, nrow = N2, ncol = ( numIter * numFiles ) );
	t.w1.e.i = matrix ( 0, nrow = N2, ncol = ( numIter * numFiles ) );
	t.w1.i.e = matrix ( 0, nrow = N2, ncol = ( numIter * numFiles ) );

	ik = 1;
	for ( iRefinement in iwStart:iwEnd ) {

		#wfileName = paste ( "Run.", iwBase, ".", iRefinement, sep="" );
		load ( file=paste ( "Run.", iwBase, ".", iRefinement, sep="" ));
		startIter = (ik - 1) * numIter + 1; endIter = ik * numIter;
		t.w1.e.0[,startIter:endIter] = w1.e.0[,twWhichCell,];
		t.w1.e.e[,startIter:endIter] = w1.e.e[,twWhichCell,];
		t.w1.e.i[,startIter:endIter] = w1.e.i[,twWhichCell,];
		t.w1.i.e[,startIter:endIter] = w1.i.e[,twWhichCell,];
		ik = ik + 1;

	} # for ( iRefinement in iStart:iEnd ) {

	x11(); par(mfrow=c(2,2));
	PlotAllWeights ( t.w1.e.0, twWhichCell, "W1.E.0" );
	PlotAllWeights ( t.w1.e.e, twWhichCell, "W1.E.E" );
	PlotAllWeights ( t.w1.e.i, twWhichCell, "W1.E.I" );
	PlotAllWeights ( t.w1.i.e, twWhichCell, "W1.I.E" );

} # WeightEvolAll = function ( N2, numIter, iwBase, iwStart, iwEnd, twWhichCell ) {


iwBase = 5;
iwStart = 0;
iwEnd = 5;
twWhichCell = 113;
N2 = 225;
numIter = 150;
WeightEvolAll ( N2, numIter, iwBase, iwStart, iwEnd, twWhichCell ); 

