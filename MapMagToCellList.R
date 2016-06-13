######################################################################################################
######################################################################################################
#
#	MapMagToCellList.R
#
#	Input:
#		xMag - magnificationMap (boolean matrix);
#				(i, j) entry tells whether layer N+1 cell j has within its receptive field
#						the ith cell on layer N; assume that layer N projects to layer N+1
#
#		fromList - the list of layer N cells of interest
#
#	Output:
#		cellList - a list (no duplicates) of layer N+1 cells that contain within their receptive
#				field any of the cells of interest of layer N
#
######################################################################################################
######################################################################################################

MapMagToCellList = function ( xMag, fromList ) {

	cellList = NULL;
	for ( k in 1:length(fromList) ) {
			cellList = c ( cellList, which( xMag[ , fromList[k] ]) );
	} # for ( i in 1:length(sLayerNodes.ref)
	
	return ( sort ( unique ( cellList ) ) );

} # MapMagToCellList = function ( xMag, fromList )


#
#	This variant incorporates the use of a stimCount matrix and cutoff.
#

MapMagToCellList1 = function ( xMag, fromList, stimCount, cutOff ) {

	cellList = NULL;
	for ( k in 1:length(fromList) ) {
		if ( stimCount[fromList[k]] == cutOff) {
			cellList = c ( cellList, which( xMag[ , fromList[k] ]) );
		} # ()
	} # for ( i in 1:length(sLayerNodes.ref)
	
	return ( sort ( unique ( cellList ) ) );

} # MapMagToCellList1 = function ( xMag, fromList 

MapMagToCellList1GE = function ( xMag, fromList, stimCount, cutOff ) {

	cellList = NULL;
	for ( k in 1:length(fromList) ) {
		if ( stimCount[fromList[k]] >= cutOff) {
			cellList = c ( cellList, which( xMag[ , fromList[k] ]) );
		} # ()
	} # for ( i in 1:length(sLayerNodes.ref)
	
	return ( sort ( unique ( cellList ) ) );

} # MapMagToCellList1 = function ( xMag, fromList 

#
#	Trim edge cells from the list of cells.
#		A value of lenEdgeTrim = 0 removes nothing;
#		A value of lenEdgeTrim = N > 0 removes the outer-most N rings;
#

TrimEdgesFromCellList = function ( iCellList, N, lenEdgeTrim ) {

	iLen = length ( iCellList );
	iRow = GetRow ( iCellList, N );
	iCol = GetCol ( iCellList, N );

	minRow = lenEdgeTrim + 1;
	maxRow = N - lenEdgeTrim;

	minCol = lenEdgeTrim + 1;
	maxCol = N - lenEdgeTrim;

	iWhich = (iRow >= minRow) & (iRow <= maxRow) & (iCol >= minCol) & (iCol <= maxCol);

	return ( GetLin ( iRow[iWhich], iCol[iWhich], N ) );	
} # TrimEdgesFromCellList = function ( iCellList, N, lenEdgeTrim ) {

	#	Return list of cells to trim.
EdgeTrimCellList = function ( iCellList, N, lenEdgeTrim ) {

	iLen = length ( iCellList );
	iRow = GetRow ( iCellList, N );
	iCol = GetCol ( iCellList, N );

	minRow = lenEdgeTrim;
	maxRow = N - lenEdgeTrim;

	minCol = lenEdgeTrim;
	maxCol = N - lenEdgeTrim;

	iWhich = (iRow <= minRow) | (iRow > maxRow) | (iCol <= minCol) | (iCol > maxCol);

	return ( GetLin ( iRow[iWhich], iCol[iWhich], N ) );

} # TrimEdgesFromCellList = function ( iCellList, N, lenEdgeTrim ) {


GetLongAxisCorticalMag = function ( magMap.ref, N, iDigit, iLongAxis, trimEnds ) {

	mapMag = rep ( 0, N*N );
	tmp = MapMagToCellList ( magMap.ref,  D3LongAxisLocs2 ( N, iDigit, iLongAxis, trimEnds ) );
	mapMag[tmp] = 1.0;
	return ( mapMag );

} # PlotLongAxisCorticalMag = function ( N, iDigit, iLongAxis, trimEnds ) {





