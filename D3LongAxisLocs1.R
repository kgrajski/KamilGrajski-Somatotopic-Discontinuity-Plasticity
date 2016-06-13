######################################################################################################
######################################################################################################
#
#	D3LongAxisLocs1: return the list of locations (in linear representation) corresponding to
#					the given digit and longitudonal track.
#
#		Input:
#			N - length (or width) of the input array
#			xDigit - digit of interest
#			xTrack - track of interest in the range 1...digitWidth
#
#		Output:
#			sectLocs - an array where ith COLUMN lists the (N/3)^2 linear representation
#					locations of the ith sector
#
######################################################################################################
######################################################################################################

D3LongAxisLocs1 = function ( N, xDigit, xTrack ) {

	return ( GetLin ( seq ( 1, N ), ( ( xDigit - 1 ) * ( N / 3 ) + xTrack), N ) );

} # D3LongAxisLocs1= function ( ) {


	#
	#	This variant will trim the ends off the list of the longitudonal axis node list.
	#
D3LongAxisLocs2 = function ( N, xDigit, xTrack, trimEnds ) {

	return ( GetLin ( seq ( 1 + trimEnds, N - trimEnds ), ( ( xDigit - 1 ) * ( N / 3 ) + xTrack), N ) );

} # D3LongAxisLocs1= function ( ) {


