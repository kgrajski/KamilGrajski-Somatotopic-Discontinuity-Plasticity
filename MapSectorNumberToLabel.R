######################################################################################################
######################################################################################################
#
#	MapSectorNumberToLabel.R
#
#	Do some scatterplots: RF extent vs response magnitude.
#	Do some plots to look at cortical magnification, inverse magnification.#
#
######################################################################################################
######################################################################################################

MapSectorNumberToLabel = function ( ix ) {

	xLabel = "d";
	if ( ix == 2 ) { xLabel = "m"; }
	if ( ix == 3 ) { xLabel = "p"; }

	return ( xLabel );

} # MapSectorNumberToLabel = function ( ix ) {