######################################################################################################
######################################################################################################
#
#	Dig3SectorLocs1: return the list of locations (in linear representation) corresponding to each
#				of the possible sectors in a three digit layer S partition, e.g., 9.
#
#		Input:
#			N - length (or width) of the input array
#
#		Output:
#			sectLocs - an array where ith COLUMN lists the (N/3)^2 linear representation
#					locations of the ith sector
#
######################################################################################################
######################################################################################################

Dig3SectorLocs1 = function ( N ) {

	numSectors = 9;
	sectorLength = as.integer ( (N/3) );
	sLocs = matrix ( 0, nrow=(sectorLength * sectorLength), ncol=numSectors );
	kSector = 1;

	for ( iSec in 1:3 ) {

		iRowStart = ( iSec - 1 ) * sectorLength + 1;
		iRowEnd = ( iSec - 1 ) * sectorLength + sectorLength;

		for ( jSec in 1:3 ) {

			iColStart = ( jSec - 1 ) * sectorLength + 1;
			iColEnd = ( jSec - 1 ) * sectorLength + sectorLength;			
			tmp = 0;
			for ( kCol in iColStart:iColEnd ) {
				tmp = c ( tmp, GetLin ( iRowStart:iRowEnd, kCol, N ) );
			} # for ( kCol in iColStart:iColEnd )
			sLocs[,kSector] = (tmp[2:length(tmp)]);
			kSector = kSector + 1;

		} # for ( jSec in 1:3)

	} # for ( iSec in 1:3 )

	return ( sLocs )

} # Dig3SectorLocs1 = function ( ) {



#xtmp = Dig3SectorLocs1 ( 15 );



