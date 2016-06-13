#	D3RFDataBoxPlots.R


######################################################################################################
######################################################################################################
#
#	D3PairLongAxisRFDataDistr:
#
#		Do a boxplot for a single pair of longitudonal axes where one is .ref and one is .exp.
#
######################################################################################################
######################################################################################################

D3PairLongAxisRFDataDistr= function ( rfData.ref, iDigit.ref, iSector.ref,
							rfData.exp, iDigit.exp, iSector.exp, dataType, titleText.ref, titleText.exp ) {

	x11();

	xMat = cbind ( rfData.ref[ D3LongAxisLocs1(N, iDigit.ref, iSector.ref) ],
					rfData.exp[ D3LongAxisLocs1(N, iDigit.exp, iSector.exp) ] );
	ylim = c ( 0, max ( xMat ) ) * 1.1; 	#	Add 10% buffer
	xMat = as.data.frame ( xMat );
	names(xMat)[1] = paste ( "Ref: Dig",iDigit.ref,"Track",iSector.ref );
	names(xMat)[2] = paste ( "Exp: Dig",iDigit.exp,"Track",iSector.exp );

	boxplot ( xMat, ylim=ylim, main=paste( paste("Distr. ",dataType," in Longitudonal Axes",sep=""), paste(titleText.ref, " : ", titleText.exp, sep=""), sep="\n" ),
				ylab=dataType );

} # D3PairLongAxisRFDataDistr= function ( rfData, titleText ) {

######################################################################################################
######################################################################################################
#
#	D3SingleDigitLongAxisRFDataDistr:
#
#		Do a boxplot for all the longitudonal axes on a single digit.
#
######################################################################################################
######################################################################################################

D3SingleDigitLongAxisRFDataDistr= function ( rfData, xDigit, dataType, titleText, newX11Flag ) {

	N = as.integer ( sqrt ( length ( rfData ) ) );
	numDigits = 3;
	numTracks = as.integer ( N / 3 );
	
	xMat = matrix ( 0, nrow = N, ncol = numTracks );
	for ( i in 1:numTracks ) {
		xMat[,i] = rfData[ D3LongAxisLocs1(N, xDigit, i) ];
	} # for ( i in 1:numTracks ) {

	ylim = c ( 0, max(xMat) ) * 1.1;	#	Add 10% buffer. 
	xMat = as.data.frame(xMat);
	for ( i in 1:numTracks ) {
		names(xMat)[i] = paste ( "D", xDigit, ".", i, sep="" );
	} # for ( i in 1:numTracks ) {

	if ( newX11Flag ) { x11(); }
	boxplot ( as.data.frame(xMat), main=paste(paste("Distr. ",dataType," in Longitudonal Axes",sep=""),titleText,sep="\n"),
			ylab=dataType);

} # D3SingleDigitLongAxisRFDataDistr= function ( rfData, titleText ) {

######################################################################################################
######################################################################################################
#
#	D3AllSectorsRFDataDistr:
#
#		Do a boxplot for each of the sectors in the model.
#
######################################################################################################
######################################################################################################

D3AllSectorsRFDataDistr = function ( rfData, dataType, titleText, newX11Flag, ylimMax=-1 ) {

	sLocs = Dig3SectorLocs1 ( as.integer ( sqrt ( length ( rfData ) ) ) );
	numSectors = dim(sLocs)[2];
	xMat = matrix ( 0, nrow = dim(sLocs)[1], ncol = numSectors );
	for ( i in 1:numSectors ) {
		xMat[,i] = rfData[ sLocs[,i] ];
	} # for ( i in 1:numSectors ) {

	if ( ylimMax > 0 ) {
		ylim = c ( 0, ylimMax );
	} else {
		ylim = c ( 0, max(xMat) ) * 1.1;	#	Add 10% buffer
	} # if ( ylimMax > 0 )
	xMat = as.data.frame(xMat);
	i = 1;
	for ( iSecLoc in 1:3 ) {
		secLabel = "p";
		if ( iSecLoc == 1 ) { secLabel = "d"; } 
		if ( iSecLoc == 2 ) { secLabel = "m"; }

		for ( iDigit in 1:3 ) {
			names(xMat)[(iSecLoc-1)*3 + iDigit] = paste ( iDigit, secLabel, sep="" );
		} # for ( i in 1:numTracks ) {
		i = i + 1;
	} # for ( iSecLoc in 1:3 )

	if ( newX11Flag ) { x11(); }
	boxplot ( as.data.frame(xMat), ylim=ylim, main=paste( paste("Distr. ",dataType," by Sector",sep=""), titleText, sep="\n" ),
			ylab=dataType);

} # D3AllSectorsRFDataDistr= function ( rfData, titleText, ... ) {

######################################################################################################
######################################################################################################
#
#	D3PairSectorsRFDataDistr:
#
#		Do a boxplot for a pair of sectors in the model.
#
######################################################################################################
######################################################################################################

D3PairSectorsRFDataDistr = function ( rfData.ref, iSector.ref, rfData.exp, iSector.exp, dataType, titleText.ref, titleText.exp, newX11Flag ) {

	numSectors = 2;

	sLocs = Dig3SectorLocs1 ( as.integer ( sqrt ( length ( rfData.ref ) ) ) );

	xMat = matrix ( 0, nrow = dim(sLocs)[1], ncol = numSectors );
	xMat[,1] = rfData.ref[ sLocs[,iSector.ref] ];
	xMat[,2] = rfData.exp[ sLocs[,iSector.exp] ];
	ylim = c ( 0, max ( xMat ) ) * 1.1;	# Add 10%
	xMat = as.data.frame(xMat);

	names(xMat)[1] = "Ref Sector";
	names(xMat)[2] = "Exp Sector";

	if ( newX11Flag ) { x11(); }
	boxplot ( as.data.frame(xMat), ylim=ylim,
			main=paste( paste("Distr. ",dataType," by Sector",sep=""), paste(titleText.ref,"v",titleText.exp,sep=" "), sep="\n" ),
			ylab=dataType );

} # D3PairSectorsRFDataDistr= function ( rfData, titleText ) {