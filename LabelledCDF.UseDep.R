######################################################################################################
######################################################################################################
#
#	CorticalMagCDF.UseDep 
#
#		Generate histograms of CORTICAL MAGNIFICATION variable partitioned by the 
#		set of unique values in the label.  Both the measured variable and
#		labels are a linear representation of a 2D map.
#
#
#
#	Input:
#		xID (e.g., relative frequency of stimulation)
#		xValue (e.g., cortical magnification)
#		titleText, xLabelText, xValueText,
#		xIDSubset 
#		
#
#	Output:
#		ggplot type histograms.
#
######################################################################################################
#####################################################################################################

CorticalMagCDF.UseDep = function ( xID, xValue1, xValue2, titleText, subText, xLabelText, yLabelText, xIDSubset=NULL, newX11Flag ) {

		#	Determine the unique labels in xLabel
	xIDList = sort ( unique ( xID ) );
	xIDListLength = length ( xIDList );

		#	Package up the data.  Make use of subset if specified.
	breaks = seq ( 0, max ( xValue1, xValue2 ), 1 );

		#	Gather some stats to make sure visual and numerical results conform.
	dStats = matrix ( 0, nrow=xIDListLength, ncol=7);

	xIDSubsetLength = length ( xIDSubset );
	cdfs1 = matrix ( 0, nrow = length(breaks)-1, ncol = xIDSubsetLength );
	cdfs2 = matrix ( 0, nrow = length(breaks)-1, ncol = xIDSubsetLength );

	for ( i in 1:xIDSubsetLength ) {

		if ( i == xIDSubsetLength ) {
			tmp = xValue1 [ which ( xID > xIDSubset[i] ) ];
		} else {
			tmp = xValue1 [ which ( xID == xIDSubset[i] ) ];
		}

		tmp = tmp [ tmp != 0 ];
		dStats[i,1] = xIDSubset[i]; dStats[i,2] = length ( tmp );
		dStats[i,3] = mean ( tmp ); dStats[i,4] = sqrt ( var ( tmp ) );
		tmp = hist ( tmp, plot=FALSE, breaks=breaks );
		cdfs1[,i] = cumsum ( tmp$counts / sum ( tmp$counts ) );

		if ( i == xIDSubsetLength ) {
			tmp = xValue2 [ which ( xID > xIDSubset[i] ) ];
		} else {
			tmp = xValue2 [ which ( xID == xIDSubset[i] ) ];
		}
		tmp = tmp [ tmp != 0 ];
		dStats[i,5] = length ( tmp ); dStats[i,6] = mean ( tmp ); dStats[i,7] = sqrt ( var ( tmp ) );
		tmp = hist ( tmp, plot=FALSE, breaks=breaks );
		cdfs2[,i] = cumsum ( tmp$counts / sum ( tmp$counts ) );
		
	} # for ( i in 1:xIDSubsetLength )
		
	if ( newX11Flag ) { x11(); }

	plot ( tmp$mids, cdfs1[,1], ylim=c(0,1), type="l", col=1, lwd=1, lty=2, main=titleText, sub=subText, xlab=xLabelText, ylab="CDF" );
	lines ( tmp$mids, cdfs2[,1], col=1, lwd=2, lty=1 );

	k=2;
	for ( i in 2:xIDSubsetLength ) {
		lines ( tmp$mids, cdfs1[,i], col=k, lwd=1, lty=2 );
		lines ( tmp$mids, cdfs2[,i], col=k, lwd=2, lty=1 );
		k = k + 1;
	} # for ( i in 2:xIDSubsetLength )	

	abline(h=0.5, col=5, lty=4);

	return ( dStats );

} # CorticalMagCDF.UseDep = function ( xLabel, xValue, titleText, xLabelText,... ) {


######################################################################################################
######################################################################################################
#
#	InvCorticalMagCDF 
#
#		Generate histograms of some measured variable partitioned by the 
#		set of unique values in the label.  Both the measured variable and
#		labels are a linear representation of a 2D map.
#
#
#
#	Input:
#		xID (e.g., relative frequency of stimulation)
#		xValue (e.g., cortical magnification)
#		titleText, xLabelText, xValueText,
#		xIDSubset 
#		
#
#	Output:
#		ggplot type histograms.
#
######################################################################################################
#####################################################################################################

InvCorticalMagCDF.UseDep = function ( xID, cMagMap.ref, rfData.ref, cMagMap.exp, rfData.exp,
				titleText, subText, xLabelText, yLabelText, xIDSubset=NULL, newX11Flag ) {

		#	Determine the unique labels in xLabel

	xIDList = sort ( unique ( xID ) );
	xIDListLength = length ( xIDList );

		#	Package up the data.  Make use of subset if specified.
	breaks = seq ( 0, max ( rfData.ref, rfData.exp ), 1 );

		#	Gather some stats to make sure visual and numerical results conform.
	dStats = matrix ( 0, nrow=xIDListLength, ncol=7);

	xIDSubsetLength = length ( xIDSubset );
	cdf.ref = matrix ( 0, nrow = length(breaks)-1, ncol = xIDSubsetLength );
	cdf.exp = matrix ( 0, nrow = length(breaks)-1, ncol = xIDSubsetLength );

	for ( i in 1:xIDSubsetLength ) {

		if ( i == xIDSubsetLength ) {
			tmp = rfData.ref [ MapMagToCellList ( cMagMap.ref, which ( xID > xIDSubset[i] ) ) ];
		} else {
			tmp = rfData.ref [ MapMagToCellList ( cMagMap.ref, which ( xID == xIDSubset[i] ) ) ];
		}

		dStats[i,1] = xIDSubset[i]; dStats[i,2] = length ( tmp );
		dStats[i,3] = mean ( tmp ); dStats[i,4] = sqrt ( var ( tmp ) );
		tmp = hist ( tmp [ tmp != 0 ], plot=FALSE, breaks=breaks );
		cdf.ref[,i] = cumsum ( tmp$counts / sum ( tmp$counts ) );

		if ( i == xIDSubsetLength ) {
			tmp = rfData.exp [ MapMagToCellList ( cMagMap.exp, which ( xID > xIDSubset[i] ) ) ];
		} else {
			tmp = rfData.exp [ MapMagToCellList ( cMagMap.exp, which ( xID == xIDSubset[i] ) ) ];
		}

		dStats[i,5] = length ( tmp ); dStats[i,6] = mean ( tmp ); dStats[i,7] = sqrt ( var ( tmp ) );
		tmp = hist ( tmp [ tmp != 0 ], plot=FALSE, breaks=breaks );
		cdf.exp[,i] = cumsum ( tmp$counts / sum ( tmp$counts ) );
		
	} # for ( i in 1:xIDSubsetLength )
		
	if ( newX11Flag ) { x11(); }

	plot ( tmp$mids, cdf.ref[,1], ylim=c(0,1), type="l", col=1, lwd=1, lty=2, main=titleText, sub=subText, xlab=xLabelText, ylab="CDF" );
	lines ( tmp$mids, cdf.exp[,1], col=1, lwd=2, lty=1 );

	k=2;
	for ( i in 2:xIDSubsetLength ) {
		lines ( tmp$mids, cdf.ref[,i], col=k, lwd=1, lty=2 );
		lines ( tmp$mids, cdf.exp[,i], col=k, lwd=2, lty=1 );
		k = k + 1;
	} # for ( i in 2:xIDSubsetLength )	

	abline(h=0.5, col=5, lty=4);

	return ( dStats );

} # InvCorticalMagCDF.UseDep = function ( xLabel, xValue, titleText, xLabelText,... ) {

