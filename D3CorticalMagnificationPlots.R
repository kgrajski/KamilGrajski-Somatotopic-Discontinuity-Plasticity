
######################################################################################################
######################################################################################################
#
#	D3CorticalMagnificationPlots: Do a variety of plots to visualize cortical magnification.
#		Input:
#			rfAreas - linear representation of an N x N grid of RF extents.
#
#		Output:
#			Various and sundry plots.
#			Scatterplots.
#			Magnification by Sector.
#			Magnification by Tracks parallel to longitudonal axis of symmetry.
#
######################################################################################################
######################################################################################################

D3CorticalMagnificationPlots = function ( rfAreas, titleText ) {

		#
		# 	Visualizations in the exploration of inverse cortical magnification rule.
		# 	BY SECTOR: For each of the 9 sectors, compute its histogram of receptive field extents
		#

		#	Set histogram parameters.
	stepSize = 5;
	xmaxVal = (as.integer(max(rfAreas)/stepSize) + 1) * stepSize; 

		#	Generate the histograms of RF extent for each sector.
	hSect = list();
	sLocs = Dig3SectorLocs1 ( N );
	for ( i in 1:9 ) {
		hSect[[i]] = hist ( rfAreas[sLocs[,i]], breaks=seq(from=0, to=xmaxVal, by=stepSize), plot=FALSE );
	} # for ( i in 1:9 )

		#	Determine plot parameters so that all histogram results can be overlayed on one plot.

	ymaxVal = -1;
	for ( i in 1:9 ) {
		ymaxVal = max ( ymaxVal, hSect[[i]]$counts )
	} # for ( i in 1:9 )
	ymaxVal = (as.integer(ymaxVal/stepSize) + 1) * stepSize;
	ylim = c ( 0, ymaxVal );
	xlim = c ( 0, xmaxVal );

		#	Plot the histograms for each sector.
	x11();
	par(mfrow=c(2,1))
	plot( hSect[[1]]$mids, hSect[[1]]$counts, type="l",   lty=3, col=1, xlim = xlim, ylim=ylim, 
		main=paste("Receptive Field Extent Histogram", titleText, sep="\n"),
		xlab="Receptive Field Extent", ylab="Count",
		sub="Digits: 3(Solid), 2(Dash), 1(Dot); Sector: Distal(Blk), Mid(Red), Prox(Green))");
	lines ( hSect[[2]]$mids, hSect[[2]]$counts, type="l", lty=2, col=1 );
	lines ( hSect[[3]]$mids, hSect[[3]]$counts, type="l", lty=1, col=1 );

	lines ( hSect[[4]]$mids, hSect[[4]]$counts, type="l", lty=3, col=2 );
	lines ( hSect[[5]]$mids, hSect[[5]]$counts, type="l", lty=2, col=2 );
	lines ( hSect[[6]]$mids, hSect[[6]]$counts, type="l",   lty=1, col=2 );

	lines ( hSect[[7]]$mids, hSect[[7]]$counts, type="l", lty=3, col=3 );
	lines ( hSect[[8]]$mids, hSect[[8]]$counts, type="l", lty=2, col=3 );
	lines ( hSect[[9]]$mids, hSect[[9]]$counts, type="l", lty=1, col=3 );

		#	Plot the average across digits by sector.
	tmp.distal = apply ( rbind ( hSect[[1]]$counts, hSect[[2]]$counts, hSect[[3]]$counts ), 2, mean );
	tmp.mid = apply ( rbind ( hSect[[4]]$counts, hSect[[5]]$counts, hSect[[6]]$counts ), 2, mean );
	tmp.prox = apply ( rbind ( hSect[[7]]$counts, hSect[[8]]$counts, hSect[[9]]$counts ), 2, mean );

	plot( hSect[[1]]$mids, tmp.distal, type="l",   lty=1, col=1, xlim = xlim, ylim=ylim, 
		main=paste("Receptive Field Extent Histogram", "Averaged by Sector Type (Dist, Mid, Prox)", titleText, sep="\n"),
		xlab="Receptive Field Extent", ylab="Count",
		sub="Sector: Distal(Blk), Mid(Red), Prox(Green))");
	lines ( hSect[[6]]$mids, tmp.mid, type="l",   lty=1, col=2 );
	lines ( hSect[[9]]$mids, tmp.prox, type="l", lty=1, col=3 );

} # D3CorticalMagnificationPlots = function ( rfAreas, titleText ) {

######################################################################################################
######################################################################################################
#
#	D3PairedCorticalMagnificationPlots: Do a variety of plots to visualize cortical magnification.
#		Input:
#			rfAreas - the reference data: linear representation of an N x N grid of RF extents.
#			rfAreas.exp - the data to be compared with the reference data.
#
#		Output:
#			Various and sundry plots.
#			Scatterplots.
#			Magnification by Sector.
#			Magnification by Tracks parallel to longitudonal axis of symmetry.
#
######################################################################################################
######################################################################################################

D3PairedCorticalMagnificationPlots = function ( rfAreas, rfAreas.exp, titleText, titleText.exp ) {

		#
		# 	Visualizations in the exploration of inverse cortical magnification rule.
		# 	BY SECTOR: For each of the 9 sectors, compute its histogram of receptive field extents
		#
	numDigits = 3;
	numSectors = 3;

		#	Set histogram parameters.
	stepSize = 5;
	xmaxVal = (as.integer(max(rfAreas)/stepSize) + 1) * stepSize;
	xmaxVal = max ( xmaxVal,  (as.integer(max(rfAreas.exp)/stepSize) + 1) * stepSize );

		#	Generate the histograms of RF extent for each sector.
	hSect = list();
	hSect.exp = list();
	sLocs = Dig3SectorLocs1 ( N );
	for ( i in 1:9 ) {
		hSect[[i]] = hist ( rfAreas[sLocs[,i]], breaks=seq(from=0, to=xmaxVal, by=stepSize), plot=FALSE );
		hSect.exp[[i]] = hist ( rfAreas.exp[sLocs[,i]], breaks=seq(from=0, to=xmaxVal, by=stepSize), plot=FALSE );
	} # for ( i in 1:9 )

		#	Determine plot parameters so that all histogram results can be overlayed on one plot.

	ymaxVal = -1;
	for ( i in 1:9 ) {
		ymaxVal = max ( ymaxVal, hSect[[i]]$counts, hSect.exp[[i]]$counts )
	} # for ( i in 1:9 )
	ymaxVal = (as.integer(ymaxVal/stepSize) + 1) * stepSize;
	ylim = c ( 0, ymaxVal );
	xlim = c ( 0, xmaxVal );

		#	Plot the histograms for each sector.
	x11(); par(mfrow=c(3,3));
	for ( iDigit in seq ( numDigits, 1, -1 ) ) {
		iLoc = iDigit;
		for ( iSector in seq ( 3, 1, -1 ) ) {
			plot( hSect[[iLoc]]$mids, hSect[[iLoc]]$counts, type="l", lty=3, col=1, xlim = xlim, ylim=ylim, 
					main=paste ( paste(titleText, " : ", titleText.exp, sep=""), 
					paste ("Digit", iDigit, "Sector", iSector, sep=" " ), sep="\n" ),
					xlab="Receptive Field Extent", ylab="Histogram Count",
					sub="Ref: (Dotted); Exp: (Solid)" );
			lines ( hSect.exp[[iLoc]]$mids, hSect.exp[[iLoc]]$counts, lty=1, col=1 );
			iLoc = iLoc + numSectors;
		} # for ( iSector in seq ( numSectors, 1, -1 ) )
	} # for ( iDigit in seq ( numDigits, 1, -1 ) )

		#	Plot the average across digits by sector.
	tmp.distal = apply ( rbind ( hSect[[1]]$counts, hSect[[2]]$counts, hSect[[3]]$counts ), 2, mean );
	tmp.mid = apply ( rbind ( hSect[[4]]$counts, hSect[[5]]$counts, hSect[[6]]$counts ), 2, mean );
	tmp.prox = apply ( rbind ( hSect[[7]]$counts, hSect[[8]]$counts, hSect[[9]]$counts ), 2, mean );

	tmp.distal.exp = apply ( rbind ( hSect.exp[[1]]$counts, hSect.exp[[2]]$counts, hSect.exp[[3]]$counts ), 2, mean );
	tmp.mid.exp = apply ( rbind ( hSect.exp[[4]]$counts, hSect.exp[[5]]$counts, hSect.exp[[6]]$counts ), 2, mean );
	tmp.prox.exp = apply ( rbind ( hSect.exp[[7]]$counts, hSect.exp[[8]]$counts, hSect.exp[[9]]$counts ), 2, mean );
	
	ymaxVal = max ( tmp.distal, tmp.distal.exp, tmp.mid, tmp.mid.exp, tmp.prox, tmp.prox.exp );
	ymaxVal = (as.integer(ymaxVal/stepSize) + 1) * stepSize;
	ylim = c ( 0, ymaxVal );

	x11(); par(mfrow=c(2,1));
	plot( hSect[[1]]$mids, tmp.distal, type="l",   lty=3, col=1, xlim = xlim, ylim=ylim, 
		main=paste ( "Averaged By Sector Histograms", paste(titleText, " : ", titleText.exp, sep=""), sep="\n" ),
					xlab="Receptive Field Extent", ylab="Histogram Count",
					sub="Ref: (Dotted); Exp: (Solid)" );
	lines ( hSect[[1]]$mids, tmp.mid, type="l",   lty=3, col=2 );
	lines ( hSect[[1]]$mids, tmp.prox, type="l", lty=3, col=3 );
	lines ( hSect[[1]]$mids, tmp.distal.exp, type="l",  lty=1, col=1 );
	lines ( hSect[[1]]$mids, tmp.mid.exp, type="l", lty=1, col=2 );
	lines ( hSect[[1]]$mids, tmp.prox.exp, type="l", lty=1, col=3 );

} # D3PairedCorticalMagnificationPlots = function ( rfAreas, titleText ) {

######################################################################################################
######################################################################################################
#
#	D3PairLongAxisRFExtentDistr:
#
#		Do a boxplot for a single pair of longitudonal axes where one is .ref and one is .exp.
#
######################################################################################################
######################################################################################################

D3PairLongAxisRFExtentDistr= function ( rfAreas.ref, iDigit.ref, iSector.ref,
							rfAreas.exp, iDigit.exp, iSector.exp, titleText.ref, titleText.exp ) {

	x11();

	xMat = cbind ( rfAreas.ref[ D3LongAxisLocs1(N, iDigit.ref, iSector.ref) ],
					rfAreas.exp[ D3LongAxisLocs1(N, iDigit.exp, iSector.exp) ] );
	ylim = c ( 0, max ( xMat ) ) * 1.1; 	#	Add 10% buffer
	xMat = as.data.frame ( xMat );
	names(xMat)[1] = paste ( "Ref: Dig",iDigit.ref,"Track",iSector.ref );
	names(xMat)[2] = paste ( "Exp: Dig",iDigit.exp,"Track",iSector.exp );

	boxplot ( xMat, ylim=ylim, main=paste( "Distr. RF Extent in Longitudonal Axes", paste(titleText, " : ", titleText.exp, sep=""), sep="\n" ),
				ylab="Receptive Field Extent" );

} # D3PairLongAxisRFExtentDistr= function ( rfAreas, titleText ) {


######################################################################################################
######################################################################################################
#
#	D3SingleDigitLongAxisRFExtentDistr:
#
#		Do a boxplot for all the longitudonal axes on a single digit.
#
######################################################################################################
######################################################################################################

D3SingleDigitLongAxisRFExtentDistr= function ( rfAreas, xDigit, titleText, newX11Flag ) {

	N = as.integer ( sqrt ( length ( rfAreas ) ) );
	numDigits = 3;
	numTracks = as.integer ( N / 3 );
	
	xMat = matrix ( 0, nrow = N, ncol = numTracks );
	for ( i in 1:numTracks ) {
		xMat[,i] = rfAreas[ D3LongAxisLocs1(N, xDigit, i) ];
	} # for ( i in 1:numTracks ) {

	ylim = c ( 0, max(xMat) ) * 1.1;	#	Add 10% buffer. 
	xMat = as.data.frame(xMat);
	for ( i in 1:numTracks ) {
		names(xMat)[i] = paste ( "D", xDigit, "Trk", i );
	} # for ( i in 1:numTracks ) {

	if ( newX11Flag ) { x11(); }
	boxplot ( as.data.frame(xMat), main=paste("Distributions RF Extent in Longitudonal Axes",titleText,sep="\n"),
			ylab="Receptive Field Extent" );

} # D3SingleDigitLongAxisRFExtentDistr= function ( rfAreas, titleText ) {

######################################################################################################
######################################################################################################
#
#	D3AllSectorsRFExtentDistr:
#
#		Do a boxplot for each of the sectors in the model.
#
######################################################################################################
######################################################################################################

D3AllSectorsRFExtentDistr = function ( rfAreas, titleText, newX11Flag ) {

	sLocs = Dig3SectorLocs1 ( as.integer ( sqrt ( length ( rfAreas ) ) ) );
	numSectors = dim(sLocs)[2];
	xMat = matrix ( 0, nrow = dim(sLocs)[1], ncol = numSectors );
	for ( i in 1:numSectors ) {
		xMat[,i] = rfAreas[ sLocs[,i] ];
	} # for ( i in 1:numSectors ) {

	ylim = c ( 0, max(xMat) ) * 1.1;	#	Add 10% buffer
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
	boxplot ( as.data.frame(xMat), ylim=ylim, main=paste( "Distributions RF Extent in Sectors", titleText, sep="\n" ),
			ylab="Receptive Field Extent" );

} # D3AllSectorsRFExtentDistr= function ( rfAreas, titleText ) {

######################################################################################################
######################################################################################################
#
#	D3PairSectorsRFExtentDistr:
#
#		Do a boxplot for a pair of sectors in the model.
#
######################################################################################################
######################################################################################################

D3PairSectorsRFExtentDistr = function ( rfAreas.ref, iSector.ref, rfAreas.exp, iSector.exp, titleText.ref, titleText.exp, newX11Flag ) {

	numSectors = 2;

	sLocs = Dig3SectorLocs1 ( as.integer ( sqrt ( length ( rfAreas.ref ) ) ) );

	xMat = matrix ( 0, nrow = dim(sLocs)[1], ncol = numSectors );
	xMat[,1] = rfAreas.ref[ sLocs[,iSector.ref] ];
	xMat[,2] = rfAreas.exp[ sLocs[,iSector.exp] ];
	ylim = c ( 0, max ( xMat ) ) * 1.1;	# Add 10%
	xMat = as.data.frame(xMat);

	names(xMat)[1] = "Ref Sector";
	names(xMat)[2] = "Exp Sector";

	if ( newX11Flag ) { x11(); }
	boxplot ( as.data.frame(xMat), ylim=ylim,
			main=paste( "Distributions RF Extent by Sector", paste(titleText.ref,"v",titleText.exp,sep=" "), sep="\n" ),
			ylab="Receptive Field Extent" );

} # D3PairSectorsRFExtentDistr= function ( rfAreas, titleText ) {

