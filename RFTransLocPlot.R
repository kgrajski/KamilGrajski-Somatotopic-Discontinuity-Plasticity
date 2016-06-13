######################################################################################################
######################################################################################################
#
#	RFTransLocPlot
#
#	Input
#		expPos - N^2 x 2 array of the "experimental" RF centroids
#		ctlPos - N^2 x 2 array of the "control" RF centroids
#
#	Output
#		Contour plot of the RF translocation vector magnitude
#		Contour plot of the RF translocation vector angle (in radians).
#
######################################################################################################
######################################################################################################

RFTransLocPlot = function ( expPos, ctlPos, mainText, xlabText, ylabText, moveMinDist=0.5 ) {

	iTrimRings = 1;
	N = dim(expPos)[1];

	C = expPos - ctlPos;
	C2 = sqrt ( apply ( ( C^2 ), 1, sum ) );
	iMoved = which(C2 >= moveMinDist);
	iDidNotMove = which ( C2 < moveMinDist );
	C2[iDidNotMove] = 0;
	thetaC = rep ( 0, N );
	thetaC[iMoved] = sign ( C[iMoved,2] ) * acos( (C[iMoved,1]) / C2[iMoved] );

	iTrimmed = TrimEdgesFromCellList ( seq(1,N), sqrt(N), iTrimRings );
	ShowVecAsMapContour ( C2, paste(mainText," Displ Vec Mag", sep=""), xlab, ylab );
	abline ( h = c ( 0.5, 11.5, 22.5 ), lty=3,col=605 );
	ShowVecAsMapContour ( thetaC, paste(mainText, " Displ Vec Angle (Radians)", sep=""), xlab, ylab );
	abline ( h = c ( 0.5, 11.5, 22.5 ), lty=3,col=1 );

} # RFTransLocPlot = function ( expPos, ctlPox, mainText, xlabText, ylabText ) {



#	Testing Code

#x11(); par( mfcol=c(2,2) );
#mainText = "E"; mainText="I";
#xlabText="Distal -> Proximal"; ylabText="Digit 1 -> Digit 3";
#expPos = selStim$rfTrackData.e[[iLast]]$rfCentroid;
#ctlPos = selStim.ctl$rfTrackData.e[[iLast]]$rfCentroid;

#expPos = selStim$rfTrackData.i[[iLast]]$rfCentroid;
#ctlPos = selStim.ctl$rfTrackData.i[[iLast]]$rfCentroid;

#expPos = rbind ( c (2, 2), c (0, 2), c (0, 0), c (2, 0) );
#ctlPos = rbind ( c (1, 1), c (1, 1), c (1, 1), c (1, 1) );