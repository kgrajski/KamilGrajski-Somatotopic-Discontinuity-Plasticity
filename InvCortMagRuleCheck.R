######################################################################################################
######################################################################################################
#
#	InvCortMagRuleCheck.R
#
#		Do some calculations and plots to check for inverse cortical
#		magnification rule conformance.
#
#	Input:
#		stimCount - spatial map of relative frequency of stimulation
#		cMap - map produced by QuantCorticalAmp
#				(magnification and RF extent can be derived)
#		xVal - spatial map of some reference quantitative value under investigation
#				(if x == 0 then just work on receptive field extent)
#				(another example might be receptive field maximum response
#
#		Both inputs are linear representations of an N x N map.
#
#	Output:
#		List structure containing variables that identify the partition values,
#		the linear representation of 
#
######################################################################################################
######################################################################################################

InvCortMagRuleCheck= function ( mapStimCount, mapRFData, xVal = 0) {

		#	Compute cortical magnification values.
	cMag = apply ( mapRFData, 2, sum );

		#	Compute receptive field extents.
	mapRFExtent = apply ( mapRFData, 1, sum );

		#	Determine number of relative stimulation frequency values.
	stimCountList = sort ( unique ( mapStimCount ) );
	stimCountListLength = length ( stimCountList );

		#	Determine the highest numerical count over all the frequency values.
	stimCountMax = -1;
	for ( i in 1:stimCountListLength  ) {
		stimCountMax = max ( stimCountMax, sum( mapStimCount == stimCountList[i] ) );
	} # for ( i in stimCountList )
		
		#	Initialize the mean and SD arrays.
	xMean.Mag = rep ( 0, stimCountListLength  );
	xSD.Mag =  rep ( 0, stimCountListLength  );

		#	For each of the nodes of a given stimulation frequency value
		#	find its magnification.  Report the mean and SD.
	for ( i in 1:stimCountListLength  ) {
		xWhich = which( mapStimCount == stimCountList[i] );
		numWhich = length ( xWhich );

		xTmp = cMag[xWhich];
		xMean.Mag[i] = mean( xTmp[ xTmp>0 ] );
		xSD.Mag[i] = sqrt ( var ( xTmp[ xTmp>0 ] ) );
	} # for ( i in stimCountList )

		#	For each of the nodes of a given stimulation frequency value
		#	determine the list of cortical cells that contribute to its
		#	magnification.  Report the mean and SD of the receptive field
		#	extents of those cortical cells.
	xMean.rfExtentByRelStim = rep ( 0, stimCountListLength  );
	xSD.rfExtentByRelStim = rep ( 0, stimCountListLength  );
	for ( k in 1:stimCountListLength  ) {

		xWhich = which( mapStimCount == stimCountList[k] );
		numWhich = length ( xWhich );
		xMean.tmp = rep ( 0, numWhich );
		xSD.tmp = rep ( 0, numWhich );
		xMin.tmp = rep ( 0, numWhich );
		xMax.tmp = rep ( 0, numWhich );
		for ( j in 1:numWhich ) {
			whichCorticalCells = which(mapRFData[,xWhich[j]]);
			xMean.tmp[j] = mean ( mapRFExtent[whichCorticalCells], na.rm=TRUE );
			xSD.tmp[j] = sqrt ( var ( mapRFExtent[whichCorticalCells], na.rm=TRUE ) );
			xMin.tmp[j] = min ( mapRFExtent[whichCorticalCells] );
			xMax.tmp[j] = max ( mapRFExtent[whichCorticalCells] );
		} # for ( j in 1:numWhich )

		xMean.rfExtentByRelStim[k] = mean ( xMean.tmp, na.rm=TRUE );
		xSD.rfExtentByRelStim[k] = mean ( xSD.tmp, na.rm=TRUE );

	} # for ( k in 1:stimCountListLength  )

		#	Return the list of results.

	rtmp = list ( stimCountListLength=stimCountListLength, stimCountList=stimCountList,
				xMean.Mag=xMean.Mag, xSD.Mag=xSD.Mag,
				xMean.rfExtentByRelStim=xMean.rfExtentByRelStim,
				xSD.rfExtentByRelStim=xSD.rfExtentByRelStim );

	return ( rtmp );

} # InvCortMagRuleCheck= function ( ) {