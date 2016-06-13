######################################################################################################
######################################################################################################
#
#	MannWhitneyTest.R
#
#	Do some scatterplots: RF extent vs response magnitude.
#	Do some plots to look at cortical magnification, inverse magnification.#
#
######################################################################################################
######################################################################################################

MannWhitneyTest = function ( x, y, pLevel ) {

	mwRes = list();
	pValues = matrix ( 0, 1, ncol=3 );
	mwRes[[1]] = wilcox.test ( x, y, exact = FALSE, correct=FALSE );
			pValues[1] = 1; pValues[2] = 1; pValues[3] = mwRes[[1]]$p.value;
	acceptAltHyp = pValues[,3] < pLevel;

	return ( list ( mwRes=mwRes, pValues=pValues, acceptAltHyp=acceptAltHyp ) );
		
} # MannWhitneyTest= function ( ) {

######################################################################################################
######################################################################################################
#
#	MWTestForCorticalMagnification .R
#
#	Do some scatterplots: RF extent vs response magnitude.
#	Do some plots to look at cortical magnification, inverse magnification.#
#
######################################################################################################
######################################################################################################

MWTestForCorticalMagnification = function ( mapStimCount, rfData1, rfData2, pLevel ) {

		#	Determine number of relative stimulation frequency values.
	stimCountList = sort ( unique ( mapStimCount ) );
	stimCountListLength = length ( stimCountList );
	
	mwRes = list();
	pValues = matrix ( 0, nrow=(stimCountListLength * ( stimCountListLength + 1 ) / 2), ncol=9 );

	k = 1;
	for ( i in 1:stimCountListLength ) {
		iWhich = which ( mapStimCount == stimCountList[i] );
		tmp.i = rfData1[iWhich]; tmp.i = tmp.i [ tmp.i != 0 ];
		if ( length ( tmp.i ) ) {
			for ( j in i:stimCountListLength ) {
				jWhich = which ( mapStimCount == stimCountList[j] );
				tmp.j = rfData2[jWhich]; tmp.j = tmp.j [ tmp.j != 0 ];
				if ( length ( tmp.j ) ) {
					mwRes[[k]] = wilcox.test ( tmp.i, tmp.j, exact = FALSE, correct=FALSE );
					pValues[k,1] = stimCountList[i]; pValues[k,2] = stimCountList[j]; pValues[k,3] = mwRes[[k]]$p.value;
					pValues[k,4] = length ( tmp.i );
					pValues[k,5] = length ( tmp.j );
					pValues[k,6] = mean ( tmp.i ); pValues[k,7] = sqrt ( var ( tmp.i ) );
					pValues[k,8] = mean ( tmp.j ); pValues[k,9] = sqrt ( var ( tmp.j ) );
				} else {
					pValues[k,3] = 1.0; 	# To make sure this never passes significance tests
				} # if ( length ( tmp.j ) )
				k = k + 1;
			} # for ( j in i:numSectors )
		} else {
			for ( j in i:stimCountListLength ) {
				pValues[k,3] = 1.0;	# Just to make sure it'll never pass a stat sig threshold test.
				k = k + 1;
			} # for ( j in i:numSectors )
		} # if ( length ( tmp.i ) )
	} # for ( i in 1:numSectors )
	acceptAltHyp = pValues[,3] < pLevel;

	return ( list ( mwRes=mwRes, pValues=pValues, acceptAltHyp=acceptAltHyp ) );
		
} # MWTestForCorticalMagnification = function ( ) {


######################################################################################################
######################################################################################################
#
#	MWTestForInvCorticalMagnification .R
#
#	Do some scatterplots: RF extent vs response magnitude.
#	Do some plots to look at cortical magnification, inverse magnification.#
#
######################################################################################################
######################################################################################################

MWTestForInvCorticalMagnification = function ( mapStimCount, cMagMap.ref, rfData.ref, cMagMap.exp, rfData.exp, pLevel ) {

		#	Determine number of relative stimulation frequency values.
	stimCountList = sort ( unique ( mapStimCount ) );
	stimCountListLength = length ( stimCountList );
	
	mwRes = list();
	pValues = matrix ( 0, nrow=(stimCountListLength * ( stimCountListLength + 1 ) / 2), ncol=9 );

	k = 1;
	for ( i in 1:stimCountListLength ) {

		iWhich = which ( mapStimCount == stimCountList[i] );
		tmp.ref = rfData.ref [ MapMagToCellList ( cMagMap.ref, iWhich ) ];

		if ( length ( tmp.ref ) ) {
			for ( j in i:stimCountListLength ) {

				jWhich = which ( mapStimCount == stimCountList[j] );
				tmp.exp = rfData.exp [ MapMagToCellList ( cMagMap.exp, jWhich ) ];
				if ( length ( tmp.exp ) ) {
					mwRes[[k]] = wilcox.test ( tmp.ref, tmp.exp, exact = FALSE, correct=FALSE );
					pValues[k,1] = stimCountList[i];
					pValues[k,2] = stimCountList[j];
					pValues[k,3] = mwRes[[k]]$p.value;
					pValues[k,4] = length ( iWhich );
					pValues[k,5] = length ( jWhich );
					pValues[k,6] = mean ( tmp.ref ); pValues[k,7] = sqrt ( var ( tmp.ref ) );
					pValues[k,8] = mean ( tmp.exp ); pValues[k,9] = sqrt ( var ( tmp.exp ) );
				} else {
					pValues[k,3] = 1.0;	# Just to make sure it'll never pass a stat sig threshold test.
				} # if ( length ( tmp.exp ) )
				k = k + 1;
			} # for ( j in i:numSectors )
		} else {
			for ( j in i:stimCountListLength ) {
				pValues[k,3] = 1.0;	# Just to make sure it'll never pass a stat sig threshold test.
				k = k + 1;
			} # for ( j in i:numSectors )
		} # if ( length ( temp.ref ) )

	} # for ( i in 1:numSectors )

	acceptAltHyp = pValues[,3] < pLevel;

	return ( list ( mwRes=mwRes, pValues=pValues, acceptAltHyp=acceptAltHyp ) );

		
} # D3AllSectorsMannWhitneyTest = function ( ) {

######################################################################################################
######################################################################################################
#
#	D3AllSectorsMannWhitneyTest.R
#
#	Do some scatterplots: RF extent vs response magnitude.
#	Do some plots to look at cortical magnification, inverse magnification.#
#
######################################################################################################
######################################################################################################

D3AllSectorsMannWhitneyTest = function ( cMagMap.ref, rfData.ref, cMagMap.exp, rfData.exp, pLevel ) {

		#	Identify which D3 layer S nodes make up each sector
	sLocs = Dig3SectorLocs1 ( as.integer ( sqrt ( length ( rfData.ref ) ) ) );
	numSectors = dim(sLocs)[2];
	
		#	Preallocate memory for the data to be returned.
	mwRes = list();
	pValues = matrix ( 0, nrow=(numSectors * ( numSectors + 1 ) / 2), ncol=9 );

		#	Cycle through all possible pairs bearing in mind symmetry.
	kmsRes = 1;
	for ( i in 1:numSectors ) {
		
			#	Find the (no duplicates) list of C layer cells whose receptive field includes any of the given sector S layer nodes.
		tmp.ref = rfData.ref [ MapMagToCellList ( cMagMap.ref, sLocs[,i] ) ];
		
		for ( j in i:numSectors ) {
		
			tmp.exp = rfData.exp [ MapMagToCellList ( cMagMap.exp, sLocs[,j] ) ];

			mwRes[[kmsRes ]] = wilcox.test ( tmp.ref, tmp.exp, exact = FALSE, correct=FALSE );
			pValues[kmsRes,1] = i; pValues[kmsRes,2] = j; pValues[kmsRes,3] = mwRes[[kmsRes]]$p.value;
			pValues[kmsRes,4] = length ( tmp.ref );
			pValues[kmsRes,5] = length ( tmp.exp );
			pValues[kmsRes,6] = mean ( tmp.ref ); pValues[kmsRes,7] = sqrt ( var ( tmp.ref ) );
			pValues[kmsRes,8] = mean ( tmp.exp ); pValues[kmsRes,9] = sqrt ( var ( tmp.exp ) );
			kmsRes = kmsRes  + 1;
		} # for ( j in i:numSectors )

	} # for ( i in 1:numSectors )
	acceptAltHyp = pValues[,3] < pLevel;

	return ( list ( mwRes=mwRes, pValues=pValues, acceptAltHyp=acceptAltHyp ) );
		
} # D3AllSectorsMannWhitneyTest = function ( ) {

######################################################################################################
######################################################################################################
#
#	D3AllLongAxesMannWhitneyTest.R
#
#
######################################################################################################
######################################################################################################

D3AllLongAxesMannWhitneyTest = function ( cMagMap.ref, rfData.ref, cMagMap.exp, rfData.exp, pLevel ) {

		#	Identify the number of longitudonal axes
	numLongAxes = as.integer ( sqrt ( length ( rfData.ref ) ) );
	digitWidth = numLongAxes / 3;
	
		#	Preallocate memory for the data to be returned.
	mwRes = list();
	pValues = matrix ( 0, nrow=(numLongAxes * ( numLongAxes + 1 ) / 2), ncol=9 );

		#	Cycle through all possible pairs bearing in mind symmetry.
	kmsRes = 1;
	for ( i in 1:numLongAxes ) {
		
		iDigit = as.integer ( (i-1)/digitWidth ) + 1;
		iLongAxis = i - as.integer ( (i-1)/digitWidth ) * digitWidth;
		tmp.ref = rfData.ref [ MapMagToCellList ( cMagMap.ref,  D3LongAxisLocs1 ( N, iDigit, iLongAxis ) ) ];
		
		for ( j in i:numLongAxes ) {
		
			jDigit = as.integer ( (j-1)/digitWidth ) + 1;
			jLongAxis = j - as.integer ( (j-1)/digitWidth ) * digitWidth;
			tmp.exp = rfData.exp [ MapMagToCellList ( cMagMap.exp,  D3LongAxisLocs1 ( N, jDigit, jLongAxis ) ) ];

			mwRes[[kmsRes ]] = wilcox.test ( tmp.ref, tmp.exp, exact = FALSE, correct=FALSE );
			pValues[kmsRes,1] = i; pValues[kmsRes,2] = j; pValues[kmsRes,3] = mwRes[[kmsRes]]$p.value;
			pValues[kmsRes,4] = length ( tmp.ref );
			pValues[kmsRes,5] = length ( tmp.exp );
			pValues[kmsRes,6] = mean ( tmp.ref ); pValues[kmsRes,7] = sqrt ( var ( tmp.ref ) );
			pValues[kmsRes,8] = mean ( tmp.exp ); pValues[kmsRes,9] = sqrt ( var ( tmp.exp ) );
			kmsRes = kmsRes  + 1;
		} # for ( j in i:numSectors )

	} # for ( i in 1:numSectors )
	acceptAltHyp = pValues[,3] < pLevel;

	return ( list ( mwRes=mwRes, pValues=pValues, acceptAltHyp=acceptAltHyp ) );
		
} # D3AllLongAxesMannWhitneyTest = function ( ) {


######################################################################################################
######################################################################################################
#
#	D3AllLongAxesMannWhitneyTest1.R
#
#	Difference from previous is to apply a test of stimulation count in order to include.
#
#
######################################################################################################
######################################################################################################

D3AllLongAxesMannWhitneyTest1 = function ( mapStimCount.ref, cutOff.ref, cMagMap.ref, rfData.ref, 
								mapStimCount.exp, cutOff.exp, cMagMap.exp, rfData.exp, pLevel, trimEnds ) {

		#	Identify the number of longitudonal axes
	numLongAxes = as.integer ( sqrt ( length ( rfData.ref ) ) );
	digitWidth = numLongAxes / 3;
	
		#	Preallocate memory for the data to be returned.
	mwRes = list();
	pValues = matrix ( 0, nrow=(numLongAxes * ( numLongAxes + 1 ) / 2), ncol=9 );

		#	Cycle through all possible pairs bearing in mind symmetry.
	kmsRes = 1;
	for ( i in 1:numLongAxes ) {
		
		iDigit = as.integer ( (i-1)/digitWidth ) + 1;
		iLongAxis = i - as.integer ( (i-1)/digitWidth ) * digitWidth;
		tmp.ref = rfData.ref [ MapMagToCellList1 ( cMagMap.ref,  D3LongAxisLocs2 ( N, iDigit, iLongAxis, trimEnds ),mapStimCount.ref, cutOff.ref ) ];
 
		if ( length ( tmp.ref ) > 2 ) {
		
			for ( j in i:numLongAxes ) {
		
				jDigit = as.integer ( (j-1)/digitWidth ) + 1;
				jLongAxis = j - as.integer ( (j-1)/digitWidth ) * digitWidth;
				tmp.exp = rfData.exp [ MapMagToCellList1 ( cMagMap.exp,  D3LongAxisLocs2 ( N, jDigit, jLongAxis, trimEnds ), mapStimCount.exp, cutOff.exp ) ];

				if ( length ( tmp.exp ) > 2 ) {
					mwRes[[kmsRes ]] = wilcox.test ( tmp.ref, tmp.exp, exact = FALSE, correct=FALSE );
					pValues[kmsRes,1] = i; pValues[kmsRes,2] = j; pValues[kmsRes,3] = mwRes[[kmsRes]]$p.value;
					pValues[kmsRes,4] = length ( tmp.ref );
					pValues[kmsRes,5] = length ( tmp.exp );
					pValues[kmsRes,6] = mean ( tmp.ref ); pValues[kmsRes,7] = sqrt ( var ( tmp.ref ) );
					pValues[kmsRes,8] = mean ( tmp.exp ); pValues[kmsRes,9] = sqrt ( var ( tmp.exp ) );
					kmsRes = kmsRes  + 1;
				} # if ( length ( tmp.exp ) > 2 )
			} # for ( j in i:numSectors )

		} # if ( length ( tmp.ref ) > 2 )

	} # for ( i in 1:numSectors )
	acceptAltHyp = (pValues[,3] < pLevel) & ( pValues[,3] > 0);

	return ( list ( mwRes=mwRes, pValues=pValues, acceptAltHyp=acceptAltHyp ) );
		
} # D3AllLongAxesMannWhitneyTest1 = function ( ) {

######################################################################################################
######################################################################################################
#
#	D3AllLongAxesMannWhitneyTest2.R
#
#	Differences from above: does NOT do the stimulus count check; does trim the ends of a longitudonal 
#	axis (gets rid of edge effects).
#
######################################################################################################
######################################################################################################

D3AllLongAxesMannWhitneyTest2 = function ( cMagMap.ref, rfData.ref, cMagMap.exp, rfData.exp, pLevel, trimEnds ) {

		#	Identify the number of longitudonal axes
	N = numLongAxes = as.integer ( sqrt ( length ( rfData.ref ) ) );
	digitWidth = numLongAxes / 3;
	
		#	Preallocate memory for the data to be returned.
	mwRes = list();
	pValues = matrix ( 0, nrow=(numLongAxes * ( numLongAxes + 1 ) / 2), ncol=9 );

		#	Cycle through all possible pairs bearing in mind symmetry.
	kmsRes = 1;
	for ( i in 1:numLongAxes ) {
		
		iDigit = as.integer ( (i-1)/digitWidth ) + 1;
		iLongAxis = i - as.integer ( (i-1)/digitWidth ) * digitWidth;
		tmp.ref.cellIDs = MapMagToCellList ( cMagMap.ref,  D3LongAxisLocs1 ( N, iDigit, iLongAxis ) );
		tmp.ref.cellIDs = TrimEdgesFromCellList ( tmp.ref.cellIDs, N, trimEnds );
		tmp.ref = rfData.ref [ tmp.ref.cellIDs ];

		
		for ( j in i:numLongAxes ) {
		
			jDigit = as.integer ( (j-1)/digitWidth ) + 1;
			jLongAxis = j - as.integer ( (j-1)/digitWidth ) * digitWidth;
			tmp.exp.cellIDs = MapMagToCellList ( cMagMap.exp,  D3LongAxisLocs1 ( N, jDigit, jLongAxis ) );
			tmp.exp.cellIDs = TrimEdgesFromCellList ( tmp.exp.cellIDs, N, trimEnds );
			tmp.exp = rfData.exp [ tmp.exp.cellIDs];


			mwRes[[kmsRes ]] = wilcox.test ( tmp.ref, tmp.exp, exact = FALSE, correct=FALSE );
			pValues[kmsRes,1] = i; pValues[kmsRes,2] = j; pValues[kmsRes,3] = mwRes[[kmsRes]]$p.value;
			pValues[kmsRes,4] = length ( tmp.ref );
			pValues[kmsRes,5] = length ( tmp.exp );
			pValues[kmsRes,6] = mean ( tmp.ref ); pValues[kmsRes,7] = sqrt ( var ( tmp.ref ) );
			pValues[kmsRes,8] = mean ( tmp.exp ); pValues[kmsRes,9] = sqrt ( var ( tmp.exp ) );
			kmsRes = kmsRes  + 1;
		} # for ( j in i:numSectors )

	} # for ( i in 1:numSectors )
	acceptAltHyp = pValues[,3] < pLevel;

	return ( list ( mwRes=mwRes, pValues=pValues, acceptAltHyp=acceptAltHyp ) );
		
} # D3AllLongAxesMannWhitneyTest2 = function ( ) {




