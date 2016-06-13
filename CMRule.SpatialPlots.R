######################################################################################################
######################################################################################################
#
#	CMRule.SpatialPlots.R
#
#		Input:
#			ctl.CMRaw.e - E cell control Boolean N2 x N2 matrix from which CM data is derived
#			exp.CMRaw.e - E cell experimental Boolean N2 x N2 matrix from which CM data is derived
#			ctl.CMRaw.i - I cell control Boolean N2 x N2 matrix from which CM data is derived
#			exp.CMRaw.i - I cell experimental Boolean N2 x N2 matrix from which CM data is derived
#
#			ctl.RFArea.e - E cell control RFArea for each of the N2 cells in the cortical network
#			exp.RFArea.e - E cell experimental RFArea for each of the N2 cells in the cortical network
#			ctl.RFArea.i - I cell control RFArea for each of the N2 cells in the cortical network
#			exp.RFArea.i - I cell experimental RFArea for each of the N2 cells in the cortical network
#
#			iTrimRings - remove the outer-most iTrimRings when displaying data.			
#
######################################################################################################
######################################################################################################

CMPointByPoint = function ( ctl.CMRaw.e, exp.CMRaw.e, ctl.CMRaw.i, exp.CMRaw.i,
					ctl.RFArea.e, exp.RFArea.e, ctl.RFArea.i, exp.RFArea.i,
					iTrimRings, delta.CM.min, delta.MeanRFArea.min, mainTextIn ) {

	N2 = ( dim ( ctl.CMRaw.e )[1] );
	N = as.integer ( sqrt ( N2 ) );

	ctl.CM.e = rep ( 0, N2 );
	ctl.CM.i = rep ( 0, N2 );
	exp.CM.e = rep ( 0, N2 );
	exp.CM.i = rep ( 0, N2 );

	ctl.MeanRFArea.e = rep ( 0, N2 );
	ctl.MeanRFArea.i = rep ( 0, N2 );

	exp.MeanRFArea.e = rep ( 0, N2 );
	exp.MeanRFArea.i = rep ( 0, N2 );

	whichHaveCM.e = rep ( FALSE, N2 );
	whichHaveCM.i = rep ( FALSE, N2 );

	for ( i in 1:N2 ) {

		tmp1.List = (MapMagToCellList( ctl.CMRaw.e, i ));
		tmp2.List = (MapMagToCellList( exp.CMRaw.e, i ));
		tmp1 = length(tmp1.List);
		tmp2 = length(tmp2.List);
		if ( length(tmp1)>0 & length(tmp2)>0 ) {
			ctl.CM.e[i] = tmp1;
			exp.CM.e[i] = tmp2;

			ctl.MeanRFArea.e[i] = median(ctl.RFArea.e[TrimEdgesFromCellList ( tmp1.List, N, iTrimRings )]);
			exp.MeanRFArea.e[i] = median(exp.RFArea.e[TrimEdgesFromCellList ( tmp2.List, N, iTrimRings )]);

			whichHaveCM.e[i] = TRUE;
		} # 	if ( length(tmp1)>0 & length(tmp2)>0 ) {

		tmp1.List = (MapMagToCellList( ctl.CMRaw.i, i ));
		tmp2.List = (MapMagToCellList( exp.CMRaw.i, i ));
		tmp1 = length(tmp1.List);
		tmp2 = length(tmp2.List);
		if ( length(tmp1)>0 & length(tmp2)>0 ) {
			ctl.CM.i[i] = tmp1;
			exp.CM.i[i] = tmp2;

			ctl.MeanRFArea.i[i] = median(ctl.RFArea.i[TrimEdgesFromCellList ( tmp1.List, N, iTrimRings )]);
			exp.MeanRFArea.i[i] = median(exp.RFArea.i[TrimEdgesFromCellList ( tmp2.List, N, iTrimRings )]);

			whichHaveCM.i[i] = TRUE;
		} # 	if ( length(tmp1)>0 & length(tmp2)>0 ) {

	} # for ( i in 1:N2 )

	toShow.e = TrimEdgesFromCellList ( which ( whichHaveCM.e == TRUE ), N, iTrimRings );
	toShow.i = TrimEdgesFromCellList ( which ( whichHaveCM.i == TRUE ), N, iTrimRings );

		#	Scatterplot: CM
	x11(); par(mfrow = c(2,2));
	mainText = paste ( "Per Input Node CM", mainTextIn, "Cortical E-Cells",sep="\n" );
	plot ( ctl.CM.e[toShow.e], exp.CM.e[toShow.e], type="p", pch=".", col=1, cex=4,
		main=mainText, xlab="Initial Baseline CM", ylab="Final Focal Stim CM" );
	lines ( c(0,max(ctl.CM.e[toShow.e])), c(0,max(ctl.CM.e[toShow.e])) );

	mainText = paste ( "Per Input Node CM", mainTextIn, "Cortical I-Cells",sep="\n" );
	plot ( ctl.CM.i[toShow.i], exp.CM.i[toShow.i], type="p", pch=".", col=1, cex=4,
		main=mainText, xlab="Initial Baseline CM", ylab="Final Focal Stim CM" );
	lines ( c(0,max(ctl.CM.i[toShow.i])), c(0,max(ctl.CM.i[toShow.i])) );

		#	Scatterplot: RFArea
	mainText = paste ( "Mean RF Area Per CM Per Input Node", mainTextIn, "Cortical E-Cells",sep="\n" );
	plot ( ctl.MeanRFArea.e[toShow.e], exp.MeanRFArea.e[toShow.e], type="p", pch=".", col=1, cex=4,
		main=mainText, xlab="Initial Baseline CM", ylab="Final Focal Stim CM" );
	lines ( c(0,max(ctl.MeanRFArea.e[toShow.e])), c(0,max(ctl.MeanRFArea.e[toShow.e])) );

	mainText = paste ( "Mean RF Area Per CM Per Input Node", mainTextIn, "Cortical I-Cells",sep="\n" );
	plot ( ctl.MeanRFArea.i[toShow.i], exp.MeanRFArea.i[toShow.i], type="p", pch=".", col=1, cex=4,
		main=mainText, xlab="Initial Baseline CM", ylab="Final Focal Stim CM" );
	lines ( c(0,max(ctl.MeanRFArea.i[toShow.i])), c(0,max(ctl.MeanRFArea.i[toShow.i])) );

		#	Spatial distribution: delta CM and delta RFAreas.

	iTrimmed = TrimEdgesFromCellList ( seq(1,N2), N, iTrimRings );

	delta.CM.e = ( exp.CM.e - ctl.CM.e ) / ctl.CM.e;
	delta.CM.i = ( exp.CM.i - ctl.CM.i ) / ctl.CM.i;

	delta.MeanRFArea.e = ( exp.MeanRFArea.e - ctl.MeanRFArea.e ) / ctl.MeanRFArea.e;
	delta.MeanRFArea.i = ( exp.MeanRFArea.i - ctl.MeanRFArea.i ) / ctl.MeanRFArea.i;

	x11(); par ( mfrow=c(2,2) );
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";
	mainText = paste ( "% Delta CM Per Input Node", mainTextIn, "Cortical E-Cells",sep="\n" );
	ShowVecAsMapContour ( delta.CM.e[iTrimmed], mainText, xlab, ylab );

	mainText = paste ( "% Delta CM Per Input Node", mainTextIn, "Cortical I-Cells",sep="\n" );
	ShowVecAsMapContour( delta.CM.i[iTrimmed], mainText, xlab, ylab );

	mainText = paste ( "% Delta RF Area Per Input Node", mainTextIn, "Cortical E-Cells",sep="\n" );
	ShowVecAsMapContour( delta.MeanRFArea.e[iTrimmed], mainText, xlab, ylab );

	mainText = paste ( "% Delta RF Area Per Input Node", mainTextIn, "Cortical I-Cells",sep="\n" );
	ShowVecAsMapContour( delta.MeanRFArea.i[iTrimmed], mainText, xlab, ylab );

	#x11(); par ( mfrow=c(2,2) );
	#tmp = rep ( 0, N2 ); tmp [ which ( delta.CM.e >= delta.CM.min ) ] = 1; ShowVecAsMap ( tmp[iTrimmed], "delta.CM.e >= 0" );
	#tmp = rep ( 0, N2 ); tmp [ which ( delta.CM.i >= delta.CM.min ) ] = 1; ShowVecAsMap ( tmp[iTrimmed], "delta.CM.i >= 0" );
	#tmp = rep ( 0, N2 ); tmp [ which ( delta.CM.e < -delta.CM.min ) ] = 1; ShowVecAsMap ( tmp[iTrimmed], "delta.CM.e < 0" );
	#tmp = rep ( 0, N2 ); tmp [ which ( delta.CM.i < -delta.CM.min ) ] = 1; ShowVecAsMap ( tmp[iTrimmed], "delta.CM.i < 0" );

	#x11(); par ( mfrow=c(2,2) );
	#tmp = rep ( 0, N2 ); tmp [ which ( delta.MeanRFArea.e >= delta.MeanRFArea.min ) ] = 1; ShowVecAsMap ( tmp[iTrimmed], "delta.RFArea.e >= 0" );
	#tmp = rep ( 0, N2 ); tmp [ which ( delta.MeanRFArea.i >= delta.MeanRFArea.min ) ] = 1; ShowVecAsMap ( tmp[iTrimmed], "delta.RFArea.i >= 0" );
	#tmp = rep ( 0, N2 ); tmp [ which ( delta.MeanRFArea.e < -delta.MeanRFArea.min ) ] = 1; ShowVecAsMap ( tmp[iTrimmed], "delta.RFArea.e < 0" );
	#tmp = rep ( 0, N2 ); tmp [ which ( delta.MeanRFArea.i < -delta.MeanRFArea.min ) ] = 1; ShowVecAsMap ( tmp[iTrimmed], "delta.CM.i < 0" );

	return ( list ( exp.CM.e=exp.CM.e, ctl.CM.e=ctl.CM.e, exp.CM.i=exp.CM.i, ctl.CM.i=ctl.CM.i,
				exp.MeanRFArea.e=exp.MeanRFArea.e, ctl.MeanRFArea.e=ctl.MeanRFArea.e,
				exp.MeanRFArea.i=exp.MeanRFArea.i, ctl.MeanRFArea.e=ctl.MeanRFArea.e,
				delta.CM.e=delta.CM.e, delta.CM.i=delta.CM.i,
				delta.MeanRFArea.e=delta.MeanRFArea.e, delta.MeanRFArea.i=delta.MeanRFArea.i ) );
	
} # CMPointByPoint = function ( ctl.CMRaw.e, exp.CMRaw.e, ctl.CMRaw.i, exp.CMRaw.i, ...


	#
	#	START: This section contains code that is useful for debugging the routine.
	#
#iTrimRings = 1;
#delta.CM.min = 0.1;
#delta.MeanRFArea.min = 0.1;
#CMPointByPoint ( init.cortical.Amp.e, final.cortical.Amp.e, init.cortical.Amp.i, final.cortical.Amp.i,
			#init.rfMap.rfAreas.e, final.rfMap.rfAreas.e, init.rfMap.rfAreas.i, final.rfMap.rfAreas.i,
			#iTrimRings, delta.CM.min, delta.MeanRFArea.min );

#ctl.CMRaw.e = init.cortical.Amp.e; exp.CMRaw.e = final.cortical.Amp.e;ctl.CMRaw.i = init.cortical.Amp.i; exp.CMRaw.i = final.cortical.Amp.i;
#ctl.RFArea.e = init.rfMap.rfAreas.e; exp.RFArea.e = final.rfMap.rfAreas.e; ctl.RFArea.i = init.rfMap.rfAreas.i; exp.RFArea.i =  final.rfMap.rfAreas.i;
#ShowVecAsMapContour ( delta.MeanRFArea.e[iTrimmed], "Title", "x", "y" );

	#
	#	END: This section contains code that is useful for debugging the routine.
	#

