		#	Initial Conditions

	iWhich = 1;
	titleText = paste ( "Col. # ", iCell, " E<-S Wghts\nBaseline Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( init.w1.e0, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 1, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 1;
	titleText = paste ( "Col. # ", iCell, " E<-E Wghts\nBaseline Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( init.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );

	iWhich = 1;
	titleText = paste ( "Col. # ", iCell, " E<-I Wghts\nBaseline Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( init.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );

	iWhich = 1;
	titleText = paste ( "Col. # ", iCell, " I<-E Wghts\nBaseline Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( init.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );