#
#	KnockoutDeepRFAnalysis.R
#

#
#	Use in conjunction with: KnockoutFigure family of routines.
#

KnockoutRFDeepDive = function ( rnet, base, test.ctl.I, test.ctl.E, test.oneside.I, test.oneside.E, test.twoside.E, test.twoside.I,
						x, y.ctl, y.oneside, iCellList, tiffFlag, iKnockOutLength ) {

	y.ctl[1] = y.ctl[1] - 0.5; y.ctl[2] = y.ctl[2] + 0.5;
	y.oneside[1] = y.oneside[1] - 0.5; y.oneside[2] = y.oneside[2] + 0.5;
	x[1] = x[1] - 0.5; x[2] = x[2] + 0.5;
	tmp.x = rbind(c(x[1], x[2]), c(x[1], x[2]), c(x[1], x[1]), c(x[2], x[2]));
	tmp.y.ctl = rbind(c(y.ctl[1], y.ctl[1]), c(y.ctl[2], y.ctl[2]), c(y.ctl[1], y.ctl[2]), c(y.ctl[1], y.ctl[2]));
	tmp.y.oneside = rbind( c ( y.oneside[1], y.oneside[1]), c(y.oneside[2], y.oneside[2]),
					c(y.oneside[1], y.oneside[2]), c(y.oneside[1], y.oneside[2]) );

	for ( iCell in iCellList ) {


		zmin = min ( 	log10(rnet$r1.e.rfMap[iCell,]), log10(base$r1.e.rfMap[iCell,]),
				log10(test.ctl.I$r1.e.rfMap[iCell,]), log10(test.ctl.E$r1.e.rfMap[iCell,]),
				log10(test.oneside.I$r1.e.rfMap[iCell,]), log10(test.oneside.E$r1.e.rfMap[iCell,]),
				log10(rnet$r1.i.rfMap[iCell,]), log10(base$r1.i.rfMap[iCell,]),
				log10(test.ctl.I$r1.i.rfMap[iCell,]), log10(test.ctl.E$r1.i.rfMap[iCell,]),
				log10(test.oneside.I$r1.i.rfMap[iCell,]), log10(test.oneside.E$r1.i.rfMap[iCell,])
			);


		zmax = max ( 	log10(rnet$r1.e.rfMap[iCell,]), log10(base$r1.e.rfMap[iCell,]),
				log10(test.ctl.I$r1.e.rfMap[iCell,]), log10(test.ctl.E$r1.e.rfMap[iCell,]),
				log10(test.oneside.I$r1.e.rfMap[iCell,]), log10(test.oneside.E$r1.e.rfMap[iCell,]),
				log10(rnet$r1.i.rfMap[iCell,]), log10(base$r1.i.rfMap[iCell,]),
				log10(test.ctl.I$r1.i.rfMap[iCell,]), log10(test.ctl.E$r1.i.rfMap[iCell,]),
				log10(test.oneside.I$r1.i.rfMap[iCell,]), log10(test.oneside.E$r1.i.rfMap[iCell,])
			);

		xlabText = "Distal -> Proximal"; ylabText = "Digit 1 -> Digit 3";
		if ( tiffFlag ) {
			tiff ( paste("Knock", "DeepRF", iCell, iKnockOutLength, "E", "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
		} else {
			x11();
		} # if ( tiffFlag )
		par(mfrow=c(2,3));
		titleText = paste("E Cell", iCell, "\nRandom Network", sep=" ");
		ShowVecAsMap2 ( log10(rnet$r1.e.rfMap[iCell,]), titleText, xlabText, ylabText, zmin, zmax );
		GenOutline1X ( tmp.x, tmp.y.ctl, "white", 1, 0.5 );
		GenOutline1X ( tmp.x, tmp.y.oneside, "white", 1, 0.5 );
		abline ( h = 15.5, lty=3, col=5 );

		titleText = paste("E Cell", iCell, "\nCore Knockout I Cells", sep=" ");
		ShowVecAsMap2 ( log10(test.ctl.I$r1.e.rfMap[iCell,]), titleText, xlabText, ylabText, zmin, zmax );
		GenOutline1X ( tmp.x, tmp.y.ctl, "white", 1, 0.5 );
		abline ( h = 15.5, lty=3, col=5 );

		titleText = paste("E Cell", iCell, "\nCore Knockout E Cells", sep=" ");
		ShowVecAsMap2 ( log10(test.ctl.E$r1.e.rfMap[iCell,]), titleText, xlabText, ylabText, zmin, zmax );
		GenOutline1X ( tmp.x, tmp.y.ctl, "white", 1, 0.5 );
		abline ( h = 15.5, lty=3, col=5 );

		titleText = paste("E Cell", iCell, "\nBaseline Refined Network", sep=" ");
		ShowVecAsMap2 ( log10(base$r1.e.rfMap[iCell,]), titleText, xlabText, ylabText, zmin, zmax );
		GenOutline1X ( tmp.x, tmp.y.ctl, "white", 1, 0.5 );
		GenOutline1X ( tmp.x, tmp.y.oneside, "white", 1, 0.5 );
		abline ( h = 15.5, lty=3, col=5 );

		titleText = paste("E Cell", iCell, "\nBorder Knockout I Cells", sep=" ");
		ShowVecAsMap2 ( log10(test.oneside.I$r1.e.rfMap[iCell,]), titleText, xlabText, ylabText, zmin, zmax );
		GenOutline1X ( tmp.x, tmp.y.oneside, "white", 1, 0.5 );
		abline ( h = 15.5, lty=3, col=5 );

		titleText = paste("E Cell", iCell, "\nBorder Knockout E Cells", sep=" ");
		ShowVecAsMap2 ( log10(test.oneside.E$r1.e.rfMap[iCell,]), titleText, xlabText, ylabText, zmin, zmax );
		GenOutline1X ( tmp.x, tmp.y.oneside, "white", 1, 0.5 );
		abline ( h = 15.5, lty=3, col=5 );

		if ( tiffFlag ) {
			dev.off();
		} # if ( tiffFlag )


		if ( tiffFlag ) {
			tiff ( paste("Knock", "DeepRF", iCell, iKnockOutLength, "I", "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
		} else {
			x11();
		} # if ( tiffFlag )
		par(mfrow=c(2,3));
		titleText = paste("I Cell", iCell, "\nRandom Network", sep=" ");
		ShowVecAsMap2 ( log10(rnet$r1.i.rfMap[iCell,]), titleText, xlabText, ylabText, zmin, zmax );
		GenOutline1X ( tmp.x, tmp.y.ctl, "white", 1, 0.5 );
		GenOutline1X ( tmp.x, tmp.y.oneside, "white", 1, 0.5 );
		abline ( h = 15.5, lty=3, col=5 );

		titleText = paste("I Cell", iCell, "\nCore Knockout I Cells", sep=" ");
		ShowVecAsMap2 ( log10(test.ctl.I$r1.i.rfMap[iCell,]), titleText, xlabText, ylabText, zmin, zmax );
		GenOutline1X ( tmp.x, tmp.y.ctl, "white", 1, 0.5 );
		abline ( h = 15.5, lty=3, col=5 );

		titleText = paste("I Cell", iCell, "\nCore Knockout E Cells", sep=" ");
		ShowVecAsMap2 ( log10(test.ctl.E$r1.i.rfMap[iCell,]), titleText, xlabText, ylabText, zmin, zmax );
		GenOutline1X ( tmp.x, tmp.y.ctl, "white", 1, 0.5 );
		abline ( h = 15.5, lty=3, col=5 );

		titleText = paste("I Cell", iCell, "\nBaseline Refined Network", sep=" ");
		ShowVecAsMap2 ( log10(base$r1.i.rfMap[iCell,]), titleText, xlabText, ylabText, zmin, zmax );
		GenOutline1X ( tmp.x, tmp.y.ctl, "white", 1, 0.5 );
		GenOutline1X ( tmp.x, tmp.y.oneside, "white", 1, 0.5 );
		abline ( h = 15.5, lty=3, col=5 );

		titleText = paste("I Cell", iCell, "\nBorder Knockout I Cells", sep=" ");
		ShowVecAsMap2 ( log10(test.oneside.I$r1.i.rfMap[iCell,]), titleText, xlabText, ylabText, zmin, zmax );
		GenOutline1X ( tmp.x, tmp.y.oneside, "white", 1, 0.5 );
		abline ( h = 15.5, lty=3, col=5 );

		titleText = paste("I Cell", iCell, "\nBorder Knockout E Cells", sep=" ");
		ShowVecAsMap2 ( log10(test.oneside.E$r1.i.rfMap[iCell,]), titleText, xlabText, ylabText, zmin, zmax );
		GenOutline1X ( tmp.x, tmp.y.oneside, "white", 1, 0.5 );
		abline ( h = 15.5, lty=3, col=5 );

		if ( tiffFlag ) {
			dev.off();
		} # if ( tiffFlag )

	} # for ( i in iCellList ) {

} # KnockoutRFDeepDive = function ( rnet, base, test.ctl.I, test.ctl.E, test.oneside.I, test.oneside.E, test.twoside.E, test.twoside.I,



#
#	Production Area
#

x = rfExpZone.x;
y.ctl = rfExpZone.y.ctl;
y.oneside = rfExpZone.y.oneside;
tiffFlag = FALSE;

iFarSideFirstRow = 45 * 15 + seq ( 16, 25, 1 );
iFarSideSecondRow = 45 * 16 + seq ( 16, 25, 1 );
iCellList = c ( iFarSideFirstRow, iFarSideSecondRow );
iCellList = c ( 697, 742 );
KnockoutRFDeepDive ( rnet, base, test.ctl.I, test.ctl.E, test.oneside.I, test.oneside.E, test.twoside.E, test.twoside.I,
						rfExpZone.x, rfExpZone.y.ctl, rfExpZone.y.oneside, iCellList, tiffFlag, iKnockLength );

iNearSideFirstRow = 45 * 14 + seq ( 16, 25, 1 );
iNearSideSecondRow = 45 * 13 + seq ( 16, 25, 1 );
iCellList = c ( iNearSideFirstRow, iNearSideSecondRow );
KnockoutRFDeepDive ( rnet, base, test.ctl.I, test.ctl.E, test.oneside.I, test.oneside.E, test.twoside.E, test.twoside.I,
						rfExpZone.x, rfExpZone.y.ctl, rfExpZone.y.oneside, iCellList, tiffFlag, iKnockLength );












