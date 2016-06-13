#
#	Do some spatial maps of RF response differences.
#

#
#	Plot a bunch of RFs side by side.
#

PlotRFPanel = function ( exp, ref, exp.titleText, ref.titleText, iCellList ) {

	boundaryMarks = c((N/3)+0.5, (2*N/3)+0.5);
	for ( iCell in iCellList ) {

		x11(); par(mfcol=c(2,2));
		xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

		titleText = paste ( "Column # ", iCell, " E Cell RF Extent\n", ref.titleText, sep="" );
		ShowVecAsMap1 ( ref$r1.e.rfMap[iCell,], titleText, xlab, ylab );
		abline ( h = boundaryMarks, lty=3, col=3 );

		iWhich = 1;
		titleText = paste ( "Column # ", iCell, " I Cell RF Extent\n", ref.titleText, sep="" );
		ShowVecAsMap1 ( ref$r1.i.rfMap[iCell,], titleText, xlab, ylab );
		abline ( h = boundaryMarks, lty=3, col=3 );

		iWhich = iBase;
		titleText = paste ( "Column # ", iCell, " E Cell RF Extent\n", exp.titleText, sep="" );
		ShowVecAsMap1 ( exp$r1.e.rfMap[iCell,], titleText, xlab, ylab );
		abline ( h = boundaryMarks, lty=3, col=3 );

		iWhich = iBase;
		titleText = paste ( "Column # ", iCell, " I Cell RF Extent\n", exp.titleText, sep="" );
		ShowVecAsMap1 ( exp$r1.i.rfMap[iCell,], titleText, xlab, ylab );
		abline ( h = boundaryMarks, lty=3, col=3 );

	} # 	for ( iCell in iCellList ) {

} # for ( iCell in iCellList ) {



ScanForDiffs1 = function ( refMap, expMap, kCell, iLookDirection, iStart, iEnd, iWindow ) {

	refMap = ref$r1.e.rfMap;
	expMap = exp$r1.e.rfMap;
	#LookDirection = 1;
	kCell = 651;

	N2 = as.integer ( sqrt ( length ( refMap ) ) );
	N = as.integer ( sqrt ( N2 ) );

	kCell.row = GetRow ( kCell, N );
	kCell.col = GetCol ( kCell, N );
	
	jRowList = seq ( kCell.row - iWindow, kCell.row + iWindow );
	jCol = kCell.col + iLookDirection;

	tList = GetLin ( jRowList, jCol, N );

	tmp = wilcox.test ( expMap[tList,kCell], refMap[tList,kCell], correct=FALSE, exact = FALSE, paired=TRUE );
	return ( tmp$p.value );

} # ScanForDiffs1 = function ( ref, exp ) {

	#
	#
ScanForDiffsFrontEnd = function ( ref, exp, iTrim, iLookDirection, iWindow ) {

	N2 = as.integer ( sqrt ( length ( ref$r1.e.rfMap ) ) );
	N = as.integer ( sqrt ( N2 ) );
	digitWidth = as.integer ( N / 3 );
	if ( iLookDirection > 0 ) {
		iStart = ( digitWidth - 1 ) * N + iTrim + 1 + iWindow;
		iEnd = digitWidth * N - iTrim - iWindow;
	} else {
		iStart = ( digitWidth - 1 + 1 ) * N + iTrim + 1 + iWindow;
		iEnd = ( digitWidth + 1) * N - iTrim - iWindow;
	} # if ( iLookDirection > 0 )

	iCellList = seq ( iStart, iEnd );
	iLen = length ( iCellList );
	res.e = matrix ( 0, nrow=iLen, ncol=2 );
	res.i = matrix ( 0, nrow=iLen, ncol=2 );

	for ( iCell in 1:iLen ) {
		res.e[iCell,1] = iCellList[iCell];
		res.i[iCell,1] = iCellList[iCell];
		res.e[iCell,2]= ScanForDiffs1 ( ref$r1.e.rfMap, exp$r1.e.rfMap, iCellList[iCell], iLookDirection, iStart, iEnd, iWindow );
		res.i[iCell,2]= ScanForDiffs1 ( ref$r1.i.rfMap, exp$r1.i.rfMap, iCellList[iCell], iLookDirection, iStart, iEnd, iWindow );
	} # for ( iCell in iStart:iEnd )

	return ( list ( res.e=res.e, res.i=res.i ) );

} # ScanForDiffsFrontEnd = function ( ref, exp, numValsPerIter, iTrim, iLookDirection, iWindow )...


	#


ref = rnet;
exp = base;
iTrim = 2;
iWindow = 3;
iLookDirection = -1;

rnet.base = ScanForDiffsFrontEnd ( rnet, base, iTrim, iLookDirection, iWindow );
base.test.ctl.I = ScanForDiffsFrontEnd ( base, test.ctl.I, iTrim, iLookDirection, iWindow );
base.test.ctl.E = ScanForDiffsFrontEnd ( base, test.ctl.E, iTrim, iLookDirection, iWindow );
base.test.oneside.I = ScanForDiffsFrontEnd ( base, test.oneside.I, iTrim, iLookDirection, iWindow );
base.test.oneside.E = ScanForDiffsFrontEnd ( base, test.oneside.E, iTrim, iLookDirection, iWindow );

tmp.e = cbind ( rnet.base$res.e[,2], base.test.ctl.I$res.e[,2], base.test.ctl.E$res.e[,2],
			base.test.oneside.I$res.e[,2], base.test.oneside.E$res.e[,2] );

tmp.i = cbind ( rnet.base$res.i[,2], base.test.ctl.I$res.i[,2], base.test.ctl.E$res.i[,2],
			base.test.oneside.I$res.i[,2], base.test.oneside.E$res.i[,2] );

x11(); par(mfrow=c(2,1));
boxplot(tmp.e);
boxplot(tmp.i);

pValue = 0.05;
tmp.summary.e = cbind ( rnet.base$res.i[,1], tmp.e<pValue );
tmp.summary.i = cbind ( rnet.base$res.i[,1], tmp.i<pValue );

for ( iCell in seq ( 680, 715 ) ) {


	zmin = min ( log10(rnet$r1.e.rfMap[,iCell]), log10(base$r1.e.rfMap[,iCell]), log10(test.ctl.I$r1.e.rfMap[,iCell]),
				log10(test.ctl.E$r1.e.rfMap[,iCell]),
				log10(test.oneside.I$r1.e.rfMap[,iCell]), log10(test.oneside.E$r1.e.rfMap[,iCell]));


	zmax = max ( log10(rnet$r1.e.rfMap[,iCell]), log10(base$r1.e.rfMap[,iCell]), log10(test.ctl.I$r1.e.rfMap[,iCell]),
				log10(test.ctl.E$r1.e.rfMap[,iCell]),
				log10(test.oneside.I$r1.e.rfMap[,iCell]), log10(test.oneside.E$r1.e.rfMap[,iCell]));

	xlabText = "TBD"; ylabText = "TBD";
	x11(); par(mfrow=c(2,3));
	ShowVecAsMap2 ( log10(rnet$r1.e.rfMap[,iCell]), paste("RNET", iCell, sep=" "), xlabText, ylabText, zmin, zmax );
	abline ( h = 15.5, lty=3, col=5 );

	ShowVecAsMap2 ( log10(test.ctl.I$r1.e.rfMap[,iCell]), "CTL_I", xlabText, ylabText, zmin, zmax );
	abline ( h = 15.5, lty=3, col=5 );

	ShowVecAsMap2 ( log10(test.ctl.E$r1.e.rfMap[,iCell]), "CTL_E", xlabText, ylabText, zmin, zmax );
	abline ( h = 15.5, lty=3, col=5 );

	ShowVecAsMap2 ( log10(base$r1.e.rfMap[,iCell]), "BASE", xlabText, ylabText, zmin, zmax );
	abline ( h = 15.5, lty=3, col=5 );

	ShowVecAsMap2 ( log10(test.oneside.I$r1.e.rfMap[,iCell]), "ONE_SIDE_I", xlabText, ylabText, zmin, zmax );
	abline ( h = 15.5, lty=3, col=5 );

	ShowVecAsMap2 ( log10(test.oneside.E$r1.e.rfMap[,iCell]), "ONE_SIDE_E", xlabText, ylabText, zmin, zmax );
	abline ( h = 15.5, lty=3, col=5 );

} # for ( i in seq ( 680, 715 ) ) {


zmin = min ( (test.ctl.I$r1.e.rfMap[,iCell]/base$r1.e.rfMap[,iCell])-1,
		(test.ctl.E$r1.e.rfMap[,iCell]/base$r1.e.rfMap[,iCell])-1,
		(test.oneside.I$r1.e.rfMap[,iCell]/base$r1.e.rfMap[,iCell])-1,
		(test.oneside.E$r1.e.rfMap[,iCell]/base$r1.e.rfMap[,iCell])-1 );

zmax = max (  (test.ctl.I$r1.e.rfMap[,iCell]/base$r1.e.rfMap[,iCell])-1,
		(test.ctl.E$r1.e.rfMap[,iCell]/base$r1.e.rfMap[,iCell])-1,
		(test.oneside.I$r1.e.rfMap[,iCell]/base$r1.e.rfMap[,iCell])-1,
		(test.oneside.E$r1.e.rfMap[,iCell]/base$r1.e.rfMap[,iCell])-1 );


xlabText = "TBD"; ylabText = "TBD";
x11(); par(mfrow=c(2,3));
ShowVecAsMap1Log ( (rnet$r1.e.rfMap[,iCell]), "RNET", xlabText, ylabText);
abline ( h = 15.5, lty=3, col=5 );

ShowVecAsMap1Log ( (test.ctl.I$r1.e.rfMap[,iCell]/base$r1.e.rfMap[,iCell]), "CTL_I", xlabText, ylabText );
abline ( h = 15.5, lty=3, col=5 );

ShowVecAsMap1Log ( (test.ctl.E$r1.e.rfMap[,iCell]/base$r1.e.rfMap[,iCell]), "CTL_E", xlabText, ylabText );
abline ( h = 15.5, lty=3, col=5 );

ShowVecAsMap1Log ( (base$r1.e.rfMap[,iCell]), "BASE", xlabText, ylabText );
abline ( h = 15.5, lty=3, col=5 );

ShowVecAsMap1Log ( (test.oneside.I$r1.e.rfMap[,iCell]/base$r1.e.rfMap[,iCell]), "ONE_SIDE_I", xlabText, ylabText);
abline ( h = 15.5, lty=3, col=5 );

ShowVecAsMap1Log ( (test.oneside.E$r1.e.rfMap[,iCell]/base$r1.e.rfMap[,iCell]), "ONE_SIDE_E", xlabText, ylabText );
abline ( h = 15.5, lty=3, col=5 );

#
#	Just plot a lot of RFs
#

iCellList = 650 + seq ( -10, 10, 1 ) * N;
exp.titleText = "EXP"
ref.titleText = "REF"

exp = test.oneside.I;
ref = base;
PlotRFPanel ( test.twoside.I, base, exp.titleText, ref.titleText, iCellList );



	#
	#	The following is an attempt to generate a summary slide.
	#

for ( iCell in seq ( 680, 715 ) ) {


	zmin = min ( 	log10(rnet$r1.e.rfMap[,iCell]), log10(base$r1.e.rfMap[,iCell]),
				log10(test.ctl.I$r1.e.rfMap[,iCell]), log10(test.ctl.E$r1.e.rfMap[,iCell]),
				log10(test.oneside.I$r1.e.rfMap[,iCell]), log10(test.oneside.E$r1.e.rfMap[,iCell]),
				log10(rnet$r1.i.rfMap[,iCell]), log10(base$r1.i.rfMap[,iCell]),
				log10(test.ctl.I$r1.i.rfMap[,iCell]), log10(test.ctl.E$r1.i.rfMap[,iCell]),
				log10(test.oneside.I$r1.i.rfMap[,iCell]), log10(test.oneside.E$r1.i.rfMap[,iCell])
			);


	zmax = max ( 	log10(rnet$r1.e.rfMap[,iCell]), log10(base$r1.e.rfMap[,iCell]),
				log10(test.ctl.I$r1.e.rfMap[,iCell]), log10(test.ctl.E$r1.e.rfMap[,iCell]),
				log10(test.oneside.I$r1.e.rfMap[,iCell]), log10(test.oneside.E$r1.e.rfMap[,iCell]),
				log10(rnet$r1.i.rfMap[,iCell]), log10(base$r1.i.rfMap[,iCell]),
				log10(test.ctl.I$r1.i.rfMap[,iCell]), log10(test.ctl.E$r1.i.rfMap[,iCell]),
				log10(test.oneside.I$r1.i.rfMap[,iCell]), log10(test.oneside.E$r1.i.rfMap[,iCell])
			);

	xlabText = "TBD"; ylabText = "TBD";
	x11(); par(mfrow=c(2,3));
	ShowVecAsMap2 ( log10(rnet$r1.e.rfMap[,iCell]), paste("E Cell RNET", iCell, sep=" "), xlabText, ylabText, zmin, zmax );
	abline ( h = 15.5, lty=3, col=5 );

	ShowVecAsMap2 ( log10(test.ctl.I$r1.e.rfMap[,iCell]), "E Type CTL_I", xlabText, ylabText, zmin, zmax );
	abline ( h = 15.5, lty=3, col=5 );

	ShowVecAsMap2 ( log10(test.ctl.E$r1.e.rfMap[,iCell]), "E Type CTL_E", xlabText, ylabText, zmin, zmax );
	abline ( h = 15.5, lty=3, col=5 );

	ShowVecAsMap2 ( log10(base$r1.e.rfMap[,iCell]), "E Type BASE", xlabText, ylabText, zmin, zmax );
	abline ( h = 15.5, lty=3, col=5 );

	ShowVecAsMap2 ( log10(test.oneside.I$r1.e.rfMap[,iCell]), "E Type ONE_SIDE_I", xlabText, ylabText, zmin, zmax );
	abline ( h = 15.5, lty=3, col=5 );

	ShowVecAsMap2 ( log10(test.oneside.E$r1.e.rfMap[,iCell]), "E Type ONE_SIDE_E", xlabText, ylabText, zmin, zmax );
	abline ( h = 15.5, lty=3, col=5 );


} # for ( i in seq ( 680, 715 ) ) {














