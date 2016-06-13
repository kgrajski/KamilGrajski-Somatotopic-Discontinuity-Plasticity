#
#	Extract the vt values for the given iCell from the CPP generated RF Raw data file.
#	Compute the rt values.
#

TSResponse = function ( N2, numValsPerRFTrial, iCell, xt, eScale ) {

	iCellList = seq ( iCell, numValsPerRFTrial, N2 );
	startOffset.e = 1; startOffset.i = numValsPerRFTrial + 1;
	startOffset.0 = 2 * numValsPerRFTrial + 1;

	vt.e = ( xt [ (startOffset.e):(startOffset.e + numValsPerRFTrial - 1) ] ) / eScale;
	vt.e = v1e.base[iCellList];

	vt.i = xt [ (startOffset.i):(startOffset.i + numValsPerRFTrial - 1) ];
	vt.i = v1i.base[iCellList];

	vt.0 = xt [ (startOffset.0):(startOffset.0 + numValsPerRFTrial - 1) ];
	vt.0 = v0.base[iCellList];

	rt.e = sigmoid ( vt.e, 4 );
	rt.i = sigmoid ( vt.i, 4 );
	rt.0 = sigmoid ( vt.0, 4 );

	vt.min = min ( vt.e, vt.i, vt.0 );
	vt.max = max ( vt.e, vt.i, vt.0 );
	rt.min = min ( rt.e, rt.i, rt.0 );
	rt.max = max ( rt.e, rt.i, rt.0 );

	return ( list ( vt.e=vt.e, vt.i=vt.i, vt.0=vt.0, vt.min=vt.min, vt.max=vt.max,
			rt.e=rt.e, rt.i=rt.i, rt.0=rt.0, rt.min=rt.min, rt.max=rt.max ) );

} # TSResponse = function ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale ) {


#
#	Prepare a baseline vs experimental dataset for the given iCell
#

ExpRefTSResponse = function ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale ) {

	xt.base = TSResponse ( N2, numValsPerRFTrial, iCell, alldata.base, eScale );
	xt.exp = TSResponse ( N2, numValsPerRFTrial, iCell, alldata.exp, eScale );
	
	v1e.delta = xt.exp$vt.e - xt.base$vt.e;
	v1i.delta = xt.exp$vt.i - xt.base$vt.i;
	v0.delta = xt.exp$vt.0 - xt.base$vt.0;
	
	vt.min = min ( xt.base$vt.min, xt.exp$vt.min ); vt.max = max ( xt.base$vt.max, xt.exp$vt.max );
	delta.min = min ( v1e.delta, v1i.delta, v0.delta ); delta.max = max ( v1e.delta, v1i.delta, v0.delta );

	return ( list ( v1e.base=xt.base$vt.e, v1i.base=xt.base$vt.i, v0.base=xt.base$vt.0,
				v1e.exp=xt.exp$vt.e, v1i.exp=xt.exp$vt.i, v0.exp=xt.exp$vt.0,
				v1e.delta=v1e.delta, v1i.delta=v1i.delta, v0.delta=v0.delta,
				vt.min=vt.min, vt.max=vt.max, delta.min=delta.min, delta.max=delta.max ) );

} # TSResponse = function ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale ) {


#
#	Handle the plotting of the time series.
#

KnockoutTimeSeriesPlot = function ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleTextRoot, subTitleTextRoot, xlabTextRoot, ylabTextRoot ) {

	x11(); par ( mfrow=c(3,3) );

	iRecordCellsList = c ( iProbeCell + N - 1, iProbeCell, iProbeCell + N + 1 );

	for ( iCell in iRecordCellsList ) {

		tmp = ExpRefTSResponse ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale );

		titleText = titleTextRoot;
		xlabText = xlabTextRoot;
		ylabText = ylabTextRoot;

		plot ( tmp$v1e.exp, type="l", ylim=c(tmp$vt.min, tmp$vt.max), col=3, lty=1,
			main=titleText, xlab=xlabText, ylab=ylabText );
		lines ( tmp$v1i.exp, col=2, lty=1 );
		lines ( tmp$v0.exp, col=1, lty=1 );

		plot ( tmp$v1e.base, type="l", ylim=c(tmp$vt.min, tmp$vt.max), col=3, lty=1,
			main=titleText, xlab=xlabText, ylab=ylabText );
		lines ( tmp$v1i.base, col=2, lty=1 );
		lines ( tmp$v0.base, col=1, lty=1 );

		plot ( tmp$v1e.delta, type="l", ylim=c(tmp$delta.min, tmp$delta.max), col=3, lty=1,
			main=titleText, xlab=xlabText, ylab=ylabText );
		lines ( tmp$v1i.delta, col=2, lty=1 );
		lines ( tmp$v0.delta, col=1, lty=1 );

	} # for ( iCell in iRecordCellsList ) {

} # KnockoutTimeSeriesPlot = function ( N, numValsPerRFTrial, iProbeCell...
