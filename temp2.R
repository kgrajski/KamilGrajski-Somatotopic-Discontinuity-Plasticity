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
	
	vt.e.delta = xt.exp$vt.e - xt.base$vt.e;
	vt.i.delta = xt.exp$vt.i - xt.base$vt.i;
	vt.0.delta = xt.exp$vt.0 - xt.base$vt.0;

} # TSResponse = function ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale ) {










