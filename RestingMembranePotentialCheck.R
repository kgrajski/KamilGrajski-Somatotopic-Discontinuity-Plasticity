#
#	RFProbeResponsesCheck.Time3.R
#

#	Clear the workspace.
rm(list = ls());

#

source ( "NMHelperFunctions.R" );


RestingMPSpatialPlot = function ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleTextRoot, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fName ) {

	if ( tiffFlag ) {
		tiff ( paste ( "RestingMP", fName, "tiff", sep="." ), compression="lzw", units="in", width=7.0, height=7.0, res=300 );
	} else {
		x11(); par ( mfcol=c(3,2) );
	} # if ( tiffFlag )

	par ( mfrow = c ( 3,2 ) );
	xLabText  = "Distal -> Proximal";
	yLabText = "Digit 1 -> Digit 3";
	boundaryMarks = c ( as.integer(N/3)+0.5, as.integer(N/3)*2+0.5 );

	restingMP.base.e = rep ( 0, N2 ); restingMP.base.i = rep ( 0, N2 );
	restingMP.exp.e = rep ( 0, N2 ); restingMP.exp.i = rep ( 0, N2 );

	iPreStim.length = 100;
	iPtrList.e = seq ( 1, N2*iPreStim.length, N2 );
	iPtrList.i = numValsPerRFTrial + iPtrList.e;

	for ( iCell in 1:N2 ) {

			#	Using knowledge of data structure details to speed things up.
		restingMP.base.e[iCell] = mean ( alldata.base[iPtrList.e] );
		restingMP.base.i[iCell] = mean ( alldata.base[iPtrList.i] );
		restingMP.exp.e[iCell] = mean ( alldata.exp[iPtrList.e] );
		restingMP.exp.i[iCell] = mean ( alldata.exp[iPtrList.i] );

		iPtrList.e = iPtrList.e + 1;
		iPtrList.i = iPtrList.i + 1;

	} # for ( iCell in 1:N2 )

	delta.restingMP.e = restingMP.exp.e - restingMP.base.e;
	delta.restingMP.i = restingMP.exp.i - restingMP.base.i;

	titleText = paste ( "E Cell Mean Prestim Resting M.P.", titleTextRoot, sep="\n" );
	ShowVecAsMap2 ( restingMP.exp.e, titleText, xLabText, yLabText, min(restingMP.exp.e), max(restingMP.exp.e) );
	abline ( h = boundaryMarks, lty=3, col=3 );

	titleText = paste ( "I Cell Mean Prestim Resting M.P.", titleTextRoot, sep="\n" );
	ShowVecAsMap2 ( restingMP.exp.i, titleText, xLabText, yLabText, min(restingMP.exp.i), max(restingMP.exp.i) );
	abline ( h = boundaryMarks, lty=3, col=3 );

	titleText = paste ( "E Cell Mean Prestim Resting M.P.", "Baseline Refined Network", sep="\n" );
	ShowVecAsMap2 ( restingMP.base.e, titleText, xLabText, yLabText, min(restingMP.base.e), max(restingMP.base.e) );
	abline ( h = boundaryMarks, lty=3, col=3 );

	titleText = paste ( "I Cell Mean Prestim Resting M.P.", "Baseline Refined Network", sep="\n" );
	ShowVecAsMap2 ( restingMP.base.i, titleText, xLabText, yLabText, min(restingMP.base.i), max(restingMP.base.i) );
	abline ( h = boundaryMarks, lty=3, col=3 );

	titleText = paste ( "E Cell Mean Prestim Delta Resting M.P.", "Knockout - Baseline", sep="\n" );
	ShowVecAsMap2 ( delta.restingMP.e, titleText, xLabText, yLabText, min(delta.restingMP.e), max(delta.restingMP.e) );
	abline ( h = boundaryMarks, lty=3, col=3 );

	titleText = paste ( "I Cell Mean Prestim Delta Resting M.P.", "Knockout - Baseline", sep="\n" );
	ShowVecAsMap2 ( delta.restingMP.i, titleText, xLabText, yLabText, min(delta.restingMP.i), max(delta.restingMP.i) );
	abline ( h = boundaryMarks, lty=3, col=3 );

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

} # RestingMPSpatialPlot = function ( N, numValsPerRFTrial, iProbeCell...




	#
	#	Additional helper functions.
	#

	#
	#	Set some constants
	#
N = 45;
N2 = N * N;
deltaT = 0.001;
trialDuration = 0.500;
numItersPerTrial = 500;
numValsPerRFTrial = numItersPerTrial * N * N;
tiffFlag = TRUE;
initialPass = FALSE;
iShortPlot = FALSE;

	#
	#	Identify the pool of recording sites of potential interest.
	#
iShowTheseRows = seq ( 18, 32 );
iCellList = c ( 13*N+iShowTheseRows, 14*N+iShowTheseRows,
			15*N+iShowTheseRows, 16*N+iShowTheseRows,
			787, 832, 877, 791, 836, 881 );

	#
	#	Identify the speciufic RF probes of interest.
	#
iSetProbeList = c ( 652 );

	#
	#	Set up titles
	#
subTitleTextRoot = NULL; xlabTextRoot = "Msec."; ylabTextRoot = "Vs(t), Vi(t)";
eScale = 1.0;

	#
	#	Identify the reference and the experimental datasets.
	#
iKnockLengthList = c ( 4, 8, 12 );
for ( iKnockLength in iKnockLengthList ) { 

	fDir = paste  ("E:/NMLab/Working/S.45.7B.Fix2.Knockout.", iKnockLength, "/", sep="" );
	dirRData = fDir;

		#
		#	Set up final final names and read the data.
		#
	fName.base = paste ( "BorderKnockout_Control_Placebo.Raw", iKnockLength, 0, "bin", sep="." );
	fName.ctl.I = paste ( "BorderKnockout_Control_I.RFMap.Raw", iKnockLength, 1, "bin", sep="." );
	fName.ctl.E = paste ( "BorderKnockout_Control_E.RFMap.Raw", iKnockLength, 2, "bin", sep="." );
	fName.exp.I = paste ( "BorderKnockout_OneSide_I.RFMap.Raw", iKnockLength, 4, "bin", sep="." );
	fName.exp.E = paste ( "BorderKnockout_OneSide_E.RFMap.Raw", iKnockLength, 5, "bin", sep="." );

	fName.base = paste ( fDir, fName.base, sep="" );
	fName.ctl.E = paste ( fDir, fName.ctl.E, sep="" );
	fName.ctl.I = paste ( fDir, fName.ctl.I, sep="" );
	fName.exp.E = paste ( fDir, fName.exp.E, sep="" );
	fName.exp.I = paste ( fDir, fName.exp.I, sep="" );

		#
		#	Get the subset of data of interest and save it.
		#	The files can be quite large so have to split doing the control and exps modes.
		#
	if ( initialPass ) {
		base = ExtractRFProbeRawData ( fName.base, N2, numValsPerRFTrial, iCellList );
		ctl.E = ExtractRFProbeRawData ( fName.ctl.E, N2, numValsPerRFTrial, iCellList );
		ctl.I = ExtractRFProbeRawData ( fName.ctl.I, N2, numValsPerRFTrial, iCellList );
		save( base, file = paste ( dirRData, paste ( "Base", iKnockLength, , "bin", sep="." ) ), sep="/" );
		save( ctl.E, file = paste ( dirRData, paste ( "Control", "E", iKnockLength, "bin", sep="." ) ), sep="/" );
		save( ctl.I, file = paste ( dirRData, paste ( "Control", "I", iKnockLength, "bin", sep="." ) ), sep="/" );
	} else {
		load( file = paste ( dirRData, paste ( "Base", iKnockLength, "bin", sep="." ), sep="" ) );
		load( file = paste ( dirRData, paste ( "Control", "E", iKnockLength, "bin", sep="." ), sep="/" ) );
		load( file = paste ( dirRData, paste ( "Control", "I", iKnockLength, "bin", sep="." ), sep="/" ) );
	} # if ( initialPass )

		#
		#	Time series difference plots.
		#
	for ( iProbe in 1:length(iSetProbeList) ) {
	
		iProbeCell = iSetProbeList[iProbe];
		indexToRFProbeRawData = which ( iCellList == iProbeCell );
		alldata.base = c ( base$v1e[,indexToRFProbeRawData], base$v1i[,indexToRFProbeRawData], base$v0[,indexToRFProbeRawData] );
		
		titleLineOne = "Control E-Only Knockout";
		titleText = paste ( titleLineOne, iKnockLength, sep=":" );
		fNameBase = paste ( titleLineOne, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( ctl.E$v1e[,indexToRFProbeRawData], ctl.E$v1i[,indexToRFProbeRawData], ctl.E$v0[,indexToRFProbeRawData] );
		RestingMPSpatialPlot ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase );

		titleLineOne = "Control I-Only Knockout";
		titleText = paste ( titleLineOne , iKnockLength, sep=":" );
		fNameBase = paste ( titleLineOne , iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( ctl.I$v1e[,indexToRFProbeRawData], ctl.I$v1i[,indexToRFProbeRawData], ctl.I$v0[,indexToRFProbeRawData] );
		RestingMPSpatialPlot ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase );

	} # for ( iProbeCell in ...

	rm ( ctl.E );
	rm ( ctl.I );

	if ( initialPass ) {
		exp.E = ExtractRFProbeRawData ( fName.exp.E, N2, numValsPerRFTrial, iCellList );
		exp.I = ExtractRFProbeRawData ( fName.exp.I, N2, numValsPerRFTrial, iCellList );
		save( exp.E, file = paste ( dirRData, paste ( "OneSide", "E", iKnockLength, "bin", sep="." ) ), sep="/" );
		save( exp.I, file = paste ( dirRData, paste ( "OneSide", "I", iKnockLength, "bin", sep="." ) ), sep="/" );
	} else {
		load( file = paste ( dirRData, paste ( "OneSide", "E", iKnockLength, "bin", sep="." ), sep="/" ) );
		load( file = paste ( dirRData, paste ( "OneSide", "I", iKnockLength, "bin", sep="." ), sep="/" ) );
	} # 	if ( initialPass ) {

		#
		#	Time series difference plots.
		#
	for ( iProbe in 1:length(iSetProbeList) ) {
	
		iProbeCell = iSetProbeList[iProbe];
		indexToRFProbeRawData = which ( iCellList == iProbeCell );
		alldata.base = c ( base$v1e[,indexToRFProbeRawData], base$v1i[,indexToRFProbeRawData], base$v0[,indexToRFProbeRawData] );

		titleLineOne = "One-Sided E-Only Knockout";
		titleText = paste ( titleLineOne, iKnockLength, sep=":" );
		fNameBase = paste ( titleLineOne, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( exp.E$v1e[,indexToRFProbeRawData], exp.E$v1i[,indexToRFProbeRawData], exp.E$v0[,indexToRFProbeRawData] );
		RestingMPSpatialPlot ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase );

		titleLineOne = "One-Sided I-Only Knockout";
		titleText = paste ( titleLineOne, iKnockLength, sep=":" );
		fNameBase = paste ( titleLineOne, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( exp.I$v1e[,indexToRFProbeRawData], exp.I$v1i[,indexToRFProbeRawData], exp.I$v0[,indexToRFProbeRawData] );
		RestingMPSpatialPlot ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase );

	} # for ( iProbeCell in ...

	rm ( base );
	rm ( exp.E );
	rm ( exp.I );

} # for ( iKnockLength in iKnockLengthlist ) {



