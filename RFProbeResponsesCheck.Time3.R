#
#	RFProbeResponsesCheck.Time3.R
#

#	Clear the workspace.
rm(list = ls());

#

source ( "NMHelperFunctions.R" );

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
tiffFlag = FALSE;
initialPass = FALSE;
iShortPlot = FALSE;

	#
	#	Identify the pool of recording sites of potential interest.
	#
iShowTheseRows = seq ( 18, 32 );
iCellList = c ( 13*N+iShowTheseRows, 14*N+iShowTheseRows,
			15*N+iShowTheseRows, 16*N+iShowTheseRows,
			787, 832, 877, 791, 836, 881 );

iCellList = c ( 6*N+iShowTheseRows, 7*N+iShowTheseRows,
			8*N+iShowTheseRows, 9*N+iShowTheseRows,
			787, 832, 877, 791, 836, 881 );

	#
	#	Identify the speciufic RF probes of interest.
	#
iSetProbeList = c ( 652, 657 );

	#
	#	Set up titles
	#
subTitleTextRoot = NULL; xlabTextRoot = "Msec."; ylabTextRoot = "Vs(t), Vi(t)";
eScale = 1.0;

	#
	#	Identify the reference and the experimental datasets.
	#
iKnockLengthList = c ( 8 );
for ( iKnockLength in iKnockLengthList ) { 

	fDir = paste  ("E:/NMLab/Working/S.45.7B.Fix2.Knockout.", iKnockLength, "/", sep="" );
	#dirRData = paste ( getwd(), "/", sep="" );
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
		save( base, file = paste ( "Base", iKnockLength, "bin", sep="." ) );
		save( ctl.E, file = paste ( "Control", "E", iKnockLength, "bin", sep="." ) );
		save( ctl.I, file = paste ( "Control", "I", iKnockLength, "bin", sep="." ) );
	} else {
		load( file = paste ( dirRData, paste ( "Base", iKnockLength, "bin", sep="." ), sep="/" ) );
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

		titleText = paste ( "Control E-Only Knockout: ", iKnockLength, sep="" );
		fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( ctl.E$v1e[,indexToRFProbeRawData], ctl.E$v1i[,indexToRFProbeRawData], ctl.E$v0[,indexToRFProbeRawData] );
		KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );

		titleText = paste ( "Control I-Only Knockout: ", iKnockLength, sep="" );
		fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( ctl.I$v1e[,indexToRFProbeRawData], ctl.I$v1i[,indexToRFProbeRawData], ctl.I$v0[,indexToRFProbeRawData] );
		KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );

	} # for ( iProbeCell in ...

	rm ( ctl.E );
	rm ( ctl.I );

	if ( initialPass ) {
		exp.E = ExtractRFProbeRawData ( fName.exp.E, N2, numValsPerRFTrial, iCellList );
		exp.I = ExtractRFProbeRawData ( fName.exp.I, N2, numValsPerRFTrial, iCellList );
		save( exp.E, file = paste ( "OneSide", "E", iKnockLength, "bin", sep="." ) );
		save( exp.I, file = paste ( "OneSide", "I", iKnockLength, "bin", sep="." ) );
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

		titleText = paste ( "One-Sided E-Only Knockout: ", iKnockLength, sep="" );
		fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( exp.E$v1e[,indexToRFProbeRawData], exp.E$v1i[,indexToRFProbeRawData], exp.E$v0[,indexToRFProbeRawData] );
		KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );

		titleText = paste ( "One-Sided I-Only Knockout: ", iKnockLength, sep="" );
		fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( exp.I$v1e[,indexToRFProbeRawData], exp.I$v1i[,indexToRFProbeRawData], exp.I$v0[,indexToRFProbeRawData] );
		KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );

	} # for ( iProbeCell in ...

	rm ( base );
	rm ( exp.E );
	rm ( exp.I );

} # for ( iKnockLength in iKnockLengthlist ) {



