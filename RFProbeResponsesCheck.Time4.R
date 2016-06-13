#
#	RFProbeResponsesCheck.Time4.R
#
#	Difference from RFProbeResponsesCheck.Time3.R is purely additive: add spatial maps of
#		statistics of membrane potential response.
#

#	Clear the workspace.
rm(list = ls());

#

source ( "NMHelperFunctions.R" );
source ( "ScratchPad.R")

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
preStimIters = 100;
stimDurIters = 50;
tiffFlag = FALSE;
build.base = build.ctl.E = build.ctl.I = build.ctl.EI = TRUE;
build.exp.E = build.exp.I = build.exp.EI = TRUE;

iShortPlot = TRUE;

	#
	#	Identify the pool of recording sites of potential interest.
	#
iShowTheseRows = seq ( 18, 32 );
iCellList = c ( 	7*N+iShowTheseRows, 8*N+iShowTheseRows,
			9*N+iShowTheseRows, 10*N+iShowTheseRows,
			13*N+iShowTheseRows, 14*N+iShowTheseRows,
			15*N+iShowTheseRows, 16*N+iShowTheseRows );

	#
	#	Identify the specific RF probes of interest.
	#
iSetProbeList = c ( 382, 652 );

	#
	#	Set up titles
	#
subTitleTextRoot = NULL; xlabTextRoot = "Msec."; ylabTextRoot = "Vs(t), Vi(t)";
eScale = 1.0;

	#
	#	Identify the reference and the experimental datasets.
	#
iKnockOffset = round( ( N / 2 ), 0 ) - 4;
iKnockLengthList = c ( 4 );

for ( iKnockLength in iKnockLengthList ) { 

		#	Derive the coordinates to draw the knockout zone
	xKnockBox = c ( iKnockOffset + 1, iKnockOffset + iKnockLength );

	#fDir = paste  ("D:/NMLab/Working/S.45.7.Lesion.", iKnockLength, "/", sep="" );
	fDir = paste  ("F:/NMLab/Working/S.45.7.Lesion.x0.2.", iKnockLength, "/", sep="" );
	dirRData = fDir;

		#
		#	Set up final final names and read the data.
		#
	fName.base = "D:/NMLab/S.45.7.Knockout.Placebo/BorderKnockout_Control_Placebo.Raw.4.0.bin";

	fName.ctl.I = paste ( "BorderKnockout_Control_I.RFMap.Raw", iKnockLength, 1, "bin", sep="." );
	fName.ctl.E = paste ( "BorderKnockout_Control_E.RFMap.Raw", iKnockLength, 2, "bin", sep="." );
	fName.ctl.EI = paste ( "BorderKnockout_Control_EI.RFMap.Raw", iKnockLength, 3, "bin", sep="." );

	fName.exp.I = paste ( "BorderKnockout_OneSide_I.RFMap.Raw", iKnockLength, 4, "bin", sep="." );
	fName.exp.E = paste ( "BorderKnockout_OneSide_E.RFMap.Raw", iKnockLength, 5, "bin", sep="." );
	fName.exp.EI = paste ( "BorderKnockout_OneSide_EI.RFMap.Raw", iKnockLength, 6, "bin", sep="." );

	fName.ctl.E = paste ( fDir, fName.ctl.E, sep="" );
	fName.ctl.I = paste ( fDir, fName.ctl.I, sep="" );
	fName.ctl.EI = paste ( fDir, fName.ctl.EI, sep="" );

	fName.exp.E = paste ( fDir, fName.exp.E, sep="" );
	fName.exp.I = paste ( fDir, fName.exp.I, sep="" );
	fName.exp.EI = paste ( fDir, fName.exp.EI, sep="" );

	if ( build.base ) {
		base = ExtractRFProbeRawData ( fName.base, N2, numValsPerRFTrial, iCellList );
		save( base, file = paste( dirRData, paste ( "Base", iKnockLength, "bin", sep="." ), sep="/" ) );
		build.base = FALSE;
	} else {
		load( file = paste( dirRData, paste ( "Base", iKnockLength, "bin", sep="." ), sep="/" ) );
	} # if ( build.base)

		#
		#	CONTROL (WITHIN REPRESENTATION) Time series and spatial difference plots.
		#
	yKnockBox = c ( round((N/6),0), round((N/6),0)+1 );		# Aids in drawing knockout zone.
	if ( build.ctl.E ) {
		ctl.E = ExtractRFProbeRawData ( fName.ctl.E, N2, numValsPerRFTrial, iCellList );
		save( ctl.E, file = paste ( dirRData, paste ( "Control", "E", iKnockLength, "bin", sep="." ), sep="/" ) );
		build.ctl.E = FALSE;
	} else {
		load( file = paste ( dirRData, paste ( "Control", "E", iKnockLength, "bin", sep="." ), sep="/" ) );
	} # if ( !initialPass )
	for ( iProbe in 1:length(iSetProbeList) ) {

		iProbeCell = iSetProbeList[iProbe];
		indexToRFProbeRawData = which ( iCellList == iProbeCell );
		alldata.base = c ( base$v1e[,indexToRFProbeRawData], base$v1i[,indexToRFProbeRawData], base$v0[,indexToRFProbeRawData] );
		titleText = paste ( "Control E-Only Knockout ", iKnockLength, sep="" );
		fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( ctl.E$v1e[,indexToRFProbeRawData], ctl.E$v1i[,indexToRFProbeRawData], ctl.E$v0[,indexToRFProbeRawData] );
		KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );
		TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase,
					xKnockBox, yKnockBox );
	} # for ( iProbe in 1:length(iSetProbeList) ) {
	rm ( ctl.E );
	
	if ( build.ctl.I ) {
		ctl.I = ExtractRFProbeRawData ( fName.ctl.I, N2, numValsPerRFTrial, iCellList );
		save( ctl.I, file = paste ( dirRData, paste ( "Control", "I", iKnockLength, "bin", sep="." ), sep="/" ) );
		build.ctl.I = FALSE;
	} else {
		load( file = paste ( dirRData, paste ( "Control", "I", iKnockLength, "bin", sep="." ), sep="/" ) );
	} # if ( !initialPass )
	for ( iProbe in 1:length(iSetProbeList) ) {
		iProbeCell = iSetProbeList[iProbe];
		indexToRFProbeRawData = which ( iCellList == iProbeCell );
		alldata.base = c ( base$v1e[,indexToRFProbeRawData], base$v1i[,indexToRFProbeRawData], base$v0[,indexToRFProbeRawData] );
		titleText = paste ( "Control I-Only Knockout ", iKnockLength, sep="" );
		fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( ctl.I$v1e[,indexToRFProbeRawData], ctl.I$v1i[,indexToRFProbeRawData], ctl.I$v0[,indexToRFProbeRawData] );
		KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );
		TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase,
					xKnockBox, yKnockBox );
	} # 	for ( iProbe in 1:length(iSetProbeList) ) {	
	rm ( ctl.I );

	if ( build.ctl.EI ) {
		ctl.EI = ExtractRFProbeRawData ( fName.ctl.EI, N2, numValsPerRFTrial, iCellList );
		save( ctl.EI, file = paste ( dirRData, paste ( "Control", "EI", iKnockLength, "bin", sep="." ), sep="/" ) );
		build.ctl.EI = FALSE;
	} else {
		load( file = paste ( dirRData, paste ( "Control", "EI", iKnockLength, "bin", sep="." ), sep="/" ) );
	} # if ( !initialPass )
	for ( iProbe in 1:length(iSetProbeList) ) {	
		iProbeCell = iSetProbeList[iProbe];
		indexToRFProbeRawData = which ( iCellList == iProbeCell );
		alldata.base = c ( base$v1e[,indexToRFProbeRawData], base$v1i[,indexToRFProbeRawData], base$v0[,indexToRFProbeRawData] );
		titleText = paste ( "Control E+I Knockout ", iKnockLength, sep="" );
		fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( ctl.EI$v1e[,indexToRFProbeRawData], ctl.EI$v1i[,indexToRFProbeRawData], ctl.EI$v0[,indexToRFProbeRawData] );
		KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );
		TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase,
					xKnockBox, yKnockBox );
	} # 	for ( iProbe in 1:length(iSetProbeList) ) {	
	rm ( ctl.EI );

		#
		#	EXPERIMENTAL (@ BOUNDARY REPRESENTATION) Time series and spatial difference plots.
		#
	yKnockBox = c ( round((N/3),0)-1, round((N/3),0)-0 );
	if ( build.exp.E ) {
		exp.E = ExtractRFProbeRawData ( fName.exp.E, N2, numValsPerRFTrial, iCellList );
		save( exp.E, file = paste ( dirRData, paste ( "OneSide", "E", iKnockLength, "bin", sep="." ), sep="/" ) );
		build.exp.E = FALSE;
	} else {
		load( file = paste ( dirRData, paste ( "OneSide", "E", iKnockLength, "bin", sep="." ), sep="/" ) );
	} # 	if ( initialPass ) {
	for ( iProbe in 1:length(iSetProbeList) ) {
	
		iProbeCell = iSetProbeList[iProbe];
		indexToRFProbeRawData = which ( iCellList == iProbeCell );
		alldata.base = c ( base$v1e[,indexToRFProbeRawData], base$v1i[,indexToRFProbeRawData], base$v0[,indexToRFProbeRawData] );
		titleText = paste ( "One-Sided E-Only Knockout ", iKnockLength, sep="" );
		fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( exp.E$v1e[,indexToRFProbeRawData], exp.E$v1i[,indexToRFProbeRawData], exp.E$v0[,indexToRFProbeRawData] );
		KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );		
		TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase,
					xKnockBox, yKnockBox );
	} #for ( iProbe in 1:length(iSetProbeList) ) {
	rm ( exp.E );

	if ( build.exp.I ) {
		exp.I = ExtractRFProbeRawData ( fName.exp.I, N2, numValsPerRFTrial, iCellList );
		save( exp.I, file = paste ( dirRData, paste ( "OneSide", "I", iKnockLength, "bin", sep="." ), sep="/" ) );
		build.exp.I = FALSE;
	} else {
		load( file = paste ( dirRData, paste ( "OneSide", "I", iKnockLength, "bin", sep="." ), sep="/" ) );
	} # 	if ( initialPass ) {
	for ( iProbe in 1:length(iSetProbeList) ) {
		iProbeCell = iSetProbeList[iProbe];
		indexToRFProbeRawData = which ( iCellList == iProbeCell );
		alldata.base = c ( base$v1e[,indexToRFProbeRawData], base$v1i[,indexToRFProbeRawData], base$v0[,indexToRFProbeRawData] );
		titleText = paste ( "One-Sided I-Only Knockout ", iKnockLength, sep="" );
		fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( exp.I$v1e[,indexToRFProbeRawData], exp.I$v1i[,indexToRFProbeRawData], exp.I$v0[,indexToRFProbeRawData] );
		KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );
		TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase,
					xKnockBox, yKnockBox );
	} # 	for ( iProbe in 1:length(iSetProbeList) ) {
	rm ( exp.I );

	if ( build.exp.EI ) {
		exp.EI = ExtractRFProbeRawData ( fName.exp.EI, N2, numValsPerRFTrial, iCellList );
		save( exp.EI, file = paste ( dirRData, paste ( "OneSide", "EI", iKnockLength, "bin", sep="." ), sep="/" ) );
		build.exp.EI = FALSE;
	} else {
		load( file = paste ( dirRData, paste ( "OneSide", "EI", iKnockLength, "bin", sep="." ), sep="/" ) );
	} # 	if ( initialPass ) {
	for ( iProbe in 1:length(iSetProbeList) ) {
		iProbeCell = iSetProbeList[iProbe];
		indexToRFProbeRawData = which ( iCellList == iProbeCell );
		alldata.base = c ( base$v1e[,indexToRFProbeRawData], base$v1i[,indexToRFProbeRawData], base$v0[,indexToRFProbeRawData] );
		titleText = paste ( "One-Sided E+I Knockout ", iKnockLength, sep="" );
		fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
		alldata.exp = c ( exp.EI$v1e[,indexToRFProbeRawData], exp.EI$v1i[,indexToRFProbeRawData], exp.EI$v0[,indexToRFProbeRawData] );
		KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );
		TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase,
					xKnockBox, yKnockBox );
	} # 	for ( iProbe in 1:length(iSetProbeList) ) {
	rm ( exp.EI );
	rm ( base );

	build.base = build.ctl.E = build.ctl.I = build.ctl.EI = TRUE;
	build.exp.E = build.exp.I = build.exp.EI = TRUE;

} # for ( iKnockLength in iKnockLengthlist ) {



