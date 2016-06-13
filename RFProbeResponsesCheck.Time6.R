#
#	RFProbeResponsesCheck.Time6.R
#
#	Difference from .5 is to add a new type of figure that consolidates
#	the variations in m.p. responses across conditions.
#
#	Difference from RFProbeResponsesCheck.Time4.R is that the source data is coming from files
#		that contain one RF probe trial only.
#
#
#	Difference from RFProbeResponsesCheck.Time3.R is purely additive: add spatial maps of
#		statistics of membrane potential response.
#

#	Clear the workspace.
rm(list = ls());

#

source ( "NMHelperFunctions.R" );
source ( "ScratchPad2.R")

	#
	#	Additional helper functions.
	#

	#
	#	Set some constants
	#
N = 75;
N2 = N * N;
deltaT = 0.001;
deltaTInMsec = deltaT*1000;
trialDuration = 0.500;
numItersPerTrial = ceiling ( trialDuration / deltaT );
numValsPerRFTrial = numItersPerTrial * N * N;
preStimIters = ceiling ( 0.1/deltaT );
stimDurIters = ceiling ( 0.05/deltaT );
tiffFlag = FALSE;
iShortPlot = FALSE;

if ( N == 75 ) {
	rowSeedShowMap = 25;
	colSeedShowMap = 5;
	subNShowMap = 25;
} else {
	rowSeedShowMap = 5;
	colSeedShowMap = 5;
	subNShowMap = 35;
} # if ( N == 75 )

	#
	#	Identify the specific RF probes of interest.
	#
iProbeCell = 1011;
iProbeCell = 1088-75;

	#
	#	Set up titles
	#
subTitleTextRoot = NULL; xlabTextRoot = "Msec."; ylabTextRoot = "V(t)";
eScale = 1.0;

	#
	#	Identify the reference and the experimental datasets.
	#
iKnockOffset = round( ( N / 2 ), 0 ) - 4;
iKnockLength = c ( 8 );

		#	Derive the coordinates to draw the knockout zone
	xKnockBox = c ( iKnockOffset + 1, iKnockOffset + iKnockLength );

	fDir = paste  ("D:/NMLab/Working/T.75.1.7.Lesion.x-2.0.", iKnockLength, "/", "/BR5/CellDeepDive", "/", sep="" );
	dirRData = fDir;

		#
		#	Set up final final names and read the data.
		#
	fName.base = paste ( "BorderKnockout_Control_Placebo.Raw", iKnockLength, 0, iProbeCell-1, "bin", sep="." );

	fName.ctl.I = paste ( "BorderKnockout_Control_I.RFMap.Raw", iKnockLength, 1, iProbeCell-1, "bin", sep="." );
	fName.ctl.E = paste ( "BorderKnockout_Control_E.RFMap.Raw", iKnockLength, 2, iProbeCell-1, "bin", sep="." );
	fName.ctl.EI = paste ( "BorderKnockout_Control_EI.RFMap.Raw", iKnockLength, 3, iProbeCell-1, "bin", sep="." );
	fName.exp.I = paste ( "BorderKnockout_OneSide_I.RFMap.Raw", iKnockLength, 4, iProbeCell-1, "bin", sep="." );
	fName.exp.E = paste ( "BorderKnockout_OneSide_E.RFMap.Raw", iKnockLength, 5, iProbeCell-1, "bin", sep="." );
	fName.exp.EI = paste ( "BorderKnockout_OneSide_EI.RFMap.Raw", iKnockLength, 6, iProbeCell-1, "bin", sep="." );
	fName.exp.2.I = paste ( "BorderKnockout_TwoSide_I.RFMap.Raw", iKnockLength, 7, iProbeCell-1, "bin", sep="." );
	fName.exp.2.E = paste ( "BorderKnockout_TwoSide_E.RFMap.Raw", iKnockLength, 8, iProbeCell-1, "bin", sep="." );
	fName.exp.2.EI = paste ( "BorderKnockout_TwoSide_EI.RFMap.Raw", iKnockLength, 9, iProbeCell-1, "bin", sep="." );

	fName.base = paste ( fDir, fName.base, sep="" );
	fName.ctl.E = paste ( fDir, fName.ctl.E, sep="" );
	fName.ctl.I = paste ( fDir, fName.ctl.I, sep="" );
	fName.ctl.EI = paste ( fDir, fName.ctl.EI, sep="" );
	fName.exp.E = paste ( fDir, fName.exp.E, sep="" );
	fName.exp.I = paste ( fDir, fName.exp.I, sep="" );
	fName.exp.EI = paste ( fDir, fName.exp.EI, sep="" );
	fName.exp.2.E = paste ( fDir, fName.exp.2.E, sep="" );
	fName.exp.2.I = paste ( fDir, fName.exp.2.I, sep="" );
	fName.exp.2.EI = paste ( fDir, fName.exp.2.EI, sep="" );

	base = ExtractRFProbeRawDataSingleShot ( fName.base, N2, numValsPerRFTrial, iProbeCell );
	alldata.base = c ( base$v1e, base$v1i, base$v0 );

		#
		#	CONTROL (WITHIN REPRESENTATION) Time series and spatial difference plots.
		#
	yKnockBox = c ( ceiling(N/6), ceiling(N/6)+1 );		# Aids in drawing knockout zone.
	ctl.E = ExtractRFProbeRawDataSingleShot ( fName.ctl.E, N2, numValsPerRFTrial, iProbeCell );
	titleText.exp = paste ( " W REP E-Only SIL ", iKnockLength, sep="" );
	fNameBase.exp = paste ( titleText.exp, iProbeCell, iKnockLength, "bin", sep="." );
	alldata.exp = c ( ctl.E$v1e, ctl.E$v1i, ctl.E$v0 );

	ctl.I = ExtractRFProbeRawDataSingleShot ( fName.ctl.I, N2, numValsPerRFTrial, iProbeCell );
	titleText.exp2 = paste ( "W REP I-Only SIL ", iKnockLength, sep="" );
	fNameBase.exp2 = paste ( titleText.exp2, iProbeCell, iKnockLength, "bin", sep="." );
	alldata.exp2 = c ( ctl.I$v1e, ctl.I$v1i, ctl.I$v0 );

	ctl.EI = ExtractRFProbeRawDataSingleShot ( fName.ctl.EI, N2, numValsPerRFTrial, iProbeCell );
	titleText.exp3 = paste ( "W REP E+I SIL ", iKnockLength, sep="" );
	fNameBase.exp3 = paste ( titleText.exp3, iProbeCell, iKnockLength, "bin", sep="." );
	alldata.exp3 = c ( ctl.EI$v1e, ctl.EI$v1i, ctl.EI$v0 );

	KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText.exp, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase.exp, iShortPlot );

	KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp2, eScale,
					titleText.exp2, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase.exp2, iShortPlot );

	KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp3, eScale,
					titleText.exp3, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase.exp3, iShortPlot );

	titleText = "Differential Response to WR Silencing";
	ConsolidatedKnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell,
					alldata.base, alldata.exp, alldata.exp2, alldata.exp3, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, deltaTInMsec );	

	TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText.exp, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase.exp,
					xKnockBox, yKnockBox, rowSeedShowMap, colSeedShowMap, subNShowMap );

	TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp2,
					titleText.exp2, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase.exp2,
					xKnockBox, yKnockBox, rowSeedShowMap, colSeedShowMap, subNShowMap );

	TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp3,
					titleText.exp3, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase.exp3,
					xKnockBox, yKnockBox, rowSeedShowMap, colSeedShowMap, subNShowMap );

		#
		#	EXPERIMENTAL (@ BOUNDARY REPRESENTATION) Time series and spatial difference plots.
		#
	yKnockBox = c ( round((N/3),0)-1, round((N/3),0)+0 );
	exp.E = ExtractRFProbeRawDataSingleShot ( fName.exp.E, N2, numValsPerRFTrial, iProbeCell );
	titleText = paste ( "BRDR 1-Side E-Only SIL ", iKnockLength, sep="" );
	fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
	alldata.exp = c ( exp.E$v1e, exp.E$v1i, exp.E$v0 );
	KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );
	TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase,
					xKnockBox, yKnockBox, rowSeedShowMap, colSeedShowMap, subNShowMap );

	exp.I = ExtractRFProbeRawDataSingleShot ( fName.exp.I, N2, numValsPerRFTrial, iProbeCell );
	titleText = paste ( "BRDR 1-Side I-Only SIL ", iKnockLength, sep="" );
	fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
	alldata.exp = c ( exp.I$v1e, exp.I$v1i, exp.I$v0 );
	KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );
	TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase,
					xKnockBox, yKnockBox, rowSeedShowMap, colSeedShowMap, subNShowMap );

	exp.EI = ExtractRFProbeRawDataSingleShot ( fName.exp.EI, N2, numValsPerRFTrial, iProbeCell );
	titleText = paste ( "BRDR 1-Side E+I SIL ", iKnockLength, sep="" );
	fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
	alldata.exp = c ( exp.EI$v1e, exp.EI$v1i, exp.EI$v0 );
	KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );
	TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase,
					xKnockBox, yKnockBox, rowSeedShowMap, colSeedShowMap, subNShowMap );

		#
		#	EXPERIMENTAL (2-Side @ BOUNDARY REPRESENTATION) Time series and spatial difference plots.
		#
	yKnockBox = c ( round((N/3),0)-1, round((N/3),0)+2 );
	exp.2.E = ExtractRFProbeRawDataSingleShot ( fName.exp.2.E, N2, numValsPerRFTrial, iProbeCell );
	titleText = paste ( "BRDR 2-Side E-Only SIL ", iKnockLength, sep="" );
	fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
	alldata.exp = c ( exp.2.E$v1e, exp.2.E$v1i, exp.2.E$v0 );
	KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );
	TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase,
					xKnockBox, yKnockBox );

	exp.2.I = ExtractRFProbeRawDataSingleShot ( fName.exp.2.I, N2, numValsPerRFTrial, iProbeCell );
	titleText = paste ( "BRDR 2-Side I-Only SIL ", iKnockLength, sep="" );
	fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
	alldata.exp = c ( exp.2.I$v1e, exp.2.I$v1i, exp.2.I$v0 );
	KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );
	TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase,
					xKnockBox, yKnockBox );

	exp.2.EI = ExtractRFProbeRawDataSingleShot ( fName.exp.2.EI, N2, numValsPerRFTrial, iProbeCell );
	titleText = paste ( "BRDR 2-Side E+I SIL ", iKnockLength, sep="" );
	fNameBase = paste ( titleText, iProbeCell, iKnockLength, "bin", sep="." );
	alldata.exp = c ( exp.2.EI$v1e, exp.2.EI$v1i, exp.2.EI$v0 );
	KnockoutTimeSeriesPlotDriver ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot );
	TSSTStatsDriver ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp,
					titleText, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase,
					xKnockBox, yKnockBox );
