######################################################################################################
######################################################################################################
#
#	Long-term state topo map evolution and various quantitative measure displays.
#
#	Run once (which is slow) to gather up all of the topo map data into a single list.
#
#	That way can rerun or do various other experiments faster.
#
#######################################################################################################
#######################################################################################################

#
#
#	Optimized for SYNDACTYLY.
#
#	There are 5 stages:
#
#	Base
#	SyndactCtl
#	SyndactReleaseCtl
#	SyndactExp
#	SyndactReleaseExp
#

#	Clear the workspace.
rm(list = ls());

#	Load up some libraries and helper functions.
require(ggplot2);
require(reshape2);
require(ellipse);
require(stats4);
source ( "NMHelperFunctions.R" );

GetPlasticityTitleText = function ( ix ) {

	if ( ix == 0 ) {
		txt = "Adaptation ON: E0 EE EI IE";
	} else if ( ix == 1 ) {
		txt = "Adaptation OFF: EI";
	} else if ( ix == 2 ) {
		txt = "Adaptation OFF: IE";
	} else if ( ix == 3 ) {
		txt = "Adaptation OFF: EI IE";
	} else if ( ix == 4 ) {
		txt = "Adaptation OFF: E0";
	} else if ( ix == 5 ) {
		txt = "Adaptation OFF: EE";
	} else {
		txt = "Adaptation OFF: E0 EE";
	} # if ( ix == 0 ) {

	return ( txt );

} # GetPlasticityTitleText = function ( ix ) {

#
#	MAIN
#

		#
		#	Some additional variables needed if/when NMExpZone isn't the source.
kRFPeakToEdgeDetect = 0.5;

		#	Set the directory where the experimental data is sitting.
N = 45;
g0 = 3;

		#	Set some additional values.
N2 = N * N;
iBase = 15;
iStart = 0;
iEnd = 15;
iStepSize = 15;
lenEdgeTrim = 0;

xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

tiffFlag = TRUE;

iIndex = 0;
fDir = paste ( "E:/NMLab/Working/T.45.7", iIndex, sep="." );
fRoot = "Base.RFMap";
fileRootName = paste(fDir, fRoot, sep="//");
base = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iTop.base = length ( base$rfTrackData.e );

rDiv = matrix ( 0, nrow=N2, ncol=8 );
rDiv[,1] = IntraColumnarRFDivergence ( base$rfTrackData.e[[1]], base$rfTrackData.i[[1]] );
rDiv[,2] = IntraColumnarRFDivergence ( base$rfTrackData.e[[iTop.base]], base$rfTrackData.i[[iTop.base]] );

rfSize.e = matrix ( 0, nrow=N2, ncol=8 );
rfSize.i = matrix ( 0, nrow=N2, ncol=8 );

rfSize.e[,2] = QuantRFSize ( base$r1.e.rfMap, kRFPeakToEdgeDetect );
rfSize.i[,2] = QuantRFSize ( base$r1.i.rfMap, kRFPeakToEdgeDetect );

	#	Have to do a special get for the initial random RF size plot.
tmp = GetRFMapData2A ( fileRootName, iBase, 0, 0, iStepSize, kRFPeakToEdgeDetect, N2 );
rfSize.e[,1] = QuantRFSize ( tmp$r1.e.rfMap, kRFPeakToEdgeDetect );
rfSize.i[,1] = QuantRFSize ( tmp$r1.i.rfMap, kRFPeakToEdgeDetect );

for ( iIndex in seq ( 1, 6, 1 ) ) {

	fDir = paste ( "E:/NMLab/Working/T.45.7", iIndex, sep="." );

		#
		#	Get the Experimental Refinement data.
		#
	fRoot = "Base.RFMap";
	fileRootName = paste(fDir, fRoot, sep="//");
	exp = GetRFMapData2A ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
	iTop.exp = length ( exp$rfTrackData.e );

	rDiv[,2+iIndex] = IntraColumnarRFDivergence ( exp$rfTrackData.e[[iTop.exp]], exp$rfTrackData.i[[iTop.exp]] );

	rfSize.e[,2+iIndex] = QuantRFSize ( exp$r1.e.rfMap, kRFPeakToEdgeDetect );
	rfSize.i[,2+iIndex] = QuantRFSize ( exp$r1.i.rfMap, kRFPeakToEdgeDetect );

		#
		#	Figure 1a.  Plot of E-Cell 
		#

	if ( tiffFlag ) {
		tiff ( paste("WghtAdapt", "E", iIndex, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )

		#	Plot spatial distribution of divergence.
	zmin = min ( rDiv[,2], rDiv[,2+iIndex] ); zmax = max ( rDiv[,2], rDiv[,2+iIndex] );
	par ( mfrow = c ( 3, 2 ) );
	titleText = paste ( "IntraCol RF Div", paste ( paste ( " Baseline Iter", iBase, sep=" " ) ), sep=" " );
	subTitleText = GetPlasticityTitleText  ( 0 );
	textText = paste ( titleText, subTitleText, sep="\n" );
	ShowVecAsMap2 ( rDiv[,2], textText, xlab, ylab, zmin, zmax );

	titleText = paste ( "IntraCol RF Div", paste ( paste ( "Iter", iBase, sep=" " ) ), sep=" " );
	subTitleText = GetPlasticityTitleText  ( iIndex );
	textText = paste ( titleText, subTitleText, sep="\n" );
	ShowVecAsMap2 ( rDiv[,2+iIndex], textText, xlab, ylab, zmin, zmax );
	
		#	Plot E RF size distributions.
	zmin = min ( rfSize.e[,2], rfSize.e[,2+iIndex] );
	zmax = max ( rfSize.e[,2], rfSize.e[,2+iIndex] );
	titleText = paste ( "RF Extent", paste ( paste ( "E-Type", " Baseline Iter", iBase, sep=" " ) ), sep=" " );
	subTitleText = GetPlasticityTitleText  ( 0 );
	textText = paste ( titleText, subTitleText, sep="\n" );
	ShowVecAsMap2 ( rfSize.e[,2], textText, xlab, ylab, zmin, zmax );

	titleText = paste ( "RF Extent", paste ( paste ( "E-Type", " Baseline Iter", iBase, sep=" " ) ), sep=" " );
	subTitleText = GetPlasticityTitleText  ( iIndex );
	textText = paste ( titleText, subTitleText, sep="\n" );
	textText = paste ( titleText, subTitleText, sep="\n" );
	ShowVecAsMap2 ( rfSize.e[,2+iIndex], textText, xlab, ylab, zmin, zmax );

		#	Plot I RF size distributions.
	zmin = min ( rfSize.i[,2], rfSize.i[,2+iIndex] );
	zmax = max ( rfSize.i[,2], rfSize.i[,2+iIndex] );
	titleText = paste ( "RF Extent", paste ( paste ( "I-Type", " Baseline Iter", iBase, sep=" " ) ), sep=" " );
	subTitleText = GetPlasticityTitleText  ( 0 );
	textText = paste ( titleText, subTitleText, sep="\n" );
	ShowVecAsMap2 ( rfSize.i[,2], textText, xlab, ylab, zmin, zmax );

	titleText = paste ( "RF Extent", paste ( paste ( "I-Type", " Baseline Iter", iBase, sep=" " ) ), sep=" " );
	subTitleText = GetPlasticityTitleText  ( iIndex );
	textText = paste ( titleText, subTitleText, sep="\n" );
	textText = paste ( titleText, subTitleText, sep="\n" );
	ShowVecAsMap2 ( rfSize.i[,2+iIndex], textText, xlab, ylab, zmin, zmax );
	
	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

} # for ( iIndex in seq ( 0, 6, 1 )

	#
	#	Do a boxplot for the rDiv distributions.
	#
x11(); par(mfrow=c(3,1));
boxplot ( rDiv );
boxplot ( rfSize.e );
boxplot ( rfSize.i );


















