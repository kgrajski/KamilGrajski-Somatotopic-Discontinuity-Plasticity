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
		txt = "Connection Type Adaptation ON: E0 EE EI IE";
	} else if ( ix == 1 ) {
		txt = "Connection Type Adaptation OFF: EI";
	} else if ( ix == 2 ) {
		txt = "Connection Type Adaptation OFF: IE";
	} else if ( ix == 3 ) {
		txt = "Connection Type Adaptation OFF: EI IE";
	} else if ( ix == 4 ) {
		txt = "Connection Type Adaptation OFF: E0";
	} else if ( ix == 5 ) {
		txt = "Connection Type Adaptation OFF: EE";
	} else {
		txt = "Connection Type Adaptation OFF: E0 EE";
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

tiffFlag = FALSE;

iIndex = 0;
fDir = paste ( "E:/NMLab/Working/T.45.7", iIndex, sep="." );
fRoot = "Base.RFMap";
fileRootName = paste(fDir, fRoot, sep="//");
base = GetRFMapData2 ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iTop.base = length ( base$rfTrackData.e );

rTopo.e = matrix ( 0, nrow=N2, ncol=8 );
rTopo.i = matrix ( 0, nrow=N2, ncol=8 );

rTopo.e[,1] = TopoQuant3 ( N, base$rfTrackData.e[[1]], 0 );
rTopo.i[,1] = TopoQuant3 ( N, base$rfTrackData.i[[1]], 0 );

rTopo.e[,2] = TopoQuant3 ( N, base$rfTrackData.e[[iTop.base]], 0 );
rTopo.i[,2] = TopoQuant3 ( N, base$rfTrackData.i[[iTop.base]], 0 );

base.q.e = c ( mean ( rTopo.e[,2] ), sqrt ( var ( rTopo.e[,2] ) ) );
base.q.i = c ( mean ( rTopo.i[,2] ), sqrt ( var ( rTopo.i[,2] ) ) );

for ( iIndex in seq ( 1, 6, 1 ) ) {

	fDir = paste ( "E:/NMLab/Working/T.45.7", iIndex, sep="." );

		#
		#	Get the Experimental Refinement data.
		#
	fRoot = "Base.RFMap";
	fileRootName = paste(fDir, fRoot, sep="//");
	exp = GetRFMapData2 ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
	iTop.exp = length ( exp$rfTrackData.e );

	rTopo.e[,2+iIndex] = TopoQuant3 ( N, exp$rfTrackData.e[[iTop.exp]], 0 );
	rTopo.i[,2+iIndex] = TopoQuant3 ( N, exp$rfTrackData.i[[iTop.exp]], 0 );

	exp.q.e = c ( mean ( rTopo.e[,2+iIndex] ), sqrt ( var ( rTopo.e[,2+iIndex] ) ) );
	exp.q.i = c ( mean ( rTopo.i[,2+iIndex] ), sqrt ( var ( rTopo.i[,2+iIndex] ) ) );

		#
		#	Figure 1a.  Plot of E-Cell Centroids
		#

	if ( tiffFlag ) {
		tiff ( paste("WghtAdapt", "E", iIndex, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )

	par ( mfrow = c ( 2, 1 ) );
	titleText = paste ( "Receptive Field Centroids", paste ( paste ( "E-Type", " Baseline Iter", iBase, sep=" " ) ), sep=" " );
	subTitleText = GetPlasticityTitleText  ( 0 );
	topoMeasureText = paste ( "Topog. Measure (Mean & SD)=(", format ( base.q.e[1], digits=3 ), ", ", format ( base.q.e[2], digits=3 ), ")", sep="" );
	textText = paste ( titleText, subTitleText, topoMeasureText, sep="\n" );
	ShowTopoMap1 ( base$rfTrackData.e[[iTop.base]], textText, FALSE, 0.5, 0 );

	titleText = paste ( "Receptive Field Centroids", paste ( paste ( "E-Type", "Iter", iBase, sep=" " ) ), sep=" " );
	subTitleText = GetPlasticityTitleText  ( iIndex );
	topoMeasureText = paste ( "Topog. Measure (Mean & SD)=(", format ( exp.q.e[1], digits=3 ), ", ", format ( exp.q.e[2], digits=3 ), ")", sep="" );
	textText = paste ( titleText, subTitleText, topoMeasureText, sep="\n" );
	ShowTopoMap1 ( exp$rfTrackData.e[[iTop.base]], textText, FALSE, 0.5, 0 );

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

		#
		#	Figure 1a.  Plot of I-Cell Centroids
		#

	if ( tiffFlag ) {
		tiff ( paste("WghtAdapt", "I", iIndex, "tiff", sep="."), compression="lzw", units="in", width=7.0, height=7.0, res=300);
	} else {
		x11();
	} # if ( tiffFlag )

	par ( mfrow = c ( 2, 1 ) );
	titleText = paste ( "Receptive Field Centroids", paste ( paste ( "I-Type", " Baseline Iter", iBase, sep=" " ) ), sep=" " );
	subTitleText = GetPlasticityTitleText  ( 0 );
	topoMeasureText = paste ( "Topog. Measure (Mean & SD)=(", format ( base.q.i[1], digits=3 ), ", ", format ( base.q.i[2], digits=3 ), ")", sep="" );
	textText = paste ( titleText, subTitleText, topoMeasureText, sep="\n" );
	ShowTopoMap1 ( base$rfTrackData.i[[iTop.base]], textText, FALSE, 0.5, 0 );

	titleText = paste ( "Receptive Field Centroids", paste ( paste ( "I-Type", "Iter", iBase, sep=" " ) ), sep=" " );
	subTitleText = GetPlasticityTitleText  ( iIndex );
	topoMeasureText = paste ( "Topog. Measure (Mean & SD)=(", format ( exp.q.i[1], digits=3 ), ", ", format ( exp.q.i[2], digits=3 ), ")", sep="" );
	textText = paste ( titleText, subTitleText, topoMeasureText, sep="\n" );
	ShowTopoMap1 ( exp$rfTrackData.i[[iTop.base]], textText, FALSE, 0.5, 0 );

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

} # for ( iIndex in seq ( 0, 6, 1 )

	#
	#	Do a boxplot for the rTopo distributions.
	#
x11(); par(mfrow=c(2,1));
boxplot ( rTopo.e );
boxplot ( rTopo.i );


















