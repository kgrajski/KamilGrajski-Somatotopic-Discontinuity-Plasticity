#
#	KnockoutFigure.R
#

#	Clear the workspace.
rm(list = ls());

#	Load up some libraries and helper functions.
require(ggplot2);
require(reshape2);
require(ellipse);
require(stats4);
source ( "NMHelperFunctions.R" );

#
#	MAIN
#

		#
		#	Some additional variables needed if/when NMExpZone isn't the source.
kRFPeakToEdgeDetect = 0.5;
N = 33;
N2 = N * N;


		#	Set the directory where the experimental data is sitting.
fDir = "C:/Users/KAG/Documents/NMLab/Simulations/S.2_Knockout_I_Only/"
fDir = "C:/Users/KAG/Documents/NMLab/Simulations/S2_Knockout_I_Only_One_Side/"

		#
		#	Get the Baseline Refinement RF Map.
		#
fRoot = "Base.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iBase = 10; iStart = 10; iEnd = 10; iStepSize = 10;
base = GetRFMapData2 ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iTop.base = length ( base$rfTrackData.e );

		#
		#	Get the Knockout RF Map.
		#

fRoot = "BorderKnockout.RFMap";
fileRootName = paste(fDir, fRoot, sep="\\");
iBase = 10; iStart = 0; iEnd = 0; iStepSize = 10;
test = GetRFMapData2 ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 );
iTop.test = length ( test$rfTrackData.e );

		#
		#	Figure 1: Compare Baseline & Knockout E and I RF Centroids
		#
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.base; ShowTopoMap1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

subTitleText = "\nBorder-Adjacent I-Type Knockout ";
iWhich = iTop.test; ShowTopoMap1 ( test$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );
iWhich = iTop.test; ShowTopoMap1 ( test$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), FALSE, 0.5, 0 );

		#
		#	Figure 2: Compare Baseline & Knockout E and I RF Sizes for Border Adjacent E-Cells Longitudonal Tracks.
		#
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrack1 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.base; ShowThreeDigitRFTrack1 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = "\nBorder-Adjacent I-Type Knockout";
iWhich = iTop.test; ShowThreeDigitRFTrack1 ( test$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.test; ShowThreeDigitRFTrack1 ( test$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
		
		#
		#	Figure 3: Compare Baseline & Knockout E and I RF Sizes for Border Adjacent E-Cells Longitudonal Tracks.
		#
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrack1A ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.base; ShowThreeDigitRFTrack1A ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = "\nBorder-Adjacent I-Type Knockout";
iWhich = iTop.test; ShowThreeDigitRFTrack1A ( test$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.test; ShowThreeDigitRFTrack1A ( test$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

		#
		#	Figure 4: Compare Baseline & Knockout E and I RF Sizes for Border Adjacent E-Cells Longitudonal Tracks.
		#
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrack2 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.base; ShowThreeDigitRFTrack2 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = "\nBorder-Adjacent I-Type Knockout";
iWhich = iTop.test; ShowThreeDigitRFTrack2 ( test$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.test; ShowThreeDigitRFTrack2 ( test$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

		#
		#	Figure 5: Compare Baseline & Knockout E and I RF Sizes for Border Adjacent E-Cells Longitudonal Tracks.
		#
x11(); par(mfcol=c(2,2));
subTitleText = "\nBaseline Refinement";
iWhich = iTop.base; ShowThreeDigitRFTrack3 ( base$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.base; ShowThreeDigitRFTrack3 ( base$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );

subTitleText = "\nBorder-Adjacent I-Type Knockout";
iWhich = iTop.test; ShowThreeDigitRFTrack3 ( test$rfTrackData.e[[iWhich]], paste(paste("E-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );
iWhich = iTop.test; ShowThreeDigitRFTrack3 ( test$rfTrackData.i[[iWhich]], paste(paste("I-Type","Iter",iEnd,sep=" "),subTitleText, sep=""), TRUE, 0.5, 0, 1 );














