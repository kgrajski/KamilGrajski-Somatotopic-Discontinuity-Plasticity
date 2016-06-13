
#
#
#	The definitive generator of Selective Stimulation figures.
#	(Not necessarily in the order that they appear in the paper.)
#
#

#	Clear the workspace.
rm(list = ls());
alreadyCleared = TRUE;
source ( "NMHelperFunctions.R" );

	#
	#	Set up some parameters and directory path name.
fDir = "D:\\NMLab\\";
base1.Tmp = "FS";
fDir = paste ( fDir, base1.Tmp, sep="" );

N = 33;
N2 = N * N;
expID = 65;
expSectorType = "SN";

expSectorID = 5;
expFactor = 45;
refinementPatchSize = 3;
expGridFactor = "G3";
fDir = paste ( fDir, expSectorType, refinementPatchSize, expGridFactor , sep="." );

wName.exp = paste("SelStimExp.25.25", expSectorID, expFactor, sep="." );
wName.ctl = paste("SelStimCTL.25.25" );

source("LTTopoMapLook.CPP.SelStim.R");

source("CorticalMagnificationAnalysis.CPP.SelStim.R");

source("SpecialLTTopoCorticalMag.CPP.SelStim.R");

save.image (file="N5.G3.P3.X4.Rdat");



