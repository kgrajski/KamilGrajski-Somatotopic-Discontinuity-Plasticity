#
#	Test Script to confirm reading of a binary file written by a C++ program.
#

#
#	This file has evolved into more of a grab bag.  Testing this and that.
#

fDir = "D:/NMLab/S.21.4.Control/";
fName = "SysInit.25.0.bin";
fName = paste ( fDir, fName, sep="" );

N = 15;
N2 = N^2;
lenPerTrial = N2 * 1000;

finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );

x = matrix ( alldata, nrow=lenPerTrial, ncol=6, byrow=FALSE );
x11(); par(mfrow=c(3,2));
for ( i in 1:6 ) {
	plot ( x[seq(lenPerTrial-2*N2, lenPerTrial),i], type="l" );
	abline(h=c(-0.01, 0.01), lty=3, col=2)
}


fName = "Base.25.0.bin";
lenPerTrial = 225 * 150;

finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );

x = matrix ( alldata, nrow=lenPerTrial, ncol=225, byrow=FALSE );
x11(); par(mfrow=c(3,2));
iList = c ( 1, 40, 100, 160, 200, 225 );
for ( i in iList ) {
	plot ( x[i,], type="l" );
}

	#	Convert the C++ rfMapRes to R format

iN = 225
rfMapRaw = list();
for ( i in 1:225 ) {

	tmp = matrix ( x[,i], nrow=150, ncol=225, byrow=TRUE );
	rfMapRaw[[i]] = list ( r0=tmp, r1.e=tmp, r1.i=tmp );	
}
durRFProbeInIters = 10;
oneSecondNumIter = 100;
startTimeInIterToRFProbeInIters = 25;

	####################
	#
	#	Testing bigger network save (includes weight matrices).
	#
	####################

fName = "Base.25.0.bin";
finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );
numItersPerTrial = 150;
numValsPerTrial = 225 * numItersPerTrial;
numWeightsPerTrial = 225 * 225 * 150;
iP1Len = numValsPerTrial * 6;
ixN = 15; ixN2 = ixN * ixN;
x = matrix ( alldata[1:iP1Len], nrow=numValsPerTrial, ncol=6, byrow=FALSE );
x11(); par(mfrow=c(3,2));
for ( i in 1:6 ) {
	plot ( x[,i], type="l" );
}
w1.e.0.start = iP1Len;
w1.e.0.end = w1.e.0.start + numWeightsPerTrial;
t0.w1.e.0 = matrix ( alldata[(w1.e.0.start+1):w1.e.0.end], nrow=225, ncol=225, byrow=FALSE );
w1.e.e.start = w1.e.0.end;
w1.e.e.end = w1.e.e.start + numWeightsPerTrial;
t0.w1.e.e = matrix ( alldata[(w1.e.e.start+1):w1.e.e.end], nrow=225, ncol=225, byrow=FALSE );
w1.e.i.start = w1.e.e.end;
w1.e.i.end = w1.e.i.start + numWeightsPerTrial;
t0.w1.e.i = matrix ( alldata[(w1.e.i.start+1):w1.e.i.end], nrow=225, ncol=225, byrow=FALSE );
w1.i.e.start = w1.e.i.end;
w1.i.e.end = w1.i.e.start + numWeightsPerTrial;
t0.w1.i.e = matrix ( alldata[(w1.i.e.start+1):w1.i.e.end], nrow=225, ncol=225, byrow=FALSE );

#t0.w1.e.0 = w1.e.0; t0.w1.e.e = w1.e.e; t0.w1.i.e = w1.i.e; t0.w1.e.i = w1.e.i;

fName = "Dump2.25.0.bin";
finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );
numItersPerTrial = 150;
numValsPerTrial = 225 * numItersPerTrial;
numWeightsPerTrial = 225 * 225 * 150;
iP1Len = numValsPerTrial * 6;
ixN = 15; ixN2 = ixN * ixN;
x = matrix ( alldata[1:iP1Len], nrow=numValsPerTrial, ncol=6, byrow=FALSE );
x11(); par(mfrow=c(3,2));
for ( i in 1:6 ) {
	plot ( x[,i], type="l" );
}
w1.e.0.start = iP1Len;
w1.e.0.end = w1.e.0.start + numWeightsPerTrial;
tn.w1.e.0 = matrix ( alldata[(w1.e.0.start+1):w1.e.0.end], nrow=225, ncol=225, byrow=FALSE );
w1.e.e.start = w1.e.0.end;
w1.e.e.end = w1.e.e.start + numWeightsPerTrial;
tn.w1.e.e = matrix ( alldata[(w1.e.e.start+1):w1.e.e.end], nrow=225, ncol=225, byrow=FALSE );
w1.e.i.start = w1.e.e.end;
w1.e.i.end = w1.e.i.start + numWeightsPerTrial;
tn.w1.e.i = matrix ( alldata[(w1.e.i.start+1):w1.e.i.end], nrow=225, ncol=225, byrow=FALSE );
w1.i.e.start = w1.e.i.end;
w1.i.e.end = w1.i.e.start + numWeightsPerTrial;
tn.w1.i.e = matrix ( alldata[(w1.i.e.start+1):w1.i.e.end], nrow=225, ncol=225, byrow=FALSE );

#tn.w1.e.0 = w1.e.0; tn.w1.e.e = w1.e.e; tn.w1.i.e = w1.i.e; tn.w1.e.i = w1.e.i;

tmp = tn.w1.e.0 - t0.w1.e.0; apply( tmp, 2, sum )
tmp = tn.w1.e.i - t0.w1.e.i; apply( tmp, 2, sum )
tmp = tn.w1.i.e - t0.w1.i.e; apply( tmp, 2, sum )
tmp = tn.w1.e.e - t0.w1.e.e; apply( tmp, 2, sum )


	#	Checking weight matrix.


fName = "w1e0.beforenorm";
finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );
ixN = 15; ixN2 = ixN * ixN; ixIters = 150;
w1.e.0 = matrix ( alldata[1:(ixN2*ixN2)], nrow=225, ncol=225, byrow=FALSE );
apply(w1.e.0,2,sum)
apply(w1.e.0,1,sum)

fName = "w1e0.afternorm";
finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );
ixN = 15; ixN2 = ixN * ixN; ixIters = 150;
w1.e.0 = matrix ( alldata[1:(ixN2*ixN2)], nrow=225, ncol=225, byrow=FALSE );
apply(w1.e.0,2,sum)
apply(w1.e.0,1,sum)


#	Checking the Di3 Refinement Stimulation
source("NMHelperFunctions.R");

#	Random number generator check.
fName = "RandCheck";
finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );

#	Checking the Di3 Refinement Stimulation
source("NMHelperFunctions.R");


#	Check Syndactyly Stimulation Patterns
fName = "CheckInputStim.Baseline";
finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );
N = 30;
N2 = N * N;
numIter = 150;
numValsPerIter = N2
numValsPerTrial = N2 * numIter;
offsetPatchInIters = 30;
durationPatchInIters = 25;

numStim = length(alldata) / (N2 + numValsPerTrial);
x11();
for ( i in c(25,26,27,28) ) {

	iStart = (i - 1) * numValsPerTrial + (i - 1) * numValsPerIter + 1;
	iEnd = iStart - 1 + numValsPerTrial;
	tmp = matrix ( alldata[iStart:iEnd], nrow = numIter, ncol=numValsPerIter, byrow=TRUE );

	# if ( testVerbose ) {
		txtTitle = " ";
		par ( mfrow = c ( 2, 2 ) );
		ViewSTDynamics ( tmp, ((offsetPatchInIters)), (offsetPatchInIters+1), txtTitle );
		ViewSTDynamics ( tmp, ((offsetPatchInIters + durationPatchInIters) ),
			((offsetPatchInIters + durationPatchInIters + 1)), txtTitle );
	# } # if ( testVerbose ) {

	ix1 = iEnd + 1;
	ix2 = iEnd + numValsPerIter;
	tmp.stimcount = alldata[ix1:ix2];


	for ( j in seq(1,10000000) ) { k = 0; }

} # for ( i in 1:numStim ) {
stimCount.Baseline = tmp.stimcount;


#	Check Syndactyly Stimulation Patterns
fName = "CheckInputStim.Syndact";
finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );
numStim = 132;
N = 15;
N2 = N * N;
numIter = 150;
numValsPerIter = N2
numValsPerTrial = N2 * numIter;
offsetPatchInIters = 30;
durationPatchInIters = 25;

for ( i in 1:numStim ) {

	iStart = (i - 1) * numValsPerTrial + (i - 1) * numValsPerIter + 1;
	iEnd = iStart - 1 + numValsPerTrial;
	tmp = matrix ( alldata[iStart:iEnd], nrow = numIter, ncol=numValsPerIter, byrow=TRUE );

	# if ( testVerbose ) {
		txtTitle = " ";
		par ( mfrow = c ( 2, 2 ) );
		ViewSTDynamics ( tmp, ((offsetPatchInIters)), (offsetPatchInIters+1), txtTitle );
		ViewSTDynamics ( tmp, ((offsetPatchInIters + durationPatchInIters) ),
			((offsetPatchInIters + durationPatchInIters + 1)), txtTitle );
	# } # if ( testVerbose ) {

	ix1 = iEnd + 1;
	ix2 = iEnd + numValsPerIter;
	tmp.stimcount = alldata[ix1:ix2];

} # for ( i in 1:numStim ) {
stimCount.Syndact = tmp.stimcount;

#	Check Syndactyly Stimulation Patterns
fName = "CheckInputStim.SyndactCtl";
finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );

tmp = GetStimCountData ( fName, N2 )
numStim = 144;
N = 15;
N2 = N * N;
numIter = 150;
numValsPerIter = N2
numValsPerTrial = N2 * numIter;
offsetPatchInIters = 30;
durationPatchInIters = 25;

for ( i in 1:numStim ) {

	iStart = (i - 1) * numValsPerTrial + (i - 1) * numValsPerIter + 1;
	iEnd = iStart - 1 + numValsPerTrial;
	tmp = matrix ( alldata[iStart:iEnd], nrow = numIter, ncol=numValsPerIter, byrow=TRUE );

	# if ( testVerbose ) {
		txtTitle = " ";
		par ( mfrow = c ( 2, 2 ) );
		ViewSTDynamics ( tmp, ((offsetPatchInIters)), (offsetPatchInIters+1), txtTitle );
		ViewSTDynamics ( tmp, ((offsetPatchInIters + durationPatchInIters) ),
			((offsetPatchInIters + durationPatchInIters + 1)), txtTitle );
	# } # if ( testVerbose ) {

	ix1 = iEnd + 1;
	ix2 = iEnd + numValsPerIter;
	tmp.stimcount = alldata[ix1:ix2];

} # for ( i in 1:numStim ) {
stimCount.SyndactCtl = tmp.stimcount;


stimCount.Baseline

stimCount.Syndact

stimCount.SyndactCtl

stimCount.Control.Total = 25 * 2 * stimCount.Baseline + 10 * stimCount.SyndactCtl

stimCount.Syndact.Total = 25 * 2 * stimCount.Baseline + 10 * stimCount.Syndact

stimCount.Syndact.Total - stimCount.Control.Total

#	Checking write/read sparse matrix.

fName = "Test.25.0.w1E0";
finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );


#	Check Syndactyly Stimulation Patterns
fName = "CheckInputStim.SyndactCtl";
finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );

tmp = GetStimCountData ( "CheckInputStim.Baseline", N2 )





#	Check SelStim Stimulation Patterns
fName = "CheckInputStim.SelStim";

stimCount.Baseline = GetStimCountData ( "CheckInputStim.Baseline", N2 )
stimCount.SelStim = GetStimCountData ( "CheckInputStim.SelStim", N2 )

finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );
N = 15;
N2 = N * N;
numIter = 150;
numValsPerIter = N2
numValsPerTrial = N2 * numIter;
offsetPatchInIters = 30;
durationPatchInIters = 25;
for ( i in 1:numStim ) {

	iStart = (i - 1) * numValsPerTrial + (i - 1) * numValsPerIter + 1;
	iEnd = iStart - 1 + numValsPerTrial;
	tmp = matrix ( alldata[iStart:iEnd], nrow = numIter, ncol=numValsPerIter, byrow=TRUE );

	# if ( testVerbose ) {
		txtTitle = " ";
		par ( mfrow = c ( 2, 2 ) );
		ViewSTDynamics ( tmp, ((offsetPatchInIters)), (offsetPatchInIters+1), txtTitle );
		ViewSTDynamics ( tmp, ((offsetPatchInIters + durationPatchInIters) ),
			((offsetPatchInIters + durationPatchInIters + 1)), txtTitle );
	# } # if ( testVerbose ) {

	ix1 = iEnd + 1;
	ix2 = iEnd + numValsPerIter;
	tmp.stimcount = alldata[ix1:ix2];

} # for ( i in 1:numStim ) {
stimCount.SelStim = tmp.stimcount;





#	Check Sparse Input Pattern Representation
source("NMHelperFunctions.R");

fDir = "D:/NMLab/S.21.4.Control/";
fName = "CheckInputStim.Baseline";
fName = paste ( fDir, fName, sep="" );

finfo = file.info ( fName );
toread = file ( fName, "rb" );
alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
close ( toread );
N = 15;
patchSize = 3;
N2 = N * N;
numIter = 1000;
patchFootPrint = (patchSize + 1)^2;
numStim = length(alldata) / (patchFootPrint + N2);
stimCount.Baseline = GetStimCountData ( fName, N2 );
x11();
for ( i in 1:numStim ) {
	iStart = (i - 1) * (patchFootPrint + N2) + 1;
	iEnd = iStart + (patchFootPrint + N2) - 1;
	tmp = alldata[iStart:iEnd];
	iPatchLocs = 1 + tmp[1:patchFootPrint];
	tmp = rep ( 0, N2 );
	tmp[iPatchLocs] = 1.0;
	ShowVecAsMap ( tmp, paste("Stim=",i,sep=" ") );
	Sys.sleep ( 1.0 );
} # for ( i in 1:numStim ) {















