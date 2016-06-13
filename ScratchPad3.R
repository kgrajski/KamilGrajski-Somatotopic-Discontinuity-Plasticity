#
#	ScracthPad3
#

RFProbeTrialMovie = function ( alldata, N2, iProbeCell, numValsPerTrial,
					titleText.rfProbeTrialMovie, xLabText, yLabText,
					movieFileName ) {

layerMarks = c(N+0.5, 2*N+0.5);
minZ = 0; maxZ = 1;
pCut = 0.01;
nCount = 0;
numColors = 128;

x11();

saveHTML ( { 
	startOffset.e = 1; startOffset.i = numValsPerRFTrial + 1;
	startOffset.0 = 2 * numValsPerRFTrial + 1;
	v1e = alldata [ (startOffset.e):(startOffset.e + numValsPerRFTrial - 1) ];
	v1i = alldata [ (startOffset.i):(startOffset.i + numValsPerRFTrial - 1) ];
	v0 = alldata [ (startOffset.0):(startOffset.0 + numValsPerRFTrial - 1) ];
	r1e = sigmoid ( v1e, 4 ); r1i = sigmoid ( v1i, 4 ); r0 = sigmoid ( v0, 4 );
	for ( iTime in seq(1,numItersPerTrial) ) {
		titleText = paste ( titleText.rfProbeTrialMovie, iProbeCell, "\nTime=", iTime*deltaT, sep="" );
		w = c ( r0[(iTime - 1) * N2 + seq ( 1:N2 )], r1i[(iTime - 1) * N2 + seq ( 1:N2 )],
				r1e[(iTime - 1) * N2 + seq ( 1:N2 )] );
		image.plot ( c(1:N), c(1:(3*N)), matrix ( w, nrow=N, ncol=3*N, byrow=FALSE ),
				col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
		zlim = c ( minZ, maxZ ), main=titleText, xlab=xLabText, ylab=yLabText );
		abline(h=layerMarks, col=1, lty=1 );
		abline(h=c(), col=1, lty=1 );
	} # for ( iTime in seq(1,numItersPerTrial) ) {
}, img.name="FR.PNG", interval = deltaT, htmlfile=paste(movieFileName, "AvgSpiking",sep="."), autoplay=FALSE );  # saveGIF


saveHTML ( { 
	startOffset.e = 1; startOffset.i = numValsPerRFTrial + 1;
	startOffset.0 = 2 * numValsPerRFTrial + 1;
	v1e = alldata [ (startOffset.e):(startOffset.e + numValsPerRFTrial - 1) ];
	v1i = alldata [ (startOffset.i):(startOffset.i + numValsPerRFTrial - 1) ];
	v0 = alldata [ (startOffset.0):(startOffset.0 + numValsPerRFTrial - 1) ];
	#minZ = min ( v1e, v1i, v0 ); maxZ = max ( v1e, v1i, v0 );
	for ( iTime in seq(1,numItersPerTrial) ) {
		titleText = paste ( titleText.rfProbeTrialMovie, iProbeCell, "\nTime=", iTime*deltaT, sep="" );
		w = c ( v0[(iTime - 1) * N2 + seq ( 1:N2 )], v1i[(iTime - 1) * N2 + seq ( 1:N2 )],
				v1e[(iTime - 1) * N2 + seq ( 1:N2 )] );
		minZ = min(w); maxZ = max(w);
		image.plot ( c(1:N), c(1:(3*N)), matrix ( w, nrow=N, ncol=3*N, byrow=FALSE ),
				col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
		zlim = c ( minZ, maxZ ), main=titleText, xlab=xLabText, ylab=yLabText );
		abline(h=layerMarks, col=1, lty=1 );
		abline(h=c(), col=1, lty=1 );
	} # for ( iTime in seq(1,numItersPerTrial) ) {
}, img.name="MP.PNG", interval = deltaT, htmlfile=paste(movieFileName, "MP",sep="."), autoplay=FALSE );  # saveGIF

} # RFProbeTrialMovie = function ( ) {

titleText.rfProbeTrialMovie = "Response at Each Layer to RF Probe # ";
yLabText="Bot: Skin Layer.  Mid: Cortical I.  Top: Cortical E.";
xLabText="Distal -> Proximal";
movieFileName = "test";
RFProbeTrialMovie ( alldata, N2, iProbeCell, numValsPerTrial, titleText.rfProbeTrialMovie, xLabText, yLabText,
					movieFileName );

