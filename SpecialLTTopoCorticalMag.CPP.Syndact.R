#
#	Special scatterplot and related cell deep-dive figures
#

#
#	Run in conjunction with CorticalMagnificationAnalysis and LTTopoMap scripts.
#	These routines do lots of the set up of the data used in the plots below.
#

#	Clear the workspace.
if ( !alreadyCleared ) { rm(list = ls()); }

	################
	#
	#	Draw a series of scatterplots across several conditions.  Pinpoint cells of interest.
	#
	################

initText = "Initial";
finalText = "Refined";
finalText2 = "Syndactyly";

x11(); par(mfrow=c(2,2));
	#
	#	Compare Initial Random with Refined
	#		E-Type
xlim.e = c ( min(init.rfMap.maxResp.e,final.rfMap.maxResp.e), max(init.rfMap.maxResp.e,final.rfMap.maxResp.e) );
ylim.e = c ( min(init.rfMap.rfAreas.e,final.rfMap.rfAreas.e), max(init.rfMap.rfAreas.e,final.rfMap.rfAreas.e) );
plot ( init.rfMap.maxResp.e, init.rfMap.rfAreas.e, xlim=xlim.e, ylim=ylim.e, pch=".", col=1, cex=2,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer E-Cells\n",initText,"(.); ",finalText,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.e, final.rfMap.rfAreas.e, pch="+", col=2, cex=1 );
for ( i in 1:length(div.base.final.cellList) ) {
	points ( final.rfMap.maxResp.e[div.base.final.cellList[i]], final.rfMap.rfAreas.e[div.base.final.cellList[i]],
		pch=symbolChoice[i], col=3, cex=2, lwd=2 );
} # for ( i in 1:length(div.base.final.cellList) ) {

	#		I-Type
xlim.i = c ( min(init.rfMap.maxResp.i,final.rfMap.maxResp.i), max(init.rfMap.maxResp.i,final.rfMap.maxResp.i) );
ylim.i = c ( min(init.rfMap.rfAreas.i,final.rfMap.rfAreas.i), max(init.rfMap.rfAreas.i,final.rfMap.rfAreas.i) );
plot ( init.rfMap.maxResp.i, init.rfMap.rfAreas.i, xlim=xlim.i, ylim=ylim.i, pch=".", col=1, cex=2,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer I-Cells\n",initText,"(.); ",finalText,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.i, final.rfMap.rfAreas.i, pch="+", col=2, cex=1 );
for ( i in 1:length(div.base.final.cellList) ) {
	points ( final.rfMap.maxResp.i[div.base.final.cellList[i]], final.rfMap.rfAreas.i[div.base.final.cellList[i]],
		pch=symbolChoice[i], col=3, cex=2, lwd=2 );
} # for ( i in 1:length(div.base.final.cellList) ) {

	#
	#	Compare Refined with Syndactyly
	#		E-Type
xlim.e = c ( min(final.rfMap.maxResp.e,final2.rfMap.maxResp.e), max(final.rfMap.maxResp.e,final2.rfMap.maxResp.e) );
ylim.e = c ( min(final.rfMap.rfAreas.e,final2.rfMap.rfAreas.e), max(final.rfMap.rfAreas.e,final2.rfMap.rfAreas.e) );
plot ( final.rfMap.maxResp.e, final.rfMap.rfAreas.e, xlim=xlim.e, ylim=ylim.e, pch=".", col=1, cex=2,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer E-Cells\n",finalText,"(.); ",finalText2,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final2.rfMap.maxResp.e, final2.rfMap.rfAreas.e, pch="+", col=2, cex=1 );
for ( i in 1:length(div.base.final.cellList) ) {
	points ( final2.rfMap.maxResp.e[div.base.final.cellList[i]], final2.rfMap.rfAreas.e[div.base.final.cellList[i]],
		pch=symbolChoice[i], col=3, cex=2, lwd=2 );
} # for ( i in 1:length(div.base.final.cellList) ) {


pchChar = 15; wtmp = which ( InGrid ( final2.rfMap.maxResp.e, final2.rfMap.rfAreas.e, 0, 22, 27, 31 ) == 1 );
xtmp = mean ( final2.rfMap.maxResp.e[wtmp] ); ytmp = mean ( final2.rfMap.rfAreas.e[wtmp]); points ( xtmp, ytmp, pch=pchChar, col=3, cex=1.5, lwd=2 );

pchChar = 16; wtmp = which ( InGrid ( final2.rfMap.maxResp.e, final2.rfMap.rfAreas.e, 20, 23, 31, 36 ) == 1 );
xtmp = mean ( final2.rfMap.maxResp.e[wtmp] ); ytmp = mean ( final2.rfMap.rfAreas.e[wtmp]); points ( xtmp, ytmp, pch=pchChar, col=3, cex=1.5, lwd=2 );

pchChar = 17; wtmp = which ( InGrid ( final2.rfMap.maxResp.e, final2.rfMap.rfAreas.e, 21, 25, 36, 40 ) == 1 );
xtmp = mean ( final2.rfMap.maxResp.e[wtmp] ); ytmp = mean ( final2.rfMap.rfAreas.e[wtmp]); points ( xtmp, ytmp, pch=pchChar, col=3, cex=1.5, lwd=2 );



	#		I-Type

xlim.i = c ( min(final.rfMap.maxResp.i,final2.rfMap.maxResp.i), max(final.rfMap.maxResp.i,final2.rfMap.maxResp.i) );
ylim.i = c ( min(final.rfMap.rfAreas.i,final2.rfMap.rfAreas.i), max(final.rfMap.rfAreas.i,final2.rfMap.rfAreas.i) );
plot ( final.rfMap.maxResp.i, final.rfMap.rfAreas.i, xlim=xlim.i, ylim=ylim.i, pch=".", col=1, cex=2,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer I-Cells\n",finalText,"(.); ",finalText2,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final2.rfMap.maxResp.i, final2.rfMap.rfAreas.i, pch="+", col=2, cex=1 );
for ( i in 1:length(div.base.final.cellList) ) {
	points ( final2.rfMap.maxResp.i[div.base.final.cellList[i]], final2.rfMap.rfAreas.i[div.base.final.cellList[i]],
		pch=symbolChoice[i], col=3, cex=2, lwd=2 );
} # for ( i in 1:length(div.base.final.cellList) ) {

	################
	#
	#	Draw a series of plots.  Deep dives on individual cell RF Extents.
	#
	################
div.base.final.cellList = 650;
for ( iCell in div.base.final.cellList ) {

	x11(); par(mfrow=c(3,2));
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

	iWhich = 1;
	titleText = paste ( "Column # ", iCell, " E Cell RF Extent\nBaseline Refinement Iter ", (iWhich-1), sep="" );
	ShowVecAsMap1 ( init.r1.e.rfMap[iCell,], titleText, xlab, ylab );

	iWhich = 1;
	titleText = paste ( "Column # ", iCell, " I Cell RF Extent\nBaseline Refinement Iter ", (iWhich-1), sep="" );
	ShowVecAsMap1 ( init.r1.i.rfMap[iCell,], titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " E Cell RF Extent\nBaseline Refinement Iter ", (iWhich-1), sep="" );
	ShowVecAsMap1 ( final.r1.e.rfMap[iCell,], titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " I Cell RF Extent\nBaseline Refinement Iter ", (iWhich-1), sep="" );
	ShowVecAsMap1 ( final.r1.i.rfMap[iCell,], titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " E Cell RF Extent\nSyndactyly Iter ", (iWhich-1), sep="" );
	ShowVecAsMap1 ( final2.r1.e.rfMap[iCell,], titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " I Cell RF Extent\nSyndactyly Iter ", (iWhich-1), sep="" );
	ShowVecAsMap1 ( final2.r1.i.rfMap[iCell,], titleText, xlab, ylab );

} # for ( iCell in div.base.final.cellList ) {

	#	Streamlined version.  Separately talking about baseline and syndactyly.
digitWidth = as.integer ( N / 3 );
div.base.final.cellList = immediate.Upper = seq ( digitWidth * N + 1 + 2, digitWidth * N + N - 2, 2 );
div.base.final.cellList = immediate.Lower = seq ( (digitWidth - 1) * N + 1 + 2, digitWidth * N - 2, 2 );

div.base.final.cellList = immediate.Upper = seq ( (digitWidth + 1) * N + 1 + 2, (digitWidth + 1) * N + N - 2, 2 );
div.base.final.cellList = immediate.Lower = seq ( (digitWidth - 2) * N + 1 + 2, (digitWidth - 1) * N - 2, 2 );

div.base.final.cellList = immediate.Upper = seq ( (digitWidth + 2) * N + 1 + 2, (digitWidth + 2) * N + N - 2, 2 );
div.base.final.cellList = immediate.Lower = seq ( (digitWidth - 3) * N + 1 + 2, (digitWidth - 2) * N - 2, 2 );

div.base.final.cellList = 688 + c ( -3, -2, -1, 0, +1, +2, +3 ) * N;

boundaryMarks = c((N/3)+0.5, (2*N/3)+0.5);

for ( iCell in div.base.final.cellList ) {

	x11(); par(mfcol=c(2,2));
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

	iWhich = 1;
	titleText = paste ( "Column # ", iCell, " E Cell RF Extent\nBaseline Refinement Iter ", (iWhich-1), sep="" );
	ShowVecAsMap1 ( init.r1.e.rfMap[iCell,], titleText, xlab, ylab );
	abline ( h = boundaryMarks, lty=3, col=3 );

	iWhich = 1;
	titleText = paste ( "Column # ", iCell, " I Cell RF Extent\nBaseline Refinement Iter ", (iWhich-1), sep="" );
	ShowVecAsMap1 ( init.r1.i.rfMap[iCell,], titleText, xlab, ylab );
	abline ( h = boundaryMarks, lty=3, col=3 );

	iWhich = iBase;
	titleText = paste ( "Column # ", iCell, " E Cell RF Extent\nBaseline Refinement Iter ", iWhich, sep="" );
	ShowVecAsMap1 ( final.r1.e.rfMap[iCell,], titleText, xlab, ylab );
	abline ( h = boundaryMarks, lty=3, col=3 );

	iWhich = iBase;
	titleText = paste ( "Column # ", iCell, " I Cell RF Extent\nBaseline Refinement Iter ", iWhich, sep="" );
	ShowVecAsMap1 ( final.r1.i.rfMap[iCell,], titleText, xlab, ylab );
	abline ( h = boundaryMarks, lty=3, col=3 );

} # for ( iCell in div.base.final.cellList ) {

	#	Streamlined version.  Separately talking about baseline and syndactyly.
digitWidth = as.integer ( N / 3 );
div.base.final.cellList = immediate.Upper = seq ( digitWidth * N + 1 + 2, digitWidth * N + N - 2, 2 );
div.base.final.cellList = immediate.Lower = seq ( (digitWidth - 1) * N + 1 + 2, digitWidth * N - 2, 2 );

div.base.final.cellList = immediate.Upper = seq ( (digitWidth + 1) * N + 1 + 2, (digitWidth + 1) * N + N - 2, 2 );
div.base.final.cellList = immediate.Lower = seq ( (digitWidth - 2) * N + 1 + 2, (digitWidth - 1) * N - 2, 2 );

div.base.final.cellList = immediate.Upper = seq ( (digitWidth + 2) * N + 1 + 2, (digitWidth + 2) * N + N - 2, 2 );
div.base.final.cellList = immediate.Lower = seq ( (digitWidth - 3) * N + 1 + 2, (digitWidth - 2) * N - 2, 2 );

div.base.final.cellList = 688 + c ( -3, -2, -1, 0, +1, +2, +3 ) * N;

boundaryMarks = c ( N/3 + 0.5, 2*N/3 + 0.5 );

for ( iCell in div.base.final.cellList ) {

	x11(); par(mfcol=c(2,2));
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

	iWhich = 1;
	titleText = paste ( "Column # ", iCell, " E Cell RF Extent\nBaseline Refinement Iter ", 15, sep="" );
	ShowVecAsMap1 ( final.r1.e.rfMap[iCell,], titleText, xlab, ylab );
	abline ( h = boundaryMarks, lty=3, col=3 );

	iWhich = 1;
	titleText = paste ( "Column # ", iCell, " I Cell RF Extent\nBaseline Refinement Iter ", 15, sep="" );
	ShowVecAsMap1 ( final.r1.i.rfMap[iCell,], titleText, xlab, ylab );
	abline ( h = boundaryMarks, lty=3, col=3 );

	iWhich = 10;
	titleText = paste ( "Column # ", iCell, " E Cell RF Extent\nSyndactyly Iter ", iWhich, sep="" );
	ShowVecAsMap1 ( final2.r1.e.rfMap[iCell,], titleText, xlab, ylab );
	abline ( h = boundaryMarks, lty=3, col=3 );

	iWhich = 10;
	titleText = paste ( "Column # ", iCell, " I Cell RF Extent\nSyndactyly Iter ", iWhich, sep="" );
	ShowVecAsMap1 ( final2.r1.i.rfMap[iCell,], titleText, xlab, ylab );
	abline ( h = boundaryMarks, lty=3, col=3 );

} # for ( iCell in div.base.final.cellList ) {


	################
	#
	#	Draw a series of plots.  Deep dive on input weight distributions.
	#
	################

	##########################
	#
	#	Zoom in.  Each adjusts to its own weight.
	#
	##########################

digitBorderline = c ( 4.5, 3.5, 0 );
iBorder = 0;

for ( iCell in div.base.final.cellList ) {

	iBorder = iBorder + 1;	# Need to show digit boundaries correctly.

	x11();
	par(mfrow=c(3,4));
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

		#	Initial Conditions

	iWhich = 1;
	titleText = paste ( "Col. # ", iCell, " E<-S Wghts\nBaseline Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( init.w1.e0, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 1, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 1;
	titleText = paste ( "Col. # ", iCell, " E<-E Wghts\nBaseline Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( init.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );

	iWhich = 1;
	titleText = paste ( "Col. # ", iCell, " E<-I Wghts\nBaseline Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( init.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );

	iWhich = 1;
	titleText = paste ( "Col. # ", iCell, " I<-E Wghts\nBaseline Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( init.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );

		#	Refined Map

	iWhich = 26;
	titleText = paste ( "Col. # ", iCell, " E<-S Wghts\nBaseline Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( refine.w1.e0, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 1, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 26;
	titleText = paste ( "Col. # ", iCell, " E<-E Wghts\nBaseline Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( refine.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Col. # ", iCell, " E<-I Wghts\nBaseline Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( refine.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Col. # ", iCell, " I<-E Wghts\nBaseline Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( refine.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );

		#	Syndactyly

	iWhich = 25;
	titleText = paste ( "Col. # ", iCell, " E<-S Wghts\nSyndactyly Iter ", (iWhich), sep="" );
	vTmp = GetInputWeights ( syndact.w1.e0, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 1, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 25;
	titleText = paste ( "Col. # ", iCell, " E<-E Wghts\nSyndactyly Iter ", (iWhich), sep="" );
	vTmp = GetInputWeights ( syndact.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );

	iWhich = 25;
	titleText = paste ( "Col. # ", iCell, " E<-I Wghts\nSyndactyly Iter ", (iWhich), sep="" );
	vTmp = GetInputWeights ( syndact.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );

	iWhich = 25;
	titleText = paste ( "Col. # ", iCell, " I<-E Wghts\nSyndactyly Iter ", (iWhich), sep="" );
	vTmp = GetInputWeights ( syndact.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ShowVecAsMap1 ( vTmp, titleText, xlab, ylab );

	for ( iDelay in seq(0, 1, 0.001) ) {
		iDelayVar = 0;
	} #	for ( iDelay in seq(0, 1, 0.0001) ) {

} # for ( iCell in div.base.final.cellList ) {



	#
	#	This section is about doing a deeper dive on some of the scatterplots.
	#
	#	This is specific to the scatterplot above for the 30.4 case.  30x30 network; (4+1)x(4+1) inputs.
	#

	#	Do some correlation analysis.

		#	E-Cell.  Init vs  Baseline
	
			#	The positively correlated cluster: cut off 

clust1 = which ( final.rfMap.maxResp.e < 25 );
clust2 = which ( final.rfMap.maxResp.e >= 25 );
final.maxResp.rfAreas.e.corr.clust1 = cor ( final.rfMap.maxResp.e[clust1], final.rfMap.rfAreas.e[clust1] );
final.maxResp.rfAreas.e.corr.clust2 = cor ( final.rfMap.maxResp.e[clust2], final.rfMap.rfAreas.e[clust2] );
final.maxResp.rfAreas.i.corr = cor ( final.rfMap.maxResp.i, final.rfMap.rfAreas.i );

clust1 = which ( final2.rfMap.maxResp.e < 25 );
clust2 = which ( final2.rfMap.maxResp.e >= 25 );
final2.maxResp.rfAreas.e.corr.clust1 = cor ( final2.rfMap.maxResp.e[clust1], final2.rfMap.rfAreas.e[clust1] );
final2.maxResp.rfAreas.e.corr.clust2 = cor ( final2.rfMap.maxResp.e[clust2], final2.rfMap.rfAreas.e[clust2] );
final2.maxResp.rfAreas.i.corr = cor ( final2.rfMap.maxResp.i, final2.rfMap.rfAreas.i );


	#	Do some mean, sd analysis and Wilcox analysis.
	
	#	E-cells maxResp and rfAreas

init.rfMap.maxResp.e.meanSD = MeanAndSD ( init.rfMap.maxResp.e );
final.rfMap.maxResp.e.meanSD = MeanAndSD ( final.rfMap.maxResp.e );
final2.rfMap.maxResp.e.meanSD = MeanAndSD ( final2.rfMap.maxResp.e );

init.rfMap.maxResp.i.meanSD = MeanAndSD ( init.rfMap.maxResp.i );
final.rfMap.maxResp.i.meanSD = MeanAndSD ( final.rfMap.maxResp.i );
final2.rfMap.maxResp.i.meanSD = MeanAndSD ( final2.rfMap.maxResp.i );

	#	I-cells maxResp and rfAreas

init.rfMap.rfAreas.e.meanSD = MeanAndSD ( init.rfMap.rfAreas.e );
final.rfMap.rfAreas.e.meanSD = MeanAndSD ( final.rfMap.rfAreas.e );
final2.rfMap.rfAreas.e.meanSD = MeanAndSD ( final2.rfMap.rfAreas.e );

init.rfMap.rfAreas.i.meanSD = MeanAndSD ( init.rfMap.rfAreas.i );
final.rfMap.rfAreas.i.meanSD = MeanAndSD ( final.rfMap.rfAreas.i );
final2.rfMap.rfAreas.i.meanSD = MeanAndSD ( final2.rfMap.rfAreas.i );

	#

baseline.maxResp.e.v.i.wilcox = wilcox.test ( final.rfMap.maxResp.e, final.rfMap.maxResp.i, exact = FALSE, correct=FALSE );
baseline.rfAreas.e.v.i.wilcox = wilcox.test ( final.rfMap.maxResp.e, final.rfMap.maxResp.i, exact = FALSE, correct=FALSE );

syndact.maxResp.e.v.i.wilcox = wilcox.test ( final2.rfMap.maxResp.e, final2.rfMap.maxResp.i, exact = FALSE, correct=FALSE );
syndact.rfAreas.e.v.i.wilcox = wilcox.test ( final2.rfMap.maxResp.e, final2.rfMap.maxResp.i, exact = FALSE, correct=FALSE );

init.v.baseline.rfAreas.e = wilcox.test ( init.rfMap.rfAreas.e, final.rfMap.rfAreas.e, exact = FALSE, correct=FALSE );
baseline.v.syndact.rfAreas.e = wilcox.test ( final.rfMap.rfAreas.e, final2.rfMap.rfAreas.e, exact = FALSE, correct=FALSE );

init.v.baseline.rfAreas.i = wilcox.test ( init.rfMap.rfAreas.i, final.rfMap.rfAreas.i, exact = FALSE, correct=FALSE );
baseline.v.syndact.rfAreas.i = wilcox.test ( final.rfMap.rfAreas.i, final2.rfMap.rfAreas.i, exact = FALSE, correct=FALSE );


	#	Analyze some of the subclusters in the scatterplots of E-cells for baseline v synactyly.
	#	What is their spatial distribution or other differentiator?

x11(); par(mfrow=c(2,2));
ShowVecAsMap ( InGrid ( final2.rfMap.maxResp.e, final2.rfMap.rfAreas.e,  0, 22, 27, 31 ), "B v S E Cluster 1" );
ShowVecAsMap ( InGrid ( final2.rfMap.maxResp.e, final2.rfMap.rfAreas.e,  20, 23, 31, 36 ), "B v S E Cluster 2" );
ShowVecAsMap ( InGrid ( final2.rfMap.maxResp.e, final2.rfMap.rfAreas.e,  21, 25, 36, 40 ), "B v S E Cluster 3" );
ShowVecAsMap ( InGrid ( final2.rfMap.maxResp.e, final2.rfMap.rfAreas.e,  25, 50, 0, 50 ), "B v S E Cluster 4" );

