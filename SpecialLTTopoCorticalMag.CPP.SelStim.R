#
#	Special scatterplot and related cell deep-dive figures
#

#
#	Run in conjunction with CorticalMagnificationAnalysis and LTTopoMap scripts.
#	These routines do lots of the set up of the data used in the plots below.
#

#	Clear the workspace.
if ( !alreadyCleared ) { rm(list = ls()); }

if ( 0 ) {

	################
	#
	#	Draw a series of scatterplots across several conditions.  Pinpoint cells of interest.
	#
	################

deepLookList = c ( 579, 580, 612 );

	#
	#	Compare Last Iter of Baseline ("Baseline") with Last Iter of Focal Stim
	#		E-Type
x11(); par(mfrow=c(2,2));
xlim.e = c ( min(init.rfMap.maxResp.e,final.rfMap.maxResp.e), max(init.rfMap.maxResp.e,final.rfMap.maxResp.e) );
ylim.e = c ( min(init.rfMap.rfAreas.e,final.rfMap.rfAreas.e), max(init.rfMap.rfAreas.e,final.rfMap.rfAreas.e) );
plot ( init.rfMap.maxResp.e, init.rfMap.rfAreas.e, xlim=xlim.e, ylim=ylim.e, pch=".", col=1, cex=2,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer E-Cells\n",initText,"(.); ",finalText,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.e, final.rfMap.rfAreas.e, pch="+", col=2, cex=1 );

for ( i in 1:length(deepLookList) ) {
	points ( final.rfMap.maxResp.e[deepLookList[i]], final.rfMap.rfAreas.e[deepLookList[i]],
		pch=symbolChoice[i], col=3, cex=2, lwd=2 );
} # for ( i in 1:length(deepLookList) ) {

	#		I-Type
xlim.i = c ( min(init.rfMap.maxResp.i,final.rfMap.maxResp.i), max(init.rfMap.maxResp.i,final.rfMap.maxResp.i) );
ylim.i = c ( min(init.rfMap.rfAreas.i,final.rfMap.rfAreas.i), max(init.rfMap.rfAreas.i,final.rfMap.rfAreas.i) );
plot ( init.rfMap.maxResp.i, init.rfMap.rfAreas.i, xlim=xlim.i, ylim=ylim.i, pch=".", col=1, cex=2,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer I-Cells\n",initText,"(.); ",finalText,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.rfMap.maxResp.i, final.rfMap.rfAreas.i, pch="+", col=2, cex=1 );
for ( i in 1:length(deepLookList) ) {
	points ( final.rfMap.maxResp.i[deepLookList[i]], final.rfMap.rfAreas.i[deepLookList[i]],
		pch=symbolChoice[i], col=3, cex=2, lwd=2 );
} # for ( i in 1:length(deepLookList) ) {

	#
	#	Compare Last Iter of Baseline ("Baseline") with Last Iter of Focal Stim Control
	#		E-Type
xlim.e = c ( min(init.rfMap.maxResp.e,final.ctl.rfMap.maxResp.e), max(init.rfMap.maxResp.e,final.ctl.rfMap.maxResp.e) );
ylim.e = c ( min(init.rfMap.rfAreas.e,final.ctl.rfMap.rfAreas.e), max(init.rfMap.rfAreas.e,final.ctl.rfMap.rfAreas.e) );
plot ( init.rfMap.maxResp.e, init.rfMap.rfAreas.e, xlim=xlim.e, ylim=ylim.e, pch=".", col=1, cex=2,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer E-Cells\n",initText,"(.); ",finalText2,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.ctl.rfMap.maxResp.e, final.ctl.rfMap.rfAreas.e, pch="+", col=2, cex=1 );
for ( i in 1:length(deepLookList) ) {
	points ( final.ctl.rfMap.maxResp.e[deepLookList[i]], final.ctl.rfMap.rfAreas.e[deepLookList[i]],
		pch=symbolChoice[i], col=3, cex=2, lwd=2 );
} # for ( i in 1:length(deepLookList) ) {

pchChar = 15; wtmp = which ( InGrid ( final.ctl.rfMap.maxResp.e, final.ctl.rfMap.rfAreas.e, 0, 22, 27, 31 ) == 1 );
xtmp = mean ( final.ctl.rfMap.maxResp.e[wtmp] ); ytmp = mean ( final.ctl.rfMap.rfAreas.e[wtmp]); points ( xtmp, ytmp, pch=pchChar, col=3, cex=1.5, lwd=2 );

pchChar = 16; wtmp = which ( InGrid ( final.ctl.rfMap.maxResp.e, final.ctl.rfMap.rfAreas.e, 20, 23, 31, 36 ) == 1 );
xtmp = mean ( final.ctl.rfMap.maxResp.e[wtmp] ); ytmp = mean ( final.ctl.rfMap.rfAreas.e[wtmp]); points ( xtmp, ytmp, pch=pchChar, col=3, cex=1.5, lwd=2 );

pchChar = 17; wtmp = which ( InGrid ( final.ctl.rfMap.maxResp.e, final.ctl.rfMap.rfAreas.e, 21, 25, 36, 40 ) == 1 );
xtmp = mean ( final.ctl.rfMap.maxResp.e[wtmp] ); ytmp = mean ( final.ctl.rfMap.rfAreas.e[wtmp]); points ( xtmp, ytmp, pch=pchChar, col=3, cex=1.5, lwd=2 );

	#		I-Type
xlim.i = c ( min(init.rfMap.maxResp.i,final.ctl.rfMap.maxResp.i), max(init.rfMap.maxResp.i,final.ctl.rfMap.maxResp.i) );
ylim.i = c ( min(init.rfMap.rfAreas.i,final.ctl.rfMap.rfAreas.i), max(init.rfMap.rfAreas.i,final.ctl.rfMap.rfAreas.i) );
plot ( init.rfMap.maxResp.i, init.rfMap.rfAreas.i, xlim=xlim.i, ylim=ylim.i, pch=".", col=1, cex=2,
	main=paste("RF Extent vs RF Peak Response\nCortical Layer I-Cells\n",initText,"(.); ",finalText2,"(*)"),
	ylab="RF Extent", xlab="RF Peak Response" );
points ( final.ctl.rfMap.maxResp.i, final.ctl.rfMap.rfAreas.i, pch="+", col=2, cex=1 );
for ( i in 1:length(deepLookList) ) {
	points ( final.ctl.rfMap.maxResp.i[deepLookList[i]], final.ctl.rfMap.rfAreas.i[deepLookList[i]],
		pch=symbolChoice[i], col=3, cex=2, lwd=2 );
} # for ( i in 1:length(deepLookList) ) {

} # if ( 0 ) {


	################
	#
	#	Draw a series of plots.  Deep dives on individual cell RF Extents.
	#	The values displayed are the R1 values from RF_Map.  This is the
	#	raw data that goes into RF size estimation routines.
	#
	################

#deepLookList = TrimEdgesFromCellList ( seq(1,N2), N, 1 );
#deepLookList = c ( 398, 399 );
#deepLookList = TrimEdgesFromCellList ( seq(9 * N +1, 24 * N, ), N, 3 );
#tmp = Dig3SectorLocs1 ( N ); deepLookList = TrimEdgesFromCellList ( tmp[,expSectorID ], N, 1 );

deepLookList = stimZone = c ( 579, 580, 612, 613 );
deepLookList = stimZone = c ( seq ( 607, 621, 1 ), seq ( 415, 811, 33 ) );

#x11(); par(mfrow=c(3,2));
for ( iCell in deepLookList ) {

	x11(); par(mfrow=c(3,2));
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " E Cell RF Extent\nBaseline Refinement Iter ", (iWhich-1), sep="" );
	ShowVecAsMapContour ( init.r1.e.rfMap[iCell,], titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " I Cell RF Extent\nBaseline Refinement Iter ", (iWhich-1), sep="" );
	ShowVecAsMapContour ( init.r1.i.rfMap[iCell,], titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " E Cell RF Extent\nFocal Stimulation Iter ", (iWhich-1), sep="" );
	ShowVecAsMapContour ( final.r1.e.rfMap[iCell,], titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " I Cell RF Extent\nFocal Stimulation Iter ", (iWhich-1), sep="" );
	ShowVecAsMapContour ( final.r1.i.rfMap[iCell,], titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " E Cell RF Extent\nFocal Stimulation Control Iter ", (iWhich-1), sep="" );
	ShowVecAsMapContour ( final.ctl.r1.e.rfMap[iCell,], titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " I Cell RF Extent\nFocal Stimulation Control Iter ", (iWhich-1), sep="" );
	ShowVecAsMapContour ( final.ctl.r1.i.rfMap[iCell,], titleText, xlab, ylab );

	Sys.sleep ( 1.0 );

} # for ( iCell in deepLookList ) {

	################
	#
	#	Draw a series of plots.  Deep dives on individual skin node CMs.
	#
	################

if ( 0 ) {

for ( iCell in deepLookList ) {

	x11(); par(mfrow=c(2,2)); hlocs=c(11.5,22.5);

	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

	base.CM.e = TrimEdgesFromCellList ((MapMagToCellList( init.cortical.Amp.e, iCell )), N, iTrimRings );
	base.CM.i = TrimEdgesFromCellList ((MapMagToCellList( init.cortical.Amp.i, iCell )), N, iTrimRings );

	final.CM.e = TrimEdgesFromCellList ((MapMagToCellList( final.cortical.Amp.e, iCell )), N, iTrimRings );
	final.CM.i = TrimEdgesFromCellList ((MapMagToCellList( final.cortical.Amp.i, iCell )), N, iTrimRings );

	final.ctl.CM.e = TrimEdgesFromCellList ((MapMagToCellList( final.ctl.cortical.Amp.e, iCell )), N, iTrimRings );
	final.ctl.CM.i = TrimEdgesFromCellList ((MapMagToCellList( final.ctl.cortical.Amp.i, iCell )), N, iTrimRings );

	#iWhich = 26;
	#titleText = paste ( "Column # ", iCell, " E Cell CM Extent\nBaseline Refinement Iter ", (iWhich-1), sep="" );
	#tmp = rep(0,N2); tmp[base.CM.e] = 1.0; ShowVecAsMap1 ( tmp, titleText, xlab, ylab ); abline(h=hlocs,lty=3,col=605);

	#iWhich = 26;
	#titleText = paste ( "Column # ", iCell, " I Cell CM Extent\nBaseline Refinement Iter ", (iWhich-1), sep="" );
	#tmp = rep(0,N2); tmp[base.CM.i] = 1.0; ShowVecAsMap1 ( tmp, titleText, xlab, ylab ); abline(h=hlocs,lty=3,col=605);

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " E Cell CM Extent\nFocal Stimulation Iter ", (iWhich-1), sep="" );
	tmp = rep(0,N2); tmp[final.CM.e] = 1.0; ShowVecAsMap1 ( tmp, titleText, xlab, ylab ); abline(h=hlocs,lty=3,col=605);

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " I Cell CM Extent\nFocal Stimulation Iter ", (iWhich-1), sep="" );
	tmp = rep(0,N2); tmp[final.CM.i] = 1.0; ShowVecAsMap1 ( tmp, titleText, xlab, ylab ); abline(h=hlocs,lty=3,col=605);

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " E Cell CM Extent\nFocal Stimulation Control Iter ", (iWhich-1), sep="" );
	tmp = rep(0,N2); tmp[final.ctl.CM.e] = 1.0; ShowVecAsMap1 ( tmp, titleText, xlab, ylab ); abline(h=hlocs,lty=3,col=605);

	iWhich = 26;
	titleText = paste ( "Column # ", iCell, " I Cell CM Extent\nFocal Stimulation Control Iter ", (iWhich-1), sep="" );
	tmp = rep(0,N2); tmp[final.ctl.CM.i] = 1.0; ShowVecAsMap1 ( tmp, titleText, xlab, ylab ); abline(h=hlocs,lty=3,col=605);

	Sys.sleep ( 1.0 );

} # for ( iCell in deepLookList ) {

} # if ( 0 )


	################
	#
	#	Draw a series of plots.  Deep dive on input weight distributions.
	#
	################

fRoot = "Base.25.0.w1E0"; init.w1.e0 = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = "Base.25.0.w1EE"; init.w1.ee = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = "Base.25.0.w1EI"; init.w1.ei = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = "Base.25.0.w1IE"; init.w1.ie = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );

fRoot = "Base.25.25.w1E0"; refine.w1.e0 = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = "Base.25.25.w1EE"; refine.w1.ee = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = "Base.25.25.w1EI"; refine.w1.ei = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = "Base.25.25.w1IE"; refine.w1.ie = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );

fRoot = paste(wName.exp, "w1E0", sep="."); selStim.exp.w1.e0 = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = paste(wName.exp, "w1EE", sep="."); selStim.exp.w1.ee = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = paste(wName.exp, "w1EI", sep="."); selStim.exp.w1.ei = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = paste(wName.exp, "w1IE", sep="."); selStim.exp.w1.ie = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );

fRoot = paste(wName.ctl, "w1E0", sep="."); selStim.ctl.w1.e0 = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = paste(wName.ctl, "w1EE", sep="."); selStim.ctl.w1.ee = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = paste(wName.ctl, "w1EI", sep="."); selStim.ctl.w1.ei = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );
fRoot = paste(wName.ctl, "w1IE", sep="."); selStim.ctl.w1.ie = GetSparseWeightMatrix ( paste(fDir, fRoot, sep="\\") );


	##########################
	#
	#	Zoom in.  Each adjusts to its own weight.
	#
	##########################

digitBorderline = c ( 4.5, 3.5, 0 );
digitBorderline = c ( 0, 0, 0, 0 );
digitBorderline = rep ( 0, length(deepLookList) );
iBorder = 0;

for ( iCell in deepLookList ) {

	iBorder = iBorder + 1;	# Need to show digit boundaries correctly.

	x11();
	par(mfrow=c(3,4));
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

		#	Refined Map

	iWhich = 26;
	titleText = paste ( "# ", iCell, " E<-S Wghts\nBase Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( refine.w1.e0, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ylab=""; ShowVecAsMapContour ( vTmp, titleText, xlab, ylab );
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 1, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 26;
	titleText = paste ( "# ", iCell, " E<-E Wghts\nBase Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( refine.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ylab=""; ShowVecAsMapContour ( vTmp, titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "# ", iCell, " E<-I Wghts\nBase Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( refine.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ylab=""; ShowVecAsMapContour ( vTmp, titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "# ", iCell, " I<-E Wghts\nBase Iter ", (iWhich-1), sep="" );
	vTmp = GetInputWeights ( refine.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ylab=""; ShowVecAsMapContour ( vTmp, titleText, xlab, ylab );

		#	SelStimExp

	iWhich = 25;
	titleText = paste ( "# ", iCell, " E<-S Wghts\nFocal EXP Iter ", (iWhich), sep="" );
	vTmp = GetInputWeights ( selStim.exp.w1.e0, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ylab=""; ShowVecAsMapContour ( vTmp, titleText, xlab, ylab );
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 1, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 25;
	titleText = paste ( "# ", iCell, " E<-E Wghts\nFocal EXP Iter ", (iWhich), sep="" );
	vTmp = GetInputWeights ( selStim.exp.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ylab=""; ShowVecAsMapContour ( vTmp, titleText, xlab, ylab );

	iWhich = 25;
	titleText = paste ( "# ", iCell, " E<-I Wghts\nFocal EXP Iter ", (iWhich), sep="" );
	vTmp = GetInputWeights ( selStim.exp.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ylab=""; ShowVecAsMapContour ( vTmp, titleText, xlab, ylab );

	iWhich = 25;
	titleText = paste ( "# ", iCell, " I<-E Wghts\nFocal EXP Iter ", (iWhich), sep="" );
	vTmp = GetInputWeights ( selStim.exp.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ylab=""; ShowVecAsMapContour ( vTmp, titleText, xlab, ylab );

		#	SelStimCTL

	iWhich = 25;
	titleText = paste ( "# ", iCell, " E<-S Wghts\nFocal CTL Iter ", (iWhich), sep="" );
	vTmp = GetInputWeights ( selStim.ctl.w1.e0, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ylab=""; ShowVecAsMapContour ( vTmp, titleText, xlab, ylab );
	if ( digitBorderline[iBorder] != 0 ) {
		abline ( h = digitBorderline[iBorder], lty = 1, col = 3 );
	} # 	if ( digitBorderline[iBorder] != 0 ) {

	iWhich = 25;
	titleText = paste ( "# ", iCell, " E<-E Wghts\nFocal CTL Iter ", (iWhich), sep="" );
	vTmp = GetInputWeights ( selStim.ctl.w1.ee, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ylab=""; ShowVecAsMapContour ( vTmp, titleText, xlab, ylab );

	iWhich = 25;
	titleText = paste ( "# ", iCell, " E<-I Wghts\nFocal CTL Iter ", (iWhich), sep="" );
	vTmp = GetInputWeights ( selStim.ctl.w1.ei, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ylab=""; ShowVecAsMapContour ( vTmp, titleText, xlab, ylab );

	iWhich = 25;
	titleText = paste ( "# ", iCell, " I<-E Wghts\nFocal CTL Iter ", (iWhich), sep="" );
	vTmp = GetInputWeights ( selStim.ctl.w1.ie, iCell ); vTmp = vTmp [ vTmp != 0 ];
	ylab=""; ShowVecAsMapContour ( vTmp, titleText, xlab, ylab );

} # for ( iCell in deepLookList ) {


if ( 0 ) {

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

clust1 = which ( final.ctl.rfMap.maxResp.e < 25 );
clust2 = which ( final.ctl.rfMap.maxResp.e >= 25 );
final.ctl.maxResp.rfAreas.e.corr.clust1 = cor ( final.ctl.rfMap.maxResp.e[clust1], final.ctl.rfMap.rfAreas.e[clust1] );
final.ctl.maxResp.rfAreas.e.corr.clust2 = cor ( final.ctl.rfMap.maxResp.e[clust2], final.ctl.rfMap.rfAreas.e[clust2] );
final.ctl.maxResp.rfAreas.i.corr = cor ( final.ctl.rfMap.maxResp.i, final.ctl.rfMap.rfAreas.i );


	#	Do some mean, sd analysis and Wilcox analysis.
	
	#	E-cells maxResp and rfAreas

init.rfMap.maxResp.e.meanSD = MeanAndSD ( init.rfMap.maxResp.e );
final.rfMap.maxResp.e.meanSD = MeanAndSD ( final.rfMap.maxResp.e );
final.ctl.rfMap.maxResp.e.meanSD = MeanAndSD ( final.ctl.rfMap.maxResp.e );

init.rfMap.maxResp.i.meanSD = MeanAndSD ( init.rfMap.maxResp.i );
final.rfMap.maxResp.i.meanSD = MeanAndSD ( final.rfMap.maxResp.i );
final.ctl.rfMap.maxResp.i.meanSD = MeanAndSD ( final.ctl.rfMap.maxResp.i );

	#	I-cells maxResp and rfAreas

init.rfMap.rfAreas.e.meanSD = MeanAndSD ( init.rfMap.rfAreas.e );
final.rfMap.rfAreas.e.meanSD = MeanAndSD ( final.rfMap.rfAreas.e );
final.ctl.rfMap.rfAreas.e.meanSD = MeanAndSD ( final.ctl.rfMap.rfAreas.e );

init.rfMap.rfAreas.i.meanSD = MeanAndSD ( init.rfMap.rfAreas.i );
final.rfMap.rfAreas.i.meanSD = MeanAndSD ( final.rfMap.rfAreas.i );
final.ctl.rfMap.rfAreas.i.meanSD = MeanAndSD ( final.ctl.rfMap.rfAreas.i );

	#

baseline.maxResp.e.v.i.wilcox = wilcox.test ( final.rfMap.maxResp.e, final.rfMap.maxResp.i, exact = FALSE, correct=FALSE );
baseline.rfAreas.e.v.i.wilcox = wilcox.test ( final.rfMap.maxResp.e, final.rfMap.maxResp.i, exact = FALSE, correct=FALSE );

syndact.maxResp.e.v.i.wilcox = wilcox.test ( final.ctl.rfMap.maxResp.e, final.ctl.rfMap.maxResp.i, exact = FALSE, correct=FALSE );
syndact.rfAreas.e.v.i.wilcox = wilcox.test ( final.ctl.rfMap.maxResp.e, final.ctl.rfMap.maxResp.i, exact = FALSE, correct=FALSE );

init.v.baseline.rfAreas.e = wilcox.test ( init.rfMap.rfAreas.e, final.rfMap.rfAreas.e, exact = FALSE, correct=FALSE );
baseline.v.syndact.rfAreas.e = wilcox.test ( final.rfMap.rfAreas.e, final.ctl.rfMap.rfAreas.e, exact = FALSE, correct=FALSE );

init.v.baseline.rfAreas.i = wilcox.test ( init.rfMap.rfAreas.i, final.rfMap.rfAreas.i, exact = FALSE, correct=FALSE );
baseline.v.syndact.rfAreas.i = wilcox.test ( final.rfMap.rfAreas.i, final.ctl.rfMap.rfAreas.i, exact = FALSE, correct=FALSE );


	#	Analyze some of the subclusters in the scatterplots of E-cells for baseline v synactyly.
	#	What is their spatial distribution or other differentiator?

x11(); par(mfrow=c(2,2));
ShowVecAsMap ( InGrid ( final.ctl.rfMap.maxResp.e, final.ctl.rfMap.rfAreas.e,  0, 22, 27, 31 ), "B v S E Cluster 1" );
ShowVecAsMap ( InGrid ( final.ctl.rfMap.maxResp.e, final.ctl.rfMap.rfAreas.e,  20, 23, 31, 36 ), "B v S E Cluster 2" );
ShowVecAsMap ( InGrid ( final.ctl.rfMap.maxResp.e, final.ctl.rfMap.rfAreas.e,  21, 25, 36, 40 ), "B v S E Cluster 3" );
ShowVecAsMap ( InGrid ( final.ctl.rfMap.maxResp.e, final.ctl.rfMap.rfAreas.e,  25, 50, 0, 50 ), "B v S E Cluster 4" );

} # if ( 0 ) {