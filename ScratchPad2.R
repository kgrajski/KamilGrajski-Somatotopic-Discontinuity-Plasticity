######################################################################################################
######################################################################################################
#
#	ConsolidatedKnockoutTimeSeriesPlotDriver
#
######################################################################################################
######################################################################################################

#
#	ConsolidatedKnockoutTimeSeriesPlotDriver...
#

ConsolidatedKnockoutTimeSeriesPlotDriver = function ( N2, numValsPerRFTrial, iProbeCell,
								alldata.base, alldata.exp, alldata.exp2, alldata.exp3, eScale,
								titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag,
								deltaTInMsec ) {
x11();
par(mfcol=c(2,2));
iStart = 1;
iEnd = as.integer ( numValsPerRFTrial / N2 );
ConsolidatedKnockoutTimeSeriesPlotter ( N2, numValsPerRFTrial, iProbeCell, iStart, iEnd,
								alldata.base, alldata.exp, alldata.exp2, alldata.exp3, eScale,
								titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag,
								deltaTInMsec );

iStart = 98;
iEnd = 130;
ConsolidatedKnockoutTimeSeriesPlotter ( N2, numValsPerRFTrial, iProbeCell, iStart, iEnd,
								alldata.base, alldata.exp, alldata.exp2, alldata.exp3, eScale,
								titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag,
								deltaTInMsec );

iStart = 90;
iEnd = 175;
RFProbeTrialTimeSeriesPlot ( N2, numValsPerRFTrial, iProbeCell, iStart, iEnd,
								alldata.base, alldata.exp, alldata.exp2, alldata.exp3, eScale,
								titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag,
								deltaTInMsec );

} # ConsolidatedKnockoutTimeSeriesPlotDriver = function ( N2, numValsPerRFTrial, iProbeCell,


######################################################################################################
######################################################################################################
#
#	ConsolidatedKnockoutTimeSeriesPlotter 
#
######################################################################################################
######################################################################################################

ConsolidatedKnockoutTimeSeriesPlotter = function ( N2, numValsPerRFTrial, iProbeCell, iStart, iEnd,
								alldata.base, alldata.exp, alldata.exp2, alldata.exp3, eScale,
								titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag,
								deltaTInMsec ) {
	
	xlabText = xlabTextRoot;
	ylabText = ylabTextRoot;

	ts.exp = ExpRefTSResponse ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale );
	ts.exp2 = ExpRefTSResponse ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp2, eScale );
	ts.exp3 = ExpRefTSResponse ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp3, eScale );

	#x11();
	titleTextToPrint = paste ( "E", titleText, "\nStim: ", iProbeCell, "Rec: ", iProbeCell );
	xVal = iStart:iEnd;
	yMin = min ( ts.exp$v1e.base[xVal], ts.exp$v1e.exp[xVal], ts.exp2$v1e.exp[xVal], ts.exp3$v1e.exp[xVal] );
	yMax = max ( ts.exp$v1e.base[xVal], ts.exp$v1e.exp[xVal], ts.exp2$v1e.exp[xVal], ts.exp3$v1e.exp[xVal] ) + max(ts.exp$v0.base[xVal]);
	ylim = c ( yMin, yMax );
	plot ( xVal*deltaTInMsec, ts.exp$v1e.base[xVal], type="l", ylim=ylim, col=1,
			main=titleTextToPrint, xlab=xlabText, ylab="" );
	lines ( xVal*deltaTInMsec, ts.exp$v1e.exp[xVal], col=2 );
	lines ( xVal*deltaTInMsec, ts.exp2$v1e.exp[xVal], col=3 );
	lines ( xVal*deltaTInMsec, ts.exp3$v1e.exp[xVal], col=4 );
	axis(1); mtext(side=2, "V.E", cex=0.75, line=2);
	abline ( h=0, lty=4, col=5 );
	par(new=TRUE);
	plot ( xVal*deltaTInMsec, ts.exp$v0.base[xVal], type="l", col=1, lty=2, xaxt="n", yaxt="n", xlab="", ylab="" );
	axis(4); mtext(side=4, "V.S", cex=0.75, line=0);

	#x11();
	titleTextToPrint = paste ( "I", titleText, "\nStim: ", iProbeCell, "Rec: ", iProbeCell );
	yMin = min ( ts.exp$v1i.base[xVal], ts.exp$v1i.exp[xVal], ts.exp2$v1i.exp[xVal], ts.exp3$v1i.exp[xVal] );
	yMax = max ( ts.exp$v1i.base[xVal], ts.exp$v1i.expv, ts.exp2$v1i.exp[xVal], ts.exp3$v1i.exp[xVal] ) + max(ts.exp$v0.base[xVal]);
	ylim = c ( yMin, yMax );
	plot ( xVal*deltaTInMsec, ts.exp$v1i.base[xVal], type="l", ylim=ylim, col=1,
		main=titleTextToPrint, xlab=xlabText, ylab="" );
	lines ( xVal*deltaTInMsec, ts.exp$v1i.exp[xVal], col=2 );
	lines ( xVal*deltaTInMsec, ts.exp2$v1i.exp[xVal], col=3 );
	lines ( xVal*deltaTInMsec, ts.exp3$v1i.exp[xVal], col=4 );
	axis(1); mtext(side=2, "V.I", cex=0.75, line=2);
	abline ( h=0, lty=4, col=5 );
	par(new=TRUE);
	plot ( xVal*deltaTInMsec, ts.exp$v0.base[xVal], type="l", col=1, lty=2, xaxt="n", yaxt="n", xlab="", ylab="" );
	axis(4); mtext(side=4, "V.S", cex=0.75, line=0);	

} # ConsolidatedKnockoutTimeSeriesPlotDriver = function ( N2, numValsPerRFTrial, iProbeCell, ...


######################################################################################################
######################################################################################################
#
#	Publication Quality Figure Showing Membrane Potential and Firing Rate 
#
######################################################################################################
######################################################################################################

RFProbeTrialTimeSeriesPlot = function ( N2, numValsPerRFTrial, iProbeCell, iStart, iEnd,
								alldata.base, alldata.exp, alldata.exp2, alldata.exp3, eScale,
								titleText, subTitleText, xlabTextRoot, ylabTextRoot, tiffFlag,
								deltaTInMsec ) {
	
	titleText = "RF Probe Trial\nAvg. Membrane Potential";
	xlabText = xlabTextRoot;
	ylabText = ylabTextRoot;

	ts.exp = ExpRefTSResponse ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale );

	x11(); par(mfrow=c(2,2));
	titleTextToPrint = paste ( titleText );
	xVal = iStart:iEnd;
	yMin = min ( ts.exp$v1e.base[xVal], ts.exp$v1i.base, ts.exp$v0.base );
	yMax = max ( ts.exp$v1e.base[xVal], ts.exp$v1i.base, ts.exp$v0.base ) + max(ts.exp$v0.base[xVal]);
	ylim = c ( yMin, yMax );
	plot ( xVal*deltaTInMsec, ts.exp$v1e.base[xVal], type="l", cex=0.25, ylim=ylim, col=1,
		main=titleTextToPrint, xlab=xlabText, ylab="" );
	lines ( xVal*deltaTInMsec, ts.exp$v1i.base[xVal], col=2, cex=0.25 );
	axis(1); mtext(side=2, "Ve, Vi", cex=0.75, line=2);
	abline ( h=0, lty=4, col=5 );
	par(new=TRUE);
	plot ( xVal*deltaTInMsec, ts.exp$v0.base[xVal], type="l", cex=0.1, col=3, lty=1, xaxt="n", yaxt="n", xlab="", ylab="" );
	#abline ( h=0, lty=4, col=5 );
	axis(4); mtext(side=4, "Vs", cex=0.75, line=0);
	legend ( 102, 0.45, c("Ve", "Vi", "Vs"), lty=c(1,1,1), lwd=c(1,1,1), col=c(1,2,3), cex=0.60);

	titleText = "RF Probe Trial\nAvg. Spiking Rate";
	xlabText = xlabTextRoot;
	ylabText = ylabTextRoot;
	titleTextToPrint = paste ( titleText );
	xVal = iStart:iEnd;
	yMin = min ( ts.exp$v1e.base[xVal], ts.exp$v1i.base, ts.exp$v0.base );
	yMax = max ( ts.exp$v1e.base[xVal], ts.exp$v1i.base, ts.exp$v0.base ) + max(ts.exp$v0.base[xVal]);
	ylim = c ( yMin, yMax );
	ylim = sigmoid ( ylim, 4.0);
	plot ( xVal*deltaTInMsec, sigmoid(ts.exp$v1e.base[xVal],4), type="l", cex=0.25, ylim=ylim, col=1,
		main=titleTextToPrint, xlab=xlabText, ylab="" );
	lines ( xVal*deltaTInMsec, sigmoid(ts.exp$v1i.base[xVal],4), col=2, cex=0.25 );
	lines ( xVal*deltaTInMsec, sigmoid(ts.exp$v0.base[xVal],4), col=3, cex=0.25 );
	axis(1); mtext(side=2, "Re, Ri, Rs", cex=0.75, line=2);

	legend ( 125, 0.45, c("Re", "Ri", "Rs"), lty=c(1,1,1), lwd=c(1,1,1), col=c(1,2,3), cex=0.60);

} # RFProbeTrialTimeSeriesPlot = function ( N2, numValsPerRFTrial, iProbeCell, ...



