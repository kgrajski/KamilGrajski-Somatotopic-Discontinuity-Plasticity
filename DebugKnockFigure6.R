#
#	DebugKnockFigures6.R
#

ref = base;
irefIter = 15;
exp = test.ctl.E;
iexpIter = 15;
exp2 = test.ctl.I;
iexpIter2 = 15;
exp3 = test.ctl.EI;
iexpIter3 = 15;
refTitleText = "\nBaseline Network";
expTitleText = "\nBoundary E Only";
expTitleText2 = "\nBoundary I Only";
expTitleText3 = "\nBoundary E and I ";
iFilter.E = test.ctl.E.iFilter.E;
iFilter.I  = test.ctl.E.iFilter.I;
iFilter.E.2 = test.ctl.I.iFilter.E;
iFilter.I.2 = test.ctl.I.iFilter.I;
iFilter.E.3 = test.ctl.EI.iFilter.E;
iFilter.I.3 = test.ctl.EI.iFilter.I;
x = rfExpZone.x;
y = rfExpZone.y;
tiffFlag = tiffFlag;
iKnockOutLength = iKnockLength;

ref = base;
irefIter = 15;
exp = test.oneside.E;
iexpIter = 15;
exp2 = test.oneside.I;
iexpIter2 = 15;
exp3 = test.oneside.EI;
iexpIter3 = 15;
refTitleText = "\nBaseline Network";
expTitleText = "\nBoundary E Only";
expTitleText2 = "\nBoundary I Only";
expTitleText3 = "\nBoundary E and I ";
iFilter.E = test.oneside.E.iFilter.E;
iFilter.I  = test.oneside.E.iFilter.I;
iFilter.E.2 = test.oneside.I.iFilter.E;
iFilter.I.2 = test.oneside.I.iFilter.I;
iFilter.E.3 = test.oneside.EI.iFilter.E;
iFilter.I.3 = test.oneside.EI.iFilter.I;
x = rfExpZone.x;
y = rfExpZone.y;
tiffFlag = tiffFlag;
iKnockOutLength = iKnockLength;

texp = QuantRFSize ( exp$r1.e.rfMap, kRFPeakToEdgeDetect );
tref = QuantRFSize ( ref$r1.e.rfMap, kRFPeakToEdgeDetect );

which(sizeDelta.e > 1.0)

x11(); par(mfrow=c(2,1));
ShowVecAsMap(exp$r1.e.rfMap[1638,],"exp")
ShowVecAsMap(ref$r1.e.rfMap[1638,],"ref")

exp1638 = exp$r1.e.rfMap[1638,];
ref1638 = ref$r1.e.rfMap[1638,];

sum ( exp1638 > (max(exp1638)*0.5) )
sum ( exp1638 > (max(exp1638)*0.1) )

sum ( ref1638 > (max(ref1638)*0.5) )
sum ( ref1638 > (max(ref1638)*0.1) )

x11(); par(mfrow=c(2,1));
hist(exp1638, breaks=100)
hist(ref1638, breaks=100)

lookat = which ( sizeDelta.e >= 1.0 );
summary(tref[lookat])

mean.tref = mean(tref)
sd.tref = sqrt(var(tref))
thresh = mean.tref - 3.0*sd.tref

sum(tref <= thresh)


rfTrackData = ref$rfTrackData.e[[length(ref$rfTrackData.e)]]
localPlotControl = FALSE;
titleBaseText = paste(paste("E-Type","Iter",irefIter,sep=" "),refTitleText, sep="")
showEllipsesFlag = TRUE;
topoMapConfLimit = 0.5;



#	Debug movie maker.
#
ref = base
exp = test.ctl.E
iRowParm = rowParmsRFFlyOver
iColParm = colParmsRFFlyOver




#
#	Easy deep-dive on an individual RF.
#
jCell = 1686;
x11();
par(mfrow=c(2,1));
ShowVecAsMap ( test.oneside.I$r1.e.rfMap[jCell,], paste ( "Cell=", jCell ) );
ShowVecAsMap ( test.oneside.I$r1.i.rfMap[jCell,], paste ( "Cell=", jCell ) );

summary( test.oneside.I$r1.e.rfMap[jCell,] )
test.oneside.I$rfTrackData.e[[1]]$rfCentroid[jCell,]
test.oneside.I$rfTrackData.e[[1]]$rfCovar[jCell,]
tmp.rfExtent.e = QuantRFSizeB ( test.oneside.I$r1.e.rfMap, kRFPeakToEdgeDetect, 2 )
tmp.rfExtent.i = QuantRFSizeB ( test.oneside.I$r1.i.rfMap, kRFPeakToEdgeDetect, 2 )

tmp.rfExtent[iTrim] = 0;
tmp.rfExtent[test.oneside.I.iFilter.E] = 0
tmp.rfExtent[jCell]
summary(tmp.rfExtent)
sum ( test.oneside.I$r1.e.rfMap[jCell,] > 0.5 * max (test.oneside.I$r1.e.rfMap[jCell,]) )


tmp = test.ctl.E$r1.e.rfMap;
tmp [ tmp < 1.0 ] = 0.0;

tmp.rfExtent.B = QuantRFSizeB ( test.ctl.E$r1.e.rfMap, kRFPeakToEdgeDetect, minMagResponse )
tmp.rfExtent.B[jCell]


#	Debug more on these cells showing large increases in size even distant
#	from knock-zone

ref.rfExtent.e = QuantRFSizeB ( ref$r1.e.rfMap, kRFPeakToEdgeDetect, minMagResponse );
ref.rfExtent.i = QuantRFSizeB ( ref$r1.i.rfMap, kRFPeakToEdgeDetect, minMagResponse );
tooSmallList = which ( ref.rfExtent.e < ( mean(ref.rfExtent.e) - 3.0 * sqrt ( var ( ref.rfExtent.e ) )  ) );
sizeDelta.e = ( QuantRFSizeB ( exp$r1.e.rfMap, kRFPeakToEdgeDetect, minMagResponse ) / ref.rfExtent.e ) - 1.0;
if ( length ( tooSmallList ) ) { sizeDelta.e[tooSmallList] = 0; }
summary ( sizeDelta.e )
which ( sizeDelta.e == max(sizeDelta.e) )
jCell = 261
x11(); par(mfrow=c(2,1));
ShowVecAsMap ( ref$r1.e.rfMap[jCell,], paste ( "Ref Cell=", jCell ) );
ShowVecAsMap ( test.ctl.E$r1.e.rfMap[jCell,], paste ( "Exp Cell=", jCell ) );
summary( test.ctl.E$r1.e.rfMap[jCell,] )
test.ctl.E$rfTrackData.e[[1]]$rfCentroid[jCell,]
test.ctl.E$rfTrackData.e[[1]]$rfCovar[jCell,]
tmp.rfExtent = QuantRFSizeB ( test.ctl.E$r1.e.rfMap, kRFPeakToEdgeDetect, 2 )
tmp.rfExtent[iTrim] = 0;
tmp.rfExtent[test.ctl.E.iFilter.E] = 0
tmp.rfExtent[jCell]
summary(tmp.rfExtent)
sum ( test.ctl.E$r1.e.rfMap[jCell,] > 0.5 * max (test.ctl.E$r1.e.rfMap[jCell,]) )




#	Looking at the phat RFs @ border manipulations.

test.oneside.I$rfTrackData.e[[1]]$rfCentroid[1686,]
test.oneside.I$rfTrackData.e[[1]]$rfCovar[1686,]

test.oneside.I$rfTrackData.i[[1]]$rfCentroid[1686,]
test.oneside.I$rfTrackData.i[[1]]$rfCovar[1686,]



#	Looking at max resp
summary ( ref.maxresp.e )
summary ( exp.maxresp.e )
summary ( exp.maxresp.e.2 )
summary ( exp.maxresp.e.3 )

ref = base;
exp = placebo;

ref.maxresp.e = apply ( ref$r1.e.rfMap, 1, max );
ref.maxresp.i = apply ( ref$r1.i.rfMap, 1, max );

exp.maxresp.e = apply ( exp$r1.e.rfMap, 1, max );
exp.maxresp.i = apply ( exp$r1.i.rfMap, 1, max );

summary ( ref.maxresp.e )
summary ( exp.maxresp.e )

x11(); par(mfrow=c(3,2));
ShowVecAsMap ( ref.maxresp.e, "Base E" );
ShowVecAsMap ( ref.maxresp.i, "Base I" );

ShowVecAsMap ( exp.maxresp.e, "Placebo E");
ShowVecAsMap ( exp.maxresp.i, "Placebo I");

ShowVecAsMap ( log10(1+(exp.maxresp.e/ref.maxresp.e)-1), "Delta E");
ShowVecAsMap ( log10(1+(exp.maxresp.i/ref.maxresp.i)-1), "Delta I");


#	Trying out sub-window display.
w = ref.rfExtent.e
titleText = paste(paste("E-Type RF Extent","Iter",iexpIter,sep=" "),refTitleText, sep="")
xLabText
yLabText
minZ = min(ref.rfExtent.e)
maxZ = max(ref.rfExtent.e)
rowSeed=20
colSeed=5
subN=35



x11(); par(mfrow=c(2,1));
ShowVecAsMap(exp2.rfExtent.e,"e2.min.2.0")
ShowVecAsMap(tmp,"e2.min.0.0")



#
#	Help analyze the RF extent, magnitude, etc., for a single column.
#


kCell = 1088

ref.rfExtent.e[kCell]
exp.rfExtent.e[kCell]
exp2.rfExtent.e[kCell]
exp3.rfExtent.e[kCell]

ref.rfExtent.i[kCell]
exp.rfExtent.i[kCell]
exp2.rfExtent.i[kCell]
exp3.rfExtent.i[kCell]

ref.maxresp.e[kCell]
exp.maxresp.e[kCell]
exp.maxresp.e.2[kCell]
exp.maxresp.e.3[kCell]

ref.maxresp.i[kCell]
exp.maxresp.i[kCell]
exp.maxresp.i.2[kCell]
exp.maxresp.i.3[kCell]

best.ref.e = which(ref$r1.e.rfMap[kCell,]==ref.maxresp.e[kCell])
best.exp.e = which(exp$r1.e.rfMap[kCell,]==exp.maxresp.e[kCell])
best.exp.e.2 = which(exp2$r1.e.rfMap[kCell,]==exp.maxresp.e.2[kCell])
best.exp.e.3 = which(exp3$r1.e.rfMap[kCell,]==exp.maxresp.e.3[kCell])

best.ref.i = which(ref$r1.i.rfMap[kCell,]==ref.maxresp.i[kCell])
best.exp.i = which(exp$r1.i.rfMap[kCell,]==exp.maxresp.i[kCell])
best.exp.i.2 = which(exp2$r1.i.rfMap[kCell,]==exp.maxresp.i.2[kCell])
best.exp.i.3 = which(exp3$r1.i.rfMap[kCell,]==exp.maxresp.i.3[kCell])

ref$r1.e.rfMap[kCell,kCell]
exp$r1.e.rfMap[kCell,kCell]
exp2$r1.e.rfMap[kCell,kCell]
exp3$r1.e.rfMap[kCell,kCell]

ref$r1.i.rfMap[kCell,kCell]
exp$r1.i.rfMap[kCell,kCell]
exp2$r1.i.rfMap[kCell,kCell]
exp3$r1.i.rfMap[kCell,kCell]


