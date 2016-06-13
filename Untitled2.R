
tx.base = base$rfTrackData.e[[1]];
tcent.base  = tx.base$rfCentroid;
tcov.base  = tx.base$rfCovar;

tx.E = test.ctl.E$rfTrackData.e[[1]];
tcent.E = tx$rfCentroid;
tcov.E = tx$rfCovar;

row1 = GetLin ( c(33:40), 11, 75 );
cbind ( tcent.base[row1,2], tcent.E[row1,2] );

row2 = GetLin ( c(33:40), 12, 75 );
cbind ( tcent.base[row2,2], tcent.E[row2,2] );

row3 = GetLin ( c(33:40), 13, 75 );
cbind ( tcent.base[row3,2], tcent.E[row3,2] );

row4 = GetLin ( c(33:40), 14, 75 );
cbind ( tcent.base[row4,2], tcent.E[row4,2] );

row5 = GetLin ( c(33:40), 15, 75 );
cbind ( tcent.base[row5,2], tcent.E[row5,2] );

row6 = GetLin ( c(33:40), 16, 75 );
cbind ( tcent.base[row5,2], tcent.E[row5,2] );

