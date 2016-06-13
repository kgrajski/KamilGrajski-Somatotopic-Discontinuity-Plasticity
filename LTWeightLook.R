#
#	Long-term state variable evolution display.
#


iBase = 5;
iStart = 0;
iEnd = 5;

fromCell.1 = 113; toCell.1 = 113;

fromCell.2 = 1; toCell.2 = 1;

lt1.w1.e.0 = 0; lt1.w1.e.e = 0; lt1.w1.i.e = 0; lt1.w1.e.i = 0;
lt2.w1.e.0 = 0; lt2.w1.e.e = 0; lt2.w1.i.e = 0; lt2.w1.e.i = 0;

for ( iRefinement in iStart:iEnd ) {

	fileName = paste ( "Run.", iBase, ".", iRefinement, sep="" );
	load ( file = fileName );

	lt1.w1.e.0 = c ( lt1.w1.e.0, w1.e.0[fromCell.1, toCell.1, ] );
	lt1.w1.e.e = c ( lt1.w1.e.e, w1.e.e[fromCell.1, toCell.1, ] );
	lt1.w1.i.e = c ( lt1.w1.i.e, w1.i.e[fromCell.1, toCell.1, ] );
	lt1.w1.e.i = c ( lt1.w1.e.i, w1.e.i[fromCell.1, toCell.1, ] );
	
	lt2.w1.e.0 = c ( lt2.w1.e.0, w1.e.0[fromCell.2, toCell.2, ] );
	lt2.w1.e.e = c ( lt2.w1.e.e, w1.e.e[fromCell.2, toCell.2, ] );
	lt2.w1.i.e = c ( lt2.w1.i.e, w1.i.e[fromCell.2, toCell.2, ] );
	lt2.w1.e.i = c ( lt2.w1.e.i, w1.e.i[fromCell.2, toCell.2, ] );	

} # for ( iRefinement in iStart:iEnd ) {

lt1.w1.e.0 = lt1.w1.e.0[-1]; lt1.w1.e.e = lt1.w1.e.e[-1];
lt1.w1.i.e = lt1.w1.i.e[-1]; lt1.w1.e.i = lt1.w1.e.i[-1];

lt2.w1.e.0 = lt2.w1.e.0[-1]; lt2.w1.e.e = lt2.w1.e.e[-1];
lt2.w1.i.e = lt2.w1.i.e[-1]; lt2.w1.e.i = lt2.w1.e.i[-1];

ylim = c ( min ( lt1.w1.e.0, lt1.w1.e.e, lt1.w1.i.e, lt1.w1.e.i ),
		max ( lt1.w1.e.0, lt1.w1.e.e, lt1.w1.i.e, lt1.w1.e.i ) );

x11(); par(mfrow=c(2,1));
plot ( lt1.w1.e.0, ylim=ylim, type="l", col=1,
	main=paste("Weight Evolution From Cell", fromCell.1, "To Cell", toCell.1),
		xlab="Iteration", ylab="Weight",
		sub="W1.E.0 (Black); W1.E.E (Green); W1.E.I (Red); W1.I.E (Blue)" );
lines ( lt1.w1.e.e, col=3 );
lines ( lt1.w1.e.i, col=2 );
lines ( lt1.w1.i.e, col=4 );

ylim = c ( min ( lt2.w1.e.0, lt2.w1.e.e, lt2.w1.i.e, lt2.w1.e.i ),
		max ( lt2.w1.e.0, lt2.w1.e.e, lt2.w1.i.e, lt2.w1.e.i ) );
plot ( lt2.w1.e.0, ylim=ylim, type="l", col=1,
	main=paste("Weight Evolution From Cell", fromCell.2, "To Cell", toCell.2),
		xlab="Iteration", ylab="Weight",
		sub="W1.E.0 (Black); W1.E.E (Green); W1.E.I (Red); W1.I.E (Blue)" );
lines ( lt2.w1.e.e, col=3 );
lines ( lt2.w1.e.i, col=2 );
lines ( lt2.w1.i.e, col=4 );






