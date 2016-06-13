
	#
	#	WeightTime.R
	#
	#	Neural Network Refinement Version 1
	#
	#	Tool to aid analysis of the time evolution of synaptic weights.
	#
	#	Typical use: Run plasticity on a patch (e.g., 10).  Do time plots
	#	for some weights coming into cell 10.  Compare to another cell
	#	(e.g., 113) where we wouldn't expect as big changes.
	#





x11();
par(mfrow=c(2,1));

fromCellID = 10;
toCellID = 10;
yaxis = c ( min(	w1.e.0[fromCellID, toCellID, ],
			w1.e.e[fromCellID, toCellID, ],
			w1.e.i[fromCellID, toCellID, ], 
			w1.i.e[fromCellID, toCellID, ] ),
		max(	w1.e.0[fromCellID, toCellID, ],
			w1.e.e[fromCellID, toCellID, ],
			w1.e.i[fromCellID, toCellID, ], 
			w1.i.e[fromCellID, toCellID, ] ) );
plot ( w1.e.0[fromCellID, toCellID, ], ylim=yaxis, type="l", col=1,
	main=paste("Weight Evolution From Cell", fromCellID, "To Cell", toCellID),
		xlab="Iteration", ylab="Weight",
		sub="W1.E.0 (Black); W1.E.E (Green); W1.E.I (Red); W1.I.E (Blue)" );
lines ( w1.e.e[fromCellID, toCellID, ], col=3 );
lines ( w1.e.i[fromCellID, toCellID, ], col=2 );
lines ( w1.i.e[fromCellID, toCellID, ], col=4 );


fromCellID = 225;
toCellID = 225;
yaxis = c ( min(	w1.e.0[fromCellID, toCellID, ],
			w1.e.e[fromCellID, toCellID, ],
			w1.e.i[fromCellID, toCellID, ], 
			w1.i.e[fromCellID, toCellID, ] ),
		max(	w1.e.0[fromCellID, toCellID, ],
			w1.e.e[fromCellID, toCellID, ],
			w1.e.i[fromCellID, toCellID, ], 
			w1.i.e[fromCellID, toCellID, ] ) );
plot ( w1.e.0[fromCellID, toCellID, ], ylim=yaxis, type="l", col=1,
	main=paste("Weight Evolution From Cell", fromCellID, "To Cell", toCellID),
		xlab="Iteration", ylab="Weight",
		sub="W1.E.0 (Black); W1.E.E (Green); W1.E.I (Red); W1.I.E (Blue)" );
lines ( w1.e.e[fromCellID, toCellID, ], col=3 );
lines ( w1.e.i[fromCellID, toCellID, ], col=2 );
lines ( w1.i.e[fromCellID, toCellID, ], col=4 );


