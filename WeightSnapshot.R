######################################################################################################
######################################################################################################
#
#	WeightSnapshot: Generate some plots and statistics on input and output connections.
#
#		Input:
#			w1.e.0, w1.e.e, w1.e.i, w1.i.e
#			cellID of the cell of interest
#			timeStamp (in units of iterations)
#
#		Output:
#			Plots and printouts.  Page one: E-cell inputs & outputs.  Page two: I-cell in/out.
#
######################################################################################################
######################################################################################################

		#	When doing certain histograms, how finely to chop up the range of values.	
	bwDivision = 10;

		# 	Process the E-cell Inputs and Outputs as Images
	x11(); par(mfrow=c(3,2));
	ViewMapImage ( w1.e.0[,,timeStamp], cellID ); title(main=paste("L0 Inputs to L1 E Cell",cellID));
	ViewMapImage ( w1.e.e[,,timeStamp], cellID ); title(main=paste("L1 E Inputs to L1 E Cell",cellID));
	ViewMapImage ( w1.e.i[,,timeStamp], cellID ); title(main=paste("L1 I Inputs to L1 E Cell",cellID));
	ViewSTDynamics ( w1.i.e[,,timeStamp], cellID, cellID, "Outputs to L1 I Cells for L1 E" );
	ViewSTDynamics ( w1.e.e[,,timeStamp], cellID, cellID, "Outputs to L1 E Cells for L1 E" );

		# 	Process the E-cell Inputs as Histograms
	E.In.w1e.0 = w1.e.0[w1.e.0[,cellID,timeStamp]>0,cellID,timeStamp];
	E.In.w1e.e = w1.e.e[w1.e.e[,cellID,timeStamp]>0,cellID,timeStamp];
	E.In.w1e.i = w1.e.i[w1.e.i[,cellID,timeStamp]>0,cellID,timeStamp];
	bw = RangeValue ( c ( E.In.w1e.0, E.In.w1e.e, E.In.w1e.i ) ) / bwDivision;
	df = data.frame ( E.In.w1e.0, E.In.w1e.e, E.In.w1e.i );
	x11(); ggplot ( melt(df), aes(value, fill=variable)) + geom_histogram(position="dodge", binwidth=bw,
		main=paste("Inputs to L1 E Cell",cellID));

		# 	Process the E-cell Outputs as Histograms
	E.Out.w1i.e = w1.i.e[cellID,w1.i.e[cellID,,timeStamp]>0,timeStamp];
	E.Out.w1e.e = w1.e.e[cellID,w1.e.e[cellID,,timeStamp]>0,timeStamp];
	bw = RangeValue ( c ( E.Out.w1i.e, E.Out.w1e.e ) ) / bwDivision;
	df = data.frame ( E.Out.w1i.e, E.Out.w1e.e );
	x11(); ggplot ( melt(df), aes(value, fill=variable)) + geom_histogram(position="dodge", binwidth=bw);

		#	Process the I-cell Inputs and Outputs as Images
	x11(); par(mfrow=c(3,2));
	ViewMapImage ( w1.i.e[,,timeStamp], cellID ); title(main=paste("L1 E Inputs to L1 I Cell",cellID));
	ViewSTDynamics ( w1.e.i[,,timeStamp], cellID, cellID, "Outputs to L1 E Cells for L1 I" );

		#	Process the I-cell Inputs and Outputs as Histograms
	tmpIn = I.In.w1i.e = w1.i.e[ w1.i.e[,cellID,timeStamp]>0, cellID, timeStamp];
	tmpOut = I.Out.w1e.i = w1.e.i[ cellID, w1.e.i[cellID,,timeStamp]>0, timeStamp];
	bw = RangeValue ( c ( I.In.w1i.e, I.Out.w1e.i ) ) / bwDivision;
	df = data.frame ( I.In.w1i.e, I.Out.w1e.i );
	x11(); ggplot ( melt(df), aes(value, fill=variable)) + geom_histogram(position="dodge", binwidth=bw);

		#	For a given set of weights, compare at two time points using histogram.
	E.In.w1e.0.ctl = w1.e.0[w1.e.0[,cellID,timeStampCtl]>0,cellID,timeStampCtl];
	bw = RangeValue ( c ( E.In.w1e.0.ctl, E.In.w1e.0 ) ) / bwDivision;
	df = data.frame ( E.In.w1e.0.ctl, E.In.w1e.0 );
	x11(); ggplot ( melt(df), aes(value, fill=variable)) + geom_histogram(position="dodge", binwidth=bw);

	E.In.w1e.e.ctl = w1.e.e[w1.e.e[,cellID,timeStampCtl]>0,cellID,timeStampCtl];
	bw = RangeValue ( c ( E.In.w1e.e.ctl, E.In.w1e.e ) ) / bwDivision;
	df = data.frame ( E.In.w1e.e.ctl, E.In.w1e.e );
	x11(); ggplot ( melt(df), aes(value, fill=variable)) + geom_histogram(position="dodge", binwidth=bw);

	E.In.w1e.i.ctl = w1.e.i[w1.e.i[,cellID,timeStampCtl]>0,cellID,timeStampCtl];
	bw = RangeValue ( c ( E.In.w1e.i.ctl, E.In.w1e.i ) ) / bwDivision;
	df = data.frame ( E.In.w1e.i.ctl, E.In.w1e.i );
	x11(); ggplot ( melt(df), aes(value, fill=variable)) + geom_histogram(position="dodge", binwidth=bw);

	E.Out.w1i.e.ctl = w1.i.e[cellID,w1.i.e[cellID,,timeStampCtl]>0,timeStampCtl];
	bw = RangeValue ( c ( E.Out.w1i.e.ctl, E.Out.w1i.e ) ) / bwDivision;
	df = data.frame ( E.Out.w1i.e.ctl, E.Out.w1i.e );
	x11(); ggplot ( melt(df), aes(value, fill=variable)) + geom_histogram(position="dodge", binwidth=bw);

	E.Out.w1e.e.ctl = w1.e.e[cellID,w1.e.e[cellID,,timeStampCtl]>0,timeStampCtl];
	bw = RangeValue ( c ( E.Out.w1e.e.ctl, E.Out.w1e.e ) ) / bwDivision;
	df = data.frame ( E.Out.w1e.e.ctl, E.Out.w1e.e );
	x11(); ggplot ( melt(df), aes(value, fill=variable)) + geom_histogram(position="dodge", binwidth=bw);

	I.In.w1i.e.ctl = w1.i.e[w1.i.e[,cellID,timeStampCtl]>0,cellID,timeStampCtl];
	bw = RangeValue ( c ( I.In.w1i.e.ctl, I.In.w1i.e ) ) / bwDivision;
	df = data.frame ( I.In.w1i.e.ctl, I.In.w1i.e );
	x11(); ggplot ( melt(df), aes(value, fill=variable)) + geom_histogram(position="dodge", binwidth=bw);

	I.Out.w1e.i.ctl = w1.e.i[cellID,w1.e.i[cellID,,timeStampCtl]>0,timeStampCtl];
	bw = RangeValue ( c ( I.Out.w1e.i.ctl, I.Out.w1e.i ) ) / bwDivision;
	df = data.frame ( I.Out.w1e.i.ctl, I.Out.w1e.i );
	x11(); ggplot ( melt(df), aes(value, fill=variable)) + geom_histogram(position="dodge", binwidth=bw);	


