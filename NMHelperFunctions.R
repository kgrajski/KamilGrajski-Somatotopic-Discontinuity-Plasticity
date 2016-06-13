#
#	July 23, 2014
#	Kamil A. Grajski
#

#
#	This is the main file containing various helper functions.
#	Over time this file is getting unwieldy.  Switch to having smaller individual function files.
#	Collect them here.
#
source ( "Dig3SectorLocs1.R" );
source ( "D3LongAxisLocs1.R" );
source ( "D3CorticalMagnificationPlots.R" );
source ( "D3RFDataBoxPlots.R" );
source ( "InvCortMagRuleCheck.R" );
source ( "MapSectorNumberToLabel.R" );
source ( "D3MannWhitneyTestFrontEnd.R" );
source ( "LabelledCDF.R" );
source ( "MapMagToCellList.R" );
source ( "CMRule.SpatialPlots.R" );
source ( "RFTransLocPlot.R");

library(ggplot2);
library(reshape2);
library(ellipse);
library(stats4);
library(animation);
library(fields);
library(graphics);
library(pracma);

#
#library(spacetime);
#library(animation);

######################################################################################################
######################################################################################################
#
#	Neural Modeling Helper Functions
#
######################################################################################################
######################################################################################################

######################################################################################################
######################################################################################################
#
#	CartesianToPolar: a convenient way to package up the mean & sd
#
######################################################################################################
######################################################################################################

CtoP = function ( x ) {
	if ( length(x) == 2 ) {
		r = sqrt(x[1]^2 + x[2]^2);
		theta = atan(x[2]/x[1]);
	} else {
		r = sqrt(apply((x^2),1,sum));
		theta = atan(x[,2]/x[,1]);
	}
	return ( cbind(r,theta) );
} # CartesianToPolar = function( x, y ) { {

######################################################################################################
######################################################################################################
#
#	MeanAndSD: a convenient way to package up the mean & sd
#
######################################################################################################
######################################################################################################

MeanAndSD = function ( x ) {
	return ( list ( mean=mean(x), sd=sqrt(var(x)) ) );
} # MeanAndSD = RangeValue ( x ) { {

######################################################################################################
######################################################################################################
#
#	InGrid: a convenient way to package up the mean & sd
#
######################################################################################################
######################################################################################################

InGrid = function ( x, y, xMin, xMax, yMin, yMax ) {

	tmp = rep ( 0, length(x) );
	tmp[ (x>=xMin) & (x<=xMax) & (y>=yMin) & (y<=yMax) ] = 1;
	return ( tmp );

} # InGrid = RangeValue ( x ) { {

######################################################################################################
######################################################################################################
#
#	VectorNormalize: return the value of the range of the input argument.
#
######################################################################################################
######################################################################################################

VectorNormalize = function ( x ) {
	return ( x / sqrt ( x %*% x ) );
} # VectorNormalize = RangeValue ( x ) { {

######################################################################################################
######################################################################################################
#
#	RangeValue: return the value of the range of the input argument.
#
######################################################################################################
######################################################################################################

RangeValue = function ( x ) {
	tmp = range ( x );
	return ( (tmp[2] - tmp[1]) );
} # rangeValue = RangeValue ( x ) { {

######################################################################################################
######################################################################################################
#
#	Sigmoid: convert "membrane potential" to "spike rate".
#
######################################################################################################
######################################################################################################

sigmoid = function ( x, beta ) {
	y = 0.5 * (1.0 + tanh(beta*(x - 0.5)));
	return ( y );
} # sigmoid = function ( x, beta ) {

######################################################################################################
######################################################################################################
#
#	logistic: convert "membrane potential" to "spike rate".
#
######################################################################################################
######################################################################################################

logistic = function ( x, beta ) {
	y = 1.0 / ( 1.0 + exp ( -2.0 * beta * x ) );
	return ( y );
} # sigmoid = function ( x, beta ) {

######################################################################################################
######################################################################################################
#
#	GetCol: support function for mapping NxN to linear representation.
#
######################################################################################################
######################################################################################################

GetCol = function ( k, N ) {
	k = 1 + trunc ( (k - 1) / N );
	return ( k );	
} # GetCol = function ( k, N )

######################################################################################################
######################################################################################################
#
#	GetRow: support function for mapping NxN to linear representation.
#
######################################################################################################
######################################################################################################

GetRow = function ( k, N ) {
	tmp = k - N * ( GetCol ( k, N ) - 1 );
	return ( tmp );	
} # GetRow = function ( k, N )

######################################################################################################
######################################################################################################
#
#	GetLin: given Col and Row, find Linear position.
#
######################################################################################################
######################################################################################################

GetLin = function ( rowVal, colVal, N ) {
	tmp = (colVal - 1) * N + rowVal;
	return ( tmp );
} # GetLin = function ( rowVal, colVal, N )

######################################################################################################
######################################################################################################
#
#	Legal: check if a number is in an interval, inclusive of the endpoints
#
######################################################################################################
######################################################################################################

Legal = function ( x, minX, maxX ) {
	tmp = 0;
	if ( (x >= minX) & (x <= maxX) ) {
		tmp = 1;
	} # if ( (x >= minX) & (x <= maxX) ) {
	return ( tmp );
} # Legal = function ( min, max )

######################################################################################################
######################################################################################################
#
#	ViewMap: use image to display a grid of data (convert from a linear representation)
#			where the relevant data is a column of a matrix.
#
######################################################################################################
######################################################################################################

ViewMapImage = function ( w, toShow ) {
	numColors = 128;
	N = sqrt ( dim(w)[1] );
	image ( c(1:N), c(1:N), matrix ( w[,toShow], nrow=N, ncol=N, byrow=FALSE ),
		col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ) );
} # ViewMap = function ( w, toShow ) {

######################################################################################################
######################################################################################################
#
#	ViewContour: use image to display a grid of data (convert from a linear representation)
#			where the relevant data is a column of a matrix.
#
######################################################################################################
######################################################################################################

ViewMapContour = function ( w, toShow ) {
	numColors = 128;
	N = sqrt ( dim(w)[1] );
	contour ( c(1:N), c(1:N), matrix ( w[,toShow], nrow=N, ncol=N, byrow=FALSE ) );
} # ViewMap = function ( w, toShow ) {

######################################################################################################
######################################################################################################
#
#	ViewSTDynamics: use image to display a grid of data (convert from a linear representation)
#			where the relevant data is a row in a matrix
#
######################################################################################################
######################################################################################################

ViewSTDynamics = function ( w, fromShow, toShow, txtTitle ) {
	numColors = 128;
	N = sqrt ( dim(w)[2] );
	for ( i in fromShow:toShow ) {
		image ( c(1:N), c(1:N), matrix (w[i,], nrow=N, ncol=N, byrow=FALSE ),
			col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67,
			end = 1.0, alpha = 1 ) );
		tmpTxtTitle = paste ( txtTitle, "Cell", i );
		title ( main=tmpTxtTitle );
	} # for ( i in fromShow:toShow ) {
} # ViewSTDynamics = function ( w, toShow ) {

######################################################################################################
######################################################################################################
#
#	ShowVecAsMap: use image to display a vector as a grid of data
#
######################################################################################################
######################################################################################################

ExtractSubNet = function ( w, N, rowSeed, colSeed, subN ) {
	w.new = rep ( 0, subN^2 );
	iTmp = 1;
	iWindow = seq ( rowSeed, rowSeed + subN - 1, 1 );
	iColList = seq ( colSeed, colSeed + subN - 1, 1 );
	for ( iCol in iColList ) {
		iWinStart = ( iCol - 1 ) * N + rowSeed;
		iWinEnd = iWinStart + subN - 1;
		w.new[iTmp:(iTmp+subN-1)] = w[iWinStart:iWinEnd];
		iTmp = iTmp + subN;
	} # for ( iCol in colSeed:(colSeed+subN-1) ) {
	return ( w.new );
} # ExtractSubNet = function ( w, N, rowSeed, colSeed, subN );

ShowVecAsMap = function ( w, titleText ) {
	numColors = 128;
	N = sqrt ( length ( w ) );
	image.plot ( c(1:N), c(1:N), matrix ( w, nrow=N, ncol=N, byrow=FALSE ),
		col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ) );
	title ( main=titleText );
} # ShowVecAsMap = function ( w, toShow ) {

ShowVecAsMapContour = function ( w, titleText, xlab, ylab ) {
	numColors = 128;
	N = sqrt ( length ( w ) );
	image.plot ( c(1:N), c(1:N), matrix ( w, nrow=N, ncol=N, byrow=FALSE ), legend.shrink=0.5,
			main=titleText, xlab=xlab, ylab=ylab );
} # ShowVecAsMapContour = function ( w, toShow ) {

ShowVecAsMap1 = function ( w, titleText, xLabText, yLabText, iHorizontal=FALSE, xaxt="s", yaxt="s" ) {
	numColors = 128;
	N = sqrt ( length ( w ) );
	image.plot ( c(1:N), c(1:N), matrix ( log10(w), nrow=N, ncol=N, byrow=FALSE ),
		col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
		main=titleText, horizontal=iHorizontal, xlab="", ylab="", xaxt=xaxt, yaxt=yaxt,
		legend.shrink=0.8, legend.cex=0.6 );
} # ShowVecAsMap1 = function ( w, toShow ) {

ShowVecAsMap1NoLog = function ( w, titleText, xLabText, yLabText, iHorizontal=FALSE, xaxt="n", yaxt="n" ) {
	numColors = 128;
	N = sqrt ( length ( w ) );
	w.range = c(min(w), max(w));
	image.plot ( c(1:N), c(1:N), matrix ( w, nrow=N, ncol=N, byrow=FALSE ),
		legend.shrink=0.80, legend.cex=0.5, legend.mar=0, legend.line=0,
		col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ),
		horizontal=iHorizontal, main=titleText, xaxt=xaxt, yaxt=yaxt, xlab="", ylab="",
		axis.args=list(at=round(seq(w.range[1], w.range[2], ((w.range[2]-w.range[1])/2)),4),
					labels=round(seq(w.range[1], w.range[2], ((w.range[2]-w.range[1])/2)),3),
					cex.axis=0.6, line=0)
		);
} # ShowVecAsMap1NoLog = function ( w, toShow ) {

ShowVecAsMap2 = function ( w, titleText, xLabText, yLabText, minZ, maxZ ) {
	numColors = 128;
	N = sqrt ( length ( w ) );

		#	The following is a kludge.  Image.plot throws an error if w is uniform.  This works around that.
	if ( min(w) == max(w) ) {	
		w[1] = 1.0;
		minZ = min(w); maxZ = max(w);
	} # 	if ( min(w) == max(w) ) {
	image.plot ( c(1:N), c(1:N), matrix ( w, nrow=N, ncol=N, byrow=FALSE ),
		col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ), zlim = c ( minZ, maxZ ),
		main=titleText, xlab=xLabText, ylab=yLabText, xaxt="n", yaxt="n", legend.shrink=0.8, legend.cex=0.75 );
} # ShowVecAsMap2 = function ( w, toShow ) {

ShowVecAsMap2ImageOnly = function ( w, titleText, xLabText, yLabText, minZ, maxZ ) {
	numColors = 128;
	N = sqrt ( length ( w ) );

		#	The following is a kludge.  Image.plot throws an error if w is uniform.  This works around that.
	if ( min(w) == max(w) ) {	
		w[1] = 1.0;
		minZ = min(w); maxZ = max(w);
	} # 	if ( min(w) == max(w) ) {
	image ( c(1:N), c(1:N), matrix ( w, nrow=N, ncol=N, byrow=FALSE ), cex=0.5,
		col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ), zlim = c ( minZ, maxZ ),
		main=titleText, xlab=xLabText, ylab=yLabText, xaxt="n", yaxt="n" );
} # ShowVecAsMap2ImageOnly = function ( w, toShow ) {

	#	Give the same cal as ShowBVecAsMap2, but print only a subset of the lattice.
ShowVecAsMap2W = function ( w, titleText, xLabText, yLabText, minZ, maxZ, rowSeed=1, colSeed=1, subN=NULL ) {
	numColors = 128;
	N = sqrt ( length ( w ) );
	if ( subN ) {
		w = ExtractSubNet ( w, N, rowSeed, colSeed, subN );
		minZ = min(w); maxZ = max(w);
	} else {
		subN = N;
	} # if ( subN )
		#	The following is a kludge.  Image.plot throws an error if w is uniform.  This works around that.
	if ( min(w) == max(w) ) {	
		w[1] = 1.0;
		minZ = min(w); maxZ = max(w);
	} # 	if ( min(w) == max(w) ) {
	image.plot ( c(rowSeed:(rowSeed+subN-1)), c(colSeed:(colSeed+subN-1)), matrix ( w, nrow=subN, ncol=subN, byrow=FALSE ),
		col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ), zlim = c ( minZ, maxZ ),
		main=titleText, xlab=xLabText, ylab=yLabText );
} # ShowVecAsMap2W = function ( w, toShow ) {

ShowVecAsMap3 = function ( w, titleText, xLabText, yLabText, minZ, maxZ ) {
	numColors = 128;
	N = sqrt ( length ( w ) );
	image.plot ( c(1:N), c(1:N), matrix ( w, nrow=N, ncol=N, byrow=FALSE ),
		col=rainbow(numColors, s = 1.0, v = 1.0, start = 0.67, end = 1.0, alpha = 1 ), zlim = c ( minZ, maxZ ),
		main=titleText, xlab=xLabText, ylab=yLabText, add=TRUE );
} # ShowVecAsMap3 = function ( w, toShow ) {


######################################################################################################
######################################################################################################
#
#	WghtInit1: weight initialization Version 1.
#	Input: 	x: N x N matrix.  Col i corresponds to ith (post-synaptic) cell.
#			gVal: length of one side of grid centered at ith cell, e.g.,
#				jth cell wants to receive inputs from a size 2*gVal+1 centered on itself.
#
#	Note:		ith Column in x is an NxN layer stored in linear form. 
#
######################################################################################################
######################################################################################################

WghtInit1 = function ( w, gVal, wghtMinValue ) {
	
	Nsq = dim(w)[1];
	N = sqrt(Nsq);
	gridLength = 2 * gVal + 1;
	for ( cellID in 1:Nsq ) {	# For each cell in the receiving layer.

		i.cellID = GetRow ( cellID, N ) - gVal;	# Find the Cartesian Coordinate
		j.cellID = GetCol ( cellID, N ) - gVal;
		
		for ( i in 1:gridLength ) {	# Run through the patch.

			for ( j in 0:(gridLength-1) ) {

				if ( Legal ( ( j + i.cellID ), 1, N ) & Legal ( j.cellID, 1, N ) ) {
					w[ GetLin ( (j + i.cellID), j.cellID, N ), cellID] = runif(1, min=wghtMinValue);
				} # if ( Legal ( k-j, 1, Nsq ) )				

			} # for ( j in 0:(gridLength-1) )
	
			j.cellID = j.cellID + 1;
			
		} # for ( k in 1:gridLength )
		
	} # for ( cellID in 1:Nsq )

	return ( w );

} # WghtInit1 ( function ( w, gVal, wghtMinValue ) ) {

######################################################################################################
#
#	WghtInit2: weight initialization Version 2.  Does the work to properly normalize weights.
#	Input: 	w: N2 x N2 matrix.  Col i corresponds to ith (post-synaptic) cell.
#			gVal: length of one side of grid centered at ith cell.
#	Output: 	return an additional vector that counts the actual number of output connections per cell.
#
#	Note:		ith Column in x is an NxN layer stored in linear form. 
#
######################################################################################################
######################################################################################################

WghtInit2 = function ( w, gVal, wghtMinValue ) {
	
	Nsq = dim(w)[1];
	N = sqrt(Nsq);
	n.Connect = rep ( 0, Nsq );

	gridLength = 2 * gVal + 1;
	for ( cellID in 1:Nsq ) {	# For each cell in the receiving layer.

		i.cellID = GetRow ( cellID, N ) - gVal;	# Find the Cartesian Coordinate
		j.cellID = GetCol ( cellID, N ) - gVal;
		
		for ( i in 1:gridLength ) {	# Run through the patch.

			for ( j in 0:(gridLength-1) ) {

				if ( Legal ( ( j + i.cellID ), 1, N ) & Legal ( j.cellID, 1, N ) ) {
					tmp = GetLin ( (j + i.cellID), j.cellID, N );
					w[ tmp, cellID] = runif(1, min=wghtMinValue);
					n.Connect[tmp] = n.Connect[tmp] + 1;
				} # if ( Legal ( k-j, 1, Nsq ) )				

			} # for ( j in 0:(gridLength-1) )
	
			j.cellID = j.cellID + 1;
			
		} # for ( k in 1:gridLength )
		
	} # for ( cellID in 1:Nsq )

	return ( list ( w=w, n.Connect=n.Connect ) );

} # WghtInit2 ( function ( w, gVal, wghtMinValue ) ) {


######################################################################################################
#
#	WghtInit3: weight initialization Version 2.  Does the work to properly normalize weights.
#	Difference from WghtInit2: thinking more carefully about the initial distribution.
#	Input: 	w: N2 x N2 matrix.  Col i corresponds to ith (post-synaptic) cell.
#			gVal: length of one side of grid centered at ith cell.
#	Output: 	return an additional vector that counts the actual number of output connections per cell.
#
#	Note:		ith Column in x is an NxN layer stored in linear form. 
#
######################################################################################################
######################################################################################################

WghtInit3 = function ( w, gVal, wghtMinValue, wghtMaxValue ) {
	
	Nsq = dim(w)[1];
	N = sqrt(Nsq);
	n.Connect = rep ( 0, Nsq );

	gridLength = 2 * gVal + 1;
	for ( cellID in 1:Nsq ) {	# For each cell in the receiving layer.

		i.cellID = GetRow ( cellID, N ) - gVal;	# Find the Cartesian Coordinate
		j.cellID = GetCol ( cellID, N ) - gVal;
		
		for ( i in 1:gridLength ) {	# Run through the patch.

			for ( j in 0:(gridLength-1) ) {

				if ( Legal ( ( j + i.cellID ), 1, N ) & Legal ( j.cellID, 1, N ) ) {
					tmp = GetLin ( (j + i.cellID), j.cellID, N );
					w[ tmp, cellID ] = runif(1, min=wghtMinValue, max=wghtMaxValue );
					n.Connect[tmp] = n.Connect[tmp] + 1;
				} # if ( Legal ( k-j, 1, Nsq ) )				

			} # for ( j in 0:(gridLength-1) )
	
			j.cellID = j.cellID + 1;
			
		} # for ( k in 1:gridLength )
		
	} # for ( cellID in 1:Nsq )

	return ( list ( w=w, n.Connect=n.Connect ) );

} # WghtInit3 ( function ( w, gVal, wghtMinValue ) ) {


######################################################################################################
######################################################################################################
#
#	NormalizeWeights1: weight normalization version 1.
#	Input: 	x: N x N matrix.
#				Row i corresponds to ith (pre-synaptic) cell.
#				Col j corresponds to the jth (post-synaptic) cell.
#			When normalizing based on inputs: normalize by column.
#			When normalizing based on outputs: normalize by row.
#
#			The input layer should be normalized on the basis of output (by row).
#			The cortical layer(s) should be based on inputs (by col).
#				
#
#	Note:		ith Column in x is an NxN layer stored in linear form. 
#
######################################################################################################
######################################################################################################

NormalizeWeights1 = function ( w, normByInput ) {

	N2 = dim(w)[1];

	if ( normByInput ) {	# Normalize by col, e.g., post-synaptic cell view to cortical layer.
		
		w = w / matrix ( apply ( w, 2, sum ), nrow=N2, ncol=N2, byrow = TRUE );

	} else {			# Normalize by row, e.g., pre-synaptic cell view from input layer.

		w = w / matrix ( apply ( w, 1, sum ), nrow=N2, ncol=N2, byrow = FALSE );

	} # if ( normByInput)
	return ( w );

} # NormalizeWeights1 = function ( w, normByInput )



######################################################################################################
######################################################################################################
#
#	NormalizeOutputWeights: output weight normalization.
#	Input: 	w: N x N matrix.
#				Row i corresponds to ith (pre-synaptic) cell.
#				Col j corresponds to the jth (post-synaptic) cell.
#			numPreSyn: jth element is number of output connections made by pre-synaptic cell j.
#			
#			When normalizing based on outputs: normalize w by row.
#				
#
######################################################################################################
######################################################################################################

NormalizeOutputWeights = function ( w, numPreSyn ) {

	N2 = dim(w)[1];

	w = w / matrix ( ( apply ( w, 1, sum ) / ( numPreSyn / max ( numPreSyn ) ) ), nrow=N2, ncol=N2, byrow = FALSE );

	return ( w );

} # NormalizeOutputWeights = function ( w, normByInput )

######################################################################################################
######################################################################################################
#
#	NormalizeOutputWeights1: output weight normalization.
#	Input: 	w: N x N matrix.
#				Row i corresponds to ith (pre-synaptic) cell.
#				Col j corresponds to the jth (post-synaptic) cell.
#			numPreSyn: jth element is number of output connections made by pre-synaptic cell j.
#			
#			When normalizing based on outputs: normalize w by row.
#
#			Difference from NormalizeOutputWeights is that here we allow the cell to
#			have a greater than unity value to which to normalize.
#				
#
######################################################################################################
######################################################################################################

NormalizeOutputWeights1 = function ( w, numPreSyn, wResource ) {

	N2 = dim(w)[1];

	w = ( w / matrix ( ( apply ( w, 1, sum ) / ( (wResource * numPreSyn) / max ( numPreSyn ) ) ), nrow=N2, ncol=N2, byrow = FALSE ) );

	return ( w );

} # NormalizeOutputWeights1 = function ( w, normByInput )


######################################################################################################
######################################################################################################
#
#	NormalizeOutputWeights2: output weight normalization.
#	Input: 	w: N x N matrix.
#				Row i corresponds to ith (pre-synaptic) cell.
#				Col j corresponds to the jth (post-synaptic) cell.
#			numPreSyn: jth element is number of output connections made by pre-synaptic cell j.
#			
#			When normalizing based on outputs: normalize w by row.
#
#			Difference from NormalizeOutputWeights is that here we allow the cell to
#			have a greater than unity value to which to normalize.
#				
#
######################################################################################################
######################################################################################################

NormalizeOutputWeights2 = function ( w, numPreSyn, wResource ) {

	N2 = dim(w)[1];
	
	w = ( w / matrix ( ( apply ( w, 1, sum ) / ( (wResource * max ( numPreSyn ) / numPreSyn) ) ),
			nrow=N2, ncol=N2, byrow = FALSE ) );

	return ( w );

} # NormalizeOutputWeights2 = function ( w, normByInput )


######################################################################################################
######################################################################################################
#
#	NormalizeInputWeights2: output weight normalization.
#	Input: 	w: N x N matrix.
#				Row i corresponds to ith (pre-synaptic) cell.
#				Col j corresponds to the jth (post-synaptic) cell.
#			
#			When normalizing based on inputs: normalize w by col.
#	This routine corrects (or tests?) a possible logic error in version 1.
#
#	NormalizeInputWeights2 fixes a possible logic error in NormalizeInputWeights1 while preserving
#	the option to go back if it turns out not to be.
#				
#
######################################################################################################
######################################################################################################

NormalizeInputWeights2 = function ( w, wResource ) {

	N2 = dim(w)[1];

	numInputWeights = apply ( (w!=0), 2, sum );
	maxNumInputWeights = max ( numInputWeights );
	w = ( w / matrix ( ( apply ( w, 2, sum ) / ( (wResource * maxNumInputWeights ) / numInputWeights ) ),
		nrow=N2, ncol=N2, byrow = TRUE ) );

	return ( w );

} # NormalizeInputWeights2 = function ( w, wResource )


######################################################################################################
######################################################################################################
#
#	NormalizeInputWeights1: output weight normalization.
#	Input: 	w: N x N matrix.
#				Row i corresponds to ith (pre-synaptic) cell.
#				Col j corresponds to the jth (post-synaptic) cell.
#			
#			When normalizing based on inputs: normalize w by col.
#
#				
#
######################################################################################################
######################################################################################################

NormalizeInputWeights1 = function ( w, wResource ) {

	N2 = dim(w)[1];

		#	Inputs for post-synaptic cell j are stored in jth column of w.  Count number of non-zero inputs.
	tmp = apply ( (w!=0), 2, sum );

		#	There is a chance that some post-synaptic cells have no inputs (e.g., amputation case.
		#	Set the count to 1.0 -- this has only the effect to prevent divide-by-zero in last step.
	tmp [ tmp == 0 ] = 1.0;
	tmp2 = ( apply ( w, 2, sum ) / ( (wResource * tmp) / max ( tmp ) ) );
	tmp2 [ tmp2 == 0 ] = 1.0;

		#	Do the normalization using vectorization.
		#	Those weights that are zero remain zero, because weight value is in numerator.
	w = ( w / matrix ( tmp2, nrow=N2, ncol=N2, byrow = TRUE ) );

	return ( w );

} # NormalizeInputWeights1 = function ( w, wResource )

######################################################################################################
######################################################################################################
#
#	NMDispatch: Call the function corresponding to nmSelect valu.
#		nmSelect = 1:
#			NO Input to Cortical E-Cell.
#			NO: Cortical I-Cell to Cortical E-Cell.
#			NO: Cortical E-Cell to Cortical I-Cell.
#			NO: Cortical E-Cell to Cortical E-Cell.
#
#		nmSelect = 2:
#			YES: Input to Cortical E-Cell.
#			NO: Cortical I-Cell to Cortical E-Cell.
#			NO: Cortical E-Cell to Cortical I-Cell.
#			NO: Cortical E-Cell to Cortical E-Cell.
#
#		nmSelect = 3:
#			YES: Input to Cortical E-Cell.
#			YES: Cortical I-Cell to Cortical E-Cell.
#			NO: Cortical E-Cell to Cortical I-Cell.
#			NO: Cortical E-Cell to Cortical E-Cell.
#
#		nmSelect = 4:
#			YES: Input to Cortical E-Cell.
#			YES: Cortical I-Cell to Cortical E-Cell.
#			YES: Cortical E-Cell to Cortical I-Cell.
#			NO: Cortical E-Cell to Cortical E-Cell.
#
#		nmSelect = 5:
#			YES: Input to Cortical E-Cell.
#			YES: Cortical I-Cell to Cortical E-Cell.
#			YES: Cortical E-Cell to Cortical I-Cell.
#			YES: Cortical E-Cell to Cortical E-Cell.
#
######################################################################################################
######################################################################################################

NMDispatchOrig = function ( nmSelect, v0, r0, v1.e, r1.e, v1.i, r1.i,
					v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta, plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ) {

	if ( nmSelect == 1 ) {
		return ( NM1 ( v0, r0, v1.e, r1.e, v1.i, r1.i, v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta, plsticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ) );
	} else if ( nmSelect == 2 ) {
		return ( NM2 ( v0, r0, v1.e, r1.e, v1.i, r1.i, v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta, plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ) );
	} else if ( nmSelect == 3 ) {
		return ( NM3 ( v0, r0, v1.e, r1.e, v1.i, r1.i, v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta, plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ) );
	} else if ( nmSelect == 4 ) {
		return ( NM4 ( v0, r0, v1.e, r1.e, v1.i, r1.i, v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta, plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ));
	} else if ( nmSelect == 5 ) {
		return ( NM5 ( v0, r0, v1.e, r1.e, v1.i, r1.i, v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta,
					plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ));
	} else {
		return ( NA );
	} # if if ( nmSelect == 1 ) {

} # NMDispatch = function ( v0, r0, v1.e, r1.e,...

NMDispatch = function ( nmSelect, v0, r0, v1.e, r1.e, v1.i, r1.i,
					v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta, plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ) {

	if ( nmSelect == 4 ) {
		return ( NM5 ( v0, r0, v1.e, r1.e, v1.i, r1.i, v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta,
				plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ));
	} else {
		return ( NM5 ( v0, r0, v1.e, r1.e, v1.i, r1.i, v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta,
				plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ));
	} # if ( nmSelect ==4 )

} # NMDispatch = function ( v0, r0, v1.e, r1.e,...

######################################################################################################
######################################################################################################
#
#	NM1: Dirt simple 2-layer model.
#		NO: Input to Cortical E-Cell.
#		NO: Cortical I-Cell to Cortical E-Cell.
#		NO: Cortical E-Cell to Cortical I-Cell.
#		NO: Cortical E-Cell to Cortical E-Cell.
#
######################################################################################################
######################################################################################################

NM1 = function ( v0, r0, v1.e, r1.e, v1.i, r1.i, v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta, plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ) {

	v0 = v0.alpha * v0 + runif(N2, min = -noiseLevel, max = noiseLevel);
	#v0 = VectorNormalize ( v0.alpha * v0 + runif(N2, min = -noiseLevel, max = noiseLevel) );

	v1.e = v1.e.alpha * v1.e + runif(N2, min = -noiseLevel, max = noiseLevel);

	v1.i = v1.i.alpha * v1.i  + runif(N2, min = -noiseLevel, max = noiseLevel);

	r0 = sigmoid ( v0, beta );

	r1.e = sigmoid ( v1.e, beta );

	r1.i = sigmoid ( v1.i, beta );

	return ( list ( v0=v0, v1.e=v1.e, v1.i=v1.i, r0=r0, r1.e=r1.e, r1.i=r1.i, w1.e.0=w1.e.0, w1.e.i=w1.e.i, w1.i.e=w1.i.e, w1.e.e=w1.e.e ) );

} # NM1 = function ( v0, v1.e, v1.i, ...

######################################################################################################
######################################################################################################
#
#	NM2: Dirt simple 2-layer model.
#		YES: Input to Cortical E-Cell.
#		NO: Cortical I-Cell to Cortical E-Cell.
#		NO: Cortical E-Cell to Cortical I-Cell.
#		NO: Cortical E-Cell to Cortical E-Cell.
#
######################################################################################################
######################################################################################################

NM2 = function ( v0, r0, v1.e, r1.e, v1.i, r1.i, v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta, plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ) {

	v0 = v0.alpha * v0 + runif(N2, min = -noiseLevel, max = noiseLevel);
	#v0 = VectorNormalize ( v0.alpha * v0 + runif(N2, min = -noiseLevel, max = noiseLevel) );
	v1.e = v1.e.alpha * v1.e + runif(N2, min = -noiseLevel, max = noiseLevel) + r0 %*% w1.e.0;
	v1.i = v1.i.alpha * v1.i  + runif(N2, min = -noiseLevel, max = noiseLevel);
	
	if ( plasticityFlag ) {
		w1.e.0 = w.tau.alpha * w1.e.0 + w.beta * cbind(r0) %*% r1.e * (w1.e.0 != 0);
	} # if ( plasticity )

	r0 = sigmoid ( v0, beta );
	r1.e = sigmoid ( v1.e, beta );
	r1.i = sigmoid ( v1.i, beta );

	return ( list ( v0=v0, v1.e=v1.e, v1.i=v1.i, r0=r0, r1.e=r1.e, r1.i=r1.i, w1.e.0=w1.e.0, w1.e.i=w1.e.i, w1.i.e=w1.i.e, w1.e.e=w1.e.e ) );

} # NM2 = function ( v0, r0, v1.e, v1.i, ...

######################################################################################################
######################################################################################################
#
#	NM3: Dirt simple 2-layer model.
#		YES: Input to Cortical E-Cell.
#		YES: Cortical I-Cell to Cortical E-Cell.
#		NO: Cortical E-Cell to Cortical I-Cell.
#		NO: Cortical E-Cell to Cortical E-Cell.
#
#
#
######################################################################################################
######################################################################################################

NM3 = function ( v0, r0, v1.e, r1.e, v1.i, r1.i, v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta, plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ) {

	v0 = v0.alpha * v0 + runif(N2, min = -noiseLevel, max = noiseLevel);
	#v0 = VectorNormalize ( v0.alpha * v0 + runif(N2, min = -noiseLevel, max = noiseLevel) );
	v1.e = v1.e.alpha * v1.e + runif(N2, min = -noiseLevel, max = noiseLevel) + r0 %*% w1.e.0 - r1.i %*% w1.e.i;
	v1.i = v1.i.alpha * v1.i  + runif(N2, min = -noiseLevel, max = noiseLevel);

	if ( plasticityFlag ) {
		w1.e.0 = w.tau.alpha * w1.e.0 + w.beta * cbind(r0) %*% r1.e * (w1.e.0 != 0);
		w1.e.i = w.tau.alpha * w1.e.i + w.beta * cbind(r1.i) %*% r1.e * (w1.e.i != 0);
	} # if ( plasticity )

	r0 = sigmoid ( v0, beta );
	r1.e = sigmoid ( v1.e, beta );
	r1.i = sigmoid ( v1.i, beta );

	return ( list ( v0=v0, v1.e=v1.e, v1.i=v1.i, r0=r0, r1.e=r1.e, r1.i=r1.i, w1.e.0=w1.e.0, w1.e.i=w1.e.i, w1.i.e=w1.i.e, w1.e.e=w1.e.e ) );

} # NM3 = function ( v0, r0, v1.e, v1.i, ...

######################################################################################################
######################################################################################################
#
#	NM4: Dirt simple 2-layer model.
#		YES: Input to Cortical E-Cell.
#		YES: Cortical I-Cell to Cortical E-Cell.
#		YES: Cortical E-Cell to Cortical I-Cell.
#		NO: Cortical E-Cell to Cortical E-Cell.
#
######################################################################################################
######################################################################################################

NM4 = function ( v0, r0, v1.e, r1.e, v1.i, r1.i, v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta, plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ) {

	v0 = v0.alpha * v0 + runif(N2, min = -noiseLevel, max = noiseLevel);
	#v0 = VectorNormalize ( v0.alpha * v0 + runif(N2, min = -noiseLevel, max = noiseLevel) );
	v1.e = v1.e.alpha * v1.e + runif(N2, min = -noiseLevel, max = noiseLevel) + r0 %*% w1.e.0 - r1.i %*% w1.e.i;
	v1.i = v1.i.alpha * v1.i  + runif(N2, min = -noiseLevel, max = noiseLevel) + r1.e %*% w1.i.e;

	if ( plasticityFlag ) {
		w1.e.0 = w.tau.alpha * w1.e.0 + w.beta * cbind(r0) %*% r1.e * (w1.e.0 != 0);
		w1.e.i = w.tau.alpha * w1.e.i + w.beta * cbind(r1.i) %*% r1.e * (w1.e.i != 0);
		w1.i.e = w.tau.alpha * w1.i.e + w.beta * cbind(r1.e) %*% r1.i * (w1.i.e != 0);
	} # if ( plasticity )

	r0 = sigmoid ( v0, beta );
	r1.e = sigmoid ( v1.e, beta );
	r1.i = sigmoid ( v1.i, beta );

	return ( list ( v0=v0, v1.e=v1.e, v1.i=v1.i, r0=r0, r1.e=r1.e, r1.i=r1.i, w1.e.0=w1.e.0, w1.e.i=w1.e.i, w1.i.e=w1.i.e, w1.e.e=w1.e.e ) );

} # NM4 = function ( v0, r0, v1.e, v1.i, ...

######################################################################################################
######################################################################################################
#
#	NM5: Dirt simple 2-layer model.
#		YES: Input to Cortical E-Cell.
#		YES: Cortical I-Cell to Cortical E-Cell.
#		YES: Cortical E-Cell to Cortical I-Cell.
#		YES: Cortical E-Cell to Cortical E-Cell.
#
######################################################################################################
######################################################################################################

NM5 = function ( v0, r0, v1.e, r1.e, v1.i, r1.i, v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta, plasticityFlag, w.tau.alpha, w.beta, w1.e.0, w1.e.i, w1.i.e, w1.e.e ) {

	v0 = v0.alpha * v0 + runif(N2, min = -noiseLevel, max = noiseLevel);
	#v0 = VectorNormalize ( v0.alpha * v0 + runif(N2, min = -noiseLevel, max = noiseLevel) );
	v1.e = v1.e.alpha * v1.e + runif(N2, min = -noiseLevel, max = noiseLevel) + r0 %*% w1.e.0 - r1.i %*% w1.e.i + r1.e %*% w1.e.e;
	v1.i = v1.i.alpha * v1.i + runif(N2, min = -noiseLevel, max = noiseLevel) + r1.e %*% w1.i.e;

	if ( plasticityFlag ) {
		w1.e.0 = w.tau.alpha * w1.e.0 + w.beta * cbind(r0) %*% r1.e * (w1.e.0 != 0);
		w1.e.i = w.tau.alpha * w1.e.i + w.beta * cbind(r1.i) %*% r1.e * (w1.e.i != 0);
		w1.i.e = w.tau.alpha * w1.i.e + w.beta * cbind(r1.e) %*% r1.i * (w1.i.e != 0);
		w1.e.e = w.tau.alpha * w1.e.e + w.beta * cbind(r1.e) %*% r1.e * (w1.e.e != 0);
	} # if ( plasticityFlag )

	r0 = sigmoid ( v0, beta );
	r1.e = sigmoid ( v1.e, beta );
	r1.i = sigmoid ( v1.i, beta );

	return ( list ( v0=v0, v1.e=v1.e, v1.i=v1.i, r0=r0, r1.e=r1.e, r1.i=r1.i, w1.e.0=w1.e.0, w1.e.i=w1.e.i, w1.i.e=w1.i.e, w1.e.e=w1.e.e ) );

} # NM5 = function ( v0, r0, v1.e, v1.i, ...

######################################################################################################
######################################################################################################
#
#	ExperimentalRFMap: Run the experiment that generates the raw data from which to build an
#					RM Map.
#
######################################################################################################
######################################################################################################

ExperimentalRFMap = function ( ) {
} # ExperimentalRFMap = function ( ) {

######################################################################################################
######################################################################################################
#
#	GenRFMap : given rfMapRaw data generate a receptive field
#			map for each cell.
#
######################################################################################################
######################################################################################################

GenRFMap1 = function ( rfMapRaw, oneSecondNumIter, startTimeInIterToRFProbeInIters, durRFProbeInIters ) {

	nExpCells = dim(rfMapRaw[[1]]$r0)[2];
	nInputCells = length(rfMapRaw);

	r1.e.rfMap = matrix ( 0, nrow=nExpCells, ncol=nInputCells );
	r1.e.rfMap.noiseEst = matrix ( 0, nrow=nExpCells , ncol=nInputCells );
	r1.e.rfMap.response = matrix ( 0, nrow=nExpCells , ncol=nInputCells );

	r1.i.rfMap = matrix ( 0, nrow=nExpCells , ncol=nInputCells );
	r1.i.rfMap.noiseEst = matrix ( 0, nrow=nExpCells , ncol=nInputCells );
	r1.i.rfMap.response = matrix ( 0, nrow=nExpCells , ncol=nInputCells );

	#preProbeSeries = seq( (startTimeInIterToRFProbeInIters - oneSecondNumIter), (startTimeInIterToRFProbeInIters - 1) );
	#expProbeSeries = seq( startTimeInIterToRFProbeInIters, (startTimeInIterToRFProbeInIters + durRFProbeInIters - 1) );

	preProbeSeries = seq( 1, startTimeInIterToRFProbeInIters );
	expProbeSeries = seq( (startTimeInIterToRFProbeInIters + 1), (startTimeInIterToRFProbeInIters + durRFProbeInIters) );

	for ( iCell in 1:nInputCells ) {

		r1.e.rfMap.noiseEst[,iCell] = apply(rfMapRaw[[iCell]]$r1.e[preProbeSeries,],2,mean);
		r1.e.rfMap.response[,iCell] = apply(rfMapRaw[[iCell]]$r1.e[expProbeSeries,],2,mean);
		r1.e.rfMap[,iCell] =  ( r1.e.rfMap.response[,iCell] / r1.e.rfMap.noiseEst[,iCell] );

		r1.i.rfMap.noiseEst[,iCell] = apply(rfMapRaw[[iCell]]$r1.i[preProbeSeries,],2,mean);
		r1.i.rfMap.response[,iCell] = apply(rfMapRaw[[iCell]]$r1.i[expProbeSeries,],2,mean);
		r1.i.rfMap[,iCell] =  ( r1.i.rfMap.response[,iCell] / r1.i.rfMap.noiseEst[,iCell] );

	} # for ( iCell in 1:nInputCells ) {
	
	return ( list ( r1.e.rfMap=r1.e.rfMap, r1.i.rfMap=r1.i.rfMap ) );

} # GenRFMap1 = function ( rfMapRaw ) {

######################################################################################################
######################################################################################################
#
#	QuantRFSize : Display in standardized form time series from each major cell class.
#		Input: rfMap: ith row of rfMap gives the Receptive Field map for ith experimental cell.
#			 startCell, endCell: sometimes want to look at one; other times the batch.
#			 kRFPeakToEdgeDetect: call a point "in" if the rfMap activity > ratio peak to edge.
#
#		Output: vector where each component is an integer count of the number of cells with
#				values > kRFPeakToEdgeDetect
#
######################################################################################################
######################################################################################################

QuantRFSize = function ( rfData, kRFPeakToEdgeDetect ) {

	N2 = dim(rfData)[1];	# Tells # cells in receptive field experiment, e.g., all Layer 1 cells.
	return ( apply ( (rfData > matrix ( apply ( rfData, 1, max ) * kRFPeakToEdgeDetect, ncol=N2, nrow=N2, byrow=FALSE )), 1, sum ) );

} # QuantRFSize = function ( rfData, kRFPeakToEdgeDetect )

	#	Modified to apply a minimum threshold to response magnitude.
QuantRFSizeB = function ( rfData, kRFPeakToEdgeDetect, minMagResponse=2.0 ) {

	N2 = dim(rfData)[1];	# Tells number celxls in receptive field experiment, e.g., all Layer 1 cells.
	threshResponse = apply ( rfData, 1, max ) * kRFPeakToEdgeDetect;
	tmp = apply ( rfData > threshResponse, 1, sum );
	tmp[ which(threshResponse < minMagResponse) ] = 0;	
	return ( tmp );

} # QuantRFSize = function ( rfData, kRFPeakToEdgeDetect )

######################################################################################################
######################################################################################################
#
#	QuantCorticalAmp : Compute cortical amplification.
#		For each input layer node, determine the extent of the output layer whose cells contain
#		within their receptive field that input layer node.
#		Input: rfData generated by GenRFMap1
#		Output:
#			Boolean matrix where each COLUMN corresponds to an input node and the rows
#			correspond to the output nodes.
#
#			Reading by row gives the receptive field for an OUTPUT layer node.
#			Reading by column gives the cortical amplification for an INPUT layer node.
#
#
######################################################################################################
######################################################################################################

QuantCorticalAmp = function ( rfData, kRFPeakToEdgeDetect ) {

	N2 = dim(rfData)[1];	# Tells # cells in receptive field experiment, e.g., all Layer 1 cells.
	return (		( (rfData >  matrix ( apply ( rfData, 1, max ) * kRFPeakToEdgeDetect , ncol=N2, nrow=N2, byrow=FALSE ))));

} # QuantCorticalAmp = function ( rfData, kRFPeakToEdgeDetect )

######################################################################################################
######################################################################################################
#
#	ShowTimeSeries : Display in standardized form time series
#				from each major cell class.
#	Recall that for membrane potential and firing rate, the
#	NxN matrices are stored in linear format.
#
######################################################################################################
######################################################################################################

ShowTimeSeries = function ( v0, r0, v1.e, v1.i, nTrials, trialDurRFProbeInIters, startOffsetInIter, cellID, knmSelect, deltaT ) {

	N = sqrt ( length ( v0[1,] ) );
	seqEnd = startOffsetInIter+nTrials*trialDurRFProbeInIters-1;
	seqStart = startOffsetInIter + 1;
	tVal = (seqStart:seqEnd) * deltaT;

	for ( i in 1:length(cellID) ) {
		x11();
		par(mfrow=c(3,2));
		iCell = cellID[i];		

		plot( tVal, v1.e[(seqStart:seqEnd),iCell], type="l", col=1, ylab="Membrane Potential", xlab="Time (sec)" );
		title(main=paste("Cortical Layer E-Cell Membrane Potential",paste("Cell ",iCell," Exp ",knmSelect),sep='\n'));

		plot( tVal, r1.e[(seqStart:seqEnd),iCell], type="l", col=1, ylab="Firing Rate", xlab="Time (sec)");
		title(main=paste("Cortical Layer E-Cell Firing Rate",paste("Cell ",iCell," Exp ",knmSelect),sep='\n'));

		plot( tVal, v1.i[(seqStart:seqEnd),iCell], type="l", col=1, ylab="Membrane Potential", xlab="Time (sec)" );
		title(main=paste("Cortical Layer I-Cell Membrane Potential",paste("Cell ",iCell," Exp ",knmSelect),sep='\n'));

		plot( tVal, r1.i[(seqStart:seqEnd),iCell], type="l", col=1, ylab="Firing Rate", xlab="Time (sec)");
		title(main=paste("Cortical Layer I-Cell Firing Rate",paste("Cell ",iCell," Exp ",knmSelect),sep='\n'));

		plot( tVal, v0[(seqStart:seqEnd),iCell], type="l", col=1, ylab="Membrane Potential", xlab="Time (sec)" );
		title(main=paste("Input Layer Membrane Potential",paste("Cell ",iCell," Exp ",knmSelect),sep='\n'));

		plot( tVal, r0[(seqStart:seqEnd),iCell], type="l", col=1, ylab="Firing Rate", xlab="Time (sec)");
		title(main=paste("Input Layer Firing Rate",paste("Cell ",iCell," Exp ",knmSelect),sep='\n'));

	} # for ( i in 1:length(cellID) ) {
	
} # ShowTimeSeries = function ( v0, r0, v1.e, v1.i, nTrials, trialDurRFProbeInIters, startTimeInIter ) {

######################################################################################################
######################################################################################################
#
#	Dig3MapRefine0: Fixed stimulus size, no overlap, but presented in random sequence.
#				This is to ensure uniformity of exposure.
#				No double-digit stimulus.
#				Ensure that each possible fixed size stimulus is generated.
#
#		Input:
#			N - Input Layer is N x N
#			trialLengthInIters - simulation runs in trial or epoch lengths
#
#		Output:
#			Batch of input stimulii.
#			v0 and r0 that will be substituted in as the simulation evolves.
#
#		Notes:
#			Burying the parameters that set details of the refinement inputs in this routine.
#
######################################################################################################
######################################################################################################

Dig3MapRefine0= function ( N, trialLengthInIters, oneSecondNumIter, patchSize, patchDuration ) {

	numDigits = 3;
	digitWidth = N / numDigits;
	N2 = N^2;

		#	Fix the patch width.
	widthPatch = patchSize;

		#	Fix the patch length.
	lengthPatch = patchSize;

		#	Fix the amplitude - so that the resulting input vector is normalized.	
	ampPatch = 1.0 / (patchSize+1);

		#	Fix the duration.
	durationPatchInIters = ceiling ( patchDuration * oneSecondNumIter );

		#	Fix the offset.
	offsetPatchInIters = trialLengthInIters / 5;

	timeMask = seq ( from = offsetPatchInIters, to = (offsetPatchInIters + durationPatchInIters - 1), by = 1 );

		#	Count the number of times each S layer unit is stimulated.
	stimCount = rep ( 0, N2 );

	inStim = list();
	iStim = 1;

	for ( iDigit in 1:numDigits ) {

		minLegalColumn = ( iDigit - 1 ) * digitWidth + 1;
		maxLegalColumn = iDigit * digitWidth - widthPatch;

		for ( iCol in minLegalColumn:maxLegalColumn ) {

			minLegalRow = 1;
			maxLegalRow = N - lengthPatch;
			colMask = c ( iCol, iCol + patchSize );

			for ( iRow in minLegalRow:maxLegalRow ) {

				rowMask = c ( iRow, iRow + lengthPatch );


						#	Construct the spatio-temporal image of the input.
						#		ith row of inputPatch is a linear representation of the NxN input layer at time i
						#		jth col of inputPatch is the time evolution of cell j on the NxN input layer
				inputPatch = matrix ( 0, nrow = trialLengthInIters, ncol = N2 );	# initialize to zero; calling routine may/may not use additive noise
				for ( iiCol in (colMask[1]:colMask[2]) ) {
					for ( iiRow in (rowMask[1]:rowMask[2]) ) {
						iiWhich = GetLin(iiRow, iiCol, N);
						inputPatch[timeMask,iiWhich] = ampPatch;
						stimCount[iiWhich] = stimCount[iiWhich] + 1;
					} # for ( iiRow in (rowMask[1]:rowMask[2]) ) {
				} # for ( iiCol in (colMask[1]:colMask[2]) ) {
				inStim[[iStim]] = list ( inputPatch=inputPatch, digit=iDigit, columnPatch=iCol, rowPatch=iRow, 
									widthPatch=widthPatch, lengthPatch=lengthPatch, ampPatch=ampPatch,
									durationPatchInIters=durationPatchInIters, offsetPatchInIters=offsetPatchInIters,
									stimCount=stimCount );
				iStim = iStim + 1;

			} # for ( iRow in minLegalRow:maxLegalRow ) {

		} # for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} # for ( iDigit in 1:numDigits ) {

	return ( inStim );			

} # Dig3MapRefine0= function ( N, trialLengthInIters, oneSecondNumIter ) {

######################################################################################################
######################################################################################################
#
#	Dig3MapRefine1: Generate random refinement inputs for a simulated three digit input layer.
#				No double-digit stimulus.
#
#		Input:
#			N - Input Layer is N x N
#			trialLengthInIters - simulation runs in trial or epoch lengths
#
#		Output:
#			v0 and r0 that will be substituted in as the simulation evolves.
#
#		Notes:
#			Burying the parameters that set details of the refinement inputs in this routine.
#
######################################################################################################
######################################################################################################

Dig3MapRefine1= function ( N, trialLengthInIters, oneSecondNumIter ) {

	numDigits = 3;

	N2 = N^2;
	digitWidth = N / numDigits;

		#	Specify the types of stimulus patches.
	#lengthPatchList = seq ( from=1, to=N, by=1 );
	#lengthPatchList = seq ( from=1, to=digitWidth, by=1 );
	lengthPatchList = seq ( from=digitWidth, to=digitWidth, by=1 );
	lenLengthPatchList = length ( lengthPatchList );

	#widthPatchList = seq ( from=1, to=digitWidth, by=1 );
	widthPatchList = seq ( from=digitWidth, to=digitWidth, by=1 );
	lenWidthPatchList = length ( widthPatchList );

	ampRangePatchList = c ( 0.9, 1.0 );		# min and max amplitude range
	durationPatchPCntTrialLengthList = ceiling( oneSecondNumIter * ( c ( 0.25, 0.75 ) ) );
	offsetPatchList = c ( oneSecondNumIter, trialLengthInIters - 2 * oneSecondNumIter );

		#	Choose a patch width.
	widthPatch = widthPatchList[ ceiling ( runif ( 1, min = 0, max = lenWidthPatchList ) ) ];

		#	Choose a patch length.
	lengthPatch = lengthPatchList[ ceiling ( runif ( 1, min = 0, max = lenLengthPatchList ) ) ];

		#	Choose a digit
	digit = ceiling ( runif ( 1, min = 0, max = numDigits ) );

		#	Identify the columns where a patch can legally (i.e., fully) land.
	minLegalColumn = ( digit - 1 ) * digitWidth + 1;
	maxLegalColumn = minLegalColumn + digitWidth - widthPatch;
	columnPatch = ceiling ( runif ( 1, min = minLegalColumn - 1, max = maxLegalColumn ) );

		#	Identify the rows where a patch can legally (i.e., fully) land.
	minLegalRow = 1;
	maxLegalRow = 1 + N - lengthPatch;
	rowPatch = ceiling ( runif ( 1, min = minLegalRow - 1, max = maxLegalRow ) );

		#	Choose an amplitude.
	ampPatch = runif ( 1, min = ampRangePatchList[1], max = ampRangePatchList[2] );

		#	Choose a duration.
	durationPatchInIters = ceiling ( runif ( 1, min = durationPatchPCntTrialLengthList[1] - 1, max = durationPatchPCntTrialLengthList[2] ) );

		#	Choose an offset.
	offsetPatchInIters = ceiling ( runif ( 1, min = offsetPatchList[1]-1, max = offsetPatchList[2] ) );

		#	Construct the spatio-temporal image of the input.
		#		ith row of inputPatch is a linear representation of the NxN input layer at time i
		#		jth col of inputPatch is the time evolution of cell j on the NxN input layer

	inputPatch = matrix ( 0, nrow = trialLengthInIters, ncol = N2 );	# initialize to zero; calling routine may/may not use additive noise
	colMask = c ( columnPatch, columnPatch + widthPatch - 1 );
	rowMask = c ( rowPatch, rowPatch + lengthPatch - 1 );
	timeMask = seq ( from = offsetPatchInIters, to = (offsetPatchInIters + durationPatchInIters - 1), by = 1 );
	for ( iCol in (colMask[1]:colMask[2]) ) {
		for ( iRow in (rowMask[1]:rowMask[2]) ) {
			inputPatch[timeMask,GetLin(iRow, iCol, N)] = ampPatch;
		} # for ( iRow in (rowMask[1]:rowMask[2]) ) {
	} # for ( iCol in (colMask[1]:colMask[2]) ) {

	return ( list ( inputPatch=inputPatch, digit=digit, columnPatch=columnPatch, rowPatch=rowPatch, 
				widthPatch=widthPatch, lengthPatch=lengthPatch, ampPatch=ampPatch,
				durationPatchInIters=durationPatchInIters, offsetPatchInIters=offsetPatchInIters ) );			

} # Dig3MapRefine1= function ( N, trialLengthInIters, oneSecondNumIter ) {


######################################################################################################
######################################################################################################
#
#	Dig3MapRefine2: Fixed stimulus size, no overlap, but presented in random sequence.
#				This is to ensure uniformity of exposure.
#				No double-digit stimulus.
#				Ensure that each possible fixed size stimulus is generated.
#
#	Derived from Dig3MapRefine0: adding some flexible patchsize handling.
#
#		Input:
#			N - Input Layer is N x N
#			trialLengthInIters - simulation runs in trial or epoch lengths
#
#		Output:
#			Batch of input stimulii.
#			v0 and r0 that will be substituted in as the simulation evolves.
#
#		Notes:
#			Burying the parameters that set details of the refinement inputs in this routine.
#
######################################################################################################
######################################################################################################

Dig3MapRefine2= function ( N, trialLengthInIters, oneSecondNumIter, patchSize, noiseLevel ) {

	numDigits = 3;
	digitWidth = N / numDigits;
	N2 = N^2;

		#	Fix the patch width.
	widthPatch = patchSize + 1;

		#	Fix the patch length.
	lengthPatch = patchSize + 1;

		#	Fix the amplitude - so that the resulting input vector is normalized.	
	ampPatch = 1.0 / (widthPatch);

		#	Fix the duration.
	durationPatchInIters = ceiling ( 0.2 * oneSecondNumIter );

		#	Fix the offset.
	offsetPatchInIters = trialLengthInIters / 5;

	numPatches = numDigits * (digitWidth - widthPatch + 1) * (N - widthPatch + 1);
	timeMask = seq ( from = offsetPatchInIters, to = (offsetPatchInIters + durationPatchInIters - 1), by = 1 );

	inStim = list();
	iStim = 1;

	for ( iDigit in 1:numDigits ) {

		minLegalColumn = ( iDigit - 1 ) * digitWidth + 1;
		maxLegalColumn = 1 + iDigit * digitWidth - widthPatch;

		for ( iCol in minLegalColumn:maxLegalColumn ) {

			minLegalRow = 1;
			maxLegalRow = N - lengthPatch + 1;
			colMask = c ( iCol, iCol + patchSize );

			for ( iRow in minLegalRow:maxLegalRow ) {

				rowMask = c ( iRow, iRow + patchSize );


						#	Construct the spatio-temporal image of the input.
						#		ith row of inputPatch is a linear representation of the NxN input layer at time i
						#		jth col of inputPatch is the time evolution of cell j on the NxN input layer
				inputPatch = matrix ( 0, nrow = trialLengthInIters, ncol = N2 );	# initialize to zero; calling routine may/may not use additive noise
				for ( iiCol in (colMask[1]:colMask[2]) ) {
					for ( iiRow in (rowMask[1]:rowMask[2]) ) {
						inputPatch[timeMask,GetLin(iiRow, iiCol, N)] = ampPatch;
					} # for ( iiRow in (rowMask[1]:rowMask[2]) ) {
				} # for ( iiCol in (colMask[1]:colMask[2]) ) {
				inStim[[iStim]] = list ( inputPatch=inputPatch, digit=iDigit, columnPatch=iCol, rowPatch=iRow, 
									widthPatch=widthPatch, lengthPatch=lengthPatch, ampPatch=ampPatch,
									durationPatchInIters=durationPatchInIters, offsetPatchInIters=offsetPatchInIters );
				iStim = iStim + 1;

			} # for ( iRow in minLegalRow:maxLegalRow ) {

		} # for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} # for ( iDigit in 1:numDigits ) {

	return ( inStim );			

} # Dig3MapRefine2= function ( N, trialLengthInIters, oneSecondNumIter ) {

######################################################################################################
######################################################################################################
#
#	Dig3MapSelStim0: Simulate an experiment in which a subset of the inputs are stimulated
#				selectively and preferentially over others.
#
#		Input:
#			N - Input Layer is N x N
#			trialLengthInIters - simulation runs in trial or epoch lengths
#			etc.
#
#		Output:
#			Batch of input stimulii.
#			v0 and r0 that will be substituted in as the simulation evolves.
#
#		Notes:
#			Burying the parameters that set details of the refinement inputs in this routine.
#
######################################################################################################
######################################################################################################

Dig3MapSelStim0= function ( N, trialLengthInIters, oneSecondNumIter, patchSize, patchDuration, stimFactor, zoneID ) {

	numDigits = 3;
	digitWidth = N / numDigits;
	N2 = N^2;

		#	Fix the patch width.
	widthPatch = patchSize;

		#	Fix the patch length.
	lengthPatch = patchSize;

		#	Fix the amplitude - so that the resulting input vector is normalized.	
	ampPatch = 1.0 / (patchSize+1);

		#	Fix the duration.
	durationPatchInIters = ceiling ( patchDuration * oneSecondNumIter );

		#	Fix the offset.
	offsetPatchInIters = trialLengthInIters / 5;

	timeMask = seq ( from = offsetPatchInIters, to = (offsetPatchInIters + durationPatchInIters - 1), by = 1 );

		#	Count the number of times each S layer unit is stimulated.
	stimCount = rep ( 0, N2 );

	inStim = list();
	iStim = 1;

		#
		#	Part I: Background stimulation over the entire input layer
		#

	for ( iDigit in 1:numDigits ) {

		minLegalColumn = ( iDigit - 1 ) * digitWidth + 1;
		maxLegalColumn = iDigit * digitWidth - widthPatch;

		for ( iCol in minLegalColumn:maxLegalColumn ) {

			minLegalRow = 1;
			maxLegalRow = N - lengthPatch;
			colMask = c ( iCol, iCol + patchSize );

			for ( iRow in minLegalRow:maxLegalRow ) {

				rowMask = c ( iRow, iRow + lengthPatch );

						#	Construct the spatio-temporal image of the input.
						#		ith row of inputPatch is a linear representation of the NxN input layer at time i
						#		jth col of inputPatch is the time evolution of cell j on the NxN input layer
				inputPatch = matrix ( 0, nrow = trialLengthInIters, ncol = N2 );	# initialize to zero; calling routine may/may not use additive noise
				for ( iiCol in (colMask[1]:colMask[2]) ) {

					for ( iiRow in (rowMask[1]:rowMask[2]) ) {
						#inputPatch[timeMask,GetLin(iiRow, iiCol, N)] = ampPatch;

						iiWhich = GetLin(iiRow, iiCol, N);
						inputPatch[timeMask,iiWhich] = ampPatch;
						stimCount[iiWhich] = stimCount[iiWhich] + 1;
					} # for ( iiRow in (rowMask[1]:rowMask[2]) ) {
				} # for ( iiCol in (colMask[1]:colMask[2]) ) {
				inStim[[iStim]] = list ( inputPatch=inputPatch, digit=iDigit, columnPatch=iCol, rowPatch=iRow, 
									widthPatch=widthPatch, lengthPatch=lengthPatch, ampPatch=ampPatch,
									durationPatchInIters=durationPatchInIters, offsetPatchInIters=offsetPatchInIters,
									stimCount=stimCount );
				iStim = iStim + 1;

			} # for ( iRow in minLegalRow:maxLegalRow ) {

		} # for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} # for ( iDigit in 1:numDigits ) {

		#
		#	Part II: Selectively stimulate a sub-area.
		#		zoneID selects the sector
		#

	if ( zoneID == 1 || zoneID == 4 || zoneID == 7 ) {
		iDigit = 1;
	} else if ( zoneID == 2 || zoneID == 5 || zoneID == 8 ) {
		iDigit = 2;
	} else {
		iDigit = 3;
	} # if ( zoneID == 1 || zoneID == 4 || zoneID == 7 )

	for ( iFactor in 1:stimFactor ) {

		minLegalColumn = ( iDigit - 1 ) * digitWidth + 1;
		maxLegalColumn = iDigit * digitWidth - widthPatch;

		for ( iCol in minLegalColumn:maxLegalColumn ) {

			colMask = c ( iCol, iCol + patchSize );

			if ( zoneID == 1 || zoneID == 2 || zoneID == 3 ) {
				minLegalRow = 1; maxLegalRow = digitWidth - lengthPatch;
			} else if ( zoneID == 4 || zoneID == 5 || zoneID == 6 ) {
				minLegalRow = digitWidth + 1; maxLegalRow = 2 * digitWidth - lengthPatch;
			} else {
				minLegalRow = 2 * digitWidth + 1; maxLegalRow = 3 * digitWidth - lengthPatch;
			} # if ( zoneID == 1 || zoneID == 4 || zoneID == 7 )
			colMask = c ( iCol, iCol + patchSize );

			for ( iRow in minLegalRow:maxLegalRow ) {

				rowMask = c ( iRow, iRow + lengthPatch );

						#	Construct the spatio-temporal image of the input.
						#		ith row of inputPatch is a linear representation of the NxN input layer at time i
						#		jth col of inputPatch is the time evolution of cell j on the NxN input layer
				inputPatch = matrix ( 0, nrow = trialLengthInIters, ncol = N2 );	# initialize to zero; calling routine may/may not use additive noise
				for ( iiCol in (colMask[1]:colMask[2]) ) {
					for ( iiRow in (rowMask[1]:rowMask[2]) ) {
						#inputPatch[timeMask,GetLin(iiRow, iiCol, N)] = ampPatch;
						iiWhich = GetLin(iiRow, iiCol, N);
						inputPatch[timeMask,iiWhich] = ampPatch;
						stimCount[iiWhich] = stimCount[iiWhich] + 1;
					} # for ( iiRow in (rowMask[1]:rowMask[2]) ) {
				} # for ( iiCol in (colMask[1]:colMask[2]) ) {
				inStim[[iStim]] = list ( inputPatch=inputPatch, digit=iDigit, columnPatch=iCol, rowPatch=iRow, 
									widthPatch=widthPatch, lengthPatch=lengthPatch, ampPatch=ampPatch,
									durationPatchInIters=durationPatchInIters, offsetPatchInIters=offsetPatchInIters,
									stimCount=stimCount );
				iStim = iStim + 1;

			} # for ( iRow in minLegalRow:maxLegalRow ) {

		} # for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} # for ( iFactor in 1:stimFactor ) {

	return ( inStim );			

} # Dig3MapSelStim0= function ( N, trialLengthInIters, oneSecondNumIter, patchSize, patchDuration, noiseLevel ) {


######################################################################################################
######################################################################################################
#
#	Dig3MapSyndact0: Generate inputs to simulate syndactyly.
#				That should look like Dig3Map type input stimulus, but with a double-wide digit
#				width.
#
#		Input:
#			N - Input Layer is N x N
#			trialLengthInIters - simulation runs in trial or epoch lengths
#			patchSize, patchDuration,
#			zoneID: if == 1 then fuse digits 1 and 2; stimulate digit 3 as normal;
#				  if == 2 then fuse digits 2 & 3; stimulate digit 1 as normal.
#
#		Output:
#			Batch of input stimulii.
#			v0 and r0 that will be substituted in as the simulation evolves.
#
#		Notes:
#			Burying the parameters that set details of the refinement inputs in this routine.
#
##############-########################################################################################
######################################################################################################

Dig3MapSyndact0= function ( N, trialLengthInIters, oneSecondNumIter, patchSize, patchDuration, zoneID ) {

	numDigits = 3;
	digitWidth = N / numDigits;
	N2 = N^2;

		#	Fix the patch width.
	widthPatch = patchSize;

		#	Fix the patch length.
	lengthPatch = patchSize;

		#	Fix the amplitude - so that the resulting input vector is normalized.	
	ampPatch = 1.0 / (patchSize+1);

		#	Fix the duration.
	durationPatchInIters = ceiling ( patchDuration * oneSecondNumIter );

		#	Fix the offset.
	offsetPatchInIters = trialLengthInIters / 5;

	timeMask = seq ( from = offsetPatchInIters, to = (offsetPatchInIters + durationPatchInIters - 1), by = 1 );

	inStim = list();
	iStim = 1;

		#	Assume zoneID == 1 (Fused digit 1 & 2; normal digit 3)
	normalDigit = 3;
	if ( zoneID == 2 ) {
		normalDigit = 1;
	} # if ( zoneID == 1 )

		#	Generate the normal digit inputs first.
	iDigit = normalDigit;
	minLegalColumn = ( iDigit - 1 ) * digitWidth + 1;
	maxLegalColumn = iDigit * digitWidth - widthPatch;
	for ( iCol in minLegalColumn:maxLegalColumn ) {

		minLegalRow = 1;
		maxLegalRow = N - lengthPatch;
		colMask = c ( iCol, iCol + patchSize );

		for ( iRow in minLegalRow:maxLegalRow ) {

			rowMask = c ( iRow, iRow + lengthPatch );
					#	Construct the spatio-temporal image of the input.
					#		ith row of inputPatch is a linear representation of the NxN input layer at time i
					#		jth col of inputPatch is the time evolution of cell j on the NxN input layer
			inputPatch = matrix ( 0, nrow = trialLengthInIters, ncol = N2 );	# initialize to zero; calling routine may/may not use additive noise
			for ( iiCol in (colMask[1]:colMask[2]) ) {
				for ( iiRow in (rowMask[1]:rowMask[2]) ) {
					inputPatch[timeMask,GetLin(iiRow, iiCol, N)] = ampPatch;
				} # for ( iiRow in (rowMask[1]:rowMask[2]) ) {
			} # for ( iiCol in (colMask[1]:colMask[2]) ) {
			inStim[[iStim]] = list ( inputPatch=inputPatch, digit=iDigit, columnPatch=iCol, rowPatch=iRow, 
								widthPatch=widthPatch, lengthPatch=lengthPatch, ampPatch=ampPatch,
								durationPatchInIters=durationPatchInIters, offsetPatchInIters=offsetPatchInIters );
			iStim = iStim + 1;

		} # for ( iRow in minLegalRow:maxLegalRow ) {

	} # for ( iCol in minLegalColumn:maxLegalColumn ) {
	
		#	Generate the fused digit inputs.
	if ( zoneID == 1 ) {
			#	Digits 1 and 2 are fused.
		minLegalColumn = 1
		maxLegalColumn = 2 * digitWidth - widthPatch;
	} else {
			#	Digits 2 and 3 are fused
		minLegalColumn = digitWidth + 1;
		maxLegalColumn = 3 * digitWidth - widthPatch;
	} # if ( zoneID == 1 )

	for ( iCol in minLegalColumn:maxLegalColumn ) {

		minLegalRow = 1;
		maxLegalRow = N - lengthPatch;
		colMask = c ( iCol, iCol + patchSize );

		for ( iRow in minLegalRow:maxLegalRow ) {

			rowMask = c ( iRow, iRow + lengthPatch );

					#	Construct the spatio-temporal image of the input.
					#		ith row of inputPatch is a linear representation of the NxN input layer at time i
					#		jth col of inputPatch is the time evolution of cell j on the NxN input layer
			inputPatch = matrix ( 0, nrow = trialLengthInIters, ncol = N2 );	# initialize to zero; calling routine may/may not use additive noise
			for ( iiCol in (colMask[1]:colMask[2]) ) {
				for ( iiRow in (rowMask[1]:rowMask[2]) ) {
					inputPatch[timeMask,GetLin(iiRow, iiCol, N)] = ampPatch;
				} # for ( iiRow in (rowMask[1]:rowMask[2]) ) {
			} # for ( iiCol in (colMask[1]:colMask[2]) ) {
			inStim[[iStim]] = list ( inputPatch=inputPatch, digit=iDigit, columnPatch=iCol, rowPatch=iRow, 
								widthPatch=widthPatch, lengthPatch=lengthPatch, ampPatch=ampPatch,
								durationPatchInIters=durationPatchInIters, offsetPatchInIters=offsetPatchInIters );
			iStim = iStim + 1;

		} # for ( iRow in minLegalRow:maxLegalRow ) {

	} # for ( iCol in minLegalColumn:maxLegalColumn ) {

	return ( inStim );			

} # Dig3MapSyndact0 = function ( N, trialLengthInIters, oneSecondNumIter ) {

	##############################################################
	##############################################################
	#
	#	Dig3Amputate1.R
	#		Modeled starting with a baseline topographic net-
	#		work and then setting input layer to layer 1 E-cell
	#		weights to zero and then following the baseline
	#		refinement protocol.
	#
	#	Inputs:
	#		w: N2 x N2 weight matrix
	#		kAmputDigit: which digit (1,2,3)
	#		kAmputZoneMax: 1=Distal; 2=Distal+Mid; 3=Distal+Mid+Prox
	#
	##############################################################
	##############################################################

Dig3Amputate1 = function ( w, kAmputDigit, kAmputZoneMax ) {

	N2 = dim( w )[1];
	N = as.integer ( sqrt ( N2 ) );

	numDigits = 3;
	digitWidth = N / numDigits;

		# Identify the cell IDs in the amputated zone.  Zero out their outputs.
		# The input weight w are a linear representation of an N x N weight matrix.
		# The ith row of weight matrix w are the output weights of pre-synaptic cell i.
		#	For example, if w is w1.e.0, then row 5 lists the output of input cell 5 to the layer 1 E-cells.

	minLegalColumn = ( kAmputDigit - 1 ) * digitWidth + 1;
	maxLegalColumn = kAmputDigit * digitWidth ;

	for ( iCol in minLegalColumn:maxLegalColumn ) {

		minLegalRow = 1;
		if ( kAmputZoneMax == 1 ) {	# The distal zone amputated.
			maxLegalRow = as.integer ( N / 3 );
		} else if ( kAmputZoneMax == 2 ) {	# The distal and mid zones amputated.
			maxLegalRow = 2 * as.integer ( N / 3 );
		} else {	# Digit amputated.
			maxLegalRow = N;
		} # if ( kAmputZoneMax == 1 ) {

		w[ GetLin ( seq ( minLegalRow, maxLegalRow, 1 ), iCol, N ), , ] = 0;
		
	} # for ( iCol in minLegalColumn:maxLegalColumn ) {

	return ( w );

} # Dig3Amputate1 = function ( w, kAmputDigit, kAmputZoneMax ) {

	##############################################################
	##############################################################
	#
	#	Dig3CLesion1
	#		Modeled starting with a baseline topographic net-
	#		work and then setting the output layer cell output
	#		weights to zero, setting their inputs to zero, etc.
	#		Adopt same convention as used for amputation.  Could
	#		evolve this routine to provide arbitrary lesion location
	#		and size rather than relying on it being a map of
	#		any particular thing.  But the present approach eases
	#		interpretation....
	#
	#	Inputs:
	#		w: N2 x N2 weight matrix
	#		kAmputDigit: which digit (1,2,3)
	#		kAmputZoneMax: 1=Distal; 2=Distal+Mid; 3=Distal+Mid+Prox
	#		zeroTheOuputs, zeroTheInputs -- whether to zero the incoming
	#		or outgoing weights.  In the case of intracortical one has
	#		to zero out both the incoming and outgoing.  For the afferent
	#		just zero the input.
	#
	##############################################################
	##############################################################

Dig3CLesion1= function ( w, kAmputDigit, kAmputZoneMax, zeroTheOutputs, zeroTheInputs ) {

	N2 = dim( w )[1];
	N = as.integer ( sqrt ( N2 ) );

	numDigits = 3;
	digitWidth = N / numDigits;

		# Identify the cell IDs in the amputated zone.  Zero out their outputs.
		# The input weight w are a linear representation of an N x N weight matrix.
		# The ith row of weight matrix w are the output weights of pre-synaptic cell i.
		#	For example, if w is w1.e.0, then row 5 lists the output of input cell 5 to the layer 1 E-cells.

	minLegalColumn = ( kAmputDigit - 1 ) * digitWidth + 1;
	maxLegalColumn = kAmputDigit * digitWidth ;

	for ( iCol in minLegalColumn:maxLegalColumn ) {

		minLegalRow = 1;
		if ( kAmputZoneMax == 1 ) {	# The distal zone amputated.
			maxLegalRow = as.integer ( N / 3 );
		} else if ( kAmputZoneMax == 2 ) {	# The distal and mid zones amputated.
			maxLegalRow = 2 * as.integer ( N / 3 );
		} else {	# Digit amputated.
			maxLegalRow = N;
		} # if ( kAmputZoneMax == 1 ) {

		if ( zeroTheOutputs ) {
			w[ GetLin ( seq ( minLegalRow, maxLegalRow, 1 ), iCol, N ), , ] = 0;
		} # if ( zeroTheOutputs ) {

		if ( zeroTheInputs ) {
			w[ , GetLin ( seq ( minLegalRow, maxLegalRow, 1 ), iCol, N ), ] = 0;
		} # if ( zeroTheInputs ) {
		
	} # for ( iCol in minLegalColumn:maxLegalColumn ) {

	return ( w );

} # Dig3CLesion1= function ( w, kAmputDigit, kAmputZoneMax ) {

######################################################################################################
######################################################################################################
#
#	Dig3MapAmput0: Generate inputs to simulate a digit that has been partially or wholly amputated.
#		Basically same protocol as normal refinement, just skip certain zones.
#
#		Input:
#			N - Input Layer is N x N
#			trialLengthInIters - simulation runs in trial or epoch lengths
#			patchSize, patchDuration, etc.
#
#		Output:
#			Batch of input stimulii.
#			v0 and r0 that will be substituted in as the simulation evolves.
#
#		Notes:
#			Burying the parameters that set details of the refinement inputs in this routine.
#
######################################################################################################
######################################################################################################

Dig3MapAmput0= function ( N, trialLengthInIters, oneSecondNumIter, patchSize, patchDuration, kDigitID, kZoneMax ) {

	numDigits = 3;
	digitWidth = N / numDigits;
	N2 = N^2;

		#	Fix the patch width.
	widthPatch = patchSize;

		#	Fix the patch length.
	lengthPatch = patchSize;

		#	Fix the amplitude - so that the resulting input vector is normalized.	
	ampPatch = 1.0 / (patchSize+1);

		#	Fix the duration.
	durationPatchInIters = ceiling ( patchDuration * oneSecondNumIter );

		#	Fix the offset.
	offsetPatchInIters = trialLengthInIters / 5;

	timeMask = seq ( from = offsetPatchInIters, to = (offsetPatchInIters + durationPatchInIters - 1), by = 1 );

	inStim = list();
	iStim = 1;

	for ( iDigit in 1:numDigits ) {

		minLegalColumn = ( iDigit - 1 ) * digitWidth + 1;
		maxLegalColumn = iDigit * digitWidth - widthPatch;

		for ( iCol in minLegalColumn:maxLegalColumn ) {

			minLegalRow = 1;
			maxLegalRow = N - lengthPatch;
			colMask = c ( iCol, iCol + patchSize );

			workToDoFlag = TRUE;
			if ( iDigit == kDigitID ) {
				if ( kZoneMax == 1 ) {	# The distal zone amputated.
					minLegalRow = as.integer ( N / 3 ) + 1;
				} else if ( kZoneMax == 2 ) {	# The distal and mid zones amputated.
					minLegalRow = as.integer ( ( 2 * N) / 3 ) + 1;
				} else {	# Digit amputated.
					workToDoFlag = FALSE;
				} # if ( )
			} # if ( iDigit == kDigitID )
			if ( workToDoFlag ) {
				for ( iRow in minLegalRow:maxLegalRow ) {

						#	Construct the spatio-temporal image of the input.
						#		ith row of inputPatch is a linear representation of the NxN input layer at time i
						#		jth col of inputPatch is the time evolution of cell j on the NxN input layer
					rowMask = c ( iRow, iRow + lengthPatch );
					inputPatch = matrix ( 0, nrow = trialLengthInIters, ncol = N2 );	# initialize to zero; calling routine may/may not use additive noise
					for ( iiCol in (colMask[1]:colMask[2]) ) {
						for ( iiRow in (rowMask[1]:rowMask[2]) ) {
							inputPatch[timeMask,GetLin(iiRow, iiCol, N)] = ampPatch;
						} # for ( iiRow in (rowMask[1]:rowMask[2]) ) {
					} # for ( iiCol in (colMask[1]:colMask[2]) ) {
					inStim[[iStim]] = list ( inputPatch=inputPatch, digit=iDigit, columnPatch=iCol, rowPatch=iRow, 
									widthPatch=widthPatch, lengthPatch=lengthPatch, ampPatch=ampPatch,
									durationPatchInIters=durationPatchInIters, offsetPatchInIters=offsetPatchInIters );
				iStim = iStim + 1;
				} # for ( iRow in minLegalRow:maxLegalRow ) {
			} # if ( workToDoFlag ) {

		} # for ( iCol in minLegalColumn:maxLegalColumn ) {
		
	} # for ( iDigit in 1:numDigits ) {

	return ( inStim );			

} # Dig3MapAmput0= function ( N, trialLengthInIters, oneSecondNumIter, patchSize, patchDuration, zoneID ) {


######################################################################################################
######################################################################################################
#
#	StimStats: Analyze the set of input stimulus patterns and prepare summary graphs.
#
#		Input:
#			inStim - list of input stimulus patterns
#
#		Output:
#			Graphs.
#
######################################################################################################
######################################################################################################

StimStats = function ( inStim ) {

	numStim = length ( inStim );
	nCells = dim(inStim[[1]]$inputPatch)[2];

	digit = rep ( 0, numStim );
	width = rep ( 0, numStim );
	length = rep ( 0, numStim );
	patchArea = rep ( 0, numStim );
	amplitude = rep ( 0, numStim );
	duration = rep ( 0, numStim );
	offset = rep ( 0, numStim );
	exposure = matrix ( 0, nrow=1, ncol=(nCells) );

	for ( i in 1:numStim ) {
		digit[i] = inStim[[i]]$digit;
		width[i] = inStim[[i]]$widthPatch;
		length[i] = inStim[[i]]$lengthPatch;
		patchArea[i] = inStim[[i]]$widthPatch * inStim[[i]]$lengthPatch;
		amplitude[i] = inStim[[i]]$ampPatch;
		duration[i] = inStim[[i]]$durationPatchInIters;
		offset[i] = inStim[[i]]$offsetPatchInIters;
		exposure = exposure + apply ( inStim[[i]]$inputPatch, 2, sum );	
	} # for ( i in 1:numStim ) {

	x11(); par ( mfrow = c ( 3,2 ) );
	hist ( digit, plot=TRUE, main=paste("Distribution of Inputs (N=",numStim,") per Zone ID"), ylab="Count", xlab="Zone ID" );
	hist ( patchArea, plot=TRUE, main=paste("Distribution of Stimulus Area"), ylab="Count", xlab="Area" );
	hist ( amplitude, plot=TRUE, main=paste("Distribution of Stimulus Amplitudes"), ylab="Count", xlab="Amplitude" );
	hist ( duration, plot=TRUE, main=paste("Distribution of Stimulus Duration"), ylab="Count", xlab="Units of Iterations" );	
	hist ( offset, plot=TRUE, main=paste("Distribution of Stimulus Offsets"), ylab="Count", xlab="Units of Iterations" );
	ViewSTDynamics ( exposure, 1, 1, paste("Distribution of Exposure by Cell") );

} # StimStats = function ( inStim ) {


######################################################################################################
######################################################################################################
#
#	ExtractRFPositionData: display RF progression along set axes of a given RF Mapping Experiment
#	Input: an rfMap, where ith ROW is the response of output cell #i to inputs
#	Output: a list of matrices
#		rfCentroid: length N^2 X 2, where column #1 is x-value, column #2 is y-value of centroid
#		rfCoVar: length N^2 X , where col #1 is x std. dev, col #2 is y std.dev, col #3 is xy correlation
#
#	This kind of output is then used by the plotting routine to make passes along various tracks
#	on the output layer. 
#
######################################################################################################
######################################################################################################

ExtractRFPositionData = function ( rfData, kRFPeakToEdgeDetect ) {

		# 	Get # cells in receptive field experiment, e.g., all Layer 1 output cells.
	N2 = dim(rfData)[1];
	N = as.integer ( sqrt ( N2 ) );

		#	Set up the output matrices.
	rfCentroid = matrix ( 0, nrow = N2, ncol = 2 );
	rfCovar = matrix ( 0, nrow = N2, ncol = 3 );

		# 	For each output cell, identify the subset of input layer cells that produced greater than kRFPeakToEdgeDetect activity.
	rfFootPrint = ( ( (rfData > matrix ( apply ( rfData, 1, max ) * kRFPeakToEdgeDetect , ncol=N2, nrow=N2, byrow=FALSE )) ) );

		#	Compute the RF centroid for each output cell
	for ( iCell in 1:N2 ) {
		rfLocs = cbind ( GetRow ( which ( rfFootPrint[iCell,] ), N ), cbind ( GetCol ( which ( rfFootPrint[iCell,] ), N ) ) );
		if ( dim(rfLocs)[1] > 1 ) {
			rfCentroid[iCell,] = apply ( rfLocs, 2, mean );
			tmp = cov ( rfLocs );
			rfCovar[iCell,1] = tmp[1,1]; rfCovar[iCell,2] = tmp[2,2];
			if ( tmp[1,1] != 0 & tmp[2,2] != 0 ) {
				rfCovar[iCell,3] = cov2cor ( tmp )[1,2];
			} # if ( tmp[1,1] != 0 & tmp[2,2] != 0 ) {
		} else {
			#rfCentroid[iCell,] = rfLocs;
			#rfCovar[iCell,] = rep ( 0, 3 );
		} # if ( dim(rfLocs)[1] > 1 ) {
	} # for ( iCell in 1:N2 )

	return ( list ( rfCentroid=rfCentroid, rfCovar=rfCovar ) );
	
} # ExtractRFPositionData= function ( k, N )


######################################################################################################
######################################################################################################
#
#	ShowThreeDigitRFTrack0: Display RFs along several set tracks on the output layer.
#		Input: List produced by ExtractRFPositionData.
#
#		Output: (Plots 1,3,5 in Column 1 of the plot page; Plots 2,4,6 in Column 2.)
#			Plot(1) - Digit 3 distal to proximal along centerline
#			Plot(2) - Digits 1-3 distal
#			Plot(3) - Digit 2 distal to proximal along centerline
#			Plot(4) - Digits 1-3 mid
#			Plot(5) - Digit 1 distal to proximal along centerline
#			Plot(6) - Digits 1-3 proximal
#
#			This layout should give the right correspondence with other image outputs and visual ease.
#
######################################################################################################
######################################################################################################

ShowThreeDigitRFTrack0= function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit  ) {


	kMaxColorID = 5;
	N = as.integer ( sqrt ( dim(rfTrackData$rfCentroid)[1] ) );
	digitWidth = N / 3;

	x11(); par(mfcol=c(3,2));

		#	Plots 1,3,5 - Digit 3,2,1, respectively, distal to proximal along centerline.

	iColList = seq(3,1,-1) * digitWidth - (floor(digitWidth/2));
	for ( iCol in iColList  ) {
		plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
		kColorID = 0;
		for ( iRow in 1:N ) {
			iTmp = GetLin ( iRow, iCol, N );
			points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch="*", col=(1+(as.integer(kColorID%%kMaxColorID))) );
			if ( showEllipsesFlag ) {
				lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale = c(rfTrackData$rfCovar[iTmp,1], rfTrackData$rfCovar[iTmp,2]),
					centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
					col=(1+(as.integer(kColorID%%kMaxColorID))) );
				kColorID = kColorID + 1;
			} # if ( showEllipsesFlag ) {
		} # for ( irow = 1:N )

			#	For ease of interpretation of the figure, draw dotted lines for digit boundaries.
		abline ( h = (digitWidth + 0.5), lty = 4, col = 4 );
		abline ( h = (2 * digitWidth + 0.5), lty = 4, col = 4 );
	} # for ( iCol in iColList  ) {

		#	Plots 2,4,6 - Distal, mid, proximal, respectively across digits 1-3
	iRowList = seq(1,3,1) * digitWidth - (floor(digitWidth/2));
	for ( iRow in iRowList  ) {
		plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
		kColorID = 0;
		for ( iCol in 1:N ) {
			iTmp = GetLin ( iRow, iCol, N );
			points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch="*", col=(1+(as.integer(kColorID%%kMaxColorID))) );
			if ( showEllipsesFlag ) {
				lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale=c(rfTrackData$rfCovar[iTmp,1],rfTrackData$rfCovar[iTmp,2]),
					centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
					col=(1+(as.integer(kColorID%%kMaxColorID))));
				kColorID = kColorID + 1;
			} # if ( showEllipsesFlag ) {
		} # for ( iCol in 1:N ) {
	} # for ( iRow in iRowList  ) {

} # ShowThreeDigitRFTrack0= function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit  ) {


ShowThreeDigitRFTrack1= function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit, localPlotControl, proxDistOnly   ) {


	kMaxColorID = 5;
	N = as.integer ( sqrt ( dim(rfTrackData$rfCentroid)[1] ) );
	digitWidth = N / 3;

	if ( localPlotControl ) {
		x11(); par(mfcol=c(3,2));
	}

		#	Plots 1,3,5 - Digit 3,2,1, respectively, distal to proximal along centerline.

	iColList = c( 3 * digitWidth - (floor(digitWidth/2)), digitWidth );
	
	plot ( 0, 0, type="n", xaxt="n", yaxt="n", xlim=c(1,N), ylim=c(1,N), xlab="", ylab="" );
	mtext(side=1, "Distal -> Proximal", cex=0.75, line=0);
	mtext(side=2, "D1 -> D3", cex=0.75, line=0);
	mtext(side=3, titleBaseText, cex=1.0, line=1, font=2);
	
	for ( iCol in iColList  ) {
		#plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
		kColorID = 0;
		for ( iRow in seq(3, N-2, 3) ) {
			iTmp = GetLin ( iRow, iCol, N );
			points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
			if ( showEllipsesFlag ) {
				lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale = c(rfTrackData$rfCovar[iTmp,1], rfTrackData$rfCovar[iTmp,2]),
					centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
					col=(1+(as.integer(kColorID%%kMaxColorID))) );
				kColorID = kColorID + 1;
			} # if ( showEllipsesFlag ) {
		} # for ( irow = 1:N )

			#	For ease of interpretation of the figure, draw dotted lines for digit boundaries.
		#abline ( h = (digitWidth + 0.5), lty = 4, col = 4 );
		#abline ( h = (2 * digitWidth + 0.5), lty = 4, col = 4 );

		abline ( h = (digitWidth + 0.5), lty = 3, col = 1, lwd=0.25 );
		abline ( h = (2 * digitWidth + 0.5), lty = 3 , col = 1, lwd=0.25 );
	} # for ( iCol in iColList  ) {

	if ( !proxDistOnly ) {
			#	Plots 2,4,6 - Distal, mid, proximal, respectively across digits 1-3
		iRowList = seq(1,3,1) * digitWidth - (floor(digitWidth/2));
		for ( iRow in iRowList  ) {
			plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
			kColorID = 0;
			for ( iCol in 1:N ) {
				iTmp = GetLin ( iRow, iCol, N );
				points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
				if ( showEllipsesFlag ) {
					lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale=c(rfTrackData$rfCovar[iTmp,1],rfTrackData$rfCovar[iTmp,2]),
						centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
						col=(1+(as.integer(kColorID%%kMaxColorID))));
					kColorID = kColorID + 1;
				} # if ( showEllipsesFlag ) {
			} # for ( iCol in 1:N ) {
		} # for ( iRow in iRowList  ) {
	} # if ( !proxDistOnly )

} # ShowThreeDigitRFTrack1= function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit  ) {

ShowThreeDigitRFTrack1A = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit, localPlotControl, proxDistOnly   ) {


	kMaxColorID = 5;
	N = as.integer ( sqrt ( dim(rfTrackData$rfCentroid)[1] ) );
	digitWidth = N / 3;

	if ( localPlotControl ) {
		x11(); par(mfcol=c(3,2));
	}

		#	Plots 1,3,5 - Digit 3,2,1, respectively, distal to proximal along centerline.

	iColList = c( 2 * digitWidth, digitWidth+1 );
	plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
	for ( iCol in iColList  ) {
		#plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
		kColorID = 0;
		for ( iRow in 1:N ) {
			iTmp = GetLin ( iRow, iCol, N );
			points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
			if ( showEllipsesFlag ) {
				lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale = c(rfTrackData$rfCovar[iTmp,1], rfTrackData$rfCovar[iTmp,2]),
					centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
					col=(1+(as.integer(kColorID%%kMaxColorID))) );
				kColorID = kColorID + 1;
			} # if ( showEllipsesFlag ) {
		} # for ( irow = 1:N )

			#	For ease of interpretation of the figure, draw dotted lines for digit boundaries.
		#abline ( h = (digitWidth + 0.5), lty = 4, col = 4 );
		#abline ( h = (2 * digitWidth + 0.5), lty = 4, col = 4 );

		abline ( h = (digitWidth + 0.5), lty = 3, col = 1, lwd=0.25 );
		abline ( h = (2 * digitWidth + 0.5), lty = 3 , col = 1, lwd=0.25 );
	} # for ( iCol in iColList  ) {

	if ( !proxDistOnly ) {
			#	Plots 2,4,6 - Distal, mid, proximal, respectively across digits 1-3
		iRowList = seq(1,3,1) * digitWidth - (floor(digitWidth/2));
		for ( iRow in iRowList  ) {
			plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
			kColorID = 0;
			for ( iCol in 1:N ) {
				iTmp = GetLin ( iRow, iCol, N );
				points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
				if ( showEllipsesFlag ) {
					lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale=c(rfTrackData$rfCovar[iTmp,1],rfTrackData$rfCovar[iTmp,2]),
						centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
						col=(1+(as.integer(kColorID%%kMaxColorID))));
					kColorID = kColorID + 1;
				} # if ( showEllipsesFlag ) {
			} # for ( iCol in 1:N ) {
		} # for ( iRow in iRowList  ) {
	} # if ( !proxDistOnly )

} # ShowThreeDigitRFTrack1A= function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit  ) {

ShowThreeDigitRFTrack2 = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit, localPlotControl, proxDistOnly   ) {


	kMaxColorID = 5;
	N = as.integer ( sqrt ( dim(rfTrackData$rfCentroid)[1] ) );
	digitWidth = N / 3;

	if ( localPlotControl ) {
		x11(); par(mfcol=c(3,2));
	}

		#	Plots 1,3,5 - Digit 3,2,1, respectively, distal to proximal along centerline.

	#iColList = c( 3 * digitWidth - (floor(digitWidth/2)), digitWidth );
	iColList = c( 2 * digitWidth + 1, (floor(digitWidth/2)));
	
	plot ( 0, 0, type="n", xaxt="n", yaxt="n", xlim=c(1,N), ylim=c(1,N), xlab="", ylab="" );
	mtext(side=1, "Distal -> Proximal", cex=0.75, line=0);
	mtext(side=2, "D1 -> D3", cex=0.75, line=0);
	mtext(side=3, titleBaseText, cex=1.0, line=1, font=2);
		
	for ( iCol in iColList  ) {
		#plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
		kColorID = 0;
		for ( iRow in 1:N ) {
			iTmp = GetLin ( iRow, iCol, N );
			points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
			if ( showEllipsesFlag ) {
				lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale = c(rfTrackData$rfCovar[iTmp,1], rfTrackData$rfCovar[iTmp,2]),
					centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
					col=(1+(as.integer(kColorID%%kMaxColorID))) );
				kColorID = kColorID + 1;
			} # if ( showEllipsesFlag ) {
		} # for ( irow = 1:N )

			#	For ease of interpretation of the figure, draw dotted lines for digit boundaries.
		#abline ( h = (digitWidth + 0.5), lty = 4, col = 4 );
		#abline ( h = (2 * digitWidth + 0.5), lty = 4, col = 4 );

		abline ( h = (digitWidth + 0.5), lty = 3, col = 1, lwd=0.25 );
		abline ( h = (2 * digitWidth + 0.5), lty = 3 , col = 1, lwd=0.25 );
	} # for ( iCol in iColList  ) {

	if ( !proxDistOnly ) {
			#	Plots 2,4,6 - Distal, mid, proximal, respectively across digits 1-3
		iRowList = seq(1,3,1) * digitWidth - (floor(digitWidth/2));
		for ( iRow in iRowList  ) {
			plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
			kColorID = 0;
			for ( iCol in 1:N ) {
				iTmp = GetLin ( iRow, iCol, N );
				points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
				if ( showEllipsesFlag ) {
					lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale=c(rfTrackData$rfCovar[iTmp,1],rfTrackData$rfCovar[iTmp,2]),
						centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
						col=(1+(as.integer(kColorID%%kMaxColorID))));
					kColorID = kColorID + 1;
				} # if ( showEllipsesFlag ) {
			} # for ( iCol in 1:N ) {
		} # for ( iRow in iRowList  ) {
	} # if ( !proxDistOnly )

} # ShowThreeDigitRFTrack2= function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit  ) {


ShowThreeDigitRFTrack3 = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit, localPlotControl, iStep=1 ) {

	kMaxColorID = 5;
	N = as.integer ( sqrt ( dim(rfTrackData$rfCentroid)[1] ) );
	digitWidth = N / 3;

	if ( localPlotControl ) {
		x11(); par(mfcol=c(3,2));
	}

				#	Plots 2,4,6 - Distal, mid, proximal, respectively across digits 1-3
		iRowList = seq(1,3,1) * digitWidth - (floor(digitWidth/2));
		plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
		abline ( h = (digitWidth + 0.5), lty = 3, col = 1, lwd=0.25 );
		abline ( h = (2 * digitWidth + 0.5), lty = 3 , col = 1, lwd=0.25 );
		for ( iRow in iRowList  ) {
			kColorID = 0;
			for ( iCol in seq(iStep,N,iStep) ) {
				iTmp = GetLin ( iRow, iCol, N );
				points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
				if ( showEllipsesFlag ) {
					lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale=c(rfTrackData$rfCovar[iTmp,1],rfTrackData$rfCovar[iTmp,2]),
						centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
						col=(1+(as.integer(kColorID%%kMaxColorID))));
					kColorID = kColorID + 1;
				} # if ( showEllipsesFlag ) {
			} # for ( iCol in 1:N ) {
		} # for ( iRow in iRowList  ) {

} # ShowThreeDigitRFTrack3 = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit  ) {

ShowThreeDigitRFTrack3A = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit, localPlotControl, iStep=1 ) {

	kMaxColorID = 5;
	N = as.integer ( sqrt ( dim(rfTrackData$rfCentroid)[1] ) );
	digitWidth = N / 3;

	if ( localPlotControl ) {
		x11(); par(mfcol=c(3,2));
	}
	iColList = sort ( c ( 10, 14, 17, 23, 29, 32, 36 ) );
				#	Plots 2,4,6 - Distal, mid, proximal, respectively across digits 1-3
		iRowList = seq(1,3,1) * digitWidth - (floor(digitWidth/2));		
	
		plot ( 0, 0, type="n", xaxt="n", yaxt="n", xlim=c(1,N), ylim=c(1,N), xlab="", ylab="" );
		mtext(side=1, "Distal -> Proximal", cex=0.75, line=0);
		mtext(side=2, "D1 -> D3", cex=0.75, line=0);
		mtext(side=3, titleBaseText, cex=1.0, line=1, font=2);
		
		abline ( h = (digitWidth + 0.5), lty = 3, col = 1, lwd=0.25 );
		abline ( h = (2 * digitWidth + 0.5), lty = 3 , col = 1, lwd=0.25 );
		for ( iRow in iRowList ) {
			kColorID = 0;
			for ( iCol in iColList ) {
				iTmp = GetLin ( iRow, iCol, N );
				points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
				if ( showEllipsesFlag ) {
					lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale=c(rfTrackData$rfCovar[iTmp,1],rfTrackData$rfCovar[iTmp,2]),
						centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
						col=(1+(as.integer(kColorID%%kMaxColorID))));
					kColorID = kColorID + 1;
				} # if ( showEllipsesFlag ) {
			} # for ( iCol in 1:N ) {
		} # for ( iRow in iRowList  ) {

} # ShowThreeDigitRFTrack3A = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit  ) {

	#
	#
	#
ShowThreeDigitRFTrackFlex = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit,
							localPlotControl, proxDistOnly, iRowList, iColList, iStep=1 ) {

	kMaxColorID = 5;
	N = as.integer ( sqrt ( dim(rfTrackData$rfCentroid)[1] ) );
	digitWidth = N / 3;

	if ( localPlotControl ) {
		x11(); par(mfcol=c(3,2));
	}

		#	Plots 1,3,5 - Digit 3,2,1, respectively, distal to proximal along centerline.
	plot ( 0, 0, type="p", pch=" ", xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
	for ( iCol in iColList  ) {
		#plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
		kColorID = 0;
		for ( iRow in iRowList ) {
			iTmp = GetLin ( iRow, iCol, N );
			points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
			if ( showEllipsesFlag ) {
				lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale = c(rfTrackData$rfCovar[iTmp,1], rfTrackData$rfCovar[iTmp,2]),
					centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
					col=(as.integer(kColorID%%kMaxColorID)+1) );
				kColorID = kColorID + 1;
			} # if ( showEllipsesFlag ) {
		} # for ( irow = 1:N )

			#	For ease of interpretation of the figure, draw dotted lines for digit boundaries.
		#abline ( h = (digitWidth + 0.5), lty = 4, col = 4 );
		#abline ( h = (2 * digitWidth + 0.5), lty = 4, col = 4 );

		abline ( h = (digitWidth + 0.5), lty = 3, col = 1, lwd=0.25 );
		abline ( h = (2 * digitWidth + 0.5), lty = 3 , col = 1, lwd=0.25 );
	} # for ( iCol in iColList  ) {

	if ( !proxDistOnly ) {
			#	Plots 2,4,6 - Distal, mid, proximal, respectively across digits 1-3
		iRowList = seq(1,3,1) * digitWidth - (floor(digitWidth/2));
		for ( iRow in iRowList  ) {
			plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
			kColorID = 0;
			for ( iCol in 1:N ) {
				iTmp = GetLin ( iRow, iCol, N );
				points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
				if ( showEllipsesFlag ) {
					lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale=c(rfTrackData$rfCovar[iTmp,1],rfTrackData$rfCovar[iTmp,2]),
						centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
						col=(1+(as.integer(kColorID%%kMaxColorID))));
					kColorID = kColorID + 1;
				} # if ( showEllipsesFlag ) {
			} # for ( iCol in 1:N ) {
		} # for ( iRow in iRowList  ) {
	} # if ( !proxDistOnly )

} # ShowThreeDigitRFTrackFlex= function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit, ...  ) {

	#
	#	Difference from ShowThreeDigitRFTrackFlex is that the output is limited to a subset of the entire map.
	#	This is needed when generating plots for publication.
	#
ShowThreeDigitRFTrackFlexW = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit,
							localPlotControl, proxDistOnly, iRowList, iColList, iStep=1 ) {

	kMaxColorID = 5;
	N = as.integer ( sqrt ( dim(rfTrackData$rfCentroid)[1] ) );
	digitWidth = N / 3;

	if ( localPlotControl ) {
		x11(); par(mfcol=c(3,2));
	}

	if ( N == 45 ) {
		maxX = 32;
		minX = 10;
	} else {
		maxX = 55;
		minX = 20;
	} # if ( N == 45)

		#	Plots 1,3,5 - Digit 3,2,1, respectively, distal to proximal along centerline.

	#iColList = c( 3 * digitWidth - (floor(digitWidth/2)), digitWidth );
	#iColList = c( 2 * digitWidth + 1, (floor(digitWidth/2)));
	#plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
	plot ( minX, 0, xlim=c(minX,maxX ), ylim=c(1,N), type="p", pch=" ",
			xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
	
	for ( iCol in iColList  ) {
		kColorID = 0;
		for ( iRow in iRowList ) {
			iTmp = GetLin ( iRow, iCol, N );
			points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
			if ( showEllipsesFlag ) {
				lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale = c(rfTrackData$rfCovar[iTmp,1], rfTrackData$rfCovar[iTmp,2]),
					centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
					col=(1+(as.integer(kColorID%%kMaxColorID))) );
				kColorID = kColorID + 1;
			} # if ( showEllipsesFlag ) {
		} # for ( irow = 1:N )
			#	For ease of interpretation of the figure, draw dotted lines for digit boundaries.
		abline ( h = (digitWidth + 0.5), lty = 3, col = 1, lwd=0.25 );
		abline ( h = (2 * digitWidth + 0.5), lty = 3 , col = 1, lwd=0.25 );
	} # for ( iCol in iColList  ) {

	if ( !proxDistOnly ) {
			#	Plots 2,4,6 - Distal, mid, proximal, respectively across digits 1-3
		iRowList = seq(1,3,1) * digitWidth - (floor(digitWidth/2));
		for ( iRow in iRowList  ) {
			plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
			kColorID = 0;
			for ( iCol in 1:N ) {
				iTmp = GetLin ( iRow, iCol, N );
				points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
				if ( showEllipsesFlag ) {
					lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale=c(rfTrackData$rfCovar[iTmp,1],rfTrackData$rfCovar[iTmp,2]),
						centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
						col=(1+(as.integer(kColorID%%kMaxColorID))));
					kColorID = kColorID + 1;
				} # if ( showEllipsesFlag ) {
			} # for ( iCol in 1:N ) {
		} # for ( iRow in iRowList  ) {
	} # if ( !proxDistOnly )

} # ShowThreeDigitRFTrackFlexW = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit, ...  ) {

######################################################################################################
######################################################################################################
#
#	ShowTopoMap: Display RFs along several set tracks on the output layer.
#		Input: List produced by ExtractRFPositionData.
#
#		Output: One plot showing centroids.
#
#			This layout should give the right correspondence with other image outputs and visual ease.
#
######################################################################################################
######################################################################################################

ShowTopoMap = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit ) {

	N = as.integer ( sqrt ( dim(rfTrackData$rfCentroid)[1] ) );
	digitWidth = N / 3;
	kMaxColorID = 4;

	x11();
	kColorID = 0;
	plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
	for ( iCol in seq(N,1,-1) ) {
		for ( iRow in 1:N ) {
			iTmp = GetLin ( iRow, iCol, N );
			points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch="*", col=(1+(as.integer(kColorID%%kMaxColorID))) );
			if ( showEllipsesFlag ) {
				lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale=c(rfTrackData$rfCovar[iTmp,1],rfTrackData$rfCovar[iTmp,2]),
					centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
					col=(1+(as.integer(kColorID%%kMaxColorID))) );
				kColorID = kColorID + 1;
			} # if ( showEllipsesFlag ) {
		} # for ( iRow in 1:N )
	} # for ( iCol in seq(N,1,-1) ) {

		#	For ease of interpretation of the figure, draw dotted lines for digit boundaries.
	abline ( h = (digitWidth + 0.5), lty = 4, col = 4 );
	abline ( h = (2 * digitWidth + 0.5), lty = 4, col = 4 );	

} # ShowTopoMap = function ( rfData, kRFPeakToEdgeDetect )

	#	Variation for different/better plotting control.  Keep the previous for backwards compatibility.
ShowTopoMap1 = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit, localPlotControl ) {

	N = as.integer ( sqrt ( dim(rfTrackData$rfCentroid)[1] ) );
	digitWidth = N / 3;
	kMaxColorID = 4;

	if ( localPlotControl ) {
		x11();
	} # if ( localPlotControl )
	kColorID = 0;
	plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), type="n", xaxt="n", yaxt="n", xlab="", ylab="" );
	mtext(side=1, "Distal -> Proximal", cex=0.75, line=0);
	mtext(side=2, "D1 -> D3", cex=0.75, line=0);
	mtext(side=3, titleBaseText, cex=1.0, line=1, font=2);

	for ( iCol in seq(N,1,-1) ) {
		for ( iRow in 1:N ) {
			iTmp = GetLin ( iRow, iCol, N );
			points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=1, col=(1+(as.integer(kColorID%%kMaxColorID))) );
			if ( showEllipsesFlag ) {
				lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale=c(rfTrackData$rfCovar[iTmp,1],rfTrackData$rfCovar[iTmp,2]),
					centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
					col=(1+(as.integer(kColorID%%kMaxColorID))) );
				kColorID = kColorID + 1;
			} # if ( showEllipsesFlag ) {
		} # for ( iRow in 1:N )
	} # for ( iCol in seq(N,1,-1) ) {

		#	For ease of interpretation of the figure, draw dotted lines for digit boundaries.
	abline ( h = (digitWidth + 0.5), lty = 3, col = 1, lwd=0.25 );
	abline ( h = (2 * digitWidth + 0.5), lty = 3 , col = 1, lwd=0.25 );	

} # ShowTopoMap1 = function ( rfData, kRFPeakToEdgeDetect )


	#	Difference from ShowTopoMap1 is to plot a "window" of the total map.
ShowTopoMap1W = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit, localPlotControl ) {

	N = as.integer ( sqrt ( dim(rfTrackData$rfCentroid)[1] ) );
	digitWidth = N / 3;
	kMaxColorID = 4;
	
	if ( N == 45 ) {
		xMaxVal = 35;
		xMinVal = 5;
		yMinVal = 0;
	} else {
		xMaxVal = 55;
		xMinVal = 20;
		yMinVal = 0;
	} # if ( N ==45 )
	if ( localPlotControl ) {
		x11();
	} # if ( localPlotControl )
	kColorID = 0;
	#plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
	plot ( xMinVal, yMinVal, xlim=c(xMinVal,xMaxVal), ylim=c(yMinVal,xMaxVal), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3",
			main=titleBaseText, type="p", pch=" " );
		#for ( iCol in seq(N,1,-1) ) {
	for ( iCol in seq(xMaxVal,yMinVal,-1) ) {
		#for ( iRow in 1:N ) {
		for ( iRow in xMinVal:xMaxVal ) {
			iTmp = GetLin ( iRow, iCol, N );
			points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=3, col=(1+(as.integer(kColorID%%kMaxColorID))) );
			if ( showEllipsesFlag ) {
				lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale=c(rfTrackData$rfCovar[iTmp,1],rfTrackData$rfCovar[iTmp,2]),
					centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
					col=(1+(as.integer(kColorID%%kMaxColorID))) );
				kColorID = kColorID + 1;
			} # if ( showEllipsesFlag ) {
		} # for ( iRow in 1:N )
	} # for ( iCol in seq(N,1,-1) ) {

		#	For ease of interpretation of the figure, draw dotted lines for digit boundaries.
	abline ( h = (digitWidth + 0.5), lty = 3, col = 1, lwd=0.25 );
	abline ( h = (2 * digitWidth + 0.5), lty = 3 , col = 1, lwd=0.25 );	

} # ShowTopoMap1W = function ( rfData, kRFPeakToEdgeDetect )


#	Variation for different/better plotting control.  Keep the previous for backwards compatibility.

ShowTopoMap2 = function ( rfTrackData, titleBaseText, showEllipsesFlag, topoMapConfLimit, localPlotControl ) {

	N = as.integer ( sqrt ( dim(rfTrackData$rfCentroid)[1] ) );
	digitWidth = N / 3;
	kMaxColorID = 4;

	if ( localPlotControl ) {
		x11();
	} # if ( localPlotControl )
	kColorID = 0;
	plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
	for ( iCol in seq(N,1,-1) ) {
		for ( iRow in 1:N ) {
			iTmp = GetLin ( iRow, iCol, N );
			points ( rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2], pch=".", cex=4, col=(1+(as.integer(kColorID%%kMaxColorID))) );
			if ( showEllipsesFlag ) {
				lines ( ellipse ( rfTrackData$rfCovar[iTmp,3], scale=c(rfTrackData$rfCovar[iTmp,1],rfTrackData$rfCovar[iTmp,2]),
					centre=c(rfTrackData$rfCentroid[iTmp,1], rfTrackData$rfCentroid[iTmp,2]), level=topoMapConfLimit ),
					col=(1+(as.integer(kColorID%%kMaxColorID))) );
				kColorID = kColorID + 1;
			} # if ( showEllipsesFlag ) {
		} # for ( iRow in 1:N )
	} # for ( iCol in seq(N,1,-1) ) {

		#	For ease of interpretation of the figure, draw dotted lines for digit boundaries.
	abline ( h = (digitWidth + 0.5), lty = 3, col = 1, lwd=0.25 );
	abline ( h = (2 * digitWidth + 0.5), lty = 3 , col = 1, lwd=0.25 );	

} # ShowTopoMap2 = function ( rfData, kRFPeakToEdgeDetect )


######################################################################################################
######################################################################################################
#
#	ShowTopoMapOverTimeForOneCell: Display RF for one cell across time.
#		Input: List produced by ExtractRFPositionData.
#
#		Output: One plot showing centroids.
#
#		This layout should give the right correspondence with other image outputs and visual ease.
#
######################################################################################################
######################################################################################################

ShowTopoMapOverTimeForOneCell = function ( rfTrackData, iCellID, iStart, iEnd, iStep, titleBaseText,
	showEllipsesFlag, topoMapConfLimit, plotAreaReadyFlag ) {

	N = as.integer ( sqrt ( dim(rfTrackData[[1]]$rfCentroid)[1] ) );
	kT = length ( rfTrackData );
	digitWidth = N / 3;
	kMaxColorID = 5;
	kColorID = 0;
	if ( !plotAreaReadyFlag ) {
		x11();
	} # if ( !plotAreaReadyFlag )

	plot ( 0, 0, xlim=c(1,N), ylim=c(1,N), xlab="Distal -> Proximal", ylab="Digit 1 -> Digit 3", main=titleBaseText );
	
	for ( iTime in seq(iStart, iEnd, iStep) ) {
		pchChar = "+"; lwdVal = 1; ltyVal = 2;
		if ( iTime == iStart ) {
			pchChar = "o"; ltyVal = 3;
		} else if( iTime == iEnd ) {
			pchChar = "*"; lwdVal = 3; ltyVal = 1;
		} # if ( iTeim == iStart )
		points ( rfTrackData[[iTime]]$rfCentroid[iCellID,1], rfTrackData[[iTime]]$rfCentroid[iCellID,2], pch=pchChar, col=(1+(as.integer(kColorID%%kMaxColorID))) );
		if ( showEllipsesFlag ) {
			lines ( ellipse ( rfTrackData[[iTime]]$rfCovar[iCellID,3], scale=c(rfTrackData[[iTime]]$rfCovar[iCellID,1], rfTrackData[[iTime]]$rfCovar[iCellID,2]),
				centre=c(rfTrackData[[iTime]]$rfCentroid[iCellID,1], rfTrackData[[iTime]]$rfCentroid[iCellID,2]), level=topoMapConfLimit ),
				col=(1+(as.integer(kColorID%%kMaxColorID))), lwd=lwdVal, lty=ltyVal );
			kColorID = kColorID + 1;
		} # if ( showEllipsesFlag ) {
	} # for ( iCol in seq(N,1,-1) ) {

} # ShowTopoMapOverTimeForOneCell = function ( rfTrackData, iCellID, ...

######################################################################################################
######################################################################################################
#
#	TopoQuant1: Return a scalar representing neighborhood preservation
#			between the "skin" input layer and an upper layer.
#
#			All the distances are 2D.
#	
#	Input:
#		N - this implies the size and positions to use for computing input layer
#		g0 - a standard deviation figure to use for defining neighborhood
#		rfTopo - this is the list returned by ExtractRFPositionData which includes
#				rfCentroid and rfCovar.
#
#	Output:
#		C - a scalar
#
######################################################################################################
######################################################################################################

TopoQuant1= function ( N, g0, rfTopo ) {

		#	Generate cell positions for the N x N input layer.

	N2 = N^2;
	inData = matrix ( 0, nrow=N2, ncol=2 );
	iTmp = 1;
	for ( iCell in 1:N2 ) {

		inData[iTmp, 1] = GetRow ( iCell, N );
		inData[iTmp, 2] = GetCol ( iCell, N );
		iTmp = iTmp + 1;

	} # for ( icol = 1:N )
	tmpF = dist ( inData, method="euclidean" );
	tmpF[(tmpF > g0 * sqrt(2.0))] = 0.

		#	Generate the pairwise distances for the input positions of RF centroids.
	tmpG = dist ( rfTopo$rfCentroid, method="euclidean" );
	tmpG[(tmpG > g0 * sqrt(2.0))] = 0.

	qVal = sum ( tmpF * tmpG ) / (N2 * (N2 - 1) );
	#qVal = sum ( tmpF * tmpG );

	return ( qVal );

} # TopoQuant1= function ( x ) { {

######################################################################################################
######################################################################################################
#
#	TopoQuant2: Return a scalar representing neighborhood preservation
#			between the "skin" input layer and an upper layer.
#
#			All the distances are 2D.
#	
#	Input:
#		N - this implies the size and positions to use for computing input layer
#		rfTopo - this is the list returned by ExtractRFPositionData which includes
#				rfCentroid and rfCovar.
#
#	Output:
#		C - a scalar representing the sum of squared difference between the RF centroid location
#			and the corresponding input layer location (were a perfect 1:1 mapping occur)
#
######################################################################################################
######################################################################################################

TopoQuant2= function ( N, rfTopo ) {

		#	Generate cell positions for the N x N input layer.
	N2 = N^2;
	inData = matrix ( 0, nrow=N2, ncol=2 );
	iTmp = 1;
	for ( iCell in 1:N2 ) {

		inData[iTmp, 1] = GetRow ( iCell, N );
		inData[iTmp, 2] = GetCol ( iCell, N );
		iTmp = iTmp + 1;

	} # for ( icol = 1:N )

		#	Calculate sum of squared differences

	qVal = sum( apply( (((inData - rfTopo$rfCentroid )^2)/rfTopo$rfCentroid[,1:2]), 2, sum ) ) / (N2 - 1);

	return ( qVal );

} # TopoQuant2 = function ( N, rfTopo ) {

######################################################################################################
######################################################################################################
#
#	ThreeDigitCorticalAmplificationPlots: assume three-digits with a distal, mid, proximal region each.
#	
#	Input:
#		Matrix (Boolean) produced by QuantCorticalAmp.
#
#	Output:
#		
#		Optional 3 x 3 plot showing the distal, mid, proximal amplification for each of three digits.
#		Read the cortical.Amp matrix by COLUMN.
#
######################################################################################################
######################################################################################################

ThreeDigitCorticalAmplification = function ( cortical.Amp, plotFlag, titleText ) {
	
	if ( plotFlag ) {
		x11(); par(mfrow=c(3,3));
	} # if ( plotFlag )

	N2 = dim(cortical.Amp)[1];
	N = as.integer ( sqrt ( N2 ) );
	numDigits = 3;
	numSectors = 3;

	digitWidth = as.integer ( N / numDigits );
	sectorWidth = N / numSectors;

	corticalAmp = rep ( 0, numDigits * numSectors );
	regionLocs = matrix ( FALSE, nrow=N2, ncol=(numDigits * 3) );
	outputMap = matrix ( 0, nrow=N, ncol=N );

	iCount = 1;
	for ( iDigit in seq ( numDigits, 1, -1 ) ) {

		iColStart = (iDigit - 1) * digitWidth + 1;
		iColEnd = iDigit * digitWidth;

		iRowStart = 1;
		for ( iSector in seq ( numSectors, 1, -1 ) ) {

			iRowEnd = iRowStart + sectorWidth - 1;
			for ( iCol in iColStart:iColEnd ) {
				regionLocs[ GetLin ( iRowStart:iRowEnd, iCol, N ), iCount ] = TRUE;
			} # for ( iCol in iColStart:iColEnd )

				#	Plot the result of logical sum.

			cAmp = apply ( cortical.Amp[, regionLocs[,iCount] ], 1, sum );
				
				#
				#	The following is a kludgey correction terms.  Depending upon the spread of afferents, after an
				#	an amputation there could be left cells that get zero external input.  During simulations, they
				#	just bubble along oblivous, but to the receptive field mapping algorithms, they have legit activity
				#	but they wind up "responding" to every single input cell.  Probably could and should take care of
				#	this in the rfMapping routines by applying an overall activity threshold of some sort.
				#

			cAmp[ cAmp == (sectorWidth^2) ] = 0;
			corticalAmp[iCount] = sum ( cAmp != 0 );
			if ( plotFlag ) {
				ShowVecAsMap ( cAmp, paste (titleText, "\n", "Digit", iDigit, "Sector", iSector, sep=" " ) );
			} # if ( plotFlag )

			iCount = iCount + 1;
			iRowStart = iRowStart + sectorWidth;

		} # for ( iSector in seq ( numSectors, 1, -1 ) ) {

	} # for ( iDigit in 1:numDigits )
	
	return ( corticalAmp );
	#return ( list ( corticalAmp = corticalAmp ) );
	
} # ThreeDigitCorticalAmplification = function ( cortical.Amp, plotFlag, titleText ) {



#	Define special-purpose functions.
GetRFMapData = function ( fileRootName, iBase, iStart, iEnd, kRFPeakToEdgeDetect, N2 ) {
	
	rfTrackData.e = list();
	rfTrackData.i = list();
	cortical.Mag.e = list();
	cortical.Mag.i = list();
	#quant.topoMap.e = rep ( 0, iEnd - iStart + 1 );
	#quant.topoMap.i = rep ( 0, iEnd - iStart + 1 );

	iCount = 1;
	iLabels = seq ( iStart, iEnd, 1 );
	for ( iRefinement in iLabels ) {

		fileName = paste ( fileRootName, ".", iBase, ".", iRefinement, ".bin", sep="" );
		finfo = file.info ( fileName );
		toread = file ( fileName, "rb" );
		alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
		close ( toread );	
	
		itmp = length ( alldata ) / 2;
		r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
		r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

		rfTrackData.e[[iCount]] = ExtractRFPositionData ( r1.e.rfMap, kRFPeakToEdgeDetect );
		rfTrackData.i[[iCount]] = ExtractRFPositionData ( r1.i.rfMap, kRFPeakToEdgeDetect );

		cortical.Mag.e[[iCount]] = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
		cortical.Mag.i[[iCount]] = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );

		iCount = iCount + 1;
	
	} # for ( iRefinement in iStart:iEnd ) {

	return ( list ( rfTrackData.e=rfTrackData.e, rfTrackData.i=rfTrackData.i,
				cortical.Mag.e=cortical.Mag.e, cortical.Mag.i=cortical.Mag.i ) );

} # GetRFMapData = function ( fileBaseName, iBase, iStart, iEnd, kRFPeakToEdgeDetect ) {

#	Define special-purpose functions.
#	This variation is a lazy way of handling a different file naming convention arising from some batch runs.
GetRFMapData1 = function ( fileRootName, iBase, iStart, iEnd, iSubZone, iFactor, kRFPeakToEdgeDetect, N2 ) {
	
	rfTrackData.e = list();
	rfTrackData.i = list();
	cortical.Mag.e = list();
	cortical.Mag.i = list();
	#quant.topoMap.e = rep ( 0, iEnd - iStart + 1 );
	#quant.topoMap.i = rep ( 0, iEnd - iStart + 1 );

	iCount = 1;
	iLabels = seq ( iStart, iEnd, 1 );
	for ( iRefinement in iLabels ) {

		fileName = paste ( fileRootName, iBase, iRefinement, iSubZone, iFactor, "bin", sep="." );
		finfo = file.info ( fileName );
		toread = file ( fileName, "rb" );
		alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
		close ( toread );	
	
		itmp = length ( alldata ) / 2;
		r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
		r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

		rfTrackData.e[[iCount]] = ExtractRFPositionData ( r1.e.rfMap, kRFPeakToEdgeDetect );
		rfTrackData.i[[iCount]] = ExtractRFPositionData ( r1.i.rfMap, kRFPeakToEdgeDetect );

		cortical.Mag.e[[iCount]] = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
		cortical.Mag.i[[iCount]] = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );

		iCount = iCount + 1;
	
	} # for ( iRefinement in iStart:iEnd ) {

	return ( list ( rfTrackData.e=rfTrackData.e, rfTrackData.i=rfTrackData.i,
				cortical.Mag.e=cortical.Mag.e, cortical.Mag.i=cortical.Mag.i ) );

} # GetRFMapData = function ( fileBaseName, iBase, iStart, iEnd, kRFPeakToEdgeDetect ) {

#	Define special-purpose functions.
#	This variation is a lazy way of handling a different file naming convention arising from some batch runs.
GetRFMapData2 = function ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 ) {
	
	rfTrackData.e = list();
	rfTrackData.i = list();
	cortical.Mag.e = list();
	cortical.Mag.i = list();
	#quant.topoMap.e = rep ( 0, iEnd - iStart + 1 );
	#quant.topoMap.i = rep ( 0, iEnd - iStart + 1 );

	iCount = 1;
	iLabels = seq ( iStart, iEnd, iStepSize );
	for ( iRefinement in iLabels ) {

		fileName = paste ( fileRootName, iBase, iRefinement, "bin", sep="." );
		finfo = file.info ( fileName );
		toread = file ( fileName, "rb" );
		alldata = readBin ( toread, "double", size=8, n=finfo$size/8, endian="little" );
		close ( toread );	
	
		itmp = length ( alldata ) / 2;
		r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
		r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

		rfTrackData.e[[iCount]] = ExtractRFPositionData ( r1.e.rfMap, kRFPeakToEdgeDetect );
		rfTrackData.i[[iCount]] = ExtractRFPositionData ( r1.i.rfMap, kRFPeakToEdgeDetect );

		cortical.Mag.e[[iCount]] = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
		cortical.Mag.i[[iCount]] = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );

		iCount = iCount + 1;
	
	} # for ( iRefinement in iStart:iEnd ) {

	return ( list ( rfTrackData.e=rfTrackData.e, rfTrackData.i=rfTrackData.i,
				cortical.Mag.e=cortical.Mag.e, cortical.Mag.i=cortical.Mag.i ) );

} # GetRFMapData2 = function ( fileBaseName, iBase, iStart, iEnd, kRFPeakToEdgeDetect ) {

#	Define special-purpose functions.
#	This variation is a lazy way of handling a different file naming convention arising from some batch runs.
GetRFMapData2A = function ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 ) {
	
	rfTrackData.e = list();
	rfTrackData.i = list();
	cortical.Mag.e = list();
	cortical.Mag.i = list();
	#quant.topoMap.e = rep ( 0, iEnd - iStart + 1 );
	#quant.topoMap.i = rep ( 0, iEnd - iStart + 1 );

	iCount = 1;
	iLabels = seq ( iStart, iEnd, iStepSize );
	for ( iRefinement in iLabels ) {

		fileName = paste ( fileRootName, iBase, iRefinement, "bin", sep="." );
		finfo = file.info ( fileName );
		toread = file ( fileName, "rb" );
		alldata = readBin ( toread, "double", size=8, n=finfo$size/8, endian="little" );
		close ( toread );	
	
		itmp = length ( alldata ) / 2;
		r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
		r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

		rfTrackData.e[[iCount]] = ExtractRFPositionData ( r1.e.rfMap, kRFPeakToEdgeDetect );
		rfTrackData.i[[iCount]] = ExtractRFPositionData ( r1.i.rfMap, kRFPeakToEdgeDetect );

		cortical.Mag.e[[iCount]] = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
		cortical.Mag.i[[iCount]] = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );

		iCount = iCount + 1;
	
	} # for ( iRefinement in iStart:iEnd ) {

	return ( list ( rfTrackData.e=rfTrackData.e, rfTrackData.i=rfTrackData.i,
				cortical.Mag.e=cortical.Mag.e, cortical.Mag.i=cortical.Mag.i, r1.e.rfMap = r1.e.rfMap, r1.i.rfMap = r1.i.rfMap ) );

} # GetRFMapData2A = function ( fileBaseName, iBase, iStart, iEnd, kRFPeakToEdgeDetect ) {

	#
	#	Difference from 2A is that in 2B the needs of "knock-in" experiments need to be handled.
	#	In "knock-in" it is observed that average cell response rates can be quite low: < 1e-10.
	#	How to handle when computing RF extent?  Under strict definition, we'd take all sensory
	#	layer points that drive response to greater than 50% peak as a legitimate RF.  But this
	#	makes no sense.  Better to set a threshold level of magnitude of response to sensor layer
	#	input; otherwise just set the values to zero.  This should propogate nicely to a variety
	#	of display and analytical routines.
	#	Impacted data structures: rfTrackData.x and r1.x.rfMap.
	#	As of 19-May-2016 focus on rfTrackData.x and r1.x.rtMap.  Return (later) for cortical.Mag.x.
	#
GetRFMapData2B = function ( fileRootName, iBase, iStart, iEnd, iStepSize, kRFPeakToEdgeDetect, N2 ) {

		#
		#	This is designed to be VERY permissive.  At the same time, it filters out
		#	what (turn out to be) nonsensical analysis of E(I) cells whose maximum response
		#	is <<1e-10.
		#
	minRespMag = 1.0;
	
	rfTrackData.e = list();
	rfTrackData.i = list();
	cortical.Mag.e = list();
	cortical.Mag.i = list();

	iCount = 1;
	iLabels = seq ( iStart, iEnd, iStepSize );
	for ( iRefinement in iLabels ) {

		fileName = paste ( fileRootName, iBase, iRefinement, "bin", sep="." );
		finfo = file.info ( fileName );
		toread = file ( fileName, "rb" );
		alldata = readBin ( toread, "double", size=8, n=finfo$size/8, endian="little" );
		close ( toread );	
	
		itmp = length ( alldata ) / 2;
		r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
		r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

		rfTrackData.e[[iCount]] = ExtractRFPositionData ( r1.e.rfMap, kRFPeakToEdgeDetect );
		rfTrackData.i[[iCount]] = ExtractRFPositionData ( r1.i.rfMap, kRFPeakToEdgeDetect );

		cortical.Mag.e[[iCount]] = QuantCorticalAmp ( r1.e.rfMap, kRFPeakToEdgeDetect );
		cortical.Mag.i[[iCount]] = QuantCorticalAmp ( r1.i.rfMap, kRFPeakToEdgeDetect );

			#
			#	Post process rf.x.rfMap
			#

		maxVals.e = apply ( r1.e.rfMap, 1, max );
		maxVals.i = apply ( r1.i.rfMap, 1, max );

		whichAreLow.e = which ( maxVals.e < minRespMag );
		whichAreLow.i = which ( maxVals.i < minRespMag );
		
		r1.e.rfMap[ whichAreLow.e ] = 0.0;
		r1.i.rfMap[ whichAreLow.i ] = 0.0;

		rfTrackData.e[[iCount]]$rfCentroid[whichAreLow.e,] = 0.0;
		rfTrackData.i[[iCount]]$rfCentroid[whichAreLow.i,] = 0.0;
		rfTrackData.e[[iCount]]$rfCovar[whichAreLow.e,] = 0.0;
		rfTrackData.i[[iCount]]$rfCovar[whichAreLow.i,] = 0.0;		

		iCount = iCount + 1;
	
	} # for ( iRefinement in iStart:iEnd ) {

	return ( list ( rfTrackData.e=rfTrackData.e, rfTrackData.i=rfTrackData.i,
				cortical.Mag.e=cortical.Mag.e, cortical.Mag.i=cortical.Mag.i, r1.e.rfMap = r1.e.rfMap, r1.i.rfMap = r1.i.rfMap ) );

} # GetRFMapData2B = function ( fileBaseName, iBase, iStart, iEnd, kRFPeakToEdgeDetect ) {

	#	Define special-purpose functions.
	#	This variation is a lazy way of handling a different file naming convention arising from some batch runs.
GetRFMapData3 = function ( fileRootName, iBase, iSequence, kRFPeakToEdgeDetect, N2 ) {
	
	rfTrackData.e = list();
	rfTrackData.i = list();
	cortical.Mag.e = list();
	cortical.Mag.i = list();
	r1.e.rfMap = list();
	r1.i.rfMap = list();
	#quant.topoMap.e = rep ( 0, iEnd - iStart + 1 );
	#quant.topoMap.i = rep ( 0, iEnd - iStart + 1 );

	iCount = 1;
	for ( iRefinement in iSequence ) {

		fileName = paste ( fileRootName, iBase, iRefinement, "bin", sep="." );
		finfo = file.info ( fileName );
		toread = file ( fileName, "rb" );
		alldata = readBin ( toread, "double", size=8, n=finfo$size/8, endian="little" );
		close ( toread );	
	
		itmp = length ( alldata ) / 2;
		r1.e.rfMap[[iCount]] = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
		r1.i.rfMap[[iCount]] = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

		rfTrackData.e[[iCount]] = ExtractRFPositionData ( r1.e.rfMap[[iCount]], kRFPeakToEdgeDetect );
		rfTrackData.i[[iCount]] = ExtractRFPositionData ( r1.i.rfMap[[iCount]], kRFPeakToEdgeDetect );

		cortical.Mag.e[[iCount]] = QuantCorticalAmp ( r1.e.rfMap[[iCount]], kRFPeakToEdgeDetect );
		cortical.Mag.i[[iCount]] = QuantCorticalAmp ( r1.i.rfMap[[iCount]], kRFPeakToEdgeDetect );

		iCount = iCount + 1;
	
	} # for ( iRefinement in iStart:iEnd ) {

	return ( list ( rfTrackData.e=rfTrackData.e, rfTrackData.i=rfTrackData.i,
				cortical.Mag.e=cortical.Mag.e, cortical.Mag.i=cortical.Mag.i,
				r1.e.rfMap=r1.e.rfMap, r1.i.rfMap=r1.i.rfMap ) );

} # GetRFMapData3 = function ( fileBaseName, iBase, iStart, iEnd, kRFPeakToEdgeDetect ) {

GetRFMapData3A = function ( fileRootName, iBase, iSequence, kRFPeakToEdgeDetect, N2 ) {
	
	rfTrackData.e = list();
	rfTrackData.i = list();
	cortical.Mag.e = list();
	cortical.Mag.i = list();
	#quant.topoMap.e = rep ( 0, iEnd - iStart + 1 );
	#quant.topoMap.i = rep ( 0, iEnd - iStart + 1 );

	iCount = 1;
	for ( iRefinement in iSequence ) {

		fileName = paste ( fileRootName, iBase, iRefinement, "bin", sep="." );
		finfo = file.info ( fileName );
		toread = file ( fileName, "rb" );
		alldata = readBin ( toread, "double", size=8, n=finfo$size/8, endian="little" );
		close ( toread );	
	
		itmp = length ( alldata ) / 2;
		r1.e.rfMap = matrix ( alldata[1:itmp], nrow=N2, ncol=N2, byrow=FALSE );
		r1.i.rfMap = matrix ( alldata[(itmp+1):length(alldata)], nrow=N2, ncol=N2, byrow=FALSE );

		rfTrackData.e[[iCount]] = ExtractRFPositionData ( r1.e.rfMap, 0.5 );
		rfTrackData.i[[iCount]] = ExtractRFPositionData ( r1.i.rfMap, 0.75 );

		cortical.Mag.e[[iCount]] = QuantCorticalAmp ( r1.e.rfMap, 0.5 );
		cortical.Mag.i[[iCount]] = QuantCorticalAmp ( r1.i.rfMap, 0.75 );

		iCount = iCount + 1;
	
	} # for ( iRefinement in iStart:iEnd ) {

	return ( list ( rfTrackData.e=rfTrackData.e, rfTrackData.i=rfTrackData.i,
				cortical.Mag.e=cortical.Mag.e, cortical.Mag.i=cortical.Mag.i ) );

} # GetRFMapData3A = function ( fileBaseName, iBase, iStart, iEnd, kRFPeakToEdgeDetect ) {



	#
	#	CPP program dumps the input patches into a checkfile.
	#		Each "entry" is the N2 * numIter input patch AND the N2 stimulus count.
	#		In this routine, we're interested mainly in the stimulus count.
	#		Return the stimulus count from the last stimulus patch.
	#		(This is owing to the way that CPP builds up the input patch list.)
	#

GetStimCountData = function ( fileRootName, N2 ) {

	finfo = file.info ( fileRootName );
	toread = file ( fileRootName, "rb" );
	alldata = readBin ( toread, "double", size=8, n=finfo$size, endian="little" );
	close ( toread );
	
	iTmp = length ( alldata );
	xTmp = alldata[((iTmp - N2 + 1):(iTmp))];
	
	return ( xTmp );

} # GetStimCountData = function ( fileRootName, N2 ) {

GetStimCountData1 = function ( fileRootName, N2 ) {

	finfo = file.info ( fileRootName );
	toread = file ( fileRootName, "rb" );

	tmp = finfo$size;
	numReads = ( tmp / 8 ) / N2;

	for ( i in 1:numReads ) {	
		alldata = readBin ( toread, "double", size=8, n=(N2), endian="little" );
	}
	close ( toread );

	return ( alldata );

} # GetStimCountData1 = function ( fileRootName, N2 ) {

######################################################################################################
######################################################################################################
#
#	IntraColumnarRFDivergence: Display RFs along several set tracks on the output layer.
#		Inputs: Lists produced by ExtractRFPositionData.
#
#		For each node, find the Euclidean distance between the E-cell and I-cell RF centroid.#
#
######################################################################################################
######################################################################################################

IntraColumnarRFDivergence = function ( rfTrackData.e, rfTrackData.i ) {

	return ( RFCentroidDelta ( rfTrackData.e, rfTrackData.i ) );

} # IntraColumnarRFDivergence = function ( rfData.e, rfData.i ) {

RFCentroidDelta = function ( rfExp, rfBase ) {

	N = as.integer ( sqrt ( dim(rfExp$rfCentroid)[1] ) );
	N2 = N * N;

	tmp = rep ( 0 , N2 );
	for ( i in 1:N2 ) {
		tmp [i] = sqrt ( (rfExp$rfCentroid[i,1] - rfBase$rfCentroid[i,1])^2 +
					(rfExp$rfCentroid[i,2] - rfBase$rfCentroid[i,2])^2 );
	} # for ( i in 1:N2 )

	return ( tmp );

} # RFCentroidDeltaMap = function ( rfData.e, rfData.i ) {


######################################################################################################
######################################################################################################
#
#	TopoQuantSeries : Utility function to manage topographic map quantitative functions.
#
######################################################################################################
######################################################################################################

TopoQuantSeries = function ( N, g0, rfTrackData ) {

	numCount = length ( rfTrackData );
	topo1 = rep ( 0, numCount );
	topo2 = rep ( 0, numCount );
	for ( iCount in 1:numCount) {

		topo1[iCount] = TopoQuant1 ( N, g0, rfTrackData[[iCount]] );	
		topo2[iCount] = TopoQuant2 ( N, rfTrackData[[iCount]] );

	} # for ( iRefinement in iStart:iEnd ) {

	return ( list ( q1 = topo1, q2 = topo2 ) );

} # TopoQuantSeries = function ( ... )

######################################################################################################
######################################################################################################
#
#	IntraColumnarRFDivergence: Display RFs along several set tracks on the output layer.
#		Inputs: Lists produced by ExtractRFPositionData.
#
#		For each node, find the Euclidean distance between the E-cell and I-cell RF centroid.#
#
######################################################################################################
######################################################################################################

IntraColumnarRFDivergence = function ( rfTrackData.e, rfTrackData.i ) {

	N = as.integer ( sqrt ( dim(rfTrackData.e$rfCentroid)[1] ) );
	N2 = N * N;

	div = rep ( 0 , N2 );
	for ( i in 1:N2 ) {
		div[i] = sqrt ( (rfTrackData.e$rfCentroid[i,1] - rfTrackData.i$rfCentroid[i,1])^2 +
					(rfTrackData.e$rfCentroid[i,2] - rfTrackData.i$rfCentroid[i,2])^2 );
	} # for ( i in 1:N2 )

	return ( div );

} # IntraColumnarRFDivergence = function ( rfData.e, rfData.i ) {


######################################################################################################
######################################################################################################
#
#	GetSparseWeightMatrix: Ingest a Sparse Weight Matrix produced by CPP.
#
#		Input: File name
#
#		Output: 	Value vector containing linear representation for all cells aggreggated.
#				Index vector pointing to the start of the information for the ith cell.
#
#
######################################################################################################
######################################################################################################

GetSparseWeightMatrix = function ( fName ) {

	finfo = file.info ( fName );
	toread = file ( fName, "rb" );
	alldata = readBin ( toread, "double", size=8, n=finfo$size/8, endian="little" );
	close ( toread );

	N2 = alldata[1] * alldata[1];
	index = rep ( 0, N2 );
	index[1] = 1;
	for ( i in 2:N2 ) {

		index[i] = index[i-1] + ( ( alldata[ index[i-1] + 2 ] * 2 ) + 4 )
		
	} # for ( i = 2:N2 )

	return ( list ( wghts=alldata, index=index ) );

} # GetSparseWeightMatrix = function ( fName ) {


######################################################################################################
######################################################################################################
#
#	GetInputWeights: Ingest a Sparse Weight Matrix produced by CPP.
#
#		Input: Weight Matrix produced by GetSparseWeightMatrix;
#			 cellID - in 1-based indexing.
#
#		Output: 	Value vector containing linear representation for all cells aggreggated.
#				Index vector pointing to the start of the information for the ith cell.
#
#
######################################################################################################
######################################################################################################

GetInputWeights = function ( w, cellID ) {

	N2 = (w$wghts[1])^2;
	wVec = rep ( 0, N2 );

	iPtrToRecordStart = w$index[cellID];

	iPtrToNumWghts = iPtrToRecordStart + 2;
	numWghts = w$wghts[iPtrToNumWghts];

	iPtrToPreSynapticCellList = iPtrToRecordStart + 4;
	iPtrToWeightValues = iPtrToPreSynapticCellList + numWghts;

		#	Subtract 1 from cellID, because CPP is 0-index based.

	if ( w$wghts[iPtrToRecordStart + 1] != ( cellID - 1 ) ) {
		wVec = rep ( 0, N2 );
	} else {

		for ( i in 0:(numWghts - 1) ) {

				#	Subtract 1 from cellID, because CPP is 0-index based.

			wVec[w$wghts[iPtrToPreSynapticCellList + i]+1] = w$wghts[iPtrToWeightValues + i];

		} # for ( i in 1:numWghts )
		
	} # if ( wghtMatrix$wghts[iStart+1] != cellID ) {

	return ( wVec );

} # GetInputWeights = function ( w, cellID ) {

######################################################################################################
######################################################################################################
#
#	GenLongAxisSummaryTable
#
######################################################################################################
######################################################################################################

GenD3LongAxisSummaryTable = function ( magMap, rfAreas, stimCountMap, trimEnds ) {

	sTable = matrix( 0, nrow=N, ncol=4 );
	k = 1;

	for ( iDigit in 1:3 ) {

		for ( iLongAxis in 1:(N/3) ) {

			tmp = TrimEdgesFromCellList ( MapMagToCellList ( magMap, D3LongAxisLocs1 ( N, iDigit, iLongAxis ) ), N, trimEnds );
			sTable[k,1] = max ( stimCountMap [ D3LongAxisLocs1 ( N, iDigit, iLongAxis ) ] );
			sTable[k,2] = length ( tmp );
			sTable[k,3] = round ( mean ( rfAreas [ tmp ] ), 2 );
			sTable[k,4] = round ( sqrt ( var ( rfAreas [ tmp ] ) ), 2 );
			k = k + 1;

		} # 	for ( iLongAxis in 1:(N/3) ) {

	} # for ( iDigit in 1:3 )

	return ( sTable );

} # GenD3LongAxisSummaryTable = function ( magMap, rfAreas, stimCountMap, trimEnds ) {

######################################################################################################
######################################################################################################
#
#	GenSectorCellList
#
######################################################################################################
######################################################################################################

GenD3SectorCellList = function ( N, zoneID, trimRing=0 ) {

	digitWidth = as.integer ( N / 3 );
	minLegalColumn = maxLegalColumn = minLegalRow = maxLegalRow = 0;

	if ( zoneID == 1 || zoneID == 4 || zoneID == 7 ) {
		iBlock = 1;
		minLegalColumn = ( iBlock - 1 ) * digitWidth + 1 + trimRing;
		maxLegalColumn = iBlock * digitWidth;
	} else if ( zoneID == 2 || zoneID == 5 || zoneID == 8 ) {
		iBlock = 2;
		minLegalColumn = ( iBlock - 1 ) * digitWidth + 1;
		maxLegalColumn = iBlock * digitWidth - trimRing;
	} else {
		iBlock = 3;
		minLegalColumn = ( iBlock - 1 ) * digitWidth + 1;
		maxLegalColumn = iBlock * digitWidth;
	} # if ( zoneID == 1 || zoneID == 4 || zoneID == 7 )

	if ( zoneID == 1 || zoneID == 2 || zoneID == 3 ) {
		iBlock = 1;
		minLegalRow = ( iBlock - 1 ) * digitWidth + 1 + trimRing;
		maxLegalRow = iBlock * digitWidth;
	} else if ( zoneID == 4 || zoneID == 5 || zoneID == 6 ) {
		iBlock = 2;
		minLegalRow = ( iBlock - 1 ) * digitWidth + 1;
		maxLegalRow = iBlock * digitWidth;
	} else {
		iBlock = 3;
		minLegalRow = ( iBlock - 1 ) * digitWidth + 1;
		maxLegalRow = iBlock * digitWidth - trimRing;
	} # if ( zoneID == 7 || zoneID == 8 || zoneID == 9 )

	iCount = (( maxLegalRow - minLegalRow ) + 1 ) * (( maxLegalColumn - minLegalColumn ) + 1 );

	iCellList = rep ( 0, iCount );
	k = 1;
	for ( iCol in minLegalColumn:maxLegalColumn ) {
		for ( iRow in minLegalRow:maxLegalRow ) {
			iCellList[k] = GetLin ( iRow, iCol, N );
			k = k + 1;
		} # for ( iRow in minLegalRow:maxLegalRow ) {
	} # 	for ( iCol in minLegalColumn:maxLegalColumn ) {

	return ( iCellList );

} # GenD3SectorCellList = function ( N, iSector ) {

######################################################################################################
######################################################################################################
#
#	Plot Receptive Field Overlap Results
#
#
#	Input:
#		ovLap - a matrix where each column is of length N-1 corresponding to the successive
#				pairwise proportion overlap of receptive fields
#
######################################################################################################
######################################################################################################

PlotRFOverlap = function ( ovLap, mainTxt, xlabTxt, ylabTxt ) {

	numPlots = dim ( ovLap )[2];

	ylimExtra = 10;
	ylim = c ( min ( ovLap ) - ylimExtra, max ( ovLap, 100 ) );
	plot ( ovLap[,1], ylim=ylim, type="l", main=mainTxt, xlab=xlabTxt, ylab=ylabTxt );

	for ( i in 2:numPlots ) {

		lines ( ovLap[,i] );

	} # for ( i in 2:numPlots )	

} # PlotRFOverlap = function ( ovLap, mainTxt, xlabTxt, ylabTxt ) {


######################################################################################################
######################################################################################################
#
#	Driver for Receptive Field Overlap Plots
#
#	Input:
#		rfMag - a N^2 x N^2 Boolean matrix; ith row of rfMag corresponds to the ith cortical unit
#		iLongAxisFlag - 	if TRUE generate overlap plots for each longitudinal axis in the net;
#					if FALSE generate overlap plots or each cross sectional axis
#
#	Output:
#		ovLap - a list of length M - 1 corresponding to the successive pairwise proportion
#				overlap of receptive fields; position i in ovLap is that value for the
#				pair (i, i+1) of nodes referenced in expList
#
######################################################################################################
######################################################################################################

RFOverlap = function ( rfMag, iLongAxisFlag ) {

	N2 = dim(rfMag)[1];
	N = sqrt ( N2 );

	ovLap = matrix ( 0, nrow = N-1, ncol = N );
	for ( i in 1:N ) {

		if ( iLongAxisFlag ) {
			expList = seq ( (i-1)*N + 1, i*N, 1);
		} else {
			expList = seq ( i, N2, N);
		} # if ( iLongAxisFlag )

		ovLap[,i] = MeasureRFOverlap ( rfMag, expList )

	} # for ( i in 1:N )	

	return ( ovLap );

} # RFOverlap = function ( rfMag, iLongAxis ) {

######################################################################################################
######################################################################################################
#
#	Receptive Field Overlap
#
#	Input:
#		rfMag - a N^2 x N^2 Boolean matrix; ith row of rfMag corresponds to the ith cortical unit
#		expList - a list of length M of row indices into rfMag, e.g., the list of cortical units
#
#	Output:
#		ovLap - a list of length M - 1 corresponding to the successive pairwise proportion
#				overlap of receptive fields; position i in ovLap is that value for the
#				pair (i, i+1) of nodes referenced in expList
#
######################################################################################################
######################################################################################################

MeasureRFOverlap = function ( rfMag, expList ) {

	ovLap = rep ( 0, length(expList) - 1 );

	for ( i in 2:length(expList) ) {

		ovLap[i-1] = 100.0 * sum ( rfMag[expList[i-1],] & rfMag [expList[i],] ) / sum ( rfMag[expList[i-1],] | rfMag [expList[i],] );

	} # for ( i in 2:length(expList) ) {

	return ( ovLap );

} # MeasureRFOverlap = function ( rfMag, expList ) {


######################################################################################################
######################################################################################################
#
#	InterleaveForDisplay
#
#	Input:
#		A set of four NxN vectors of data (each in linear form).
#
#	Output:
#		Interleaved onto a 2Nx2N grid for use in animation.
#
######################################################################################################
######################################################################################################

InterleaveForDisplay = function ( w1, w2, w3, w4 ) {

	N = as.integer ( sqrt ( length ( w1 ) ) );
	wLen = 4 * N * N;
	w = rep ( 0, wLen );

	iLower = 1;
	iUpper = as.integer ( wLen/2 + 1 );

	for ( i in 1:N ) {

		w[iLower:(iLower+N-1)] = w1[((i-1)*N+1):(i*N)];
		iLower = iLower + N;

		w[iLower:(iLower+N-1)] = w3[((i-1)*N+1):(i*N)];
		iLower = iLower + N;

		w[iUpper:(iUpper+N-1)] = w2[((i-1)*N+1):(i*N)];
		iUpper = iUpper + N;

		w[iUpper:(iUpper+N-1)] = w4[((i-1)*N+1):(i*N)];
		iUpper = iUpper + N;

	} # for ( i in (2 * N) ) {

	return ( w );
	
} # InterleaveForDisplay = function ( w1, w2, w3, w4 ) {

######################################################################################################
######################################################################################################
#
#	RF Probe Responses Check Knockout vs Baseline Helper Functions
#
#
######################################################################################################
######################################################################################################

ExtractRFProbeRawData = function ( fName, N2, numValsPerRFTrial, iCellList ) {

	v1e = matrix ( 0, nrow=numValsPerRFTrial, ncol=length(iCellList) );
	v1i = matrix ( 0, nrow=numValsPerRFTrial, ncol=length(iCellList) );
	v0 = matrix ( 0, nrow=numValsPerRFTrial, ncol=length(iCellList) );
	iCount = 1;

	finfo = file.info ( fName ); toread = file ( fName, "rb" );
	for ( iProbeCell in 1:max(iCellList) ) {
		alldata = readBin ( toread, "double", size=8, n=(3*numValsPerRFTrial), endian="little" );
		if ( sum( iCellList == iProbeCell ) > 0 ) {
			startOffset.e = 1; startOffset.i = numValsPerRFTrial + 1;
			startOffset.0 = 2 * numValsPerRFTrial + 1;
			v1e[,iCount] = alldata [ (startOffset.e):(startOffset.e + numValsPerRFTrial - 1) ];
			v1i[,iCount] = alldata [ (startOffset.i):(startOffset.i + numValsPerRFTrial - 1) ];
			v0[,iCount] = alldata [ (startOffset.0):(startOffset.0 + numValsPerRFTrial - 1) ];
			iCount = iCount + 1;
		} # if ( sum( iCellList == iProbeCell ) > 0 ) {
	} # for ( iProbeCell in 1:max(iCellList) ) {
	close ( toread );

	return ( list ( v1e=v1e, v1i=v1i, v0=v0, rfProbeList=iCellList ) );

} # ExtractRFProbeRawData = function ( fName, N2, numValsPerRFTrial ) {

ExtractRFProbeRawDataSingleShot = function ( fName, N2, numValsPerRFTrial, iProbeCellID ) {

	#fName = fName.base;
	v1e = rep ( 0, numValsPerRFTrial  );
	v1i = rep ( 0, numValsPerRFTrial  );
	v0 = rep ( 0, numValsPerRFTrial  );
	iCount = 1;

	finfo = file.info ( fName ); toread = file ( fName, "rb" );

	alldata = readBin ( toread, "double", size=8, n=(3*numValsPerRFTrial), endian="little" );
	startOffset.e = 1; startOffset.i = numValsPerRFTrial + 1;
	startOffset.0 = 2 * numValsPerRFTrial + 1;
	v1e = alldata [ (startOffset.e):(startOffset.e + numValsPerRFTrial - 1) ];
	v1i = alldata [ (startOffset.i):(startOffset.i + numValsPerRFTrial - 1) ];
	v0 = alldata [ (startOffset.0):(startOffset.0 + numValsPerRFTrial - 1) ];
	close ( toread );

	return ( list ( v1e=v1e, v1i=v1i, v0=v0, rfProbeList=iProbeCellID ) );

} # ExtractRFProbeRawDataSingleShot = function ( fName, N2, numValsPerRFTrial ) {

	#
	#	Three-panel time series plot: r1E, r1I, r0.
	#
QuickCheckTimeSeries = function ( r1E, r1I, rZ, iCell, deltaT, plotFlag, ylim ) {

	tVal = 1:(length(r1E)) * deltaT;

	if ( plotFlag ) {
		x11();
		par(mfrow=c(3,1));
	} # if ( plotFlag )		

	plot( tVal, r1E, type="p", col=1, ylim=ylim, ylab="Firing Rate", xlab="Time (sec)");
	title(main=paste("Cortical Layer E-Cell Firing Rate",paste("Cell ",iCell),sep='\n'));

	plot( tVal, r1I, type="p", col=1, ylim=ylim, ylab="Firing Rate", xlab="Time (sec)");
	title(main=paste("Cortical Layer I-Cell Firing Rate",paste("Cell ",iCell),sep='\n'));

	plot( tVal, rZ, type="p", col=1, ylim=ylim, ylab="Firing Rate", xlab="Time (sec)");
	title(main=paste("Input Layer Firing Rate",paste("Cell ",iCell),sep='\n'));

} # QuickCheckTimeSeries = function ( r1E, r1I, r0, cellID, deltaT, ... ) {

QuickCheckTimeSeries2 = function ( tRange, r1E, r1I, rZ, iCell, deltaT, plotFlag, titleText, ylabText ) {

	tVal = tRange * deltaT;

	if ( plotFlag ) {
		x11();
		par(mfrow=c(3,1));
	} # if ( plotFlag )		

	plot( tVal, r1E, type="p", col=1, ylab=ylabText, xlab="Time (sec)");
	title ( main = paste("Cortical Layer E-Cell", titleText, paste("Cell ",iCell),sep='\n'));

	plot( tVal, r1I, type="p", col=1, ylab=ylabText, xlab="Time (sec)");
	title(main=paste("Cortical Layer I-Cell", titleText, paste("Cell ",iCell),sep='\n'));

	plot( tVal, rZ, type="p", col=1, ylab=ylabText, xlab="Time (sec)");
	title(main=paste("Input Layer Cell", titleText, paste("Cell ",iCell),sep='\n'));

} # QuickCheckTimeSeries2 = function ( r1E, r1I, r0, cellID, deltaT, ... ) {

GoodToShow = function ( N, iCandidate ) {

	candidateList = TrimEdgesFromCellList ( seq ( 1, N2, 1 ), N, 3 );
	show = FALSE;
	if ( sum ( candidateList == iCandidate ) ) { show = TRUE; }

	#if ( iCandidate <= N ) { show = FALSE; }
	#if ( iCandidate > (N^2 - N) ) { show = FALSE; }
	#if ( iCandidate %% N == 0 ) { show = FALSE; }
	#if ( iCandidate %% N == 1 ) { show = FALSE; }

	return ( show );

} # GoodToShow = function ( N, iCandidate ) {

	#
	#	Extract the vt values for the given iCell from the CPP generated RF Raw data file.
	#	Compute the rt values.
	#
TSResponse = function ( N2, numValsPerRFTrial, iCell, xt, eScale ) {

	iCellList = seq ( iCell, numValsPerRFTrial, N2 );
	startOffset.e = 1; startOffset.i = numValsPerRFTrial + 1;
	startOffset.0 = 2 * numValsPerRFTrial + 1;

	vt.e = ( xt [ (startOffset.e):(startOffset.e + numValsPerRFTrial - 1) ] ) / eScale;
	vt.e = vt.e[iCellList];

	vt.i = xt [ (startOffset.i):(startOffset.i + numValsPerRFTrial - 1) ];
	vt.i = vt.i[iCellList];

	vt.0 = xt [ (startOffset.0):(startOffset.0 + numValsPerRFTrial - 1) ];
	vt.0 = vt.0[iCellList];

	rt.e = sigmoid ( vt.e, 4 );
	rt.i = sigmoid ( vt.i, 4 );
	rt.0 = sigmoid ( vt.0, 4 );

	vt.min = min ( vt.e, vt.i, vt.0 );
	vt.max = max ( vt.e, vt.i, vt.0 );
	rt.min = min ( rt.e, rt.i, rt.0 );
	rt.max = max ( rt.e, rt.i, rt.0 );

	return ( list ( vt.e=vt.e, vt.i=vt.i, vt.0=vt.0, vt.min=vt.min, vt.max=vt.max,
			rt.e=rt.e, rt.i=rt.i, rt.0=rt.0, rt.min=rt.min, rt.max=rt.max ) );

} # TSResponse = function ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale ) {

	#
	#	Prepare a baseline vs experimental dataset for the given iCell
	#
ExpRefTSResponse = function ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale ) {

	xt.base = TSResponse ( N2, numValsPerRFTrial, iCell, alldata.base, eScale );
	xt.exp = TSResponse ( N2, numValsPerRFTrial, iCell, alldata.exp, eScale );
	
	v1e.delta = xt.exp$vt.e - xt.base$vt.e;
	v1i.delta = xt.exp$vt.i - xt.base$vt.i;
	v0.delta = xt.exp$vt.0 - xt.base$vt.0;
	
	vt.min = min ( xt.base$vt.min, xt.exp$vt.min ); vt.max = max ( xt.base$vt.max, xt.exp$vt.max );
	delta.min = min ( v1e.delta, v1i.delta, v0.delta ); delta.max = max ( v1e.delta, v1i.delta, v0.delta );

	return ( list ( v1e.base=xt.base$vt.e, v1i.base=xt.base$vt.i, v0.base=xt.base$vt.0,
				v1e.exp=xt.exp$vt.e, v1i.exp=xt.exp$vt.i, v0.exp=xt.exp$vt.0,
				v1e.delta=v1e.delta, v1i.delta=v1i.delta, v0.delta=v0.delta,
				vt.min=vt.min, vt.max=vt.max, delta.min=delta.min, delta.max=delta.max ) );

} # TSResponse = function ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale ) {

	#
	#
	#
KnockoutTimeSeriesPlotDriver = function ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleTextRoot, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fNameBase, iShortPlot ) {


	iRecordCellsList = c ( iProbeCell, iProbeCell + N, iProbeCell + 2*N, iProbeCell + 3*N );
	iRecordCellsList = c ( iProbeCell, iProbeCell + 2*N, iProbeCell + 4*N, iProbeCell + 6*N );
	fName = paste ( "D2Side", fNameBase, "tiff", sep="." );
	KnockoutTimeSeriesPlot ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleTextRoot, subTitleTextRoot, xlabTextRoot, ylabTextRoot, iRecordCellsList, tiffFlag, fName, iShortPlot );

	iRecordCellsList = c ( iProbeCell, iProbeCell - N, iProbeCell - 2*N, iProbeCell - 3*N );
	iRecordCellsList = c ( iProbeCell, iProbeCell + 4*N, iProbeCell + 6*N, iProbeCell + 8*N );
	fName = paste ( "D1Side", fNameBase, "tiff", sep="." );
	KnockoutTimeSeriesPlot ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleTextRoot, subTitleTextRoot, xlabTextRoot, ylabTextRoot, iRecordCellsList, tiffFlag, fName, iShortPlot );

	if ( iShortPlot ) {

		iRecordCellsList = c ( iProbeCell, iProbeCell + N );
		fName = paste ( "2Col.D2Side", fNameBase, "tiff", sep="." );
		KnockoutTimeSeriesPlot ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleTextRoot, subTitleTextRoot, xlabTextRoot, ylabTextRoot, iRecordCellsList, tiffFlag, fName, iShortPlot );

		iRecordCellsList = c ( iProbeCell, iProbeCell - N );
		fName = paste ( "2Col.D1Side", fNameBase, "tiff", sep="." );
		KnockoutTimeSeriesPlot ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleTextRoot, subTitleTextRoot, xlabTextRoot, ylabTextRoot, iRecordCellsList, tiffFlag, fName, iShortPlot );

	} # if ( iShortPlot )


} # KnockoutTimeSeriesPlotDriver = function ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,

	#
	#
	#
KnockoutTimeSeriesPlot = function ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,
					titleTextRoot, subTitleTextRoot, xlabTextRoot, ylabTextRoot, iRecordCellsList, tiffFlag, fName, iShortPlot  ) {

	if ( tiffFlag ) {
		tiff ( paste(fName, "tiff", sep="." ), compression="lzw", units="in", width=7.0, height=7.0, res=300 );
	} else {
		x11();
	} # if ( tiffFlag )

	if ( iShortPlot ) {
		par ( mfcol=c(3,2) );
	} else {
		par ( mfcol=c(3,4) );
	} # if ( iShortPlot )

	for ( iCell in iRecordCellsList ) {

		tmp = ExpRefTSResponse ( N2, numValsPerRFTrial, iCell, alldata.base, alldata.exp, eScale );

		titleText = paste ( titleTextRoot, "\nStim: ", iProbeCell, "Rec: ", iCell );
		xlabText = xlabTextRoot;
		ylabText = ylabTextRoot;

		ylim = c(tmp$vt.min, tmp$vt.max);
		ylim = c(0, 2);
		plot ( tmp$v0.exp, type="l", ylim=ylim, col=1, lty=1,
			main=titleText, xlab=xlabText, ylab=ylabText );
		lines ( tmp$v1i.exp, col=2, lty=1 );
		par(new=TRUE);
		plot ( tmp$v1e.exp, type="l", col=3, lty=1, xaxt="n", yaxt="n", xlab="", ylab="" );
		axis(4); mtext(side=4, "Ve(t)", cex=0.67, line=2);
		abline ( h=0, col=3, lty=3 );

		titleText = paste ( "Baseline Refined", "\nStim: ", iProbeCell, "Rec: ", iCell );
		plot ( tmp$v0.base, type="l", ylim=ylim, col=1, lty=1,
			main=titleText, xlab=xlabText, ylab=ylabText );
		lines ( tmp$v1i.base, col=2, lty=1 );
		par(new=TRUE);
		plot ( tmp$v1e.base, type="l", col=3, lty=1, xaxt="n", yaxt="n", xlab="", ylab="" );
		axis(4); mtext(side=4, "Ve(t)", cex=0.67, line=2);
		abline ( h=0, col=3, lty=3 );

		tmin = min ( tmp$v0.delta, tmp$v1i.delta ); tmax = max ( tmp$v0.delta, tmp$v1i.delta );
		titleText = paste ( "Diff. Exp - Baseline", "\nStim: ", iProbeCell, "Rec: ", iCell );
		plot ( tmp$v0.delta, type="l", ylim=c(tmin, tmax), col=1, lty=1,
			main=titleText, xlab=xlabText, ylab=ylabText );
		lines ( tmp$v1i.delta, col=2, lty=1 );
		par(new=TRUE);
		plot ( tmp$v1e.delta, type="l", col=3, lty=1, xaxt="n", yaxt="n", xlab="", ylab="" );
		axis(4); mtext(side=4, "Delta Ve(t)", cex=0.67, line=2);
		abline ( h=0, col=3, lty=3 );

		ylabTextRoot = "";

	} # for ( iCell in iRecordCellsList ) {

	if ( tiffFlag ) {
		dev.off();
	} # if ( tiffFlag )

} # KnockoutTimeSeriesPlot = function ( N, numValsPerRFTrial, iProbeCell...
	#
	#
	#
QuickCheckTimeSeries = function ( r1E, r1I, rZ, iCell, deltaT, plotFlag, ylim ) {

	tVal = 1:(length(r1E)) * deltaT;

	if ( plotFlag ) {
		x11();
		par(mfrow=c(3,1));
	} # if ( plotFlag )		

	plot( tVal, r1E, type="p", col=1, ylim=ylim, ylab="Firing Rate", xlab="Time (sec)");
	title(main=paste("Cortical Layer E-Cell Firing Rate",paste("Cell ",iCell),sep='\n'));

	plot( tVal, r1I, type="p", col=1, ylim=ylim, ylab="Firing Rate", xlab="Time (sec)");
	title(main=paste("Cortical Layer I-Cell Firing Rate",paste("Cell ",iCell),sep='\n'));

	plot( tVal, rZ, type="p", col=1, ylim=ylim, ylab="Firing Rate", xlab="Time (sec)");
	title(main=paste("Input Layer Firing Rate",paste("Cell ",iCell),sep='\n'));

} # QuickCheckTimeSeries = function ( r1E, r1I, r0, cellID, deltaT, ... ) {

	#
	#
	#
QuickCheckTimeSeries2 = function ( tRange, r1E, r1I, rZ, iCell, deltaT, plotFlag, titleText, ylabText ) {

	tVal = tRange * deltaT;

	if ( plotFlag ) {
		x11();
		par(mfrow=c(3,1));
	} # if ( plotFlag )		

	plot( tVal, r1E, type="p", col=1, ylab=ylabText, xlab="Time (sec)");
	title ( main = paste("Cortical Layer E-Cell", titleText, paste("Cell ",iCell),sep='\n'));

	plot( tVal, r1I, type="p", col=1, ylab=ylabText, xlab="Time (sec)");
	title(main=paste("Cortical Layer I-Cell", titleText, paste("Cell ",iCell),sep='\n'));

	plot( tVal, rZ, type="p", col=1, ylab=ylabText, xlab="Time (sec)");
	title(main=paste("Input Layer Cell", titleText, paste("Cell ",iCell),sep='\n'));

} # QuickCheckTimeSeries2 = function ( r1E, r1I, r0, cellID, deltaT, ... ) {

	#
	#
	#
GoodToShow = function ( N, iCandidate ) {

	candidateList = TrimEdgesFromCellList ( seq ( 1, N2, 1 ), N, 3 );
	show = FALSE;
	if ( sum ( candidateList == iCandidate ) ) { show = TRUE; }

	#if ( iCandidate <= N ) { show = FALSE; }
	#if ( iCandidate > (N^2 - N) ) { show = FALSE; }
	#if ( iCandidate %% N == 0 ) { show = FALSE; }
	#if ( iCandidate %% N == 1 ) { show = FALSE; }

	return ( show );

} # GoodToShow = function ( N, iCandidate ) {

	#
	#	Aids in the display of knockout zone in spatial plots.
	#
GenOutline1X = function ( tmp.x, tmp.y, kCol, kLty, kLwd ) {
	
	lines ( tmp.x, tmp.y, col=kCol, lty=kLty, lwd=kLwd );

} # GenOutline1X = function ( tmp.x, tmp.y, kCol, kLty, kLwd ) {
GenOutline4X = function ( N, x, y, doPlot ) {
	
	y[1] = y[1] - 0.5;
	y[2] = y[2] + 0.5;
	x[1] = x[1] - 0.5;
	x[2] = x[2] + 0.5;

	tmp.x = rbind(c(x[1], x[2]), c(x[1], x[2]), c(x[1], x[1]), c(x[2], x[2]));
	tmp.y = rbind(c(y[1], y[1]), c(y[2], y[2]), c(y[1], y[2]), c(y[1], y[2]));

	if ( doPlot ) {
		GenOutline1X ( tmp.x, tmp.y, "white", 1, .5 );
		GenOutline1X ( tmp.x, N+tmp.y, "white", 1, 0.5 );
		GenOutline1X ( N+tmp.x, tmp.y, "white", 1, 0.5 );
		GenOutline1X ( N+tmp.x, N+tmp.y, "white", 1, 0.5 );
	} # if ( doPlot )

} # GenOutline4X = function ( N, x, y, doPlot ) {

	#
	#
	#	These routines help with spatio-temporal series analysis of knockout effects.
	#
	#
	#
	#	PeakStats
	#
FirstPosPeakIfAny = function ( x, iStart, iLen, numPeaks ) {

	#iStart = preStimIters;
	#iLen = stimDurIters;

	thresholdVal = 0.1;
	nups = ndowns = 2;

	pkStats = rep ( 0, 3 );

	tmp = findpeaks ( x[(iStart):(iStart+iLen-1)], npeaks=numPeaks, threshold=thresholdVal, nups=nups, ndowns=ndowns );
	if ( length(tmp) ) {
		pkStats[1] = tmp[1];
		pkStats[2] = tmp[2];
		halfHeight = 0.5 * pkStats[1];
		leftWindowLoc = max ( which( x[(iStart):(iStart+pkStats[2]-1)] <= halfHeight) );	# Left-side of peak.
		rightWindowLoc = iLen - pkStats[2] + 1;
		iTmpWindowLoc = ( which( x[(iStart+pkStats[2]):(iStart+iLen-1)] >= halfHeight) );	# Left-side of peak.
		if ( length ( iTmpWindowLoc ) ) {
			rightWindowLoc = max ( iTmpWindowLoc );
		} # if ( length ( iTmpWindowLoc ) )
		pkStats[3] = ( pkStats[2] - leftWindowLoc + 1 ) + ( rightWindowLoc );
	} # if ( length(tmp) ) 

	return ( pkStats );

} # FirstPosPeakIfAny = function ( x )

FirstNegPeakIfAny = function ( x, iStart, iLen, numPeaks ) {

	#iStart = preStimIters;
	#iLen = stimDurIters;

	x = -x;

	thresholdVal = 1.0;
	nups = ndowns = 5;

	pkStats = rep ( 0, 3 );

	tmp = findpeaks ( x[(iStart):length(x)], npeaks=numPeaks, threshold=thresholdVal, nups=nups, ndowns=ndowns );
	if ( length(tmp) ) {
		pkStats[1] = tmp[1];
		pkStats[2] = tmp[2];
		halfHeight = 0.5 * pkStats[1];
		leftWindowLoc = max ( which( x[(iStart):(iStart+pkStats[2]-1)] <= halfHeight) );	# Left-side of peak.
		rightWindowLoc = iLen - pkStats[2] + 1;
		iTmpWindowLoc = ( which( x[(iStart+pkStats[2]):(iStart+iLen-1)] >= halfHeight) );	# Left-side of peak.
		if ( length ( iTmpWindowLoc ) ) {
			rightWindowLoc = max ( iTmpWindowLoc );
		} # if ( length ( iTmpWindowLoc ) )
		pkStats[3] = ( pkStats[2] - leftWindowLoc + 1 ) + ( rightWindowLoc );
	} # if ( length(tmp) ) 
	pkStats[1] = - pkStats[1];

	return ( pkStats );

} # FirstNegPeakIfAny = function ( x )


	#
	#
	#
TSSTStatsDelta = function( base, exp ) {

	pk.e = matrix ( 0, nrow=(dim(base$pk.e)[1]), ncol=6 );
	pk.i = pk.e;

	for ( i in 1:6 ) {
		pk.e[,i] = exp$pk.e[,i] - base$pk.e[,i];
		pk.i[,i] = exp$pk.i[,i] - base$pk.i[,i];
	} # for ( i in 1:6 )
	
	return ( list ( pk.e=pk.e, pk.i=pk.i ) );
	
} # TSSTStatsDelta = function( base, exp ) {

	#
	#	TS = Time Series 
	#	ST = Spatio-temporal
	#
TSSTStatsDriver = function ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base, alldata.exp, 
					titleTextRoot, subTitleTextRoot, xlabTextRoot, ylabTextRoot, tiffFlag, fName,
					xKnockBox, yKnockBox ) {

		#	xKnockBox, yKnockBox used to generate outline of the knockout zone.
	xKnockBox[1] = xKnockBox[1] - 0.5; xKnockBox[2] = xKnockBox[2] + 0.5;
	yKnockBox[1] = yKnockBox[1] - 0.5; yKnockBox[2] = yKnockBox[2] + 0.5;
	x.KnockBox = rbind(c(xKnockBox[1], xKnockBox[2]), c(xKnockBox[1], xKnockBox[2]),
					c(xKnockBox[1], xKnockBox[1]), c(xKnockBox[2], xKnockBox[2]));
	y.KnockBox = rbind( c ( yKnockBox[1], yKnockBox[1]), c(yKnockBox[2], yKnockBox[2]),
					c(yKnockBox[1], yKnockBox[2]), c(yKnockBox[1], yKnockBox[2]) );


	stats.base = TSSTStats ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.base );
	stats.exp = TSSTStats ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, alldata.exp );
	stats.delta = TSSTStatsDelta ( stats.base, stats.exp );

	xLabText  = "Distal -> Proximal";
	yLabText = "Digit 1 -> Digit 3";
	boundaryMarks = c ( as.integer(sqrt(N2)/3)+0.5, as.integer(sqrt(N2)/3)*2+0.5 );
	skipIPlots = FALSE;

	for ( iStatType in 1:3 ) {

		if ( tiffFlag ) {
			tiff ( paste("PkAz", iStatType, fName, "tiff", sep="." ), compression="lzw", units="in", width=7.0, height=7.0, res=300 );
		} else {
			x11();
		} # if ( tiffFlag )

		par(mfrow=c(3,2) );
		
		if ( iStatType == 1 ) {
			tagText = "Init Pos Pk Magnitude";
		} else if ( iStatType == 2 ) {
			tagText = "Init Pos Pk Latency";
		} else if ( iStatType == 3 ) {
			tagText = "Init Pos Pk Half Width";
		} else {
			tagText = NULL;
		} # if ( iStatType == 1 )

		titleText = paste ( "Spatial Plot E Cell", tagText, "\n", titleTextRoot, "Stim: ", iProbeCell );	
		ShowVecAsMap2 ( stats.exp$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.exp$pk.e[,iStatType ]), max(stats.exp$pk.e[,iStatType ]) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleText = paste ( "Spatial Plot I Cell", tagText, "\n", titleTextRoot, "Stim: ", iProbeCell );	
		ShowVecAsMap2 ( stats.exp$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.exp$pk.i[,iStatType ]), max(stats.exp$pk.i[,iStatType ]) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleTextRoot.base = "Baseline Refined Network";
		titleText = paste ( "Spatial Plot E Cell", tagText, "\n", titleTextRoot.base, "Stim: ", iProbeCell );	
		ShowVecAsMap2 ( stats.base$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.base$pk.e[,iStatType ]), max(stats.base$pk.e[,iStatType ]) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );


		titleText = paste ( "Spatial Plot I Cell", tagText, "\n", titleTextRoot.base, "Stim: ", iProbeCell );	
		ShowVecAsMap2 ( stats.base$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.base$pk.i[,iStatType ]), max(stats.base$pk.i[,iStatType ]) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleTextRoot.diff = "Diff Exp - Baseline";
		titleText = paste ( "Spatial Plot E Cell", tagText, "\n", titleTextRoot.diff, "Stim: ", iProbeCell );	
		ShowVecAsMap2 ( stats.delta$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.delta$pk.e[,iStatType ]), max(stats.delta$pk.e[,iStatType ]) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleText = paste ( "Spatial Plot I Cell", tagText, "\n", titleTextRoot.diff, "Stim: ", iProbeCell );	
		ShowVecAsMap2 ( stats.delta$pk.i[,iStatType ], titleText, xLabText, yLabText, min(stats.delta$pk.i[,iStatType ]), max(stats.delta$pk.i[,iStatType ]) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		if ( tiffFlag ) {
			dev.off();
		} # if ( tiffFlag )

	} # 	for ( iStatType in 1:3 ) {

	for ( iStatType in 4:6 ) {

		if ( tiffFlag ) {
			tiff ( paste("PkAz", iStatType, fName, "tiff", sep="." ), compression="lzw", units="in", width=7.0, height=7.0, res=300 );
		} else {
			x11();
		} # if ( tiffFlag )

		par(mfcol=c(3,2) );
		
		if ( iStatType == 4 ) {
			tagText = "Neg Pk Magnitude";
			skipIPlots = TRUE;
		} else if ( iStatType == 5 ) {
			tagText = "Neg Pk Latency";
			skipIPlots = TRUE;
		} else if ( iStatType == 6 ) {
			tagText = "Neg Pk Half Width";
			skipIPlots = TRUE;
		} else {
			tagText = NULL;
		} # if ( iStatType == 1 )

		titleText = paste ( "Spatial Plot E Cell", tagText, "\n", titleTextRoot, "Stim: ", iProbeCell );	
		ShowVecAsMap2 ( stats.exp$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.exp$pk.e[,iStatType ]), max(stats.exp$pk.e[,iStatType ]) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleText = paste ( "Spatial Plot E Cell", tagText, "\n", titleTextRoot, "Stim: ", iProbeCell );	
		ShowVecAsMap2 ( stats.exp$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.exp$pk.e[,iStatType ]), max(stats.exp$pk.e[,iStatType ]) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleTextRoot.base = "Baseline Refined Network";
		titleText = paste ( "Spatial Plot E Cell", tagText, "\n", titleTextRoot.base, "Stim: ", iProbeCell );	
		ShowVecAsMap2 ( stats.base$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.base$pk.e[,iStatType ]), max(stats.base$pk.e[,iStatType ]) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleTextRoot.base = "Baseline Refined Network";
		titleText = paste ( "Spatial Plot E Cell", tagText, "\n", titleTextRoot.base, "Stim: ", iProbeCell );	
		ShowVecAsMap2 ( stats.base$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.base$pk.e[,iStatType ]), max(stats.base$pk.e[,iStatType ]) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleTextRoot.diff = "Diff Exp - Baseline";
		titleText = paste ( "Spatial Plot E Cell", tagText, "\n", titleTextRoot.diff, "Stim: ", iProbeCell );	
		ShowVecAsMap2 ( stats.delta$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.delta$pk.e[,iStatType ]), max(stats.delta$pk.e[,iStatType ]) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		titleTextRoot.diff = "Diff Exp - Baseline";
		titleText = paste ( "Spatial Plot E Cell", tagText, "\n", titleTextRoot.diff, "Stim: ", iProbeCell );	
		ShowVecAsMap2 ( stats.delta$pk.e[,iStatType ], titleText, xLabText, yLabText, min(stats.delta$pk.e[,iStatType ]), max(stats.delta$pk.e[,iStatType ]) );
		abline ( h = boundaryMarks, lty=3, col=3 );
		GenOutline1X ( x.KnockBox, y.KnockBox, "white", 1, 0.5 );

		if ( tiffFlag ) {
			dev.off();
		} # if ( tiffFlag )

	} # 	for ( iStatType in 4:6 ) {

} # MPSTStats = function ( N2, numValsPerRFTrial, iProbeCell, alldata.base, alldata.exp, eScale,

	#
	#	Compute various statistics on the raw time series
	#
TSSTStats = function ( N2, numValsPerRFTrial, numItersPerTrial, preStimIters, stimDurIters, rData ) {

	numPeaks = 2;
	peakStats.e = matrix ( 0, nrow=N2, ncol=3*numPeaks );	# Col 1, 4: peakMagnitude; Col 2, 5: peak latency; Col 3, 6: 1/2 peak width.
	peakStats.i = matrix ( 0, nrow=N2, ncol=3*numPeaks );

	tsTmp.e = rep ( 0, numItersPerTrial);
	tsTmp.i = rep ( 0, numItersPerTrial);

	iPtrList.e = seq ( 1, N2*numItersPerTrial, N2 );
	iPtrList.i = numValsPerRFTrial + iPtrList.e;

	for ( iCell in 1:N2 ) {

				#	Using knowledge of data structure details to speed things up.
		tsTmp.e = ( rData [iPtrList.e] );
		tsTmp.i = ( rData [numValsPerRFTrial+iPtrList.e] );
		iPtrList.e = iPtrList.e + 1;

		peakStats.e[iCell,1:3] = FirstPosPeakIfAny ( tsTmp.e, preStimIters, stimDurIters, numPeaks );
		peakStats.i[iCell,1:3] = FirstPosPeakIfAny ( tsTmp.i, preStimIters, stimDurIters, numPeaks );

		peakStats.e[iCell,4:6] = FirstNegPeakIfAny ( tsTmp.e, preStimIters, stimDurIters, numPeaks );
		peakStats.i[iCell,4:6] = FirstNegPeakIfAny ( tsTmp.i, preStimIters, stimDurIters, numPeaks );

		#print ( iCell ); print ( warnings() );

	} # for ( iCell in 1:N2 )

	return ( list ( pk.e=peakStats.e, pk.i=peakStats.i ) );

} # TSSTStats = function ( N2, numValsPerRFTrial, iNumItersPerTrial, preStimIters, stimDurIters, rData ) {

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

iStart = 40;
iEnd = 75
ConsolidatedKnockoutTimeSeriesPlotter ( N2, numValsPerRFTrial, iProbeCell, iStart, iEnd,
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







