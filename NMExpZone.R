#
#	July 23, 2014 - present
#	Kamil A. Grajski
#

#
#	Neural Modeling Experiment Zone
#

#
require(ggplot2);
require(reshape2);
require(stats);
require(stats4);

######################################################################################################
######################################################################################################
#
#	Visualization #1: Help find beta for exponential function.
#
######################################################################################################
######################################################################################################

rm(list = ls());
source ( "NMHelperFunctions.R" );
x = seq ( -2, 2, 0.001 );
beta = c ( 1, 2, 3, 4, 5 );
y = matrix ( 0, nrow = length(x), ncol = length(beta) );

y[,1] = sigmoid ( x, beta[1] );
x11();
plot( x, y[,1], type="l", col=1, main="Membrane Potential to Spike Rate Conversion",
	xlab="Membrane Potential (v)", ylab="Spike Rate (r)" );

for ( i in 2:length(beta) ) {
	y[,i] = sigmoid ( x, beta[i] );
	lines( x, y[,i], type="l", col=i );
} # for ( i in 1:length(beta) )

######################################################################################################
######################################################################################################
#
#	Experiment #1-#5: Start with 2-layer model: input, output.  Successively add complexit.
#				No plasticity.  Confirm time evolution, receptive field maps, distribution
#				of receptive field sizes, etc.
#				See: NMDispatch for details.
#
######################################################################################################
######################################################################################################

rm(list = ls());
source ( "NMHelperFunctions.R" );

#	Set random number seed.
set.seed(1);

#	Configuration Parameters

verbose = FALSE;				# How much output to generate.
trackWeights = TRUE;			# Flag whether to allocate memory to track time evolution of each weight in the system.
plasticityFlag = FALSE;			# Plasticity on or off.  Applies to all wghts in a given configuration.
#wghtNormalizeFlag = 0;			# No weight normalization.
#wghtNormalizeFlag = 1;			# Weight normalization pre-synaptically, meaning over the output weights of a given cell.
wghtNormalizeFlag = 2.1;		# Weight normalization post-synaptically, meaning over the input weights of a given cell.

nmSelect = 5;				# Selects model (which set of differential equations to iterate).  See NMDispatch.
N = 15;					# Each layer will be NxN.
N2 = N^2;					# Each layer will be NxN.

noiseLevel = 0.01;			# For background noise level sample uniform dist max value.
beta = 4.0;					# Sigmoid beta value.
tau = 0.050;				# Baseline time constant.
deltaT = tau/5;				# Time step size.
oneSecondNumIter = 1.0 / deltaT;	# The number of iterations to simulate one second of time.
v0.tau = 0.50 * tau;			# Input layer activity fade out is relatively fast.
v0.alpha = 1.0-(deltaT/v0.tau);	# Save some computations and make differential equations more readable.
v.e.tau = tau;				# Time constant for e cells membrane potential relative to baseline.
v1.e.alpha = 1.0-(deltaT/v.e.tau); 	# Save some computations and make differential equations more readable.
v.i.tau = tau;				# Time constant for i cells membrane potential relative to baseline.
v1.i.alpha = 1.0-(deltaT/v.i.tau); 	# Save some computations and make differential equations more readable.
w.tau = tau * 1000.0;			# Time constant for synaptic weight decay relative to baseline.
w.tau.alpha = 1.0-(deltaT/w.tau);	# Save some computation and make differential equations more readable.
w.beta = w.beta.copy = 0.0025;	# Maximum weight change step size per deltaT.
w.beta.decay = 0.95;			# Decay constant for reducing w.beta over time.

numRefinements = 25;			# Number of refinement runs to complete.
refinementPatchSize = 3;		# Refinement patch size.  n->(n+1)x(n+1).
refinementStimDuration = 0.25;	# Refinement stimulus duration in seconds.

numSelStimRefinements = 10;		# Number of refinement runs to complete.
selectiveStimPatchSize = 3;		# Refinement patch size.  n->(n+1)x(n+1).
selectiveStimDuration = 0.25;		# Refinement stimulus duration in seconds.
selectiveStimFactor = 10;		# How many more times to stimulate the selected area vs everything else.
selectiveStimZoneID = 1;		# Zone 1 = Digit 1 Distal; Zone 2 = Digit 2 Distal; and so on.

numSyndactRefinements = 10;		# Number of refinement runs to complete following syndactyly.
syndactPatchSize = 3;			# Syndactyly refinement patch size.  n->(n+1)x(n+1).
syndactStimDuration = 0.25;		# Syndactyly refinementstimulus duration in seconds.
syndactStimZoneID = 1;			# Zone 1 = Fuse Digits 1 & 2; Zone 2 = Fuse Digits 2 & 3

numAmputRefinements = 20;		# Number of post-amputation refinement runs to complete.
refinementAmputPatchSize = 3;		# Post-amputation Refinement patch size.  n->(n+1)x(n+1).
refinementAmputStimDuration = 0.25;	# Post-amputation Refinement stimulus duration in seconds.
kAmputDigitList = 1;		# Which digit will have a zone amputated.
kAmputZoneMaxList = 1;		# 1 = Distal; 2 = Distal + Mid; 3 = Distal + Mid + Proximal.

numCLesionRefinements = 5;		# Number of post-cortical lesion refinement runs to complete.
refinementCLesionPatchSize = 3;	# Post-cortical lesion Refinement patch size.  n->(n+1)x(n+1).
refinementCLesionStimDuration = 0.25;	# Post-cortical lesion Refinement stimulus duration in seconds.
kCLesionDigitList = c ( 1 );		# Which digit will have a zone lesioned (in baseline-type start point).
kCLesionZoneMaxList = c ( 1, 2, 3);	# 1 = Distal; 2 = Distal + Mid; 3 = Distal + Mid + Proximal.

e.wResource = 2.0; 			# Total resource available for E cell to distribute across outputs or inputs (see wghtNormalizeFlag ).
i.wResource = 1.0;			# Total resource available for I cell to distribute across outputs or inputs (see wghtNormalizeFlag ).
numSpinUp = 5;				# Number of spin up cycles.  Enough to let IC transients die out.

magRFProbe = 1.0;				# Magnitude of the receptive field mapping input layer probe.
durRFProbe = 0.10;			# Duration (in seconds) of the receptive field mapping input layer probe.
trialDurRFProbe = 1.5;			# Trial length (in seconds) of each probe during receptive field mapping, refinement, etc.
offsetToRFProbe = 0.250;		# Offset length (in seconds) to start of the RF probe.
wghtMinValue = 0.0;			# Minimum value to use when assigning intial interconnection weights.
wghtMaxValue = 2 * noiseLevel;	# Minimum value to use when assigning intial interconnection weights.
kRFPeakToEdgeDetect = 0.5;		# Determine extent of a receptive field by including all cells with peak to edge ratio > kRFPeakToEdgeDetect.

numIter = as.integer(trialDurRFProbe * oneSecondNumIter);				# Default number of iterations max (when needed).  Some experiments cap iteration max count differently.
trialDurRFProbeInIters = trialLengthInIters = numIter;				# The basic experimental epoch duration in this simulation in units of iterations.
offsetToRFProbeInIters = as.integer(offsetToRFProbe  * oneSecondNumIter);	# Stimulus ONSET time relative to trial start in units of iterations.
durRFProbeInIters = as.integer(durRFProbe * oneSecondNumIter);			# Stimulus ON time in units of iterations.
timeOffRFProbeInIters = offsetToRFProbeInIters + durRFProbeInIters;	# Stimulus OFFSET time relative to trial start in units of iterations.

trackCells = c ( 1, trunc(N2/2)+1, N2 );	# Use certain cells for diagnostics and monitoring, developing figures.

g0 = 3;	# Grid size (radius) for Layer 0 (input).  If value = 0 then it is a single point to point connection from input to cortical layer.
g1.e.e = 3;	# Grid size (radius) for Layer 1 e-cells to Layer 1 e-cells.  If value = 1 then grid is 3 x 3.
g1.e.i = 3;	# Grid size (radius) for Layer 1 e-cells to Layer 1 i-cells.  If value = 2 then grid is 5 x 5.
g1.i.e = 3; # Grid size (radius) for Layer 1 e-cells to Layer 1 i-cells.  If value = 3 then grid is 7 x 7.

#	Initialize some state variables: interconnection weights.  
if ( trackWeights ) {
	w1.e.0 = array ( 0, c ( N2, N2, numIter ) );	# For ith Layer 1 e cell, col i are weights for Layer 0 node inputs.
	w1.e.e = array ( 0, c ( N2, N2, numIter  ) );	# For ith Layer 1 e cell, col i are weights for Layer 1 e cells inputs.
	w1.e.i = array ( 0, c ( N2, N2, numIter  ) );	# For ith Layer 1 e cell, col i are weights for Layer 1 i cells inputs.
	w1.i.e = array ( 0, c ( N2, N2, numIter  ) );	# For ith Layer 1 i cell, col i are weights for Layer 1 e cells inputs.
} else {
	w1.e.0 = array ( 0, c ( N2, N2, 1 ) );	# For ith Layer 1 e cell, col i are weights for Layer 0 node inputs.
	w1.e.e = array ( 0, c ( N2, N2, 1  ) );	# For ith Layer 1 e cell, col i are weights for Layer 1 e cells inputs.
	w1.e.i = array ( 0, c ( N2, N2, 1  ) );	# For ith Layer 1 e cell, col i are weights for Layer 1 i cells inputs.
	w1.i.e = array ( 0, c ( N2, N2, 1  ) );	# For ith Layer 1 i cell, col i are weights for Layer 1 e cells inputs.
} # if ( trackWeights )

#	Initialize the connection matrices.  WghtInit2: Records # of output connections of each type (for normalization purposes).
#wtmp = WghtInit2 ( w1.e.0[,,1], g0, wghtMinValue  ); w1.e.0[,,1] = wtmp$w; n0.numPreSyn = wtmp$n.Connect;
#wtmp = WghtInit2 ( w1.e.e[,,1], g1.e.e, wghtMinValue  ); w1.e.e[,,1] = wtmp$w; n1.e.e.numPreSyn = wtmp$n.Connect;
#wtmp = WghtInit2 ( w1.e.i[,,1], g1.e.i, wghtMinValue  ); w1.e.i[,,1] = wtmp$w; n1.e.i.numPreSyn = wtmp$n.Connect;
#wtmp = WghtInit2 ( w1.i.e[,,1], g1.i.e, wghtMinValue  ); w1.i.e[,,1] = wtmp$w; n1.i.e.numPreSyn = wtmp$n.Connect;

#	Initialize the connection matrices.  WghtInit3: Same as WghtInit2, but uses a min and max initial wght range.
wtmp = WghtInit3 ( w1.e.0[,,1], g0, wghtMinValue, wghtMaxValue  );
w1.e.0[,,1] = wtmp$w; n0.numPreSyn = wtmp$n.Connect;
wtmp = WghtInit3 ( w1.e.e[,,1], g1.e.e, wghtMinValue, wghtMaxValue );
w1.e.e[,,1] = wtmp$w; n1.e.e.numPreSyn = wtmp$n.Connect;
wtmp = WghtInit3 ( w1.e.i[,,1], g1.e.i, wghtMinValue, wghtMaxValue );
w1.e.i[,,1] = wtmp$w; n1.e.i.numPreSyn = wtmp$n.Connect;
wtmp = WghtInit3 ( w1.i.e[,,1], g1.i.e, wghtMinValue, wghtMaxValue );
w1.i.e[,,1] = wtmp$w; n1.i.e.numPreSyn = wtmp$n.Connect;

#	Initialize some state variables: membrane potentials and firing rates.
#		ith row is the vector of cell values at time i
#		jth column is the time evolution of the jth cell
v0 = matrix ( 0, nrow = numIter, ncol = N2 );		# Layer 0 (input) time series: stimulus on/off indicator.
r0 = matrix ( 0, nrow = numIter, ncol = N2 );		# Layer 0 (input) time series: spatio-temporally filtered.
v1.e = matrix ( 0, nrow = numIter, ncol = N2 );		# Layer 1 E-cell "membrane potential" time series.
r1.e = matrix ( 0, nrow = numIter, ncol = N2 ); 	# Layer 1 E-cell "spiking rate" time series.
v1.i = matrix ( 0, nrow = numIter, ncol = N2 ); 	# Layer 1 I-cell "membrane potential" time series.
r1.i = matrix ( 0, nrow = numIter, ncol = N2 ); 	# Layer 1 I-cell "spiking rate" time series.

#	Set initial conditions.
initialConditions.v0 = v0[1,] = runif(N2, min = -noiseLevel, max = noiseLevel);
initialConditions.v0 = v0[1,] = runif(N2, min = -noiseLevel, max = noiseLevel);
initialConditions.r0 = r0[1,] = sigmoid ( v0[1,], beta );

initialConditions.v1.e = v1.e[1,] = runif(N2, min = -noiseLevel, max = noiseLevel);
initialConditions.r1.e = r1.e[1,] = sigmoid ( v1.e[1,], beta );
initialConditions.v1.i = v1.i[1,] = runif(N2, min = -noiseLevel, max = noiseLevel);
initialConditions.r1.i = r1.i[1,] = sigmoid ( v1.i[1,], beta );

if ( wghtNormalizeFlag == 1 ) {
		# Normalize pre-synaptically, i.e., the output weights normalize to wResource (adjusted for boundary conditions).
	initialConditions.w1.e.0 = w1.e.0[,,1] = NormalizeOutputWeights2 ( w1.e.0[,,1], n0.numPreSyn, e.wResource ); 
	initialConditions.w1.e.i = w1.e.i[,,1] = NormalizeOutputWeights2 ( w1.e.i[,,1], n1.e.i.numPreSyn, i.wResource );
	initialConditions.w1.i.e = w1.i.e[,,1] = NormalizeOutputWeights2 ( w1.i.e[,,1], n1.i.e.numPreSyn, e.wResource );
	initialConditions.w1.e.e = w1.e.e[,,1] = NormalizeOutputWeights2 ( w1.e.e[,,1], n1.e.e.numPreSyn, e.wResource );
} else if ( wghtNormalizeFlag == 2.1 ) {
		# Normalize pre-synaptically, i.e., the output weights normalize to wResource (adjusted for boundary conditions).
	initialConditions.w1.e.0 = w1.e.0[,,1] = NormalizeInputWeights1 ( w1.e.0[,,1], e.wResource ); 
	initialConditions.w1.e.i = w1.e.i[,,1] = NormalizeInputWeights1 ( w1.e.i[,,1], e.wResource );
	initialConditions.w1.i.e = w1.i.e[,,1] = NormalizeInputWeights1 ( w1.i.e[,,1], i.wResource );
	initialConditions.w1.e.e = w1.e.e[,,1] = NormalizeInputWeights1 ( w1.e.e[,,1], e.wResource );
} else if ( wghtNormalizeFlag == 2.2 ) {
		# Normalize pre-synaptically, i.e., the output weights normalize to wResource (adjusted for boundary conditions).
	initialConditions.w1.e.0 = w1.e.0[,,1] = NormalizeInputWeights2 ( w1.e.0[,,1], e.wResource ); 
	initialConditions.w1.e.i = w1.e.i[,,1] = NormalizeInputWeights2 ( w1.e.i[,,1], e.wResource );
	initialConditions.w1.i.e = w1.i.e[,,1] = NormalizeInputWeights2 ( w1.i.e[,,1], i.wResource );
	initialConditions.w1.e.e = w1.e.e[,,1] = NormalizeInputWeights2 ( w1.e.e[,,1], e.wResource );
} else {
	initialConditions.w1.e.0 = w1.e.0[,,1];
	initialConditions.w1.e.i = w1.e.i[,,1];
	initialConditions.w1.i.e = w1.i.e[,,1];
	initialConditions.w1.e.e = w1.e.e[,,1];
} # if ( wghtNormalizeFlag == 1 ) { 

#
#	For each of a range of experimental configurations:
#		Evolve the system to some equilibrium then do an RF Map.
#

tauShow = c ( tau, v0.tau, v.e.tau, v.i.tau, w.tau );
numTauShow = length(tauShow);
showIter = numIter/4;
yShow = matrix ( 0, nrow=showIter, ncol=numTauShow );
yShow[1,1:length(tauShow)] = rep ( 1.0, numTauShow );

for ( i in 2:showIter ) {
	yShow[i,1:numTauShow] = yShow[(i-1),1:numTauShow] * ( 1.0 - 1.0 / tauShow * deltaT );
} # for ( i in 2:showIter ) {
x11();
plot ( (1:showIter)*deltaT, yShow[,1], type="l", col=1, main="Time Constants",
	xlab="Time (seconds)", ylab="Relative Value");
for ( i in 2:numTauShow ) {
	lines ( (1:showIter)*deltaT, yShow[,i], col=i );
} # for ( i in 2:numTauShow )
abline(v=1.0,col=2,lty=2);

knmSelect = nmSelect;

		# Reload original initial conditions at the start of each experimental configuration.
	v0[1,] = initialConditions.v0; r0[1,] = initialConditions.r0;
	v1.e[1,] = initialConditions.v1.e; r1.e[1,] = initialConditions.r1.e;
	v1.i[1,] = initialConditions.v1.i; r1.i[1,] = initialConditions.r1.i;
	w1.e.0[,,1] = initialConditions.w1.e.0; w1.e.i[,,1] = initialConditions.w1.e.i;
	w1.i.e[,,1] = initialConditions.w1.i.e; w1.e.e[,,1] = initialConditions.w1.e.e;

		#########################################################
		#########################################################
		#
		#	Spin up the system.  Let it settle down.
		#	All weights fixed.
		#
		#########################################################
		#########################################################

	for ( iSpin in 1:numSpinUp ) {

		for ( iter in 2:numIter ) {
			k = iter - 1; source("NMDispatch.R");
		} # for ( iter in 2:numIter) )

			# 	Save the last iteration as it becomes the starting point
			#	for multiple experiments that may write over state variables v0, r0, wx.x.x, ...
		lastIter.v0 = v0[numIter,]; lastIter.r0 = r0[numIter,];
		lastIter.v1.e = v1.e[numIter,]; lastIter.r1.e = r1.e[numIter,];
		lastIter.v1.i = v1.i[numIter,]; lastIter.r1.i = r1.i[numIter,];

			#	Show an epoch of system evolution for a few sample cells on each layer.
		if ( verbose ) {
			nTrials = 1; startOffsetInIter = numIter - nTrials * trialDurRFProbeInIters + 1;
			ShowTimeSeries( v0, r0, v1.e, v1.i, nTrials, trialDurRFProbeInIters, startOffsetInIter, trackCells, knmSelect, deltaT );
		} # if ( verbose )
			#	Before going on to next experiment, restore the system to the last free-running iteration.
		v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
		v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
		v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;

	} # for ( iSpin in 1:numSpinUp ) {

		#	At this point the system should be free of transients and just bubbling along ready for experiments.
	#source("FreeRunWithPlasticityON.R");		#	A useful check that weight plasticity is operating.	

		#	Experimentally determine the receptive field map for cortical layer E and I cells.
	source("RFMapExp1.R");	# Use source rather than function to save memory and time.
	save.image ( file = paste ("TRun.",numRefinements,".0",sep="" ) );

		#	Before going on to next experiment, restore the system to the desired state.
	v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
	v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
	v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;

		##############################################################
		##############################################################
		#
		#	Evolve a baseline somatotopic map.
		#	Weights will be adapted.
		#
		##############################################################
		##############################################################
	
	for ( iRefinements in 1:numRefinements ) {

		source("NMRefinement1.R");
		w.beta = w.beta.decay * w.beta;
	
			# 	Save the last iteration as it becomes the controlled starting point
			#	for next stage, i.e., we don't want to start from where RF mapping ended.
		lastIter.v0 = v0[numIter,]; lastIter.r0 = r0[numIter,];
		lastIter.v1.e = v1.e[numIter,]; lastIter.r1.e = r1.e[numIter,];
		lastIter.v1.i = v1.i[numIter,]; lastIter.r1.i = r1.i[numIter,];
		lastIter.w1.e.0 = w1.e.0[,,numIter]; lastIter.w1.e.i = w1.e.i[,,numIter];
		lastIter.w1.e.e = w1.e.e[,,numIter]; lastIter.w1.i.e = w1.i.e[,,numIter];

			#	Optionally do an RF Mapping.
		source("RFMapExp1.R");
		save.image ( file = paste ("Run.",numRefinements,".",iRefinements,sep="" ) );

			# 	Restore the last iteration as it becomes the controlled starting point
			#	for next stage, i.e., we don't want to start from where RF mapping ended.
		v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
		v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
		v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
		w1.e.0[,,1] = lastIter.w1.e.0; w1.e.i[,,1] = lastIter.w1.e.i;
		w1.e.e[,,1] = lastIter.w1.e.e; w1.i.e[,,1] = lastIter.w1.i.e;

	} # for ( iRefinements in 1:numRefinements ) {

		##############################################################
		##############################################################
		#
		#	COACTIVATION OF A LIMITED INPUT SURFACE IN AN INTACT SYSTEM
		#	e.g., Driving correlated input during a behavioral task.
		#
		#	Begin experimental manipulations with baseline soma-
		#	totopic map as a starting point.  Same basic set up
		#	as the refinment phase, just different stimulus pattern.
		#
		##############################################################
		##############################################################

			#	Optionally load the last saved image.
	load( file = paste ("Run.",numRefinements,".",numRefinements,sep="" ) );
	v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
	v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
	v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
	w1.e.0[,,1] = lastIter.w1.e.0; w1.e.i[,,1] = lastIter.w1.e.i;
	w1.e.e[,,1] = lastIter.w1.e.e; w1.i.e[,,1] = lastIter.w1.i.e;

			#	Restore certain variables.
	w.beta = w.beta.copy

	for ( iRefinements in 1:numSelStimRefinements ) {

		source("NMSelectiveStimulation1.R");
		w.beta = w.beta.decay * w.beta;
	
			# 	Save the last iteration as it becomes the controlled starting point.
		lastIter.v0 = v0[numIter,]; lastIter.r0 = r0[numIter,];
		lastIter.v1.e = v1.e[numIter,]; lastIter.r1.e = r1.e[numIter,];
		lastIter.v1.i = v1.i[numIter,]; lastIter.r1.i = r1.i[numIter,];
		lastIter.w1.e.0 = w1.e.0[,,numIter]; lastIter.w1.e.i = w1.e.i[,,numIter];
		lastIter.w1.e.e = w1.e.e[,,numIter]; lastIter.w1.i.e = w1.i.e[,,numIter];

			#	Optionally do an RF Mapping.
		source("RFMapExp1.R");
		save.image ( file = paste ("SelStim.", selectiveStimZoneID, ".", numSelStimRefinements, ".", iRefinements, sep="" ) );

			# 	Restore the last iteration as it becomes the controlled starting point
		v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
		v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
		v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
		w1.e.0[,,1] = lastIter.w1.e.0; w1.e.i[,,1] = lastIter.w1.e.i;
		w1.e.e[,,1] = lastIter.w1.e.e; w1.i.e[,,1] = lastIter.w1.i.e;

	} # for ( iRefinements in 1:numSelStimRefinements ) {


		##############################################################
		##############################################################
		#
		#	SYNDACTYLY: Part I
		#
		#	Begin experimental manipulations with baseline soma-
		#	totopic map as a starting point.  Syndactyly is implemented
		#	in this model as double-digit input stimulus patterns
		#
		##############################################################
		##############################################################

			#	Load a refined baseline somatotopic network.
	load( file = paste ("Run.",numRefinements,".",numRefinements,sep="" ) );
	v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
	v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
	v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
	w1.e.0[,,1] = lastIter.w1.e.0; w1.e.i[,,1] = lastIter.w1.e.i;
	w1.e.e[,,1] = lastIter.w1.e.e; w1.i.e[,,1] = lastIter.w1.i.e;

			#	Restore certain variables.
	w.beta = w.beta.copy

	for ( iRefinements in 1:numSyndactRefinements ) {

		source("NMSyndactyly1.R");
		w.beta = w.beta.decay * w.beta;
	
			# 	Save the last iteration as it becomes the controlled starting point.
		lastIter.v0 = v0[numIter,]; lastIter.r0 = r0[numIter,];
		lastIter.v1.e = v1.e[numIter,]; lastIter.r1.e = r1.e[numIter,];
		lastIter.v1.i = v1.i[numIter,]; lastIter.r1.i = r1.i[numIter,];
		lastIter.w1.e.0 = w1.e.0[,,numIter]; lastIter.w1.e.i = w1.e.i[,,numIter];
		lastIter.w1.e.e = w1.e.e[,,numIter]; lastIter.w1.i.e = w1.i.e[,,numIter];

			#	Do an RF Mapping.
		source("RFMapExp1.R");
		save.image ( file = paste ("Syndact.", syndactStimZoneID, ".", numSyndactRefinements , ".", iRefinements, sep="" ) );

			# 	Restore the last iteration as it becomes the controlled starting point
		v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
		v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
		v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
		w1.e.0[,,1] = lastIter.w1.e.0; w1.e.i[,,1] = lastIter.w1.e.i;
		w1.e.e[,,1] = lastIter.w1.e.e; w1.i.e[,,1] = lastIter.w1.i.e;

	} # for ( iRefinements in 1:numSyndactRefinements ) {

		##############################################################
		##############################################################
		#
		#	SYNDACTYLY: Part II.  Release of syndactyly.
		#		Modeled simply as evolving a baseline type
		#		map.
		#
		##############################################################
		##############################################################

	load ( file = paste ("Syndact.", syndactStimZoneID, ".", numSyndactRefinements, ".", numSyndactRefinements, sep="" ) );
	v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
	v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
	v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
	w1.e.0[,,1] = lastIter.w1.e.0; w1.e.i[,,1] = lastIter.w1.e.i;
	w1.e.e[,,1] = lastIter.w1.e.e; w1.i.e[,,1] = lastIter.w1.i.e;

			#	Restore beta weight factor.
	w.beta = w.beta.copy;

	for ( iRefinements in 1:numRefinements ) {

		source("NMRefinement1.R");
		w.beta = w.beta.decay * w.beta;
	
			# 	Save the last iteration as it becomes the controlled starting point
			#	for next stage, i.e., we don't want to start from where RF mapping ended.
		lastIter.v0 = v0[numIter,]; lastIter.r0 = r0[numIter,];
		lastIter.v1.e = v1.e[numIter,]; lastIter.r1.e = r1.e[numIter,];
		lastIter.v1.i = v1.i[numIter,]; lastIter.r1.i = r1.i[numIter,];
		lastIter.w1.e.0 = w1.e.0[,,numIter]; lastIter.w1.e.i = w1.e.i[,,numIter];
		lastIter.w1.e.e = w1.e.e[,,numIter]; lastIter.w1.i.e = w1.i.e[,,numIter];

			#	Optionally do an RF Mapping.
		source("RFMapExp1.R");
		save.image ( file = paste ("Syndact.Release.", syndactStimZoneID, ".", numRefinements, ".", iRefinements, sep="" ) );

			# 	Restore the last iteration as it becomes the controlled starting point
			#	for next stage, i.e., we don't want to start from where RF mapping ended.
		v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
		v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
		v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
		w1.e.0[,,1] = lastIter.w1.e.0; w1.e.i[,,1] = lastIter.w1.e.i;
		w1.e.e[,,1] = lastIter.w1.e.e; w1.i.e[,,1] = lastIter.w1.i.e;

	} # for ( iRefinements in 1:numRefinements ) {

		##############################################################
		##############################################################
		#
		#	AMPUTATION
		#		Modeled starting with a baseline topographic net-
		#		work and then setting input layer to layer 1 E-cell
		#		weights to zero and then following the baseline
		#		refinement protocol.
		#
		##############################################################
		##############################################################

			#	Optionally load the last saved baseline image.
	load( file = paste ("Run.",numRefinements,".",numRefinements,sep="" ) );
	source ( "NMHelperFunctions.R" );	# The Run.x.x file may have carry old Helper Functions.
	v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
	v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
	v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
	w1.e.0[,,1] = lastIter.w1.e.0; w1.e.i[,,1] = lastIter.w1.e.i;
	w1.e.e[,,1] = lastIter.w1.e.e; w1.i.e[,,1] = lastIter.w1.i.e;

			#	Do the "amputation" and immediately do an RF map experiment.
			#	In this configuration, there is no renormalization of weights in the target cells.
			#	A variation on this experiment could be to renormalize, or phase in renormalization.
	w1.e.0 = Dig3Amputate1 ( w1.e.0, kAmputDigit, kAmputZoneMax );
	source("RFMapExp1.R");
	save.image ( file = paste ("Amput", kAmputDigit, kAmputZoneMax, numAmputRefinements, "0", sep="." ) );
	v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
	v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
	v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
	w1.e.i[,,1] = lastIter.w1.e.i; w1.e.e[,,1] = lastIter.w1.e.e; w1.i.e[,,1] = lastIter.w1.i.e;
			
			#	Restore beta weight factor.
	w.beta = w.beta.copy;

			#	Now follow the normal map refinement protocol.  Same as NMRefinement1.R,
			#	but knowing that some input zones are knocked out, can make the computation faster
			# 	by skipping those zones.
	for ( iRefinements in 1:numAmputRefinements ) {

		source("NMPostAmputRefinement1.R");
		w.beta = w.beta.decay * w.beta;
	
			# 	Save the last iteration as it becomes the controlled starting point
			#	for next stage, i.e., we don't want to start from where RF mapping ended.
		lastIter.v0 = v0[numIter,]; lastIter.r0 = r0[numIter,];
		lastIter.v1.e = v1.e[numIter,]; lastIter.r1.e = r1.e[numIter,];
		lastIter.v1.i = v1.i[numIter,]; lastIter.r1.i = r1.i[numIter,];
		lastIter.w1.e.0 = w1.e.0[,,numIter]; lastIter.w1.e.i = w1.e.i[,,numIter];
		lastIter.w1.e.e = w1.e.e[,,numIter]; lastIter.w1.i.e = w1.i.e[,,numIter];

			#	Optionally do an RF Mapping.
		source("RFMapExp1.R");
		save.image ( file = paste ("Amput", kAmputDigit, kAmputZoneMax, numAmputRefinements, iRefinements, sep="." ) );

			# 	Restore the last iteration as it becomes the controlled starting point
			#	for next stage, i.e., we don't want to start from where RF mapping ended.
		v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
		v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
		v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
		w1.e.0[,,1] = lastIter.w1.e.0; w1.e.i[,,1] = lastIter.w1.e.i;
		w1.e.e[,,1] = lastIter.w1.e.e; w1.i.e[,,1] = lastIter.w1.i.e;

	} # for ( iRefinements in 1:numAmputRefinements ) {


		##############################################################
		##############################################################
		#
		#	CORTICAL LESION
		#		Modeled starting with a baseline topographic net-
		#		work and then setting output layer E and I cell
		#		output weights to zero.  (In such case no need
		#		worry about the v and r variables.  The zero
		#		weight is enough.  Need to double check this.)
		#
		##############################################################
		##############################################################

	for ( kCLesionDigit in kCLesionDigitList ) {

		for ( kCLesionZoneMax in kCLesionZoneMaxList ) {

				#	Load the last saved baseline image.
			load( file = paste ("Run.",numRefinements,".",numRefinements,sep="" ) );
			source ( "NMHelperFunctions.R" );	# The Run.x.x file may have carry old Helper Functions.

				#	Bootstrap the circular buffers (with numIter-th outputs).
			v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
			v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
			v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
			w1.e.0[,,1] = lastIter.w1.e.0; w1.e.i[,,1] = lastIter.w1.e.i;
			w1.e.e[,,1] = lastIter.w1.e.e; w1.i.e[,,1] = lastIter.w1.i.e;

				#	Do the "lesion".
			w1.e.0 = Dig3CLesion1 ( w1.e.0, kCLesionDigit, kCLesionZoneMax, FALSE, TRUE );
			w1.e.e = Dig3CLesion1 ( w1.e.e, kCLesionDigit, kCLesionZoneMax, TRUE, TRUE );
			w1.e.i = Dig3CLesion1 ( w1.e.i, kCLesionDigit, kCLesionZoneMax, TRUE, TRUE );
			w1.i.e = Dig3CLesion1 ( w1.i.e, kCLesionDigit, kCLesionZoneMax, TRUE, TRUE );

				#	Do an RF map before any refinement.
			source("RFMapExp1.R");
			save.image ( file = paste ("CLesion1", kCLesionDigit, kCLesionZoneMax, numCLesionRefinements, "0", sep="." ) );
			v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
			v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
			v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
			
				#	Restore beta weight factor.
			w.beta = w.beta.copy;

				#	Now follow the normal map refinement protocol.  Same as NMRefinement1.R,
				#	but knowing that some input zones are knocked out, can make the computation faster
				# 	by skipping those zones.
			for ( iRefinements in 1:numCLesionRefinements ) {

				source("NMRefinement1.R");
				w.beta = w.beta.decay * w.beta;
	
					# 	Save the last iteration as it becomes the controlled starting point
					#	for next stage, i.e., we don't want to start from where RF mapping ended.
				lastIter.v0 = v0[numIter,]; lastIter.r0 = r0[numIter,];
				lastIter.v1.e = v1.e[numIter,]; lastIter.r1.e = r1.e[numIter,];
				lastIter.v1.i = v1.i[numIter,]; lastIter.r1.i = r1.i[numIter,];
				lastIter.w1.e.0 = w1.e.0[,,numIter]; lastIter.w1.e.i = w1.e.i[,,numIter];
				lastIter.w1.e.e = w1.e.e[,,numIter]; lastIter.w1.i.e = w1.i.e[,,numIter];

					#	Do RF Map.
				source("RFMapExp1.R");
				save.image ( file = paste ( "CLesion1", kCLesionDigit, kCLesionZoneMax, numCLesionRefinements, iRefinements, sep="." ) );

					# 	Restore the last iteration as it becomes the controlled starting point
					#	for next stage, i.e., we don't want to start from where RF mapping ended.
				v0[1,] = lastIter.v0; r0[1,] = lastIter.r0;
				v1.e[1,] = lastIter.v1.e; r1.e[1,] = lastIter.r1.e;
				v1.i[1,] = lastIter.v1.i; r1.i[1,] = lastIter.r1.i;
				w1.e.0[,,1] = lastIter.w1.e.0; w1.e.i[,,1] = lastIter.w1.e.i;
				w1.e.e[,,1] = lastIter.w1.e.e; w1.i.e[,,1] = lastIter.w1.i.e;

			} # for ( iRefinements in 1:numCLesionRefinements ) {

		} # for ( kCLesionZone in kCLesionZoneMaxList ) {

	} # 	for ( kCLesionDigit in kCLesionDigitList ) {





