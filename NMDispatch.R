
#
#	NMDispatch.R
#

#
#	This is the basic call made when iterating the system.  It is used in many places.  So put it into a script.
#	Use script rather than function to prevent needless matrix copying and wasted memory and time.
#
if ( trackWeights ) {
	tmp = NMDispatch ( knmSelect, v0[k,], r0[k,], v1.e[k,], r1.e[k,], v1.i[k,], r1.i[k,], 
		v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta,
		plasticityFlag, w.tau.alpha, w.beta, w1.e.0[,,k], w1.e.i[,,k], w1.i.e[,,k], w1.e.e[,,k] );
	w1.e.0[,,iter] = tmp$w1.e.0;
	w1.e.i[,,iter] = tmp$w1.e.i;
	w1.i.e[,,iter] = tmp$w1.i.e;
	w1.e.e[,,iter] = tmp$w1.e.e;
} else {
	tmp = NMDispatch ( knmSelect, v0[k,], r0[k,], v1.e[k,], r1.e[k,], v1.i[k,], r1.i[k,], 
		v0.alpha, v1.e.alpha, v1.i.alpha, noiseLevel, beta,
		plasticityFlag, w.tau.alpha, w.beta, w1.e.0[,,1], w1.e.i[,,1], w1.i.e[,,1], w1.e.e[,,1] );
	w1.e.0[,,1] = tmp$w1.e.0;
	w1.e.i[,,1] = tmp$w1.e.i;
	w1.i.e[,,1] = tmp$w1.i.e;
	w1.e.e[,,1] = tmp$w1.e.e;
} # if ( trackWeights )

v0[iter,] = tmp$v0;
v1.e[iter,] = tmp$v1.e;
v1.i[iter,] = tmp$v1.i;

r0[iter,] = tmp$r0;
r1.e[iter,] = tmp$r1.e;
r1.i[iter,] = tmp$r1.i;



