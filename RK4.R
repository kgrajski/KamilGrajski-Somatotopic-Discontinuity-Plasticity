#	RK4.R
#

rm(list = ls());

	#
	#	Refamiliarization with RK using toy problems.
	#

		#	f1: Simple exponential decay.

x11(); par(mfrow=c(2,1));

tau = 0.050;
noiseLevel = 0.1
h = 0.001;
numIter = as.integer ( 1 / h );
toShow = 1:(numIter);

I = rep ( 0, numIter );
I.start = as.integer(numIter/10);
I.end = (1/2) * I.start - 1;
I[I.start:I.end] = 0.25;
I = I + runif ( numIter, -noiseLevel, noiseLevel );

y.E = rep ( 0, numIter );
y.RK4 = rep ( 0, numIter );
y.E[1] = y.RK4[1] = I[1];

#
#	The function

f1= function ( tau, y ) {
	return ( -y/tau );
} # f = function ( i, y, tau ) {

#
#	Euler
#
for ( i in 2:numIter ) {
	y.E[i] = y.E[i-1] + h * f1(tau, y.E[i-1]) + I[i-1];
} # for ( i in 2:numIter ) {

#
#	RK4
#
for ( i in 2:numIter ) {
	k1 = f1(tau, y.E[i-1]);
	k2 = f1(tau, y.E[i-1] + h*k1/2);
	k3 = f1(tau, y.E[i-1] + h*k2/2);
	k4 = f1(tau, y.E[i-1] + h*k3);
	T4 = (1/6) * (k1 + 2*k2 + 2 * k3 + k4);
	y.RK4[i] = y.RK4[i-1] + h * T4 + I[i-1];	
} # for ( i in 2:numIter ) {


ylim = c(min(y.E, y.RK4, I), max(y.E, y.RK4, I));
plot(toShow,y.E[toShow], type="l", col=1, ylim=ylim);
lines(toShow,y.RK4[toShow],col=2);
lines(toShow, I[toShow], col=3);





#	f1: Simple exponential decay.

x11(); par(mfrow=c(2,1));

tau = 0.025;
noiseLevel = 0.0
h = 0.001;
numIter = as.integer ( 1 / h );
toShow = 1:(numIter);

I = rep ( 0, numIter );
I.start = as.integer(numIter/10);
I.end = 2 * I.start - 1;
I[I.start:I.end] = 1.0;
I = I + runif ( numIter, -noiseLevel, noiseLevel );

y.E = rep ( 0, numIter );
y.RK4 = rep ( 0, numIter );
y.E[1] = y.RK4[1] = I[1];

#
#	The function

f2= function ( tau, y, I ) {
	return ( -y/tau + I );
} # f2 = function ( i, y, tau ) {

#
#	Euler
#
for ( i in 2:numIter ) {
	y.E[i] = y.E[i-1] + h * f2(tau, y.E[i-1], I[i-1]);
} # for ( i in 2:numIter ) {

#
#	RK4
#
for ( i in 2:numIter ) {
	k1 = f2(tau, y.E[i-1], I[i-1]);
	k2 = f2(tau, y.E[i-1] + h*k1/2, I[i-1]);
	k3 = f2(tau, y.E[i-1] + h*k2/2, I[i-1]);
	k4 = f2(tau, y.E[i-1] + h*k3, I[i-1]);
	T4 = (1/6) * (k1 + 2*k2 + 2 * k3 + k4);
	y.RK4[i] = y.RK4[i-1] + h * T4 + I[i-1];	
} # for ( i in 2:numIter ) {


ylim = c(min(y.E, y.RK4, I), max(y.E, y.RK4, I));
plot(toShow,y.E[toShow], type="l", col=1, ylim=ylim);
lines(toShow,y.RK4[toShow],col=2);
lines(toShow, I[toShow], col=3);





































#
#	RK4
#

for ( i in 2:numIter ) {
	k1 = f1(tau, y.E[i-1], I[i-1]);
	k2 = f1(tau, y.E[i-1] + h*k1/2, I[i-1]);
	k3 = f1(tau, y.E[i-1] + h*k2/2, I[i-1]);
	k4 = f1(tau, y.E[i-1] + h*k3, I[i-1]);
	T4 = (1/6) * (k1 + 2*k2 + 2 * k3 + k4);
	y.RK4[i] = y.RK4[i-1] + h * T4;	
} # for ( i in 2:numIter ) {

x11();
ylim = c(min(y.E, y.RK4, I), max(y.E, y.RK4, I));
plot(toShow,y.E[toShow], type="l", col=1, ylim=ylim);
lines(toShow,y.RK4[toShow],col=2);
lines(toShow, I[toShow], col=3);






















#
#	From Video Example: OK
#

#
# 	dy/dx = (3)(x^2)(y)
#	x0 = 1; y0 = 2; h = 0.1;
#

numIter = 10;
y = rep ( 0, numIter );
x = rep ( 0, numIter );
x[1] = 1.0; y[1] = 2.0;
h = 0.1

for ( i in 2:numIter ) {
	k1 = 3.0 * (x[i-1])^2 * y[i-1];
	k2 = 3.0 * (x[i-1] + h/2)^2 * (y[i-1] + (h/2)*k1);
	k3 = 3.0 * (x[i-1] + h/2)^2 * (y[i-1] + (h/2)*k2);
	k4 = 3.0 * (x[i-1] + h)^2 * (y[i-1] + h * k3);
	T4 = (1/6) * (k1 + 2*k2 + 2 * k3 + k4);
	y[i] = y[i-1] + h * T4;	
} # for ( i in 2:numIter ) {



#
#	From Book Example: NOT OK
#

#	dx/dt = -x; x(0) = 1; Find x(1) with different h.


h = c ( 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.001, 0.0001 );
numH = length ( h );

computedPoint = rep ( 0, numH );

f = function ( x ) {
	return ( x );
} # f = function ( x )

for ( ih in 1:length(h) ) {

	deltaT = h[ih];
	numIter = 2 * as.integer ( 1.0/deltaT );
	x = rep ( 0, numIter );
	x[1] = 1.0;

	for ( j in 2:numIter ) {

		k1 = f(x[j-1]);
		k2 = f(x[j-1] + (deltaT/2)*k1);
		k3 = f(x[j-1] + (deltaT/2)*k2);
		k4 = f(x[j-1] + deltaT*k3);
		T4 = (1/6) * ( k1 + 2.0 * k2 + 2.0 * k3 + k4 );
		x[j] = x[j-1] + deltaT * T4;

	}

	computedPoint[ih] = x[numIter];
}



#
#	From Sample Code: OK
#

f = function ( t, y ) {
	return ( y - t^2 + 1 );
}

#rk = function () {
	h = 0.2;
	t = 0;
	w = 0.5;
	format('Step 0: t = %12.8f, w = %12.8f\n', t, w);
	for ( i in 1:(2*as.integer(1/h)) ) {
		k1 = h * f(t,w);
		k2 = h * f(t+h/2, w+k1/2);
		k3 = h * f(t+h/2, w+k2/2);
		k4 = h * f(t+h, w+k3);
		w = w + (k1 + 2*k2 + 2*k3 + k4)/6;
		t = t + h;
		format('Step %d: t = %6.4f, w = %18.15f\n', i, t, w);
	}
#} # rk = function ()

