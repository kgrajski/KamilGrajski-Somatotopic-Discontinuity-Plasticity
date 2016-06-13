		#	To support use of a common scale do some type-specific min/max.

for ( iCell in div.base.final.cellList ) {

	x11(); par(mfrow=c(3,3));
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

	min.w1.e0 = min ( GetInputWeights ( init.w1.e0, iCell ), GetInputWeights ( refine.w1.e0, iCell ), GetInputWeights ( syndact.w1.e0, iCell ) );
	max.w1.e0 = max ( GetInputWeights ( init.w1.e0, iCell ), GetInputWeights ( refine.w1.e0, iCell ), GetInputWeights ( syndact.w1.e0, iCell ) );

	min.w1.ee = min ( GetInputWeights ( init.w1.ee , iCell ), GetInputWeights ( refine.w1.ee , iCell ), GetInputWeights ( syndact.w1.ee , iCell ) );
	max.w1.ee = max ( GetInputWeights ( init.w1.ee , iCell ), GetInputWeights ( refine.w1.ee , iCell ), GetInputWeights ( syndact.w1.ee , iCell ) );
	
	min.w1.ei = min ( GetInputWeights ( init.w1.ei , iCell ), GetInputWeights ( refine.w1.ei , iCell ), GetInputWeights ( syndact.w1.ei , iCell ) );
	max.w1.ei = max ( GetInputWeights ( init.w1.ei , iCell ), GetInputWeights ( refine.w1.ei , iCell ), GetInputWeights ( syndact.w1.ei , iCell ) );
		
		#	Initial Conditions

	iWhich = 1;
	titleText = paste ( "Col. # ", iCell, " E<-S Weights\nBaseline Iter ", (iWhich-1), sep="" );
	ShowVecAsMap2 ( GetInputWeights ( init.w1.e0, iCell ), titleText, xlab, ylab, min.w1.e0, max.w1.e0 );

	iWhich = 1;
	titleText = paste ( "Col. # ", iCell, " E<-E Weights\nBaseline Iter ", (iWhich-1), sep="" );
	ShowVecAsMap2 ( GetInputWeights ( init.w1.ee, iCell ), titleText, xlab, ylab, min.w1.ee, max.w1.ee );

	iWhich = 1;
	titleText = paste ( "Col. # ", iCell, " E<-I Weights\nBaseline Iter ", (iWhich-1), sep="" );
	ShowVecAsMap2 ( GetInputWeights ( init.w1.ei, iCell ), titleText, xlab, ylab, min.w1.ei, max.w1.ei );

		#	Refined Map

	iWhich = 26;
	titleText = paste ( "Col. # ", iCell, " E<-S Weights\nBaseline Iter ", (iWhich-1), sep="" );
	ShowVecAsMap2 ( GetInputWeights ( refine.w1.e0, iCell ), titleText, xlab, ylab, min.w1.e0, max.w1.e0 );

	iWhich = 26;
	titleText = paste ( "Col. # ", iCell, " E<-E Weights\nBaseline Iter ", (iWhich-1), sep="" );
	ShowVecAsMap2 ( GetInputWeights ( refine.w1.ee, iCell ), titleText, xlab, ylab, min.w1.ee, max.w1.ee );

	iWhich = 26;
	titleText = paste ( "Col. # ", iCell, " E<-I Weights\nBaseline Iter ", (iWhich-1), sep="" );
	ShowVecAsMap2 ( GetInputWeights ( refine.w1.ei, iCell ), titleText, xlab, ylab, min.w1.ei, max.w1.ei );

		#	Syndactyly

	iWhich = 25;
	titleText = paste ( "Col. # ", iCell, " E<-S Weights\nBaseline Iter ", (iWhich), sep="" );
	ShowVecAsMap2 ( GetInputWeights ( syndact.w1.e0, iCell ), titleText, xlab, ylab, min.w1.e0, max.w1.e0 );

	iWhich = 25;
	titleText = paste ( "Col. # ", iCell, " E<-E Weights\nBaseline Iter ", (iWhich), sep="" );
	ShowVecAsMap2 ( GetInputWeights ( syndact.w1.ee, iCell ), titleText, xlab, ylab, min.w1.ee, max.w1.ee );

	iWhich = 25;
	titleText = paste ( "Col. # ", iCell, " E<-I Weights\nBaseline Iter ", (iWhich), sep="" );
	ShowVecAsMap2 ( GetInputWeights ( syndact.w1.ei, iCell ), titleText, xlab, ylab, min.w1.ei, max.w1.ei );

} # for ( iCell in div.base.final.cellList ) {

	#	Same as above, but dropping the initial conditions to focus on the refined & syndactyly.
	#	Also, letting each adjust to its own scale.

for ( iCell in div.base.final.cellList ) {

	x11(); par(mfrow=c(2,2));
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

		#	Initial Conditions

	#iWhich = 1;
	#titleText = paste ( "Col. # ", iCell, " E<-S Weights\nBaseline Iter ", (iWhich-1), sep="" );
	#ShowVecAsMap2 ( GetInputWeights ( init.w1.e0, iCell ), titleText, xlab, ylab, min.w1.e0, max.w1.e0 );

	#iWhich = 1;
	#titleText = paste ( "Col. # ", iCell, " E<-E Weights\nBaseline Iter ", (iWhich-1), sep="" );
	#ShowVecAsMap2 ( GetInputWeights ( init.w1.ee, iCell ), titleText, xlab, ylab, min.w1.ee, max.w1.ee );

	#iWhich = 1;
	#titleText = paste ( "Col. # ", iCell, " E<-I Weights\nBaseline Iter ", (iWhich-1), sep="" );
	#ShowVecAsMap2 ( GetInputWeights ( init.w1.ei, iCell ), titleText, xlab, ylab, min.w1.ei, max.w1.ei );

		#	Refined Map

	iWhich = 26;
	titleText = paste ( "Col. # ", iCell, " E<-S Weights\nBaseline Iter ", (iWhich-1), sep="" );
	ShowVecAsMap1 ( GetInputWeights ( refine.w1.e0, iCell ), titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Col. # ", iCell, " E<-E Weights\nBaseline Iter ", (iWhich-1), sep="" );
	ShowVecAsMap1 ( GetInputWeights ( refine.w1.ee, iCell ), titleText, xlab, ylab );

	#iWhich = 26;
	#titleText = paste ( "Col. # ", iCell, " E<-I Weights\nBaseline Iter ", (iWhich-1), sep="" );
	#ShowVecAsMap1 ( GetInputWeights ( refine.w1.ei, iCell ), titleText, xlab, ylab );

		#	Syndactyly

	iWhich = 25;
	titleText = paste ( "Col. # ", iCell, " E<-S Weights\nBaseline Iter ", (iWhich), sep="" );
	ShowVecAsMap1 ( GetInputWeights ( syndact.w1.e0, iCell ), titleText, xlab, ylab );

	iWhich = 25;
	titleText = paste ( "Col. # ", iCell, " E<-E Weights\nBaseline Iter ", (iWhich), sep="" );
	ShowVecAsMap1 ( GetInputWeights ( syndact.w1.ee, iCell ), titleText, xlab, ylab );

	#iWhich = 25;
	#titleText = paste ( "Col. # ", iCell, " E<-I Weights\nBaseline Iter ", (iWhich), sep="" );
	#ShowVecAsMap1 ( GetInputWeights ( syndact.w1.ei, iCell ), titleText, xlab, ylab );

	x11(); par(mfrow=c(2,2));
	xlab="Distal -> Proximal"; ylab="Digit 1 -> Digit 3";

		#	Initial Conditions

	#iWhich = 1;
	#titleText = paste ( "Col. # ", iCell, " E<-S Weights\nBaseline Iter ", (iWhich-1), sep="" );
	#ShowVecAsMap2 ( GetInputWeights ( init.w1.e0, iCell ), titleText, xlab, ylab, min.w1.e0, max.w1.e0 );

	#iWhich = 1;
	#titleText = paste ( "Col. # ", iCell, " E<-E Weights\nBaseline Iter ", (iWhich-1), sep="" );
	#ShowVecAsMap2 ( GetInputWeights ( init.w1.ee, iCell ), titleText, xlab, ylab, min.w1.ee, max.w1.ee );

	#iWhich = 1;
	#titleText = paste ( "Col. # ", iCell, " E<-I Weights\nBaseline Iter ", (iWhich-1), sep="" );
	#ShowVecAsMap2 ( GetInputWeights ( init.w1.ei, iCell ), titleText, xlab, ylab, min.w1.ei, max.w1.ei );

		#	Refined Map

	#iWhich = 26;
	#titleText = paste ( "Col. # ", iCell, " E<-S Weights\nBaseline Iter ", (iWhich-1), sep="" );
	#ShowVecAsMap1 ( GetInputWeights ( refine.w1.e0, iCell ), titleText, xlab, ylab );

	#iWhich = 26;
	#titleText = paste ( "Col. # ", iCell, " E<-E Weights\nBaseline Iter ", (iWhich-1), sep="" );
	#ShowVecAsMap1 ( GetInputWeights ( refine.w1.ee, iCell ), titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Col. # ", iCell, " E<-I Weights\nBaseline Iter ", (iWhich-1), sep="" );
	ShowVecAsMap1 ( GetInputWeights ( refine.w1.ei, iCell ), titleText, xlab, ylab );

	iWhich = 26;
	titleText = paste ( "Col. # ", iCell, " I<-E Weights\nBaseline Iter ", (iWhich-1), sep="" );
	ShowVecAsMap1 ( GetInputWeights ( refine.w1.ie, iCell ), titleText, xlab, ylab );

		#	Syndactyly

	#iWhich = 25;
	#titleText = paste ( "Col. # ", iCell, " E<-S Weights\nBaseline Iter ", (iWhich), sep="" );
	#ShowVecAsMap1 ( GetInputWeights ( syndact.w1.e0, iCell ), titleText, xlab, ylab );

	#iWhich = 25;
	#titleText = paste ( "Col. # ", iCell, " E<-E Weights\nBaseline Iter ", (iWhich), sep="" );
	#ShowVecAsMap1 ( GetInputWeights ( syndact.w1.ee, iCell ), titleText, xlab, ylab );

	iWhich = 25;
	titleText = paste ( "Col. # ", iCell, " E<-I Weights\nBaseline Iter ", (iWhich), sep="" );
	ShowVecAsMap1 ( GetInputWeights ( syndact.w1.ei, iCell ), titleText, xlab, ylab );

	iWhich = 25;
	titleText = paste ( "Col. # ", iCell, " E<-I Weights\nBaseline Iter ", (iWhich), sep="" );
	ShowVecAsMap1 ( GetInputWeights ( syndact.w1.ie, iCell ), titleText, xlab, ylab );

} # for ( iCell in div.base.final.cellList ) {