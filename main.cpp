//	main.cpp

//	Kamil A. Grajski
//	2-November-2014
//

//
//	Driver program for neural network modeling.
//

#include "simmgr.h"

void main() {

	string simExpFileName;

	cout << "Welcome to NMExp" << endl;
	cout << "Enter Simulation Parameters File Name: ";
	cin >> simExpFileName;
	
	SimMgr SimMgr ( simExpFileName );

		// Work the "fix."
	char theFix;
	cout << endl << endl << "Enter any character to quit." << endl;
	cin >> theFix;

};	// void main()
