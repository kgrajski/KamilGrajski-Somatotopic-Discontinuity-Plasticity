//	expmgr.h

//	Kamil A. Grajski
//	5-November-2014
//

//
//	Header file for some i/o utilities.
//

#ifndef IOUTIL_H
#define IOUTIL_H

#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

void LoadNumericVectorFromBinaryFile ( const string &, fstream &, const long &, vector<double> &, const long & );

void SaveNumericVectorToBinaryFile ( const string &, const long &, vector<double> &, const long & , const long & );

#endif