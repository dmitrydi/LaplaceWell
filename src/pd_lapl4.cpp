//============================================================================
// Name        : pd_lapl4.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include"gwell.h"

using namespace std;
using namespace Rectangular;

int main() {
	double xed = 10.;
	double xwd = 5.;
	double yed = 10.;
	double ywd = 5.;
	double Fcd = 3.14;
	Boundary boundary = Boundary::NNNN;
	Fracture frac(boundary, xwd, xed,
			 ywd, yed,
			Fcd);
	cout << frac.pwd(1000.) << " ";
	//cout << Inverse(frac, pwd, 1000.);
}
