//============================================================================
// Name        : pd_lapl4.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <future>
#include <vector>
#include <iterator>
#include <cmath>
#include"gwell.h"
#include "profile.h"
#include "auxillary.h"
#include "tests.h"

using namespace std;
using namespace Rectangular;

int main() {
	TestRunner tr;
	//RUN_TEST(tr, TestMakeGrid);
	//RUN_TEST(tr, TestNPaginate);
	//RUN_TEST(tr, TestParallelDum);
	//RUN_TEST(tr, TestPDXY);
	//RUN_TEST(tr, TestMatrix3DV);
	//RUN_TEST(tr, Test_pd_lapl_m);
	//RUN_TEST(tr, TestPDXY_old_new);
	RUN_TEST(tr, TestCincoNew);
}
