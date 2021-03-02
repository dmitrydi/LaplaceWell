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
//	double xed = 10.;
//	double xwd = 5.;
//	double yed = 10.;
//	double ywd = 5.;
//	double Fcd = 3.14;
//	Boundary boundary = Boundary::NNNN;
//	Fracture frac(boundary, xwd, xed,
//			 ywd, yed,
//			Fcd);
//	int N = 150;
//	vector<double> tds = LogSpaced(0.001, 1000. , N);
//	vector<double> ans_single(tds.size());
//	vector<double> ans_multi(tds.size());
//	{
//		LOG_DURATION("single");
//		for (int i = 0; i < N; ++i) {
//			ans_single[i] = frac.pwd(tds[i]);
//		}
//	}
//	{
//		LOG_DURATION("multi");
//		frac.pwd_parallel(tds, ans_multi, 4);
//	}
//	for (int i = 0; i < N; ++i) {
//		cout << ans_single[i] << " " << ans_multi[i] << endl;
//	}
	TestRunner tr;
	//RUN_TEST(tr, TestMakeGrid);
	//RUN_TEST(tr, TestNPaginate);
	//RUN_TEST(tr, TestParallelDum);
	RUN_TEST(tr, TestPDXY);

}
