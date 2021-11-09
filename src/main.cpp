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
	RUN_TEST(tr, TestCinco);
	RUN_TEST(tr, TestPdLaplParallel);
	RUN_TEST(tr, ShowPdXY);
}
