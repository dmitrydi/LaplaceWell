/*
 * tests.cpp
 *
 *  Created on: 2 мар. 2021 г.
 *      Author: Dmitry_Di
 */

#include "tests.h"

using namespace std;

ostream& operator<<(ostream& os, const PointXYZV& p) {
	os << "(" << p.x << ", " << p.y << ", " << p.z << ", " << p.val << ")";
	return os;
}

void TestMakeGrid() {
	size_t nx = 11;
	size_t ny = 21;
	const double tiny = numeric_limits<double>::epsilon();
	vector<double> xs = LinSpaced(0.,10.,nx);
	vector<double> ys = LinSpaced(0.,20.,ny);
	auto grid = MakeGrid(xs, ys);
	ASSERT_EQUAL(grid.size(), ny);
	ASSERT_EQUAL(grid[0].size(), nx);
	for (size_t row = 0; row < ny; ++ row) {
		for (size_t col = 0; col < nx; ++col) {
			const auto& p = grid[row][col];
			ASSERT(abs(p.x-xs[col]) <= tiny);
			ASSERT(abs(p.y-ys[row]) <= tiny);
			ASSERT(abs(p.z-0.) <= tiny);
			ASSERT(abs(p.val - 0.) <= tiny);
		}
	}
}

void TestNPaginate() {
	TestNPaginateXY(100, 1861, 4);
	TestNPaginateXY(100, 1000, 3);
	TestNPaginateXY(100, 1000, 4);
}

void TestNPaginateXY(size_t nx, size_t ny, size_t nthread) {
	const double tiny = numeric_limits<double>::epsilon();
	vector<double> xs = LinSpaced(0.,10.,nx);
	vector<double> ys = LinSpaced(0.,20.,ny);
	auto grid = MakeGrid(xs, ys);
	auto pages = NPaginate(grid, nthread);
	// Check num of pages
	ASSERT_EQUAL(pages.size(), nthread);
	//
	int cntr = 0;
	auto func = [&cntr](IteratorRange<MatrixXYZV::iterator>& page){
		for (auto row_it = page.begin(); row_it != page.end(); ++row_it) {
			for (auto c_it = row_it->begin(); c_it != row_it->end(); ++c_it) {
				++cntr;
				c_it->val  = 1.;
			}
		}
	};
	for (auto p: pages) {
		func(p);
	}
	// Chesk all elemants were visited by lambda
	ASSERT_EQUAL(cntr, nx*ny);
	for (size_t j = 0; j < ny; ++j) {
		for (size_t i = 0; i < nx; ++i) {
			ASSERT(abs(grid[j][i].val - 1.) <= tiny);
		}
	}
	//
}

void TestParallelDum() {
	vector<double> xs = LinSpaced(0.,10.,1000);
	vector<double> ys = LinSpaced(0.,10.,5000);
	auto grid = MakeGrid(xs, ys);
	{ LOG_DURATION("single");
		_TestParallelDum(grid, 1);
	}
	{LOG_DURATION("multi");
		_TestParallelDum(grid, 4);

	}
	cout << grid[0][0] << endl;
}

void _TestParallelDum(MatrixXYZV& grid, int nthread) {
	vector<future<void>> futures;
	auto pages = NPaginate(grid, nthread);
	for (auto p: pages) {
		futures.push_back(async(launch::async, [](auto p) {
			for (auto& row: p) {
				for (auto& el : row) {
					el.val = exp(cos(sin(el.x)+sin(el.y)));
				}
			}
		}, p));
	}
}


