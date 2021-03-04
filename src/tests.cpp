/*
 * tests.cpp
 *
 *  Created on: 2 мар. 2021 г.
 *      Author: Dmitry_Di
 */

#include "tests.h"

using namespace std;
using namespace Rectangular;

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
	vector<double> xs = LinSpaced(0.,10.,11);
	vector<double> ys = LinSpaced(0.,10.,11);
	auto grid = _TestParallelDum(xs, ys, 4);
	for (auto& row: grid) {
		for (auto el: row) {
			cout << el << " ";
		}
		cout << endl;
	}
}

MatrixXYZV _TestParallelDum(const vector<double>& xs, const vector<double>& ys, int nthread) {
	MatrixXYZV grid = MakeGrid(xs, ys);
	vector<future<void>> futures;
	auto pages = NPaginate(grid, nthread);
	for (auto p: pages) {
		futures.push_back(async(launch::async, [](auto p) {
			for (auto& row: p) {
				for (auto& el : row) {
					el.val = el.x+el.y;
				}
			}
		}, p));
	}
	return grid;
}

void TestPDXY() {
	double xed = 10.;
	double xwd = 5.;
	double yed = 10.;
	double ywd = 5.;
	double Fcd = 3;
	Boundary boundary = Boundary::NNNN;
	Fracture frac(boundary, xwd, xed,
			 ywd, yed,
			Fcd);
	vector<double> xds = LinSpaced(0.1, 3.9, 10);
	for (auto x: LinSpaced(4., 6., 41)) {
		xds.push_back(x);
	}
	for (auto x: LinSpaced(6.1, 9.9, 10)) {
		xds.push_back(x);
	}
	vector<double> yds = LinSpaced(4., 6., 21);
	for (int i = 0; i < xds.size(); ++i) {
		cout << xds[i] << " ";
	}
	cout << endl << endl;

	for (int j = 0; j < yds.size(); ++j) {
		cout << yds[j] << endl;
	}
	cout << endl << endl;
	LOG_DURATION("four");

	auto grid_pd = frac.pd_parallel(1., 4, xds, yds);
	for (int j = 0; j < yds.size(); ++j) {
		for (int i = 0; i < xds.size(); ++i) {
			cout << grid_pd[j][i].val << " ";
		}
		cout << endl;
	}
}

void TestPDXY_old_new() {
	double xed = 10.;
	double xwd = 5.;
	double yed = 10.;
	double ywd = 5.;
	double Fcd = 3;
	Boundary boundary = Boundary::NNNN;
	Fracture frac(boundary, xwd, xed,
			 ywd, yed,
			Fcd);
	vector<double> xds = LinSpaced(0.1, 3.9, 10);
	for (auto x: LinSpaced(4., 6., 41)) {
		xds.push_back(x);
	}
	for (auto x: LinSpaced(6.1, 9.9, 10)) {
		xds.push_back(x);
	}
	vector<double> yds = LinSpaced(4., 6., 21);
	MatrixXYZV grid_pd_old;
	Matrix3DV grid_pd_new;
	{ LOG_DURATION("old");
	grid_pd_old = frac.pd_parallel(1., 4, xds, yds);
	}
	{ LOG_DURATION("new");
	grid_pd_new = frac.pd_m_parallel(1., 4, xds, yds);
	}
	for (int i = 0; i < xds.size(); ++i) {
		for (int j =0; j <yds.size(); ++j) {
			cout << grid_pd_old[j][i] << " " << grid_pd_new(i,j,0) << " " << grid_pd_old[j][i].val - grid_pd_new(i,j,0).val << endl;
		}
	}
}

void TestMatrix3DV() {
	size_t nx = 5;
	size_t ny = 10;
	size_t nz = 2;
	Matrix3DV m(nx, ny, nz);
	for (size_t k = 0; k < nz; ++k)
	for (size_t i = 0; i < nx; ++i) {
		for (size_t j = 0; j < ny; ++j) {
			m(i,j,k) = PointXYZV(i+1,j+1,k+1,(i+1)*(j+1)*(k+1));
		}
	}
	Matrix3DV m2(nx,ny,nz, 0.1);
	m.AddVals(m2);
	m.MultVals(0.1);
	m.DivVals(1);
	cout << m << endl;
}


void ProcessPage(IteratorRange<vector<PointXYZV>::iterator> page) {
	//this_thread::sleep_for(chrono::seconds(1));
	for (auto& p: page) {
		p.val = p.x+p.y;
	}
}

void TestMatrix3DV_Pages() {
	size_t nx = 500000;
	size_t ny = 100;
	size_t nz = 1;
	Matrix3DV m(nx, ny, nz);
	for (size_t k = 0; k < nz; ++k) {
		for (size_t i = 0; i < nx; ++i) {
			for (size_t j = 0; j < ny; ++j) {
				m(i,j,k).x = i+1;
				m(i,j,k).y = j+1;
				m(i,j,k).z = k+1;
				m(i,j,k).val = 0;
			}
		}
	}

	{
		LOG_DURATION("4");
		vector<future<void>> futures;
		auto pages = NPaginate(m, 4);
		for (auto& page: pages) {
			futures.push_back(async(launch::async, ProcessPage, ref(page)));
		}
		for (int i = 0; i < 2; ++i) {
			futures[i].get();
		}
		cout << m(0,0,0) << endl;
//		for (auto p: m.get()) {
//			ASSERT(abs(p.val-(p.x+p.y)) < 1e-16);
//		}
	}
	{
		LOG_DURATION("1");
		vector<future<void>> futures;
			auto pages = NPaginate(m, 1);
			for (auto& page: pages) {
				futures.push_back(async(launch::async, ProcessPage, ref(page)));
			}
			for (int i = 0; i < 1; ++i) {
				futures[i].get();
			}
			cout << m(0,0,0) << endl;
//			for (auto p: m.get()) {
//				ASSERT(abs(p.val-(p.x+p.y)) < 1e-16);
//			}
	}
}

void Test_pd_lapl_m() {
	double xed = 10.;
	double xwd = 5.;
	double yed = 10.;
	double ywd = 5.;
	double Fcd = 3;
	double u = 1.;
	Boundary boundary = Boundary::NNNN;
	Fracture frac(boundary, xwd, xed,
			 ywd, yed,
			Fcd);
	vector<double> xds = LinSpaced(0.1, 9.9, 100);
	vector<double> yds = xds;
	vector<vector<double>> pd_lapl_old;
	{ LOG_DURATION("old");
	for (auto x: xds) {
		vector<double> row;
		for (auto y: yds) {
			row.push_back(frac.pd_lapl(u, x, y));
		}
		pd_lapl_old.push_back(row);
	}
	}
	Matrix3DV grid = MakeGrid2(xds, yds);
	Matrix3DV pd_lapl_new;
	{LOG_DURATION("new");
		pd_lapl_new = frac.pd_lapl_m(u, grid, 4);
	}
	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < 10; ++j) {
			ASSERT(abs(pd_lapl_old[i][j] - pd_lapl_new(i,j,0).val) < numeric_limits<double>::epsilon());
		}
	}
}

void TestCincoNew() {
	const double PI = 3.141592653589793;
	const string datafile = "./test_data/cinco_infinite_reservoir.csv";
	const vector<double> Fcds = {PI*0.2, PI, 2.*PI, 10.*PI, 20.*PI, 100.*PI};
	ifstream dfile(datafile);
	vector<double> tds;
	vector<vector<double>> pwd_expected(Fcds.size());
	vector<vector<double>> calculated(Fcds.size());
	if (dfile) {
		string line;
		while (getline(dfile, line)) {
			if (line[0] == '/') continue;
			istringstream is(line);
			string s;
			getline(is, s, ';');
			tds.push_back(stod(s));
			for (size_t i = 0; i < Fcds.size(); ++i) {
				getline(is, s, ';');
				pwd_expected[i].push_back(stod(s));
			}
		}
	}
	const double xed = 100.;
	const double yed = 100.;
	const double xwd = xed*0.5;
	const double ywd = yed*0.5;
	for (size_t i = 0; i < Fcds.size(); ++i) {
		double Fcd = Fcds[i];
		Boundary boundary = Boundary::NNNN;
		Fracture well(boundary, xwd, xed,
					 ywd, yed,
					Fcd);
		for (size_t j = 0; j < tds.size(); ++j) {
			double td  = tds[j];
			double pwd = well.pwd(td);
			double pwd_e = pwd_expected[i][j];
			double eps = abs(pwd - pwd_e)/pwd_e;
			cout << "Fcd: " << Fcd << " td: " << td << " pwd: " << pwd << " expected: " << pwd_e << " eps: " << eps << endl;
			calculated[i].push_back(pwd);
		}
	}

	ofstream out("./test_data/cinco_test_out_new.csv");
	for (size_t i = 0; i < tds.size(); ++i) {
		out << tds[i] << ';';
		for (size_t j = 0; j < Fcds.size(); ++j) {
			out << calculated[j][i];
			if (j < Fcds.size() - 1) out << ';';
		}
		out << '\n';
	}
}


