#include "tests.h"

using namespace std;
using namespace Rectangular;

void TestCinco() {
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

void TestPdLaplParallel() {
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
	Matrix3DV grid = MakeGrid(xds, yds);
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


void ShowPdXY() {
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
	Matrix3DV grid_pd_new;

	//{ LOG_DURATION("new");
		grid_pd_new = frac.pd_m_parallel(1., 4, xds, yds);
	//}
	VectorD xgrid = grid_pd_new.GetAxis(MatrixAxis::X);
	VectorD ygrid = grid_pd_new.GetAxis(MatrixAxis::Y);
	MatrixD vals = grid_pd_new.GetField(MatrixField::Val);
	for (auto x: xgrid) cout << x << " ";
	cout << endl;
	for (auto y: ygrid) cout << y << " ";
	cout << endl;
	cout << vals;
	cout << xgrid.size() << endl;
	cout << ygrid.size() << endl;
	ofstream of("out.txt");
	of << grid_pd_new;
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
