/*
 * matrix3dv.cpp
 *
 *  Created on: 23 апр. 2021 г.
 *      Author: Dmitry_Di
 */

#include "matrix3dv.h"

using namespace std;

MatrixDimentions Matrix3DV::GetDimentions() const {
	return {nx, ny, nz};
}

void Matrix3DV::UnsafeResize(size_t nx, size_t ny, size_t nz) {
	this->nx = nx;
	this->ny = ny;
	this->nz = nz;
	this->n = nx*ny*nz;
	this->data.resize(this->n);
}

std::ofstream& operator<<(std::ofstream& of, const Matrix3DV& m) {
	auto dims = m.GetDimentions();
	of << dims.nx << " " << dims.ny << " " << dims.nz << "\n";
	for (size_t k = 0; k < dims.nz; ++k) {
		for (size_t i = 0; i < dims.nx; ++ i) {
			for (size_t j = 0; j < dims.ny; ++j) {
				of << m(i,j,k).x << " " << m(i,j,k).y << " " << m(i,j,k).z << " " << m(i,j,k).val << "\n";
			}
		}
	}
	return of;
}

std::ifstream& operator>>(std::ifstream& ifs, Matrix3DV& m) {
	string line;
	getline(ifs, line);
	istringstream iss(line);
	size_t nx, ny, nz;
	iss >> nx >> ny >> nz;
	assert(iss.eof());
	m.UnsafeResize(nx, ny, nz);
	for (size_t k = 0; k < nz; ++k) {
		for (size_t i = 0; i < nx; ++ i) {
			for (size_t j = 0; j < ny; ++j) {
				getline(ifs, line);
				istringstream _iss(line);
				_iss >> m(i,j,k).x >> m(i,j,k).y >> m(i,j,k).z >> m(i,j,k).val;
			}
		}
	}
	return ifs;
}

ostream& operator <<(ostream& os, const PointXYZV& p) {
	os << "(" << p.x << ", " << p.y << ", " << p.z << ", " << p.val << ")";
	return os;
}

ostream& operator<<(ostream& os, const Matrix3DV& m) {
	for (size_t k = 0; k < m.nz; ++k) {
		for (size_t i = 0; i < m.nx; ++i) {
			for (size_t j = 0; j < m.ny; ++j) {
				os << m(i,j,k) << " ";
			}
			os << endl;
		}
	}
	return os;
}

Matrix3DV::Matrix3DV(size_t nx, size_t ny, size_t nz): nx(nx),
		ny(ny), nz(nz), n(nx*ny*nz), data(n) {

}

Matrix3DV::Matrix3DV(size_t nx, size_t ny, size_t nz, double val): nx(nx),
		ny(ny), nz(nz), n(nx*ny*nz), data(n, {0,0,0,val}) {

}

Matrix3DV::Matrix3DV(): nx(0), ny(0), nz(0), n(0), data(0) {};

std::vector<PointXYZV>::iterator Matrix3DV::begin() {
	return data.begin();
};

std::vector<PointXYZV>::iterator Matrix3DV::end() {
	return data.end();
};
PointXYZV& Matrix3DV::operator()(size_t i, size_t j, size_t k) {
	return data[ny*i + j + nx*ny*k];
};

const PointXYZV& Matrix3DV::operator()(size_t i, size_t j, size_t k) const {
	return data[ny*i + j + nx*ny*k];
}

void Matrix3DV::AddVals(const Matrix3DV& other) {
	assert (n == other.n);
	for (size_t i = 0; i < n; ++i)
		data[i].val += other.data[i].val;
};
void Matrix3DV::AddVals(Matrix3DV&& other) {
	assert (n == other.n);
	for (size_t i = 0; i < n; ++i)
		data[i].val += other.data[i].val;
};
void Matrix3DV::MultVals(const double d) {
	for (size_t i = 0; i < n; ++i)
		data[i].val *= d;
};
void Matrix3DV::DivVals(const double d) {
	for (size_t i = 0; i < n; ++i)
		data[i].val /= d;
};

size_t Matrix3DV::size() const {
	return n;
}

std::vector<PointXYZV> Matrix3DV::get() {
	return data;
}

FieldGetter::FieldGetter(const MatrixField f): field(f) {};
double FieldGetter::operator()(const PointXYZV& p) const {
	switch (field) {
	case MatrixField::X:
		return p.x;
		break;
	case MatrixField::Y:
		return p.y;
		break;
	case MatrixField::Z:
		return p.z;
		break;
	case MatrixField::Val:
		return p.val;
		break;
	default:
		throw std::invalid_argument("unknown field in FieldGetter");
	}
}

AxisGetter::AxisGetter(MatrixAxis a): axis(a) {};
double AxisGetter::operator ()(const PointXYZV& p) const {
	switch (axis) {
	case MatrixAxis::X:
		return p.x;
	case MatrixAxis::Y:
		return p.y;
	case MatrixAxis::Z:
		return p.z;
	default:
		throw std::invalid_argument("unknown axis in AxisGetter");
	}
}

MatrixD Matrix3DV::GetField(MatrixField f) const  {
	FieldGetter getter(f);
	size_t nrows = nx*nz;
	MatrixD ans;
	vector<double> row(ny);
	size_t cntr = 0;
	for (size_t i = 0; i < nrows; ++i) {
		for (size_t j = 0; j < ny; ++j) {
			row[j] = getter(data[cntr++]);
		}
		ans.push_back(row);
	}
	return ans;
}

size_t Matrix3DV::GetSizeOnAxis(MatrixAxis axis) const {
	switch (axis) {
	case MatrixAxis::X:
		return nx;
	case MatrixAxis::Y:
		return ny;
	case MatrixAxis::Z:
		return nz;
	default:
		throw std::invalid_argument("Matrix3DV::GetSizeOnAxis: unknown axis\n");
	}
}

VectorD Matrix3DV::GetAxis(MatrixAxis axis) const {
	VectorD ans;
	AxisGetter getter(axis);
	size_t imax, step;
	switch (axis) {
	case MatrixAxis::X:
		imax = nx*ny;
		step = ny;
		break;
	case MatrixAxis::Y:
		imax = ny;
		step = 1;
		break;
	case MatrixAxis::Z:
		imax = nx*ny*nz;
		step = nx*ny;
		break;
	default:
		throw std::invalid_argument("unknown axis in Matrix3DV::GetAxis");
	}
	for (size_t i = 0; i < imax; i+=step) {
		ans.push_back(getter(data[i]));
	}
	return ans;
}


std::ostream& operator <<(std::ostream& os, const MatrixD& m) {
	size_t nrows = m.size();
	size_t ncols = m[0].size();
	for (size_t i = 0; i < nrows; ++i) {
		for (size_t j = 0; j < ncols; ++j) {
			os << m[i][j] << " ";
		}
		os << endl;
	}
	return os;
}

Matrix3DV MakeGrid(const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs) {
	size_t nx = xs.size(), ny = ys.size(), nz = zs.size();
	Matrix3DV grid(nx, ny, nz);
	for (size_t k = 0; k < nz; ++k) {
		for (size_t i = 0; i < nx; ++i) {
			for (size_t j = 0; j < ny; ++j) {
				grid(i,j,k).x = xs[i];
				grid(i,j,k).y = ys[j];
				grid(i,j,k).z = zs[k];
				grid(i,j,k).val = 0.;
			}
		}
	}
	return grid;
}



