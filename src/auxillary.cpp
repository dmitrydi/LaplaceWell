/*
 * auxillary.cpp
 *
 *  Created on: 2 мар. 2021 г.
 *      Author: Dmitry_Di
 */


#include "auxillary.h"

using namespace std;

MatrixXYZV MakeGrid(const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs) {
	size_t nx = xs.size();
	size_t ny = ys.size();
	size_t nz = zs.size();
	MatrixXYZV m(ny*nz, VectorXYZV(nx));
	for (size_t k = 0; k < nz; ++k)
		for (size_t j = 0; j < ny; ++j)
			for(size_t i = 0; i < nx; ++i)
				m[k*nz+j][i] = PointXYZV(xs[i], ys[j], zs[k]);
	return m;
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

};

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

Matrix3DV MakeGrid2(const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs) {
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

