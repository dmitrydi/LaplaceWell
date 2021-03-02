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

