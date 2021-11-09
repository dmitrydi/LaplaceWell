/*
 * matrix3dv.h
 *
 *  Created on: 23 ���. 2021 �.
 *      Author: Dmitry_Di
 */

#ifndef MATRIX3DV_H_
#define MATRIX3DV_H_

#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <sstream>

struct PointXYZV {
	PointXYZV(): x(0), y(0), z(0), val(0) {};
	PointXYZV(double x, double y, double z): x(x), y(y), z(z), val(0) {};
	PointXYZV(double x, double y, double z, double val): x(x), y(y), z(z), val(val) {};
	double x, y, z, val;
};

typedef std::vector<PointXYZV> VectorXYZV;
typedef std::vector<std::vector<double>> MatrixD;
typedef std::vector<double> VectorD;

std::ostream& operator <<(std::ostream& os, const PointXYZV& p);

std::ostream& operator <<(std::ostream& os, const MatrixD& m);

enum class MatrixField {
	X,
	Y,
	Z,
	Val
};

enum class MatrixAxis {
	X,
	Y,
	Z
};

class FieldGetter {
public:
	FieldGetter(const MatrixField f);
	double operator()(const PointXYZV& p) const;
private:
	MatrixField field;
};

class AxisGetter {
public:
	AxisGetter(const MatrixAxis a);
	double operator ()(const PointXYZV& p) const;
private:
	MatrixAxis axis;
};

struct MatrixDimentions {
	size_t nx;
	size_t ny;
	size_t nz;
};

class Matrix3DV {
public:
	Matrix3DV(size_t nx, size_t ny, size_t nz = 1);
	Matrix3DV(size_t nx, size_t ny, size_t nz, double val);
	Matrix3DV();
	void UnsafeResize(size_t nx, size_t ny, size_t nz);
	std::vector<PointXYZV>::iterator begin();
	std::vector<PointXYZV>::iterator end();
	PointXYZV& operator()(size_t i, size_t j, size_t k=0);
	const PointXYZV& operator()(size_t i, size_t j, size_t k=0) const;
	size_t GetSizeOnAxis(MatrixAxis axis) const;
	void AddVals(const Matrix3DV& other);
	void AddVals(Matrix3DV&& other);
	void MultVals(const double d);
	void DivVals(const double d);
	size_t size() const;
	std::vector<PointXYZV> get();
	MatrixD GetField(MatrixField f) const;
	VectorD GetAxis(MatrixAxis a) const;
	MatrixDimentions GetDimentions() const;
private:
	size_t nx, ny, nz, n;
	std::vector<PointXYZV> data;
	friend std::ostream& operator<<(std::ostream& os, const Matrix3DV&);
};

std::ofstream& operator<<(std::ofstream& of, const Matrix3DV&);
std::ifstream& operator>>(std::ifstream& ifs, Matrix3DV&);

Matrix3DV MakeGrid(const std::vector<double>& xs, const std::vector<double>& ys, const std::vector<double>& zs = {0.});

#endif /* MATRIX3DV_H_ */
